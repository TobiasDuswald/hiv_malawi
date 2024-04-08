// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN and the University of Geneva for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
//
// -----------------------------------------------------------------------------

#ifndef PERSON_BEHAVIOR_H_
#define PERSON_BEHAVIOR_H_

#include "categorical-environment.h"
#include "datatypes.h"
#include "person.h"
#include "population-initialization.h"

namespace bdm {
namespace hiv_malawi {

////////////////////////////////////////////////////////////////////////////////
// BioDynaMo's Agent / Individual Behaviors
////////////////////////////////////////////////////////////////////////////////

// A behavior that allows agents to randomly migrate between the categorical
// loactions. It uses a gausion random process to determine the next location.
struct RandomMigration : public Behavior {
  BDM_BEHAVIOR_HEADER(RandomMigration, Behavior, 1);

  RandomMigration() {}

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* env = bdm_static_cast<CategoricalEnvironment*>(sim->GetEnvironment());
    auto* random = sim->GetRandom();
    auto* person = bdm_static_cast<Person*>(agent);
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();

    // Probability to migrate
    float rand_num = static_cast<float>(random->Uniform());
    // Adult men and adult single women can initiate migration
    if (rand_num <= sparam->migration_probability && person->age_ >= 15 &&
        ((person->sex_ == Sex::kMale) ||
         (person->sex_ == Sex::kFemale && !person->hasPartner()))) {
      // Randomly determine the migration location
      // AM: Sample migration location. It depends on the current year and
      // current location
      float rand_num_loc = static_cast<float>(random->Uniform());
      // Get (cumulative) probability distribution that agent relocates the
      // current year, to each location
      const auto& migration_location_distribution_ =
          env->GetMigrationLocDistribution(person->location_);

      int new_location =
          SampleLocation(rand_num_loc, migration_location_distribution_);

      person->Relocate(new_location);
    }
  }
};

// This is the mating and therefore also the infection behaviour. The Behavior
// is only executed by male agents. For each male agent, we determine the
// number of sex partners with a gaussian random process and randomly select
// from all available females at his location. For each contact, HIV (if
// present) can be transmitted from the infected to the healthy individual with
// a certain probability.
struct MatingBehaviour : public Behavior {
  BDM_BEHAVIOR_HEADER(MatingBehaviour, Behavior, 1);

  MatingBehaviour() {}

  int SampleCompoundCategory(float rand_num,
                             const std::vector<float>& category_distribution) {
    for (size_t i = 0; i < category_distribution.size(); i++) {
      if (rand_num <= category_distribution[i]) {
        return i;
      }
    }

    // This line of code should never be reached
    Log::Warning("SampleCompoundCategory()",
                 "Could not sample the category. Recieved inputs: ", rand_num,
                 ". Use location 0.");
    return 0;
  }

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* env = bdm_static_cast<CategoricalEnvironment*>(sim->GetEnvironment());
    auto* random = sim->GetRandom();
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();
    auto* person = bdm_static_cast<Person*>(agent);

    // Randomly determine the number of mates
    // AM: Mean and standard deviation of the number of mates depend on the
    // current year and socio-behavioural category of agent
    int year = static_cast<int>(
        sparam->start_year +
        sim->GetScheduler()->GetSimulatedSteps());  // Current year
    // If no transition year is higher than current year, then use last
    // transition year
    int year_index = sparam->no_mates_year_transition.size() - 1;
    for (size_t y = 0; y < sparam->no_mates_year_transition.size() - 1; y++) {
      if (year < sparam->no_mates_year_transition[y + 1]) {
        year_index = y;
        break;
      }
    }
    /*int no_mates = static_cast<int>(random->Gaus(
        sparam->no_mates_mean[year_index][person->social_behaviour_factor_],
        sparam->no_mates_sigma[year_index][person->social_behaviour_factor_]));*/

    // Poisson Distribution
    int no_mates = random->Poisson(
        sparam->no_mates_mean[year_index][person->social_behaviour_factor_]);

    // This part is only executed for male persons in a certain age group, since
    // the infection goes into both directions.
    if (no_mates > 0 && person->sex_ == Sex::kMale &&
        person->age_ >= env->GetMinAge() && person->age_ < env->GetMaxAge()) {
      // Compute male agent's age category
      size_t age_category =
          person->GetAgeCategory(env->GetMinAge(), env->GetNoAgeCategories());
      // Get (cumulative) probability distribution that the male agent selects a
      // female mate from each compound category
      const std::vector<float>& mate_compound_category_distribution =
          env->GetMateCompoundCategoryDistribution(
              person->location_, age_category,
              person->social_behaviour_factor_);
      // Reset to 0 for this year
      // person->no_casual_partners_ = 0;

      for (int i = 0; i < no_mates; i++) {
        // AM: select compound category of mate
        float rand_num = static_cast<float>(random->Uniform());

        size_t mate_compound_category = SampleCompoundCategory(
            rand_num, mate_compound_category_distribution);

        // AM: Choose a random female mate at the selected mate compound
        // category (location, age group and sociobehavioral category
        AgentPointer<Person> mate =
            env->GetRandomCasualFemaleFromIndex(mate_compound_category);

        if (mate == nullptr) {
          Log::Fatal("MatingBehaviour()",
                     "Received nullptr as AgentPointer mate.");
        }

        // Increment number of casual partners for both agents
        person->no_casual_partners_ = person->no_casual_partners_ + 1;
        mate->no_casual_partners_ = mate->no_casual_partners_ + 1;

        int no_acts = static_cast<int>(random->Gaus(
            sparam->no_acts_mean[year_index][person->social_behaviour_factor_],
            sparam
                ->no_acts_sigma[year_index][person->social_behaviour_factor_]));

        // Scenario healthy male has intercourse with infected acute female
        if (mate->state_ == GemsState::kAcute &&
            person->state_ == GemsState::kHealthy &&
            mate->newly_infected == false &&
            random->Uniform() <
                (1.0 -
                 pow(1.0 - sparam->infection_probability_acute_fm, no_acts))) {
          person->state_ = GemsState::kAcute;
          person->transmission_type_ = TransmissionType::kCasualPartner;
          person->infection_origin_state_ = mate->state_;
          person->infection_origin_sb_ = mate->social_behaviour_factor_;
          //std::cout << "C"; 
          // AM: Add MatingBehaviour only when male gets infected
          /*person->AddBehavior(new MatingBehaviour());
          std::cout << "This should not currently happen: AddBehavior(new
          MatingBehaviour()) in MatingBehaviour::Run()" << std::endl;*/
        }
        // Scenario healthy male has intercourse with infected chronic female
        else if (mate->state_ == GemsState::kChronic &&
                 person->state_ == GemsState::kHealthy &&
                 mate->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_chronic_fm,
                                no_acts))) {
          person->state_ = GemsState::kAcute;
          person->transmission_type_ = TransmissionType::kCasualPartner;
          person->infection_origin_state_ = mate->state_;
          person->infection_origin_sb_ = mate->social_behaviour_factor_;
          //std::cout << "c";
          // AM: Add MatingBehaviour only when male gets infected
          /*person->AddBehavior(new MatingBehaviour());
          std::cout << "This should not currently happen: AddBehavior(new
          MatingBehaviour()) in MatingBehaviour::Run()" << std::endl;*/
        }
        // Scenario healthy male has intercourse with infected treated female
        else if (mate->state_ == GemsState::kTreated &&
                 person->state_ == GemsState::kHealthy &&
                 mate->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_treated_fm,
                                no_acts))) {
          person->state_ = GemsState::kAcute;
          person->transmission_type_ = TransmissionType::kCasualPartner;
          person->infection_origin_state_ = mate->state_;
          person->infection_origin_sb_ = mate->social_behaviour_factor_;
          // AM: Add MatingBehaviour only when male gets infected
          /*person->AddBehavior(new MatingBehaviour());
          std::cout << "This should not currently happen: AddBehavior(new
          MatingBehaviour()) in MatingBehaviour::Run()" << std::endl;*/
        }
        // Scenario healthy male has intercourse with infected failing treatment
        // female
        else if (mate->state_ == GemsState::kFailing &&
                 person->state_ == GemsState::kHealthy &&
                 mate->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_failing_fm,
                                no_acts))) {
          person->state_ = GemsState::kAcute;
          person->transmission_type_ = TransmissionType::kCasualPartner;
          person->infection_origin_state_ = mate->state_;
          person->infection_origin_sb_ = mate->social_behaviour_factor_;
          // AM: Add MatingBehaviour only when male gets infected
          /*person->AddBehavior(new MatingBehaviour());
          std::cout << "This should not currently happen: AddBehavior(new
          MatingBehaviour()) in MatingBehaviour::Run()" << std::endl;*/
        }
        // Scenario infected acute male has intercourse with healthy female
        else if (mate->state_ == GemsState::kHealthy &&
                 person->state_ == GemsState::kAcute &&
                 person->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_acute_mf,
                                no_acts))) {
          mate->state_ = GemsState::kAcute;
          mate->transmission_type_ = TransmissionType::kCasualPartner;
          mate->infection_origin_state_ = person->state_;
          mate->infection_origin_sb_ = person->social_behaviour_factor_;
          //std::cout << "C";
        }  // Scenario infected chronic male has intercourse with healthy female
        else if (mate->state_ == GemsState::kHealthy &&
                 person->state_ == GemsState::kChronic &&
                 person->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_chronic_mf,
                                no_acts))) {
          mate->state_ = GemsState::kAcute;
          mate->transmission_type_ = TransmissionType::kCasualPartner;
          mate->infection_origin_state_ = person->state_;
          mate->infection_origin_sb_ = person->social_behaviour_factor_;
          //std::cout << "c";
        }  // Scenario infected treated male has intercourse with healthy female
        else if (mate->state_ == GemsState::kHealthy &&
                 person->state_ == GemsState::kTreated &&
                 person->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_treated_mf,
                                no_acts))) {
          mate->state_ = GemsState::kAcute;
          mate->transmission_type_ = TransmissionType::kCasualPartner;
          mate->infection_origin_state_ = person->state_;
          mate->infection_origin_sb_ = person->social_behaviour_factor_;
        }  // Scenario infected failing treatment male has intercourse with
           // healthy female
        else if (mate->state_ == GemsState::kHealthy &&
                 person->state_ == GemsState::kFailing &&
                 person->newly_infected == false &&
                 random->Uniform() <
                     (1.0 - pow(1.0 - sparam->infection_probability_failing_mf,
                                no_acts))) {
          mate->state_ = GemsState::kAcute;
          mate->transmission_type_ = TransmissionType::kCasualPartner;
          mate->infection_origin_state_ = person->state_;
          mate->infection_origin_sb_ = person->social_behaviour_factor_;
        } else {
          ;  // if both are infected or both are healthy, do nothing
        }
      }
    }
  }
};

// This is the regular partnership behaviour. The Behavior
// is only executed by male agents. For each single male agent, we determine
// whether he wants to engage in a regular partnership with a certain
// probability. For each male in regular partnership, we determine if he will
// separate from partner.
struct RegularPartnershipBehaviour : public Behavior {
  BDM_BEHAVIOR_HEADER(RegularPartnershipBehaviour, Behavior, 1);

  RegularPartnershipBehaviour() {}

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();
    auto* person = bdm_static_cast<Person*>(agent);

    // Adult men in regular partnership can break up (symmetric for female)
    if (person->IsAdult() && person->hasPartner() &&
        random->Uniform() <= sparam->break_up_probability) {
      // Set female partner to single
      person->partner_->partner_ = nullptr;
      // Set male agent to single
      person->partner_ = nullptr;
    }

    // Adult single men can decide to engage in a regular partnership
    if (person->IsAdult() && !person->hasPartner() &&
        random->Uniform() <= sparam->regular_partnership_probability) {
      person->seek_regular_partnership_ = true;
    } else {
      person->seek_regular_partnership_ = false;
    }
  }
};

// This is the regular mating behaviour. The Behavior
// is only executed by male agents. In serodiscordant regular partners,
// the infected partner transmits HIV to his healthy partner with a certain
// probability.
struct RegularMatingBehaviour : public Behavior {
  BDM_BEHAVIOR_HEADER(RegularMatingBehaviour, Behavior, 1);

  RegularMatingBehaviour() {}

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* env = bdm_static_cast<CategoricalEnvironment*>(sim->GetEnvironment());
    auto* random = sim->GetRandom();
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();
    auto* person = bdm_static_cast<Person*>(agent);

    // Determine the number of regular acts
    // AM: Number of regular acts depends on the current year
    int year = static_cast<int>(
        sparam->start_year +
        sim->GetScheduler()->GetSimulatedSteps());  // Current year
    // If no transition year is higher than current year, then use last
    // transition year
    int year_index = sparam->no_regacts_year_transition.size() - 1;
    for (int y = 0; y < sparam->no_regacts_year_transition.size() - 1; y++) {
      if (year < sparam->no_regacts_year_transition[y + 1]) {
        year_index = y;
        break;
      }
    }

    if (person->hasPartner() && person->age_ < env->GetMaxAge()) {
      // Scenario healthy male has intercourse with infected acute female
      // partner
      if (person->partner_->state_ == GemsState::kAcute &&
          person->state_ == GemsState::kHealthy &&
        //  person->partner_->newly_infected == false &&
          random->Uniform() <
              (1.0 - pow(1.0 - sparam->infection_probability_acute_fm,
                         sparam->no_regular_acts_mean[year_index]))) {
        person->state_ = GemsState::kAcute;
        person->newly_infected = true;
        person->transmission_type_ = TransmissionType::kRegularPartner;
        person->infection_origin_state_ = person->partner_->state_;
        //std::cout << "R";
        // AM: Add MatingBehaviour only when infected
        // person->AddBehavior(new MatingBehaviour());
      }
      // Scenario healthy male has intercourse with infected chronic female
      // partner
      else if (person->partner_->state_ == GemsState::kChronic &&
               person->state_ == GemsState::kHealthy &&
               //person->partner_->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_chronic_fm,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->state_ = GemsState::kAcute;
        person->newly_infected = true;
        person->transmission_type_ = TransmissionType::kRegularPartner;
        person->infection_origin_state_ = person->partner_->state_;
        //std::cout << "r";
        // AM: Add MatingBehaviour only when infected
        // person->AddBehavior(new MatingBehaviour());
      }
      // Scenario healthy male has intercourse with infected treated female
      // partner
      else if (person->partner_->state_ == GemsState::kTreated &&
               person->state_ == GemsState::kHealthy &&
               //person->partner_->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_treated_fm,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->state_ = GemsState::kAcute;
        person->newly_infected = true;
        person->transmission_type_ = TransmissionType::kRegularPartner;
        person->infection_origin_state_ = person->partner_->state_;
        // AM: Add MatingBehaviour only when infected
        // person->AddBehavior(new MatingBehaviour());
      }
      // Scenario healthy male has intercourse with infected failing treatment
      // female partner
      else if (person->partner_->state_ == GemsState::kFailing &&
               person->state_ == GemsState::kHealthy &&
               //person->partner_->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_failing_fm,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->state_ = GemsState::kAcute;
        person->newly_infected = true;
        person->transmission_type_ = TransmissionType::kRegularPartner;
        person->infection_origin_state_ = person->partner_->state_;
        // AM: Add MatingBehaviour only when infected
        // person->AddBehavior(new MatingBehaviour());
      }
      // Scenario infected acute male has intercourse with healthy female
      // partner
      else if (person->partner_->state_ == GemsState::kHealthy &&
               person->state_ == GemsState::kAcute &&
               //person->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_acute_mf,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->partner_->state_ = GemsState::kAcute;
        person->partner_->newly_infected = true;
        person->partner_->transmission_type_ =
            TransmissionType::kRegularPartner;
        person->partner_->infection_origin_state_ = person->state_;
        //std::cout << "R";
      }  // Scenario infected chronic male has intercourse with healthy female
         // partner
      else if (person->partner_->state_ == GemsState::kHealthy &&
               person->state_ == GemsState::kChronic &&
               //person->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_chronic_mf,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->partner_->state_ = GemsState::kAcute;
        person->partner_->newly_infected = true;
        person->partner_->transmission_type_ =
            TransmissionType::kRegularPartner;
        person->partner_->infection_origin_state_ = person->state_;
        //std::cout << "r";
      }  // Scenario infected treated male has intercourse with healthy female
         // partner
      else if (person->partner_->state_ == GemsState::kHealthy &&
               person->state_ == GemsState::kTreated &&
               //person->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_treated_mf,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->partner_->state_ = GemsState::kAcute;
        person->partner_->newly_infected = true;
        person->partner_->transmission_type_ =
            TransmissionType::kRegularPartner;
        person->partner_->infection_origin_state_ = person->state_;
      }  // Scenario infected failing treatment male has intercourse with
         // healthy female partner
      else if (person->partner_->state_ == GemsState::kHealthy &&
               person->state_ == GemsState::kFailing &&
               //person->newly_infected == false &&
               random->Uniform() <
                   (1.0 - pow(1.0 - sparam->infection_probability_failing_mf,
                              sparam->no_regular_acts_mean[year_index]))) {
        person->partner_->state_ = GemsState::kAcute;
        person->partner_->newly_infected = true;
        person->partner_->transmission_type_ =
            TransmissionType::kRegularPartner;
        person->partner_->infection_origin_state_ = person->state_;
      } else {
        ;  // if both are infected or both are healthy, do nothing
      }
    }
  }
};

// The GiveBirth behavior is assigned to all female agents. If a female is in a
// certain age range, she can give birth to a child that is located at the same
// place. If she is HIV positive, there is a certain chance to infect the child
// while giving birth.
struct GiveBirth : public Behavior {
  BDM_BEHAVIOR_HEADER(GiveBirth, Behavior, 1);

  GiveBirth() {}

  // Helper function to create a single child
  Person* CreateChild(Random* random_generator, Person* mother,
                      const SimParam* sparam, size_t year) {
    // Create new child
    Person* child = new Person();
    // BioDynaMo API: Add agent (child) to simulation
    Simulation::GetActive()->GetExecutionContext()->AddAgent(child);

    // Assign sex
    child->sex_ =
        SampleSex(random_generator->Uniform(), sparam->probability_male);
    // Assign age - possibly -1 ?
    child->age_ = random_generator->Uniform();
    // Assign location
    child->location_ = mother->location_;
    // Compute risk factors
    child->social_behaviour_factor_ = 0;
    child->biomedical_factor_ = 0;
    // Stores the current GemsState of the child.
    if (mother->state_ == GemsState::kHealthy) {
      child->state_ = GemsState::kHealthy;

      ///! The aguments below are currently either not used or repetitive.
      // // Store the year when the agent got infected
      // child->year_of_infection_ = std::numeric_limits<float>::max();
    }
    ///! The aguments below are currently either not used or repetitive.
    // // year of infection to present year, Question: Ask Lukas how to get
    // iter child->year_of_infection_ = 2000;
    //}
    // AM: birth infection probability depends on whether mother is treated and
    // current year
    else if (mother->state_ == GemsState::kTreated) {
      if (random_generator->Uniform() <
          sparam->birth_infection_probability_treated) {
        child->state_ = GemsState::kAcute;
        child->transmission_type_ = TransmissionType::kMotherToChild;
        child->infection_origin_state_ = mother->state_;
      } else {
        child->state_ = GemsState::kHealthy;
      }
    } else if (year < 2003 ||
               mother->state_ ==
                   GemsState::kFailing) {  // AM: Mother is not healthy and not
                                           // treated
      if (random_generator->Uniform() <
          sparam->birth_infection_probability_untreated) {
        child->state_ = GemsState::kAcute;
        child->transmission_type_ = TransmissionType::kMotherToChild;
        child->infection_origin_state_ = mother->state_;
      } else {
        child->state_ = GemsState::kHealthy;
      }
    } else {
      if (random_generator->Uniform() <
          sparam->birth_infection_probability_prophylaxis) {
        child->state_ = GemsState::kAcute;
        child->transmission_type_ = TransmissionType::kMotherToChild;
        child->infection_origin_state_ = mother->state_;
      } else {
        child->state_ = GemsState::kHealthy;
      }
    }

    // Register child with mother
    mother->AddChild(child->GetAgentPtr<Person>());

    // Assign mother to child. When, the child becomes adult, break the link.
    child->mother_ = mother->GetAgentPtr<Person>();

    ///! The aguments below are currently either not used or repetitive.
    // // NOTE: we do not assign a specific mother or partner at the moment. Use
    // // nullptr instead.
    // child->mother_id_ = nullptr;
    // child->partner_id_ = nullptr;

    // BioDynaMo API: Add the behaviors to the Agent
    child->AddBehavior(new RandomMigration());
    if (child->sex_ == Sex::kFemale) {
      child->AddBehavior(new GiveBirth());
    } else {
      child->AddBehavior(new MatingBehaviour());
      /*if (child->state_ != GemsState::kHealthy){
          child->AddBehavior(new MatingBehaviour());
      }*/
      child->AddBehavior(new RegularMatingBehaviour());
      child->AddBehavior(new RegularPartnershipBehaviour());
    }

    return child;
  }

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();
    auto* mother = bdm_static_cast<Person*>(agent);

    // Skip soon to die agents
    if (mother->will_be_removed_) {
      return;
    }

    // Each potential mother gives birth with a certain probability.
    if (random->Uniform() < sparam->give_birth_probability &&
        mother->age_ <= sparam->max_age_birth &&
        mother->age_ >= sparam->min_age) {
      // The probability of the child to be infected depends on the current year
      // (ex. prophylaxis)
      int year = static_cast<int>(
          sparam->start_year +
          sim->GetScheduler()->GetSimulatedSteps());  // Current year

      // Create a child
      auto* new_child = CreateChild(random, mother, sparam, year);

      // Protect mother from death.
      if (sparam->protect_mothers_at_birth) {
        mother->LockProtection();
      }

      // DEBUG: CHECK MOTHER AND CHILD HAVE SAME LOCATIONS
      if (mother->location_ != new_child->location_) {
        Log::Warning("\n\nGiveBirth::Run()",
                     "Mother created a child who is at a different location");
      }
      // DEBUG: CHECK MOTHER AND CHILD POINT ON EACH OTHERS
      if (!mother->IsParentOf(new_child->GetAgentPtr<Person>())) {
        Log::Warning("\n\nGiveBirth::Run()", "Mother does not point on child ");
      }
      if (new_child->mother_ != mother->GetAgentPtr<Person>()) {
        Log::Warning("\n\nGiveBirth::Run()", "Child does not point on mother ");
      }
    }
  }
};

}  // namespace hiv_malawi
}  // namespace bdm

#endif  // PERSON_BEHAVIOR_H_
