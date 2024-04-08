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

#include <limits>
#include <vector>

#include "biodynamo.h"

#include "datatypes.h"
#include "person-behavior.h"
#include "population-initialization.h"

// All hard-coded numbers are taken from Janne's work (Parameters_D1.R)

namespace bdm {
namespace hiv_malawi {

float SampleAge(float rand_num_1, float rand_num_2, int sex,
                const std::vector<float>& age_distribution) {
  for (size_t i = 0; i < age_distribution.size(); i++) {
    if (rand_num_1 <= age_distribution[i]) {
      return 5 * (i + rand_num_2);
    } else {
      continue;
    }
  }
  // This line of code should never be reached
  Log::Warning("SampleAge()",
               "Could not sample the age. Recieved inputs:", rand_num_1, ", ",
               rand_num_2, ", ", sex, ". Use age 0.");
  return 0;
}

int SampleLocation(float rand_num,
                   const std::vector<float>& location_distribution) {
  for (size_t i = 0; i < location_distribution.size(); i++) {
    if (rand_num <= location_distribution[i]) {
      return i;
    }
  }

  // This line of code should never be reached
  Log::Warning("SampleLocation()",
               "Could not sample the location. Recieved inputs: ", rand_num,
               ". Use location 0.");
  for (size_t i = 0; i < location_distribution.size(); i++) {
    std::cout << location_distribution[i] << ", ";
  }
  std::cout << std::endl;
  return 0;
}

int SampleSex(float rand_num, float probability_male) {
  if (rand_num <= probability_male) {
    return Sex::kMale;
  } else {
    return Sex::kFemale;
  }
}

int SampleState(float rand_num_1, float rand_num_2,
                float initial_healthy_probability,
                const std::vector<float>& initial_infection_probability) {
  if (rand_num_1 <= initial_healthy_probability) {
    return GemsState::kHealthy;
  } else {
    // std::cout << "SampleState ==> HIV, " << rand_num_1 << std::endl;
    for (size_t i = 0; i < initial_infection_probability.size(); i++) {
      if (rand_num_2 <= initial_infection_probability[i]) {
        // std::cout << "SampleState ==> " << i+1 << ", rand_num_2 " <<
        // rand_num_2 << std::endl;
        return i + 1;  // AM: GemsState Enum values are by default associated
                       // with int values
      } else {
        continue;
      }
    }
  }

  // This line of code should never be reached
  Log::Warning("SampleState()",
               "Could not sample the state. Recieved inputs:", rand_num_1,
               rand_num_2, ". Use state GemsState::kHealthy.");
  return GemsState::kHealthy;
}

int ComputeState(float rand_num_1, float rand_num_2, int age, int min_age,
                 int max_age, size_t location,
                 const std::vector<bool>& seed_districts,
                 const float initial_healthy_probability,
                 const std::vector<float>& initial_infection_probability) {
  // Younger than min_age, older than max_age and living outside the subset of
  // seed locations are all healthy.
  if (age < min_age || age >= max_age ||
      (seed_districts[location] ==
       false)) {                 // AM: 15 yo is not a child anymore
    return GemsState::kHealthy;  // Healthy
  }
  // Else, sample the state given the initial infection probability.
  return SampleState(rand_num_1, rand_num_2, initial_healthy_probability,
                     initial_infection_probability);
}

int ComputeSociobehavioural(float rand_num, int age,
                            float sociobehavioural_risk_probability) {
  if (age < 15) {  // AM: 15 yo is not a child anymore
    return 0;
  }
  if (rand_num <= sociobehavioural_risk_probability) {
    return 1;
  } else {
    return 0;
  }
}

int ComputeBiomedical(float rand_num, int age,
                      float biomedical_risk_probability) {
  if (age < 15) {  // / AM: 15 yo is not a child anymore
    return 0;
  }
  if (rand_num <= biomedical_risk_probability) {
    return 1;
  } else {
    return 0;
  }
}

auto CreatePerson(Random* random_generator, const SimParam* sparam) {
  // Get all random numbers for initialization
  std::vector<float> rand_num{};
  rand_num.resize(10);
  for (int i = 0; i < 10; i++) {
    rand_num[i] = static_cast<float>(random_generator->Uniform());
  }

  // Create new person
  Person* person = new Person();
  // Assign sex
  person->sex_ = SampleSex(rand_num[0], sparam->probability_male);
  // Assign age
  if (person->sex_ == Sex::kMale) {
    person->age_ = SampleAge(rand_num[1], rand_num[2], person->sex_,
                             sparam->male_age_distribution);
  } else {
    person->age_ = SampleAge(rand_num[1], rand_num[2], person->sex_,
                             sparam->female_age_distribution);
  }
  // Assign location
  person->location_ =
      SampleLocation(rand_num[3], sparam->location_distribution);
  // Assign the GemsState of the person.
  // AM: This depends on the person's age and location. HIV+ are only between
  // min_age and max_age, ans in a subset of locations.
  person->state_ =
      ComputeState(rand_num[4], rand_num[5], person->age_, sparam->min_age,
                   sparam->max_age, person->location_, sparam->seed_districts,
                   sparam->initial_healthy_probability,
                   sparam->initial_infection_probability);

  // If the person is infected at initialisation, set that it got infected
  // through casual partnership.
  if (person->state_ != GemsState::kHealthy) {
    person->transmission_type_ = TransmissionType::kCasualPartner;
  }

  // Compute risk factors.
  // AM: social_behaviour_factor_ depends on the age and health/hiv state
  person->social_behaviour_factor_ = ComputeSociobehavioural(
      rand_num[6], person->age_,
      sparam->sociobehavioural_risk_probability[0][person->state_]);
  person->biomedical_factor_ = ComputeBiomedical(
      rand_num[7], person->age_, sparam->biomedical_risk_probability);

  // DEBUG
  /*if (person->state_ == GemsState::kAcute){
    std::cout << "Right after ComputeState, state_ = " << person->state_
    << ", IsAcute = " << person->IsAcute()
    << ", IsMale = " << person->IsMale()
    << " (isFemale = "<< person->IsFemale() << ")"
    << ", HasHighRiskSocioBehav = " << person->HasHighRiskSocioBehav()
    << " (HasLowRiskSocioBehav = " << person->HasLowRiskSocioBehav() << ")"
    << std::endl;
  }*/

  ///! The aguments below are currently either not used or repetitive.
  // // Store the year when the agent got infected
  // person->year_of_infection_ = std::numeric_limits<float>::max();
  // // NOTE: we do not assign a specific mother or partner during the
  // population
  // // initialization. Use nullptr.
  // person->mother_id_ = nullptr;
  // person->partner_id_ = nullptr;

  // BioDynaMo API: Add the behaviors to the Agent
  person->AddBehavior(new RandomMigration());
  if (person->sex_ == Sex::kFemale) {
    person->AddBehavior(new GiveBirth());
  } else {
    /*if (person->state_ != GemsState::kHealthy){
      person->AddBehavior(new MatingBehaviour());
    }*/
    person->AddBehavior(new MatingBehaviour());
    //person->AddBehavior(new RegularPartnershipBehaviour());
    person->AddBehavior(new RegularMatingBehaviour());
    person->AddBehavior(new RegularPartnershipBehaviour());
  }
  return person;
};

void InitializePopulation() {
#pragma omp parallel
  {
    auto* sim = Simulation::GetActive();
    auto* ctxt = sim->GetExecutionContext();
    auto* random_generator = sim->GetRandom();
    auto* param = sim->GetParam();
    const auto* sparam = param->Get<SimParam>();

#pragma omp for
    for (uint64_t x = 0; x < sparam->initial_population_size; x++) {
      // Create a person
      auto* new_person = CreatePerson(random_generator, sparam);
      // BioDynaMo API: Add agent (person) to simulation
      ctxt->AddAgent(new_person);
    }
  }
}

}  // namespace hiv_malawi
}  // namespace bdm
