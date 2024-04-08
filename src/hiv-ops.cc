// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, UniGe. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#include "hiv-ops.h"
#include "person.h"
#include "sim-param.h"

namespace bdm {
namespace hiv_malawi {

void GetOlderOperation::operator()() {
  auto* sim = Simulation::GetActive();
  auto* rm = sim->GetResourceManager();
  auto process = L2F([&](Agent* a, AgentHandle) { ProcessAgent(a); });
  rm->ForEachAgentParallel(process);
}

float GetOlderOperation::GetMortalityRateAge(
    float age, const std::vector<int>& mortality_rate_age_transition,
    const std::vector<float>& mortality_rate_by_age) {
  size_t age_index = mortality_rate_by_age.size() - 1;
  for (size_t i = 0; i < mortality_rate_age_transition.size(); i++) {
    if (age < mortality_rate_age_transition[i]) {
      age_index = i;
      break;
    }
  }
  // std::cout << "Age " << age << " => Mortality rate " <<
  // mortality_rate_by_age[age_index] << std::endl;
  return mortality_rate_by_age[age_index];
}

float GetOlderOperation::GetMortalityRateHiv(
    int state, const std::vector<float>& hiv_mortality_rate) {
  return hiv_mortality_rate[state];
}

void GetOlderOperation::ProcessAgent(Agent* agent) {
  auto* sim = Simulation::GetActive();
  auto* random = sim->GetRandom();
  auto* param = sim->GetParam();
  const auto* sparam = param->Get<SimParam>();
  auto* person = bdm_static_cast<Person*>(agent);

  // Assign or reassign risk factors
  if (floor(person->age_) ==
      sparam->min_age) {  // Assign potentially high risk
                          // factor at first year of adulthood
    // Probability of being at high risk depends on year and HIV status
    int year = static_cast<int>(
        sparam->start_year +
        sim->GetScheduler()->GetSimulatedSteps());  // Current year
    // Check transition year
    // If no sociobehavioural risk transition year is higher than current
    // year, then use last transition year
    int year_index = sparam->sociobehavioural_risk_year_transition.size() - 1;
    for (size_t y = 0;
         y < sparam->sociobehavioural_risk_year_transition.size() - 1; y++) {
      if (year < sparam->sociobehavioural_risk_year_transition[y + 1]) {
        year_index = y;
        break;
      }
    }

    if (random->Uniform() <=
        sparam->sociobehavioural_risk_probability[year_index][person->state_]) {
      person->social_behaviour_factor_ = 1;
    } else {
      person->social_behaviour_factor_ = 0;
    }
    if (random->Uniform() <= sparam->biomedical_risk_probability) {
      person->biomedical_factor_ = 1;
    } else {
      person->biomedical_factor_ = 0;
    }
  } else if (person->age_ > sparam->min_age) {
    // Potential change in risk factor foradults (after first year of
    // adulthood)
    // Update risk factors stochastically like in initialization
    if (random->Uniform() <=
        sparam
            ->sociobehaviour_transition_matrix[person->social_behaviour_factor_]
                                              [person->sex_][0]) {
      person->social_behaviour_factor_ = 0;

    } else {
      person->social_behaviour_factor_ = 1;
    }
    if (random->Uniform() > sparam->biomedical_risk_probability) {
      person->biomedical_factor_ = 0;
    } else {
      person->biomedical_factor_ = 1;
    }
  } else {  // Low risk factor for children
    person->social_behaviour_factor_ = 0;
    person->biomedical_factor_ = 0;
  }

  // AM: HIV  transition, depending on current year and population
  // category (important for transition to treatment)
  int year_population_category = -1;

  int start_year = sim->GetParam()->Get<SimParam>()->start_year;
  int year =
      static_cast<int>(start_year + sim->GetScheduler()->GetSimulatedSteps());

  // TO DO AM: Replace code below with function :
  // ComputeYearPopulationCategory(int year, float age, int sex)
  // year_population_category =
  // sim->GetParam()->Get<SimParam>()->ComputeYearPopulationCategory(year,
  // person->age_, person->sex_);
  if (year < 2003) {  // Prior to 2003
    year_population_category =
        0;  // All (No difference in ART between people. ART not available.)
  } else if (year < 2011) {  // Between 2003 and
                                           // 2010
    if (person->sex_ == Sex::kFemale && person->age_ >= 15  and
        person->age_ <= 40) {
      year_population_category = 1;  // Female between 15 and 40
    } else if (person->age_ < 15) {
      year_population_category = 2;  // Child
    } else {
      year_population_category = 3;  // Others (Male over 15 and Female over 40)
    }
  } else {                           // After 2011
    if (person->sex_ == Sex::kFemale && person->age_ >= 15  and
        person->age_ <= 40) {
      year_population_category = 4;  // Female between 15 and 40
    } else if (person->age_ < 15 ) {
      year_population_category = 5;  // Child
    } else {
      year_population_category =
          6;  // Others (Male over 15, and Female over 40)
    }
  }
  const auto& transition_proba =
      sparam->hiv_transition_matrix[person->state_][year_population_category];
  for (size_t i = 0; i < transition_proba.size(); i++) {
    if (random->Uniform() < transition_proba[i] && person->newly_infected==false) {
      person->state_ = i;
      break;
    }
  }
  if(person->newly_infected==true){
    //std::cout << "W";
    person->newly_infected=false;
  }
  // Possibly die - if not, just get older
  bool stay_alive{true};

  // AM: Mortality
  // HIV-related mortality
  float rand_num_hiv = static_cast<float>(random->Uniform());
  if (rand_num_hiv <
      GetMortalityRateHiv(person->state_, sparam->hiv_mortality_rate)) {
    stay_alive = false;
  }
  // Age-related mortality
  float rand_num_age = static_cast<float>(random->Uniform());
  if (rand_num_age < GetMortalityRateAge(person->age_,
                                         sparam->mortality_rate_age_transition,
                                         sparam->mortality_rate_by_age)) {
    stay_alive = false;
  }

  // We protect mothers that just gave birth. This should not have a large
  // impact on the simulation. Essentially, if a mother gives birth, she
  // cannot die in this particular year. This is an additional safety
  // mechanism. By default it is not active. If the problem of children not
  // having the right AgentPtr for their mother occurs again, consider
  // activating this switch.
  if (sparam->protect_mothers_at_birth && person->IsProtected()) {
    stay_alive = true;
    // Allow death of agent in the next year.
    person->UnockProtection();
  }

  if (!stay_alive) {
    // Person dies, i.e. is removed from simulation.
    person->RemoveFromSimulation();
    person->will_be_removed_ = true;
  } else {
    // increase age
    person->age_ += 1;
  }
}

}  // namespace hiv_malawi
}  // namespace bdm
