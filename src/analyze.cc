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

#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"

#include "datatypes.h"

#include "biodynamo.h"
#include "core/util/log.h"
#include "person.h"
#include "sim-param.h"

namespace bdm {
namespace hiv_malawi {

using experimental::Counter;
using experimental::GenericReducer;

void DefineAndRegisterCollectors() {
  // Get population statistics, i.e. extract data from simulation
  // Get the pointer to the TimeSeries
  auto* ts = Simulation::GetActive()->GetTimeSeries();

  // Define how to get the time values of the TimeSeries
  auto get_year = [](Simulation* sim) {
    // AM: Starting Year Variable set in sim_params
    int start_year = sim->GetParam()->Get<SimParam>()->start_year;

    return static_cast<double>(start_year +
                               sim->GetScheduler()->GetSimulatedSteps());
  };

  // Define how to count the healthy individuals
  auto healthy = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsHealthy();
  };
  ts->AddCollector("healthy_agents", new Counter<double>(healthy), get_year);

  // Define how to count the infected individuals
  auto infected = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy());
  };
  ts->AddCollector("infected_agents", new Counter<double>(infected), get_year);

  // Some are  infected but not acute or chronic

  auto noac = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (!(person->IsHealthy()) && !(person->IsAcute()) && !(person->IsChronic()));
  };
  ts->AddCollector("noac", new Counter<double>(noac), get_year);




  auto totalpop  = [](Agent *a) {
    auto* person = bdm_static_cast<Person*>(a);
    return(true);
  };
  ts->AddCollector("all agents", new Counter<double>(totalpop), get_year);

  // AM: Define how to count the infected acute individuals
  auto acute = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute();
  };
  ts->AddCollector("acute_agents", new Counter<double>(acute), get_year);

  // Infected for location1 
  auto acute1 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person-> Isloc1();
  };
  ts->AddCollector("acute_agents loc1", new Counter<double>(acute1), get_year);

// Infected for location2 
  auto acute2 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person-> Isloc2();
  };
  ts->AddCollector("acute_agents loc2", new Counter<double>(acute2), get_year);

  // AM: Define how to count the infected acute male individuals
  auto acute_male_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsMale();
  };
  ts->AddCollector("acute_male_agents", new Counter<double>(acute_male_agents),
                   get_year);

  // AM: Define how to count the infected acute male individuals with low risk
  // behaviours
  auto acute_male_low_sb_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsMale() &&
           person->HasLowRiskSocioBehav();
  };
  ts->AddCollector("acute_male_low_sb_agents",
                   new Counter<double>(acute_male_low_sb_agents), get_year);

  // AM: Define how to count the infected acute male individuals with high risk
  // behaviours
  auto acute_male_high_sb_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsMale() &&
           person->HasHighRiskSocioBehav();
  };
  ts->AddCollector("acute_male_high_sb_agents",
                   new Counter<double>(acute_male_high_sb_agents), get_year);

  // AM: Define how to count the infected acute female individuals
  auto acute_female_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsFemale();
  };
  ts->AddCollector("acute_female_agents",
                   new Counter<double>(acute_female_agents), get_year);

  // AM: Define how to count the infected acute female individuals with low risk
  // sociobehaviours
  auto acute_female_low_sb_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsFemale() &&
           person->HasLowRiskSocioBehav();
  };
  ts->AddCollector("acute_female_low_sb_agents",
                   new Counter<double>(acute_female_low_sb_agents), get_year);

  // AM: Define how to count the infected acute female individuals with high
  // risk sociobehaviours
  auto acute_female_high_sb_agents = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsAcute() && person->IsFemale() &&
           person->HasHighRiskSocioBehav();
  };
  ts->AddCollector("acute_female_high_sb_agents",
                   new Counter<double>(acute_female_high_sb_agents), get_year);

  // AM: Define how to count the infected chronic individuals
  auto chronic = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsChronic();
  };
  ts->AddCollector("chronic_agents", new Counter<double>(chronic), get_year);

  // AM: Define how to count the infected treated individuals
  auto treated = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsTreated();
  };
  ts->AddCollector("treated_agents", new Counter<double>(treated), get_year);

  // AM: Define how to count the infected failing individuals
  auto failing = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsFailing();
  };
  ts->AddCollector("failing_agents", new Counter<double>(failing), get_year);

  // AM: Define how to count the individuals infected at birth
  auto mtct = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->MTCTransmission();
  };
  ts->AddCollector("mtct_agents", new Counter<double>(mtct), get_year);

  // AM: Define how to count the male individuals infected at birth
  auto mtct_transmission_to_male = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->MTCTransmission() && person->IsMale();
  };
  ts->AddCollector("mtct_transmission_to_male",
                   new Counter<double>(mtct_transmission_to_male), get_year);

  // AM: Define how to count the female individuals infected at birth
  auto mtct_transmission_to_female = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->MTCTransmission() && person->IsFemale();
  };
  ts->AddCollector("mtct_transmission_to_female",
                   new Counter<double>(mtct_transmission_to_female), get_year);

  // AM: Define how to count the individuals infected through casual mating
  auto casual = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->CasualTransmission();
  };
  ts->AddCollector("casual_transmission_agents", new Counter<double>(casual),
                   get_year);

  // AM: Define how to count the male individuals infected through casual mating
  auto casual_transmission_to_male = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->CasualTransmission() && person->IsMale();
  };
  ts->AddCollector("casual_transmission_to_male",
                   new Counter<double>(casual_transmission_to_male), get_year);
  // AM: Define how to count the female individuals infected through casual
  // mating
  auto casual_transmission_to_female = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->CasualTransmission() && person->IsFemale();
  };
  ts->AddCollector("casual_transmission_to_female",
                   new Counter<double>(casual_transmission_to_female),
                   get_year);

  // AM: Define how to count the individuals infected through regular mating
  auto regular = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->RegularTransmission();
  };
  ts->AddCollector("regular_transmission_agents", new Counter<double>(regular),
                   get_year);
  // AM: Define how to count the male individuals infected through regular
  // mating
  auto regular_transmission_to_male = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->RegularTransmission() && person->IsMale();
  };
  ts->AddCollector("regular_transmission_to_male",
                   new Counter<double>(regular_transmission_to_male), get_year);
  // AM: Define how to count the female individuals infected through regular
  // mating
  auto regular_transmission_to_female = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->RegularTransmission() && person->IsFemale();
  };
  ts->AddCollector("regular_transmission_to_female",
                   new Counter<double>(regular_transmission_to_female),
                   get_year);

  // AM: Define how to count the individuals that were infected by an Acute HIV
  // partner/Mother
  auto acute_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->AcuteTransmission();
  };
  ts->AddCollector("acute_transmission",
                   new Counter<double>(acute_transmission), get_year);

  // AM: Define how to count the individuals that were infected by an Chronic
  // HIV partner/Mother
  auto chronic_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->ChronicTransmission();
  };
  ts->AddCollector("chronic_transmission",
                   new Counter<double>(chronic_transmission), get_year);

  // AM: Define how to count the individuals that were infected by an Treated
  // HIV partner/Mother
  auto treated_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->TreatedTransmission();
  };
  ts->AddCollector("treated_transmission",
                   new Counter<double>(treated_transmission), get_year);

  // AM: Define how to count the individuals that were infected by an Failing
  // HIV partner/Mother
  auto failing_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->FailingTransmission();
  };
  ts->AddCollector("failing_transmission",
                   new Counter<double>(failing_transmission), get_year);

  // AM: Define how to count the individuals that were infected by an low risk
  // HIV partner
  auto low_sb_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->LowRiskTransmission();
  };
  ts->AddCollector("low_sb_transmission",
                   new Counter<double>(low_sb_transmission), get_year);

  // AM: Define how to count the individuals that were infected by a high risk
  // HIV partner
  auto high_sb_transmission = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HighRiskTransmission();
  };
  ts->AddCollector("high_sb_transmission",
                   new Counter<double>(high_sb_transmission), get_year);

  // Define how to compute mean number of casual partners for males with
  // low-risk sociobehaviours
  //
  // Define how to count adult males younger than 50 with low risk social
  // behavior
  auto adult_male_age_lt50_low_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (person->IsMale() && person->IsAdult() && person->age_ < 50 &&
            person->HasLowRiskSocioBehav());
  };
  ts->AddCollector("adult_male_age_lt50_low_sb",
                   new Counter<double>(adult_male_age_lt50_low_sb), get_year);

  // Sum all casual partners for adult_male_age_lt50_low_sb
  auto sum_casual_partners = [](Agent* agent, uint64_t* tl_result) {
    *tl_result += bdm_static_cast<Person*>(agent)->no_casual_partners_;
  };
  auto sum_tl_results = [](const SharedData<uint64_t>& tl_results) {
    uint64_t result = 0;
    for (auto& el : tl_results) {
      result += el;
    }
    return result;
  };
  ts->AddCollector(
      "total_nocas_men_low_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_male_age_lt50_low_sb),
      get_year);

  auto mean_nocas_men_low_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_men_low_sb").back();
    auto agents = ts->GetYValues("adult_male_age_lt50_low_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_men_low_sb", mean_nocas_men_low_sb, get_year);

  // Define how to compute mean number of casual partners for males with
  // high-risk sociobehaviours
  //
  // Define how to count adult males younger than 50 with high risk social
  // behavior
  auto adult_male_age_lt50_high_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (person->IsMale() && person->IsAdult() && person->age_ < 50 &&
            person->HasHighRiskSocioBehav());
  };
  ts->AddCollector("adult_male_age_lt50_high_sb",
                   new Counter<double>(adult_male_age_lt50_high_sb), get_year);

  // Sum all casual partners for adult_male_age_lt50_high_sb
  ts->AddCollector(
      "total_nocas_men_high_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_male_age_lt50_high_sb),
      get_year);

  auto mean_nocas_men_high_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_men_high_sb").back();
    auto agents = ts->GetYValues("adult_male_age_lt50_high_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_men_high_sb", mean_nocas_men_high_sb, get_year);

  // Define how to compute mean number of casual partners for females with
  // low-risk sociobehaviours
  //
  // Define how to count adult females younger than 50 with low risk social
  // behavior
  auto adult_female_age_lt50_low_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (person->IsFemale() && person->IsAdult() && person->age_ < 50 &&
            person->HasLowRiskSocioBehav());
  };
  ts->AddCollector("adult_female_age_lt50_low_sb",
                   new Counter<double>(adult_female_age_lt50_low_sb), get_year);

  // Sum all casual partners for adult_female_age_lt50_low_sb
  ts->AddCollector(
      "total_nocas_women_low_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_female_age_lt50_low_sb),
      get_year);

  auto mean_nocas_women_low_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_women_low_sb").back();
    auto agents = ts->GetYValues("adult_female_age_lt50_low_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_women_low_sb", mean_nocas_women_low_sb,
                   get_year);

  // Define how to compute mean number of casual partners for females with
  // high-risk sociobehaviours
  //
  // Define how to count adult females younger than 50 with high risk social
  // behavior
  auto adult_female_age_lt50_high_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (person->IsFemale() && person->IsAdult() && person->age_ < 50 &&
            person->HasHighRiskSocioBehav());
  };
  ts->AddCollector("adult_female_age_lt50_high_sb",
                   new Counter<double>(adult_female_age_lt50_high_sb),
                   get_year);

  // Sum all casual partners for adult_female_age_lt50_high_sb
  ts->AddCollector(
      "total_nocas_women_high_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_female_age_lt50_high_sb),
      get_year);

  auto mean_nocas_women_high_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_women_high_sb").back();
    auto agents = ts->GetYValues("adult_female_age_lt50_high_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_women_high_sb", mean_nocas_women_high_sb,
                   get_year);

  // AM: Define how to compute mean number of casual partners for HIV infected
  // females with high-risk sociobehaviours
  //
  // Define how to count adult HIV infected females younger than 50 with high
  // risk social behavior
  auto adult_hiv_female_age_lt50_high_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (!person->IsHealthy() && person->IsFemale() && person->IsAdult() &&
            person->age_ < 50 && person->HasHighRiskSocioBehav());
  };
  ts->AddCollector("adult_hiv_female_age_lt50_high_sb",
                   new Counter<double>(adult_hiv_female_age_lt50_high_sb),
                   get_year);

  // Sum all casual partners for adult_hiv_female_age_lt50_high_sb
  ts->AddCollector(
      "total_nocas_hiv_women_high_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_hiv_female_age_lt50_high_sb),
      get_year);

  auto mean_nocas_hiv_women_high_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_hiv_women_high_sb").back();
    auto agents = ts->GetYValues("adult_hiv_female_age_lt50_high_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_hiv_women_high_sb", mean_nocas_hiv_women_high_sb,
                   get_year);

  // AM: Define how to compute mean number of casual partners for HIV infected
  // females with high-risk sociobehaviours
  //
  // Define how to count adult HIV infected females younger than 50 with high
  // risk social behavior
  auto adult_hiv_female_age_lt50_low_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (!person->IsHealthy() && person->IsFemale() && person->IsAdult() &&
            person->age_ < 50 && person->HasLowRiskSocioBehav());
  };
  ts->AddCollector("adult_hiv_female_age_lt50_low_sb",
                   new Counter<double>(adult_hiv_female_age_lt50_low_sb),
                   get_year);

  // Sum all casual partners for adult_hiv_female_age_lt50_low_sb
  ts->AddCollector(
      "total_nocas_hiv_women_low_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_hiv_female_age_lt50_low_sb),
      get_year);

  auto mean_nocas_hiv_women_low_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_hiv_women_low_sb").back();
    auto agents = ts->GetYValues("adult_hiv_female_age_lt50_low_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_hiv_women_low_sb", mean_nocas_hiv_women_low_sb,
                   get_year);

  // AM: Define how to compute mean number of casual partners for HIV infected
  // males with high-risk sociobehaviours
  //
  // Define how to count adult HIV infected males younger than 50 with high risk
  // social behavior
  auto adult_hiv_male_age_lt50_high_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (!person->IsHealthy() && person->IsMale() && person->IsAdult() &&
            person->age_ < 50 && person->HasHighRiskSocioBehav());
  };
  ts->AddCollector("adult_hiv_male_age_lt50_high_sb",
                   new Counter<double>(adult_hiv_male_age_lt50_high_sb),
                   get_year);

  // Sum all casual partners for adult_hiv_male_age_lt50_high_sb
  ts->AddCollector(
      "total_nocas_hiv_men_high_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_hiv_male_age_lt50_high_sb),
      get_year);

  auto mean_nocas_hiv_men_high_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_hiv_men_high_sb").back();
    auto agents = ts->GetYValues("adult_hiv_male_age_lt50_high_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_hiv_men_high_sb", mean_nocas_hiv_men_high_sb,
                   get_year);

  // AM: Define how to compute mean number of casual partners for HIV infected
  // males with low-risk sociobehaviours
  //
  // Define how to count adult HIV infected males younger than 50 with high risk
  // social behavior
  auto adult_hiv_male_age_lt50_low_sb = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return (!person->IsHealthy() && person->IsMale() && person->IsAdult() &&
            person->age_ < 50 && person->HasLowRiskSocioBehav());
  };
  ts->AddCollector("adult_hiv_male_age_lt50_low_sb",
                   new Counter<double>(adult_hiv_male_age_lt50_low_sb),
                   get_year);

  // Sum all casual partners for adult_hiv_male_age_lt50_low_sb
  ts->AddCollector(
      "total_nocas_hiv_men_low_sb",
      new GenericReducer<uint64_t, double>(sum_casual_partners, sum_tl_results,
                                           adult_hiv_male_age_lt50_low_sb),
      get_year);

  auto mean_nocas_hiv_men_low_sb = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto num_casual = ts->GetYValues("total_nocas_hiv_men_low_sb").back();
    auto agents = ts->GetYValues("adult_hiv_male_age_lt50_low_sb").back();
    return num_casual / agents;
  };
  ts->AddCollector("mean_nocas_hiv_men_low_sb", mean_nocas_hiv_men_low_sb,
                   get_year);

  // AM: Define how to compute general prevalence
  auto pct_prevalence = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_agents").back();
    return infected / sim->GetResourceManager()->GetNumAgents();
  };
  ts->AddCollector("prevalence", pct_prevalence, get_year);

  // AM: Define how to compute prevalence between 15 and 49 year olds
  auto infected_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) && person->age_ >= 15 && person->age_ < 50;
  };
  ts->AddCollector("infected_15_49", new Counter<double>(infected_15_49),
                   get_year);

  auto all_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->age_ >= 15 && person->age_ < 50;
  };
  ts->AddCollector("all_15_49", new Counter<double>(all_15_49), get_year);

  auto pct_prevalence_15_49 = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_15_49").back();
    auto all = ts->GetYValues("all_15_49").back();
    return infected / all;
  };
  ts->AddCollector("prevalence_15_49", pct_prevalence_15_49, get_year);

  // AM: Define how to compute prevalence among women
  auto infected_females = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) && person->IsFemale();
  };
  ts->AddCollector("infected_females", new Counter<double>(infected_females),
                   get_year);

  auto females = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsFemale();
  };
  ts->AddCollector("females", new Counter<double>(females), get_year);

  auto pct_prevalence_females = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_females").back();
    auto all = ts->GetYValues("females").back();
    return infected / all;
  };
  ts->AddCollector("prevalence_females", pct_prevalence_females, get_year);

  // AM: Define how to compute prevalence among women between 15 and 49
  auto infected_women_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) && person->IsFemale() && person->age_ >= 15 &&
           person->age_ < 50;
  };
  ts->AddCollector("infected_women_15_49",
                   new Counter<double>(infected_women_15_49), get_year);

  auto women_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsFemale() && person->age_ >= 15 && person->age_ < 50;
  };
  ts->AddCollector("women_15_49", new Counter<double>(women_15_49), get_year);

  auto pct_prevalence_women_15_49 = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_women_15_49").back();
    auto all = ts->GetYValues("women_15_49").back();
    return infected / all;
  };
  ts->AddCollector("prevalence_women_15_49", pct_prevalence_women_15_49,
                   get_year);

  // AM: Define how to compute prevalence among men
  auto infected_males = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) && person->IsMale();
  };
  ts->AddCollector("infected_males", new Counter<double>(infected_males),
                   get_year);

  auto males = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsMale();
  };
  ts->AddCollector("males", new Counter<double>(males), get_year);

  auto pct_prevalence_males = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_males").back();
    auto all = ts->GetYValues("males").back();
    return infected / all;
  };
  ts->AddCollector("prevalence_males", pct_prevalence_males, get_year);

  // AM: Define how to compute prevalence among men between 15 and 49
  auto infected_men_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) && person->IsMale() && person->age_ >= 15 &&
           person->age_ < 50;
  };
  ts->AddCollector("infected_men_15_49",
                   new Counter<double>(infected_men_15_49), get_year);

  auto men_15_49 = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsMale() && person->age_ >= 15 && person->age_ < 50;
  };
  ts->AddCollector("men_15_49", new Counter<double>(men_15_49), get_year);

  auto pct_prevalence_men_15_49 = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto infected = ts->GetYValues("infected_men_15_49").back();
    auto all = ts->GetYValues("men_15_49").back();
    return infected / all;
  };
  ts->AddCollector("prevalence_men_15_49", pct_prevalence_men_15_49, get_year);

  // AM: Define how to compute general incidence
  auto pct_incidence = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto acute = ts->GetYValues("acute_agents").back();
    return acute / sim->GetResourceManager()->GetNumAgents();
  };
  ts->AddCollector("incidence", pct_incidence, get_year);

  // AM: Define how to compute proportion of people with high-risk
  // socio-beahviours among hiv+
  auto high_risk_hiv = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and !(person->IsHealthy());
  };
  ts->AddCollector("high_risk_hiv", new Counter<double>(high_risk_hiv),
                   get_year);

  auto pct_high_risk_hiv = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto high_risk_hiv = ts->GetYValues("high_risk_hiv").back();
    auto infected = ts->GetYValues("infected_agents").back();
    return high_risk_hiv / infected;
  };
  ts->AddCollector("high_risk_sb_hiv", pct_high_risk_hiv, get_year);

  // AM: Define how to compute proportion of people with low-risk
  // socio-beahviours among hiv+
  auto low_risk_hiv = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and !(person->IsHealthy());
  };
  ts->AddCollector("low_risk_hiv", new Counter<double>(low_risk_hiv), get_year);

  auto pct_low_risk_hiv = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto low_risk_hiv = ts->GetYValues("low_risk_hiv").back();
    auto infected = ts->GetYValues("infected_agents").back();
    return low_risk_hiv / infected;
  };
  ts->AddCollector("low_risk_sb_hiv", pct_low_risk_hiv, get_year);

  // AM: Define how to compute proportion of people with high-risk
  // socio-beahviours among healthy
  auto high_risk_healthy = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and person->IsHealthy();
  };
  ts->AddCollector("high_risk_healthy", new Counter<double>(high_risk_healthy),
                   get_year);

  auto pct_high_risk_healthy = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto high_risk_healthy = ts->GetYValues("high_risk_healthy").back();
    auto healthy = ts->GetYValues("healthy_agents").back();
    return high_risk_healthy / healthy;
  };
  ts->AddCollector("high_risk_sb_healthy", pct_high_risk_healthy, get_year);

  // AM: Define how to compute proportion of people with low-risk
  // socio-beahviours among healthy
  auto low_risk_healthy = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and person->IsHealthy();
  };
  ts->AddCollector("low_risk_healthy", new Counter<double>(low_risk_healthy),
                   get_year);

  auto pct_low_risk_healthy = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    auto low_risk_healthy = ts->GetYValues("low_risk_healthy").back();
    auto healthy = ts->GetYValues("healthy_agents").back();
    return low_risk_healthy / healthy;
  };
  ts->AddCollector("low_risk_sb_healthy", pct_low_risk_healthy, get_year);

  // AM: Define how to compute proportion of high-risk socio-beahviours among
  // hiv adult women
  auto high_risk_hiv_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and !(person->IsHealthy()) and
           person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("high_risk_hiv_women",
                   new Counter<double>(high_risk_hiv_women), get_year);

  auto hiv_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) and person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("hiv_women", new Counter<double>(hiv_women), get_year);

  auto pct_high_risk_hiv_women = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("high_risk_hiv_women").back() /
           ts->GetYValues("hiv_women").back();
  };
  ts->AddCollector("high_risk_sb_hiv_women", pct_high_risk_hiv_women, get_year);

  // AM: Define how to compute proportion of low-risk socio-beahviours among hiv
  // adult women
  auto low_risk_hiv_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and !(person->IsHealthy()) and
           person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("low_risk_hiv_women",
                   new Counter<double>(low_risk_hiv_women), get_year);

  auto pct_low_risk_hiv_women = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("low_risk_hiv_women").back() /
           ts->GetYValues("hiv_women").back();
  };
  ts->AddCollector("low_risk_sb_hiv_women", pct_low_risk_hiv_women, get_year);

  // AM: Define how to compute proportion of high-risk socio-beahviours among
  // hiv adult men
  auto high_risk_hiv_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and !(person->IsHealthy()) and
           person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("high_risk_hiv_men", new Counter<double>(high_risk_hiv_men),
                   get_year);

  auto hiv_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return !(person->IsHealthy()) and person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("hiv_men", new Counter<double>(hiv_men), get_year);

  auto pct_high_risk_hiv_men = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("high_risk_hiv_men").back() /
           ts->GetYValues("hiv_men").back();
  };
  ts->AddCollector("high_risk_sb_hiv_men", pct_high_risk_hiv_men, get_year);

  // AM: Define how to compute proportion of low-risk socio-beahviours among hiv
  // adult men
  auto low_risk_hiv_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and !(person->IsHealthy()) and
           person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("low_risk_hiv_men", new Counter<double>(low_risk_hiv_men),
                   get_year);

  auto pct_low_risk_hiv_men = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("low_risk_hiv_men").back() /
           ts->GetYValues("hiv_men").back();
  };
  ts->AddCollector("low_risk_sb_hiv_men", pct_low_risk_hiv_men, get_year);

  // AM: Define how to compute proportion of high-risk socio-beahviours among
  // healthy adult women
  auto high_risk_healthy_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and person->IsHealthy() and
           person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("high_risk_healthy_women",
                   new Counter<double>(high_risk_healthy_women), get_year);

  auto healthy_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsHealthy() and person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("healthy_women", new Counter<double>(healthy_women),
                   get_year);

  auto pct_high_risk_healthy_women = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("high_risk_healthy_women").back() /
           ts->GetYValues("healthy_women").back();
  };
  ts->AddCollector("high_risk_sb_healthy_women", pct_high_risk_healthy_women,
                   get_year);

  // AM: Define how to compute proportion of low-risk socio-beahviours among
  // healthy adult women
  auto low_risk_healthy_women = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and person->IsHealthy() and
           person->IsAdult() and person->IsFemale();
  };
  ts->AddCollector("low_risk_healthy_women",
                   new Counter<double>(low_risk_healthy_women), get_year);

  auto pct_low_risk_healthy_women = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("low_risk_healthy_women").back() /
           ts->GetYValues("healthy_women").back();
  };
  ts->AddCollector("low_risk_sb_healthy_women", pct_low_risk_healthy_women,
                   get_year);

  // AM: Define how to compute proportion of high-risk socio-beahviours among
  // healthy adult men
  auto high_risk_healthy_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasHighRiskSocioBehav() and person->IsHealthy() and
           person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("high_risk_healthy_men",
                   new Counter<double>(high_risk_healthy_men), get_year);

  auto healthy_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->IsHealthy() and person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("healthy_men", new Counter<double>(healthy_men), get_year);

  auto pct_high_risk_healthy_men = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("high_risk_healthy_men").back() /
           ts->GetYValues("healthy_men").back();
  };
  ts->AddCollector("high_risk_sb_healthy_men", pct_high_risk_healthy_men,
                   get_year);

  // AM: Define how to compute proportion of low-risk socio-beahviours among
  // healthy adult men
  auto low_risk_healthy_men = [](Agent* a) {
    auto* person = bdm_static_cast<Person*>(a);
    return person->HasLowRiskSocioBehav() and person->IsHealthy() and
           person->IsAdult() and person->IsMale();
  };
  ts->AddCollector("low_risk_healthy_men",
                   new Counter<double>(low_risk_healthy_men), get_year);

  auto pct_low_risk_healthy_men = [](Simulation* sim) {
    auto* ts = sim->GetTimeSeries();
    return ts->GetYValues("low_risk_healthy_men").back() /
           ts->GetYValues("healthy_men").back();
  };
  ts->AddCollector("low_risk_sb_healthy_men", pct_low_risk_healthy_men,
                   get_year);
}

// -----------------------------------------------------------------------------
int PlotAndSaveTimeseries() {
  // Get pointers for simulation and TimeSeries data
  auto sim = Simulation::GetActive();
  auto* ts = sim->GetTimeSeries();

  // Save the TimeSeries Data as JSON to the folder <date_time>
  ts->SaveJson(Concat(sim->GetOutputDir(), "/data.json"));

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g(ts, "Population - Healthy/Infected", "Time",
                                 "Number of agents", true);
  g.Add("healthy_agents", "Healthy", "L", kBlue, 1.0);
  g.Add("infected_agents", "HIV", "L", kRed, 1.0);
  g.Draw();
  g.SaveAs(Concat(sim->GetOutputDir(), "/simulation_hiv"), {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2(ts, "HIV stages", "Time", "Number of agents",
                                  true);
  // g2.Add("healthy_agents", "Healthy", "L", kBlue, 1.0, 1);
  g2.Add("infected_agents", "HIV", "L", kOrange, 1.0, 1);
  g2.Add("acute_agents", "Acute", "L", kRed, 1.0, 10);
  g2.Add("chronic_agents", "Chronic", "L", kMagenta, 1.0, 10);
  g2.Add("treated_agents", "Treated", "L", kGreen, 1.0, 10);
  g2.Add("failing_agents", "Failing", "L", kGray, 1.0, 10);
  g2.Draw();
  g2.SaveAs(Concat(sim->GetOutputDir(), "/simulation_hiv_with_states"),
            {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_1(ts, "Acute HIV by sex", "Time",
                                    "Number of agents", true);
  g2_1.Add("acute_male_agents", "Male Acute", "L", kBlue, 1.0, 1);
  g2_1.Add("acute_female_agents", "Female Acute", "L", kMagenta, 1.0, 1);
  g2_1.Draw();
  g2_1.SaveAs(Concat(sim->GetOutputDir(), "/simulation_hiv_acute_sex"),
              {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_1_1(ts, "Acute HIV by sex and risk", "Time",
                                      "Number of agents", true);
  g2_1_1.Add("acute_male_low_sb_agents", "Male Acute - Low risk", "L", kBlue,
             1.0, 2);
  g2_1_1.Add("acute_male_high_sb_agents", "Male Acute - High risk", "L", kBlue,
             1.0, 1);
  g2_1_1.Add("acute_female_low_sb_agents", "Female Acute - Low risk", "L",
             kMagenta, 1.0, 2);
  g2_1_1.Add("acute_female_high_sb_agents", "Female Acute - High risk", "L",
             kMagenta, 1.0, 1);

  g2_1_1.Draw();
  g2_1_1.SaveAs(Concat(sim->GetOutputDir(), "/simulation_hiv_acute_sex_sb"),
                {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_2(ts, "Transmission", "Time",
                                    "Number of agents", true);
  g2_2.Add("mtct_agents", "MTCT", "L", kGreen, 1.0, 3);
  g2_2.Add("casual_transmission_agents", "Casual Transmission", "L", kRed, 1.0,
           3);
  g2_2.Add("regular_transmission_agents", "Regular Transmission", "L", kBlue,
           1.0, 3);
  g2_2.Draw();
  g2_2.SaveAs(Concat(sim->GetOutputDir(), "/simulation_transmission_types"),
              {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_2_1(ts, "Transmission", "Time",
                                      "Number of agents", true);
  g2_2_1.Add("casual_transmission_to_male", "Casual Transmission - to Male",
             "L", kBlue, 1.0, 1);
  g2_2_1.Add("casual_transmission_to_female", "Casual Transmission - to Female",
             "L", kMagenta, 1.0, 1);
  g2_2_1.Add("regular_transmission_to_male", "Regular Transmission - to Male",
             "L", kBlue, 1.0, 2);
  g2_2_1.Add("regular_transmission_to_male", "Regular Transmission - to Female",
             "L", kMagenta, 1.0, 2);
  g2_2_1.Add("mtct_transmission_to_male", "MTCT - to Male", "L", kBlue, 1.0, 3);
  g2_2_1.Add("mtct_transmission_to_female", "MTCT Transmission - to Female",
             "L", kMagenta, 1.0, 3);
  g2_2_1.Draw();
  g2_2_1.SaveAs(
      Concat(sim->GetOutputDir(), "/simulation_transmission_types_by_sex"),
      {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_3(ts, "Source of infection - HIV stage",
                                    "Time", "Number of agents", true);
  g2_3.Add("acute_transmission", "Infected by Acute", "L", kRed, 1.0, 10);
  g2_3.Add("chronic_transmission", "Infected by Chronic", "L", kMagenta, 1.0,
           10);
  g2_3.Add("treated_transmission", "Infected by Treated", "L", kGreen, 1.0, 10);
  g2_3.Add("failing_transmission", "Infected by Failing", "L", kGray, 1.0, 10);
  g2_3.Draw();
  g2_3.SaveAs(
      Concat(sim->GetOutputDir(), "/simulation_transmission_sources_state"),
      {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2_4(ts, "Source of infection - Risk level",
                                    "Time", "Number of agents", true);
  g2_4.Add("low_sb_transmission", "Infected by Low Risk", "L", kRed, 1.0, 10);
  g2_4.Add("high_sb_transmission", "Infected by High Risk", "L", kMagenta, 1.0,
           10);
  g2_4.Draw();
  g2_4.SaveAs(
      Concat(sim->GetOutputDir(), "/simulation_transmission_sources_sb"),
      {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g3(ts, "HIV", "Time", "", true);
  g3.Add("prevalence", "Prevalence", "L", kOrange, 1.0, 3, 1, kOrange, 1.0, 5);
  g3.Add("prevalence_females", "Prevalence - Females", "L", kRed, 1.0, 3, 1,
         kRed, 1.0, 10);
  g3.Add("prevalence_males", "Prevalence - Males", "L", kBlue, 1.0, 3, 1, kBlue,
         1.0, 10);

  g3.Add("prevalence_15_49", "Prevalence (15-49)", "L", kOrange, 1.0, 1, 1);
  g3.Add("prevalence_women_15_49", "Prevalence - Women (15-49)", "L", kRed, 1.0,
         1, 1);
  g3.Add("prevalence_men_15_49", "Prevalence - Men (15-49)", "L", kBlue, 1.0, 1,
         1);

  g3.Add("incidence", "Incidence", "L", kRed, 1.0, 3, 1, kRed, 1.0, 5);

  g3.Draw();
  g3.SaveAs(Concat(sim->GetOutputDir(), "/simulation_hiv_prevalence_incidence"),
            {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g4(ts, "my result", "Time", "Proportion", true);
  g4.Add("high_risk_sb_hiv", "High Risk SB - HIV", "L", kRed, 1.0, 1);
  g4.Add("low_risk_sb_hiv", "Low Risk SB - HIV", "L", kBlue, 1.0, 1);
  g4.Add("high_risk_sb_healthy", "High Risk SB - Healthy", "L", kOrange, 1.0,
         1);
  g4.Add("low_risk_sb_healthy", "Low Risk SB - Healthy", "L", kGreen, 1.0, 1);
  g4.Add("high_risk_sb_hiv_women", "High Risk SB - HIV Women", "L", kRed, 1.0,
         10);
  g4.Add("low_risk_sb_hiv_women", "Low Risk SB - HIV Women", "L", kBlue, 1.0,
         10);

  g4.Add("high_risk_sb_hiv_men", "High Risk SB - HIV Men", "L", kRed, 1.0, 2);
  g4.Add("low_risk_sb_hiv_men", "Low Risk SB - HIV Men", "L", kBlue, 1.0, 2);

  g4.Add("high_risk_sb_healthy_women", "High Risk SB - Healthy Women", "L",
         kOrange, 1.0, 10);
  g4.Add("low_risk_sb_healthy_women", "Low Risk SB - Healthy Women", "L",
         kGreen, 1.0, 10);

  g4.Add("high_risk_sb_healthy_men", "High Risk SB - Healthy Men", "L", kOrange,
         1.0, 5);
  g4.Add("low_risk_sb_healthy_men", "Low Risk SB - Healthy Men", "L", kGreen,
         1.0, 5);
  g4.Draw();
  g4.SaveAs(Concat(sim->GetOutputDir(), "/simulation_sociobehaviours"),
            {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g5(ts, "Casual sex partners", "Time", "Number",
                                  true);

  g5.Add("mean_nocas_men_low_sb", "Mean - Men w/ Low Risk SB", "L", kGreen, 1.0,
         2);
  g5.Add("mean_nocas_men_high_sb", "Mean - Men w/ High Risk SB", "L", kGreen,
         1.0, 1);
  g5.Add("mean_nocas_women_low_sb", "Mean - Women w/ Low Risk SB", "L", kRed,
         1.0, 2);
  g5.Add("mean_nocas_women_high_sb", "Mean - Women w/ High Risk SB", "L", kRed,
         1.0, 1);
  g5.Add("mean_nocas_hiv_men_low_sb", "Mean - Men w/ HIV & Low Risk SB", "L",
         kBlue, 1.0, 2);
  g5.Add("mean_nocas_hiv_men_high_sb", "Mean - Men w/ HIV & High Risk SB", "L",
         kBlue, 1.0, 1);
  g5.Add("mean_nocas_hiv_women_low_sb", "Mean - Women w/ HIV & Low Risk SB",
         "L", kMagenta, 1.0, 2);
  g5.Add("mean_nocas_hiv_women_high_sb", "Mean - Women w/ HIV & High Risk SB",
         "L", kMagenta, 1.0, 1);

  g5.Draw();
  g5.SaveAs(Concat(sim->GetOutputDir(), "/simulation_casual_mating_mean"),
            {".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g6(ts, "Casual sex partners", "Time", "Number",
                                  true);

  g6.Add("total_nocas_men_low_sb", "Total - Men w/ Low Risk SB", "L", kGreen,
         1.0, 2);
  g6.Add("total_nocas_men_high_sb", "Total - Men w/ High Risk SB", "L", kRed,
         1.0, 2);
  g6.Add("total_nocas_women_low_sb", "Total - Women w/ Low Risk SB", "L",
         kGreen, 1.0, 1);
  g6.Add("total_nocas_women_high_sb", "Total - Women w/ High Risk SB", "L",
         kRed, 1.0, 1);

  g6.Draw();
  g6.SaveAs(Concat(sim->GetOutputDir(), "/simulation_casual_mating_total"),
            {".png"});

  // Print info for user to let him/her know where to find simulation results
  std::string info =
      Concat("<PlotAndSaveTimeseries> ", "Results of simulation were saved to ",
             sim->GetOutputDir(), "/");
  std::cout << "Info: " << info << std::endl;

  return 0;
}

}  // namespace hiv_malawi
}  // namespace bdm
