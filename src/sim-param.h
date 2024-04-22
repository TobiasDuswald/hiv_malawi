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

#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include <vector>
#include "biodynamo.h"
#include "datatypes.h"  //AM: Added to access GemState Enum

namespace bdm {
namespace hiv_malawi {

/// This class defines parameters that are specific to this simulation.
class SimParam : public ParamGroup {
 private:
  BDM_CLASS_DEF_OVERRIDE(SimParam, 1);

 public:
  static const ParamGroupUid kUid;
  SimParam() { Initialize(); }
  virtual ~SimParam() {}
  ParamGroup* NewCopy() const override { return new SimParam(*this); }
  ParamGroupUid GetUid() const override { return kUid; }

  // Starting year
  int start_year = 1975;  // 1975;//1960;

  // The number of iterations that BioDynaMo simulates. (#iterations = #years)
  uint64_t number_of_iterations = 55;  // 20; //45;//5;// 60;// (1960-2020)

  // Number of agents that are present at the first iteration of the simulation
  uint64_t initial_population_size = 53020;  // 3600000;//5302000;

  // Activate an additional safety mechanism: protect mothers from death in the
  // year in which they give birth
  bool protect_mothers_at_birth = false;

  // Age when agents start to engage in sexual activities, e.g. possibly give
  // birth, infect, or get infected
  int min_age = 15;

  // Age when agents stop to engage in casual sexual activities
  int max_age = 50;

  // Age when agents stop to give birth
  int max_age_birth = 50;  // 40;

  // Maximum age that an agent can reach. Once the simulation is calibrated and
  // the parameters are fitted, one may consider removing all healthy agents
  // from the simulation that are older than max_age.
  int age_of_death = 90;

  // Mortality rate by age. Size(mortality_rate_by_age) must be equal to
  // 1+size(mortality_rate_age_transition).
  std::vector<int> mortality_rate_age_transition{15, 50, 90};
  std::vector<float> mortality_rate_by_age{0.01, 0.005, 0.05, 1.0};
  // Test - No death
  // std::vector<float> mortality_rate_by_age{0.0, 0.0, 0.0, 0.0};

  // HIV-related mortality. For Healthy, Acute, Chronic, Treated, Failing states
  std::vector<float> hiv_mortality_rate{0.0, 0.0, 0.05, 0.01, 0.1};
  // Test - No death
  //std::vector<float> hiv_mortality_rate{0.0, 0.0, 0.0, 0.0, 0.0};

  // AM: Probability to migrate
  // TO DO AM: Make this probability dependent on the origin location?
  float migration_probability =  0.01; // 0.0; //  No Mogration //
  // AM: Migration year index
  const std::vector<int> migration_year_transition{1960};
  // AM: Migration Matrix. Year index x Location x Location
  std::vector<std::vector<std::vector<float>>> migration_matrix;

  // AM: Probability that a single man wants to engage in regular partnership
  float regular_partnership_probability = 1.0;

  // AM: Probability that a couple in regular partnership separate
  float break_up_probability = 1.0;

  // Years where number of mates per socio-behavioural factors changes
  const std::vector<int> no_mates_year_transition  //{1960, 1990, 2000};
      {1960, 1991, 1992, 1993, 1994, 1995, 1996,
       1997, 1998, 1999, 2000, 2001, 2002};

  // The mating behaviour is modeled with a random process. For each male agent,
  // we sample the number of female sex partners per year from a Gaussian
  // distribution.
  // Gaussian distribution defining the number of casual partners per year
  // depending on year (see no_mates_year_transition) and socio-behaviour
  const std::vector<std::vector<float>> no_mates_mean
  //    {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
  //    {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
  //    {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};/*{{40.0,80.0},
  //                                                   {30.0,60.0},
  //                                                    {20.0,40.0}};*/
      {{22, 88}, {22*0.9375, 88*0.9375}, {22*0.8750, 88*0.8750},
       {22*0.8125, 88*0.8125}, {22*0.750, 88*0.750}, {22*0.6875, 88*0.6875},
        {22*0.6250, 88*0.6250}, {22*0.5625, 88*0.5625}, {22*0.500, 88*0.500},
         {22*0.4375, 88*0.4375}, {22*0.3750, 88*0.3750}, {22*0.3125, 88*0.3125},
          {22*0.25, 88*0.25}};
  //{{20.0, 70.0}, {15.0, 53.0}, {10.0, 35.0}};
  //{{2.0, 8.0}, {1.0, 4.0}, {1.0, 4.0}};

  const std::vector<std::vector<float>> no_mates_sigma /*{{100.0,100.0},
                                                    {100.0,100.0},
                                                    {100.0,100.0}};*/
      {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
       {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
       {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  //{{8.0, 10.0}, {8.0, 10.0}, {8.0, 10.0}};
  //{{10.0, 10.0}, {10.0, 10.0}, {10.0, 10.0}};

  // We sample the number of sex acts with each female sex partner per year
  // from a Gaussian distribution.
  const std::vector<std::vector<float>> no_acts_mean{
      {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
      {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
      {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
  //{{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
  //{{10.0, 10.0}, {10.0, 10.0}, {10.0, 10.0}};

  const std::vector<std::vector<float>> no_acts_sigma{
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  //{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  //{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};

  const std::vector<int> no_regacts_year_transition{
      1960, 1991, 1992, 1993, 1994, 1995, 1996,
      1997, 1998, 1999, 2000, 2001, 2002};

  // Mean number of sexual acts with regular partner per year.
  // const float no_regular_acts_mean = 50.0;  // 150;
  const std::vector<float> no_regular_acts_mean{50,
   50*0.9375, 50*0.8750, 50*0.8125, 50*0.750, 50*0.6875,
    50*0.6250, 50*0.5625, 50*0.500, 50*0.4375, 50*0.3750, 50*0.3125, 50*0.25};
                                            //50, 47, 44, 41, 38, 34, 31,
                                             //   28, 25, 22, 19, 16, 12};
  //{90, 84, 79, 73, 68, 62, 56, 51, 45, 39, 34, 28, 22};

  // AM: Probability of getting infected depends on
  // 1) disease state, 2) sex of partners Male-to-female
  float coef_infection_probability = 2;
  float infection_probability_acute_mf = 9.3e-3 * coef_infection_probability;
  float infection_probability_chronic_mf = 1.9e-3 * coef_infection_probability;
  float infection_probability_treated_mf = 1.3e-4 * coef_infection_probability;
  float infection_probability_failing_mf = 7.6e-4 * coef_infection_probability;
  // Female-to-male
  float infection_probability_acute_fm = 4.8e-3 * coef_infection_probability;
  float infection_probability_chronic_fm = 9.5e-4 * coef_infection_probability;
  float infection_probability_treated_fm = 6.5e-4 * coef_infection_probability;
  float infection_probability_failing_fm = 3.9e-4 * coef_infection_probability;
  // Male-to-male
  float infection_probability_acute_mm = 9.3e-2 * coef_infection_probability;
  float infection_probability_chronic_mm = 1.9e-2 * coef_infection_probability;
  float infection_probability_treated_mm = 1.3e-3 * coef_infection_probability;
  float infection_probability_failing_mm = 7.6e-3 * coef_infection_probability;

  // AM: Transition Matrix between HIV states.
  // GemState->Year-and-Population-category->GemsState
  std::vector<std::vector<std::vector<float>>> hiv_transition_matrix;

  // AM: Transition Matrix between socio-behavioural categories.
  // Used for yearly update of agents' socio-behaviours
  // nb_sociobehav_categories x Sex x nb_sociobehav_categories
  std::vector<std::vector<std::vector<float>>> sociobehaviour_transition_matrix{
      {
          // sb = 0 (Low risk)
          {1.0, 0.0},  // sex=0 (kMale)
          {1.0, 0.0}   // sex=1 (kFemale)
      },
      {
          // sb=1 (high risk)
         // {0.0, 1.0},  // sex=0  (kMale)
         // {0.0, 1.0}     // sex=1  (kFemale)
          {0.04, 0.96},  // sex=0  (kMale)
          {0.1, 0.9}     // sex=1  (kFemale)
      }};

  // Number of locations
  int nb_locations =  Location::kLocLast; //1; //   

  // AM: Location Mixing Matrix used for casual partner selection.
  // Location->Location
  std::vector<std::vector<float>> location_mixing_matrix;

  // Five-years age categories 15-19, 20-24, ...,65-69,70+
  int nb_age_categories = 12;

  // AM: Age Mixing Matrix used for casual partner selection. Age Category ->
  // Age Category
  std::vector<std::vector<float>> age_mixing_matrix;

  // AM: Age Mixing Matrix used for regular partner selection. Age Category ->
  // Age Category
  std::vector<std::vector<float>> reg_partner_age_mixing_matrix;

  // AM: Number of socio-behavioural categories
  int nb_sociobehav_categories = 2;

  // AM: Socio-behavioural Mixing Matrix used for casual partner selection.
  // Socio-beahioural Category -> Socio-behavioural Category
  std::vector<std::vector<float>> sociobehav_mixing_matrix{{1.0, 4.0},
                                                           {1.0, 4.0}};

  // AM: Socio-behavioural Mixing Matrix used for regular partner selection.
  // Socio-beahioural Category -> Socio-behavioural Category
  std::vector<std::vector<float>> reg_partner_sociobehav_mixing_matrix;

  // Initial prevalence among 15-50 years old.
  float initial_prevalence = 18e-4;  // 15e-3;  // 30e-4; //15e-4;

  // AM: Probability for agent infected by HIV to be at a certain HIV
  // progression state at beginning of simulation. Given in a summed
  // up/cumulative form. First vector x1 = p_Acute (|HIV+) component corresponds
  // to probability of being Acute. x2 = p_Acute(|HIV+) + p_Chronic(|HIV+), x3
  // = p_Acute(|HIV+) + p_Chronic(|HIV+) + p_Treated(|HIV+), x4 = p_Acute(|HIV+)
  // + p_Chronic(|HIV+) + p_Treated(|HIV+) + p_Failing(|HIV+) These probablities
  // involve 15-49 years old agents, located in seed districts. 1/5 of HIV
  // infected are in acute phase, others are chronic.
  std::vector<float> initial_infection_probability{1.0, 1.0, 1.0, 1.0};
  // Initial probability to be healthy for 15-49 years old in seed districts
  float initial_healthy_probability;

  // Districts where HIV infected agents are initially located
  const std::vector<bool> seed_districts{
    //true};
    //    false, true,  false, false, false, false, false, false, false, false,
    //    false,  false, false,  false, false,  false,  false,  false,  false,  false,
    //    false,  false,  false,  false,  false,  false,  false,  false,  false};


    false, true,  false, false, false, false, false, false, false, false,
    true,  false, true,  false, true,  true,  true,  true,  true,  true,
    true,  true,  true,  true,  true,  true,  true,  true,  false};

  // One Location
  //const std::vector<bool> seed_districts{
  //    true, false};

  // Parameter 0.18 is chosen because our GiveBirth Behaviour is based on a
  // Bernoulli experiment. A binomial distribuition peaks at around 6 for 25
  // tries (=40-15) and a birth probability of 0.24. This corresponds to the
  // typical birth rate in the region. We substracted 0.06 to account for child
  // motability and reach a realistic demographic development from 1960-2020.
  // Parameter 0.21 is used in Janne's R implementation.
  float give_birth_probability = 0.188;  // 0.18

  // AM : Probability for agent to be infected at birth, if its mother is
  // infected and treated
  float birth_infection_probability_treated = 0.05;
  // AM : Probability for agent to be infected at birth, if its mother is
  // infected and untreated (i.e. acute, chronic, failing)
  float birth_infection_probability_untreated = 0.35;
  // Probability to be infected at birth, if its mother is under prophylaxis
  // (after 2003)
  float birth_infection_probability_prophylaxis = 0.2;

  // Probability of creating a male agent, used in population initialization and
  // giving birth
  float probability_male = 0.499;

  // Years where probability of high socio-behavioural factor changes
  const std::vector<int> sociobehavioural_risk_year_transition{1960, 1976};

  // Probability of assigning 1 to socio-behavioural factor (high risk)
  // depending on year (see sociobehavioural_risk_year_transition) and health
  // state (Healthy, Acute, Chronic, Treated, Failing)
  const std::vector<std::vector<float>> sociobehavioural_risk_probability{
      //{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}};
      {0.05, 0.5, 0.5, 0.5, 0.5}, {0.05, 0.05, 0.05, 0.05, 0.05}};

  float biomedical_risk_probability = 0.0; //0.05;

  // Age distribution for population initialization. Given in a summed up form.
  // Each vector component corresponds to a age bin of 5 years, e.g. x1 (0-5),
  // x2 (6-10). Usually, these vectors are given in probabilities p1, p2, p3, ..
  // Here: x1 = p1, x2 = p1+p2, x3 = p1+p2+p3, ..
  const std::vector<float> male_age_distribution{
      0.156, 0.312, 0.468, 0.541, 0.614, 0.687, 0.76,  0.833, 0.906,
      0.979, 0.982, 0.985, 0.988, 0.991, 0.994, 0.997, 1};
  const std::vector<float> female_age_distribution{
      0.156, 0.312, 0.468, 0.540, 0.612, 0.684, 0.756,  0.828, 0.9,
      0.972, 0.976, 0.98, 0.984, 0.988, 0.992, 0.996, 1};

  // Location distribution for population initialization, same logic as for age
  // distribution.
  const std::vector<float> location_distribution{
    //1};
      0.012, 0.03,  0.031, 0.088, 0.104, 0.116, 0.175, 0.228, 0.273, 0.4,
      0.431, 0.453, 0.498, 0.517, 0.54,  0.569, 0.645, 0.679, 0.701, 0.736,
      0.794, 0.834, 0.842, 0.86,  0.903, 0.925, 0.995, 1,     1};
  // One Location
  // const std::vector<float> location_distribution{
  //    1.0, 1.0};

  ///////////////////////////////////////////////////////////////////////////
  // Initalizer Functions
  ///////////////////////////////////////////////////////////////////////////

  // Defines the right sizes for the matrices, for some fills it's entries, and
  // is called in constructor.
  void Initialize() {
    SetInitialInfectionProbability();

    // SetSociobehavMixingMatrix();

    SetAgeMixingMatrix();
    SetLocationMixingMatrix();

    SetRegPartnerSociobehavMixingMatrix();
    SetRegPartnerAgeMixingMatrix();

    SetHivTransitionMatrix();
    SetMigrationMatrix();
  };

  // Resizes vector to number of GemsStates and fills
  // as a function of initial_prevalence and seed districts.
  void SetInitialInfectionProbability();

  // Resizes matrix to (nb_age_categories x nb_sociobehav_categories) and fills
  // with ones.
  void SetSociobehavMixingMatrix();

  // Resizes matrix to (nb_age_categories x nb_age_categories) and fills with
  // ones.
  void SetAgeMixingMatrix();

  // Resizes matrix to (nb_locations x nb_locations) and fills with
  // ones.
  void SetLocationMixingMatrix();

  // Resizes matrix to (nb_age_categories x nb_sociobehav_categories) and fills
  // with ones.
  void SetRegPartnerSociobehavMixingMatrix();

  // Resizes matrix to (nb_age_categories x nb_age_categories) and fills with
  // ones.
  void SetRegPartnerAgeMixingMatrix();

  // Computes "index" representing the year and population category.
  // Used in hiv_transition_matrix
  // int ComputeYearPopulationCategory(int year, float age, int sex);

  // Resizes to (nb_states x year_population_category x nb_states) and fills
  // entries with appropriate hard coded values.
  void SetHivTransitionMatrix();

  // Resizes to (migration_year_transitions x nb_locations x nb_locations) and
  // fills entries with appropriate hard coded values.
  void SetMigrationMatrix();

  // Resizes to (migration_year_transitions x nb_locations x nb_locations) and
  // fills entries with normalized and cumulative probabilities
  void SetMigrationLocationProbability();
};

}  // namespace hiv_malawi
}  // namespace bdm

#endif  // SIM_PARAM_H_
