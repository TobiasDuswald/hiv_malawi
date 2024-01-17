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

#ifndef BDM_SIMULAION_H_
#define BDM_SIMULAION_H_

#include <fstream>
#include <iostream>

#include "core/operation/operation_registry.h"
#include "core/operation/reduction_op.h"

#include "analyze.h"
#include "categorical-environment.h"
#include "custom-operations.h"
#include "hiv-ops.h"
#include "population-initialization.h"
#include "sim-param.h"

namespace bdm {
namespace hiv_malawi {

////////////////////////////////////////////////////////////////////////////////
// BioDynaMo's main simulation
////////////////////////////////////////////////////////////////////////////////
inline int Simulate(int argc, const char** argv) {
  // Register the Siulation parameter
  Param::RegisterParamGroup(new SimParam());

  // Initialize the Simulation
  gAgentPointerMode = AgentPointerMode::kDirect;
  auto set_param = [&](Param* param) {
    param->show_simulation_step = 1;
    param->remove_output_dir_contents = false;
    param->statistics = true;
  };
  Simulation simulation(argc, argv, set_param);

  // Get a pointer to the param object
  auto* param = simulation.GetParam();
  // Get a pointer to an instance of SimParam
  auto* sparam = param->Get<SimParam>();

  // AM: Construct Environment with numbers of age and socio-behavioral
  // categories.
  auto* env = new CategoricalEnvironment(
      sparam->min_age, sparam->max_age, sparam->nb_age_categories,
      sparam->nb_locations, sparam->nb_sociobehav_categories);

  simulation.SetEnvironment(env);

  // Randomly initialize a population
  {
    Timing timer_init("RUNTIME POPULATION INITIALIZATION: ");
    InitializePopulation();
  }

  DefineAndRegisterCollectors();

  // Unschedule some default operations
  auto* scheduler = simulation.GetScheduler();
  // Don't compute forces
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);
  // Don't run load balancing, not working with custom environment.
  scheduler->UnscheduleOp(scheduler->GetOps("load balancing")[0]);

  // Add a operation that resets the number of casual partners at the beginning
  // of each iteration
  OperationRegistry::GetInstance()->AddOperationImpl(
      "ResetCasualPartners", OpComputeTarget::kCpu, new ResetCasualPartners());
  auto* reset_casual_partners = NewOperation("ResetCasualPartners");
  scheduler->ScheduleOp(reset_casual_partners, OpType::kPreSchedule);

  // Pre-schedule an operation that implements the previous GetOlder behavior
  OperationRegistry::GetInstance()->AddOperationImpl(
      "GetOlder", OpComputeTarget::kCpu, new GetOlderOperation());
  auto* get_older = NewOperation("GetOlder");
  scheduler->ScheduleOp(get_older, OpType::kPostSchedule);

  // Run simulation for <number_of_iterations> timesteps
  {
    Timing timer_sim("RUNTIME");
    scheduler->Simulate(sparam->number_of_iterations);
  }

  {
    Timing timer_post("RUNTIME POSTPROCESSING:            ");

    // Generate ROOT plot to visualize the number of healthy and infected
    // individuals over time.
    PlotAndSaveTimeseries();
  }

  // DEBUG - AM - TO DO: Works only when selection depended soloely on locations
  /*env->NormalizeMateLocationFrequencies();
  env->PrintMateLocationFrequencies();*/

  return 0;
}

}  // namespace hiv_malawi
}  // namespace bdm

#endif  // BDM_SIMULAION_H_
