// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN (Tobias Duswald, Lukas Breitwieser, Ahmad Hesam, Fons
// Rademakers) for the benefit of the BioDynaMo collaboration. All Rights
// Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef PERSON_H_
#define PERSON_H_

#include "biodynamo.h"

namespace bdm {
////////////////////////////////////////////////////////////////////////////////
// BioDynaMo's Agent / Individual
////////////////////////////////////////////////////////////////////////////////

class Person : public Cell {
  // BioDynaMo API
  BDM_AGENT_HEADER(Person, Cell, 1);

 public:
  Person() {}
  explicit Person(const Double3& position) : Base(position) {}
  virtual ~Person() {}

  /// Stores the current GemsState of the person.
  int state_;
  // Stores the age of the agent
  float age_;
  // Stores the sex of the agent
  int sex_;
  // Stores the location as categorical variable
  int location_;
  // Stores a factor representing the socio-behavioural risk
  int social_behaviour_factor_;
  // Stores a factor representing the biomedical risk
  int biomedical_factor_;

  ///! The aguments below are currently either not used or repetitive.
  // // Stores if an agent is infected or not
  // bool infected_;
  // // Store the year when the agent got infected
  // float year_of_infection_;
  // // Stores the ID of the mother
  // AgentPointer<Person> mother_id_;
  // // Stores the id of the partner
  // AgentPointer<Person> partner_id_;
};

}  // namespace bdm

#endif  // PERSON_H_