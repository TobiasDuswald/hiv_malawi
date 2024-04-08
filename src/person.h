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

#ifndef PERSON_H_
#define PERSON_H_

#include "biodynamo.h"
#include "core/simulation.h"
#include "datatypes.h"

namespace bdm {
namespace hiv_malawi {
////////////////////////////////////////////////////////////////////////////////
// BioDynaMo's Agent / Individual
////////////////////////////////////////////////////////////////////////////////

class Person : public Agent {
  // BioDynaMo API
  BDM_AGENT_HEADER(Person, Agent, 1);

 public:
  // Warning: many member unititialized!
  Person() {
    mother_ = nullptr;
    partner_ = nullptr;
    children_.clear();
    children_.reserve(3);
    protected_ = false;
    no_casual_partners_ = 0;
  }
  virtual ~Person() {}

  // Overwrite abstract methods from Base class
  Shape GetShape() const override { return Shape::kSphere; }
  double GetDiameter() const override { return 0; }
  const Double3& GetPosition() const override {
    static Double3 kDefault{0, 0, 0};
    return kDefault;
  }
  void SetDiameter(double diameter) override {}
  void SetPosition(const Double3& position) override {}
  Double3 CalculateDisplacement(const InteractionForce* force,
                                double squared_radius, double dt) override {
    return {0, 0, 0};
  }
  void ApplyDisplacement(const Double3& displacement) override {}
  // Overwrite end

  /// Stores the current GemsState of the person.
  int state_;
  // Stores how the agent was infected.
  int transmission_type_;
  // Stores the state of the agent who infected them.
  int infection_origin_state_;
  // Stores the socio-behavioural risk of the agent who infected them.
  int infection_origin_sb_;
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
  // Protect a person against death. Currently only used for mothers in the year
  // in which they give birth and if sparam->protect_mothers_at_birth is true.
  // The associated member functions LockProtection, UnlockProtection, and
  // IsProtected appear in the GiveBirth and GetOlder Behaviors.
  bool protected_;
  // Single adult men can seek for a regular partnership. All male seeking
  // for regular partnership are then indexed, select the compound category of
  // their partner, and are mapped to females corresponding to the selected
  // category
  bool seek_regular_partnership_ = true;
  // Prevent newly infected to move on to next state the same year of infection
  bool newly_infected = false;
  // Number of casual partners
  int no_casual_partners_;
  // Tag to track if this person will be removed at the end of the iteration
  bool will_be_removed_;

  ///! The aguments below are currently either not used or repetitive.
  // // Stores if an agent is infected or not
  // bool infected_;
  // // Store the year when the agent got infected
  // float year_of_infection_;
  // Stores the ID of the mother. Useful to unlink child from mother, when child
  // dies.
  AgentPointer<Person> mother_ = nullptr;
  // Stores the IDs of the children. Useful, when mother migrates, and takes her
  // children. Unlink mother from child, when mother dies
  std::vector<AgentPointer<Person>> children_;
  // Stores the ID of the regular partner. Useful for infection in
  // serodiscordant regular relationships, and family migration.
  AgentPointer<Person> partner_ = nullptr;

  // The following function is used to avoid simultaneous modification of
  // related agents. (Tread safety)
  void CriticalRegion(std::vector<AgentPointer<Agent>>* aptrs) override {
    aptrs->push_back(GetAgentPtr<Agent>());
    if (partner_ != nullptr) {
      aptrs->push_back(partner_);
    }
    for (auto& child : children_) {
      aptrs->push_back(child);
    }
    if (mother_ != nullptr) {
      aptrs->push_back(mother_);
    }
  }

  void RemoveFromSimulation() override {
    // If has regular partner, end partnership
    if (hasPartner()) {
      SeparateFromPartner();
    }
    // If mother dies, children have no mother anymore
    for (uint64_t c = 0; c < children_.size(); c++) {
      children_[c]->mother_ = nullptr;
    }
    // If a child dies and has a mother, remove him from mother's list of
    // children
    if (mother_ != nullptr) {
      // std::cout << "A Child dies" << std::endl;
      mother_->RemoveChild(GetAgentPtr<Person>());
      // std::cout << " ==> Removed from mother's list of children" <<
      // std::endl;
    }
    Base::RemoveFromSimulation();
  }

  // Returns True if the agent is healthy
  bool IsHealthy() { return state_ == GemsState::kHealthy; }
  // AM: Added below functions for more detailed follow up of HIV state
  // Returns True if the agent is infected in acute state
  bool IsAcute() { return state_ == GemsState::kAcute; }
  // Returns True if the agent is infected in chronic state
  bool IsChronic() { return state_ == GemsState::kChronic; }
  // Returns True if the agent is infected and treated state
  bool IsTreated() { return state_ == GemsState::kTreated; }
  // Returns True if the agent is infected in failing treatement state
  bool IsFailing() { return state_ == GemsState::kFailing; }

  // Return True if recently infected at birth
  bool MTCTransmission() {
    return IsAcute() && transmission_type_ == TransmissionType::kMotherToChild;
  }
  // Return True if recently infected during casual mating
  bool CasualTransmission() {
    return IsAcute() && transmission_type_ == TransmissionType::kCasualPartner;
  }
  // Return True if recently infected during regular mating
  bool RegularTransmission() {
    return IsAcute() && transmission_type_ == TransmissionType::kRegularPartner;
  }

  // Return True if recently infected by an acute partner/mother
  bool AcuteTransmission() {
    return IsAcute() && infection_origin_state_ == GemsState::kAcute;
  }
  // Return True if recently infected by an chronic partner/mother
  bool ChronicTransmission() {
    return IsAcute() && infection_origin_state_ == GemsState::kChronic;
  }
  // Return True if recently infected by an treated partner/mother
  bool TreatedTransmission() {
    return IsAcute() && infection_origin_state_ == GemsState::kTreated;
  }
  // Return True if recently infected by an failing partner/mother
  bool FailingTransmission() {
    return IsAcute() && infection_origin_state_ == GemsState::kFailing;
  }
  // Return True if recently infected by a low risk partner
  bool LowRiskTransmission() { return IsAcute() && infection_origin_sb_ == 0; }
  // Return True if recently infected by an high risk partner
  bool HighRiskTransmission() { return IsAcute() && infection_origin_sb_ == 1; }

  // Returns True if the agent has high-risk socio-behaviours
  bool HasHighRiskSocioBehav() { return social_behaviour_factor_ == 1; }
  // Returns True if the agent is at low-risk socio-behaviours
  bool HasLowRiskSocioBehav() { return social_behaviour_factor_ == 0; }
  // Returns True if the agent is adult, is at least 15 years old
  bool IsAdult() { return age_ >= 15; }
  // Returns True if the agent is a male
  bool IsMale() { return sex_ == Sex::kMale; }
  // Returns True if the agent is a female
  bool IsFemale() { return sex_ == Sex::kFemale; }
 // Returns True if the agent is in location 1
  bool Isloc1() { return location_ == 1; }
// Returns True if the agent is in location 2
  bool Isloc2() { return location_ == 8; }

  // AM - Get Age Category from 0 to no_age_categories. 5-years interval
  // categories from min_age.
  int GetAgeCategory(size_t min_age, size_t no_age_categories) {
    int age_category;
    if (age_ >= min_age + (no_age_categories - 1) * 5) {
      age_category = no_age_categories - 1;
    } else {
      age_category = (int)(age_ - min_age) / 5;
    }
    // DEBUG:
    // std::cout << "age " << age_ << " --> age_category " << age_category << "
    // (min_age " << min_age << ",  no_age_categories " << no_age_categories <<
    // ")" << std::endl;
    return age_category;
  }

  // Note: This does not assign the child a mother. The mother must be assigned
  // separately.
  void AddChild(AgentPointer<Person> child) {
    /*if (child->location_ != location_){
        Log::Warning("Person::AddChild()", "Adding a child who is at a different
    location");
    }*/
    children_.push_back(child);
  }

  // Note: this does not remove the mother from the child. The mother must be
  // removed separately.
  void RemoveChild(AgentPointer<Person> child) {
    bool found = false;
    for (int c = 0; c < GetNumberOfChildren(); c++) {
      if (children_[c] == child) {
        found = true;
        // std::cout << "Before, nb_children = " << GetNumberOfChildren() <<
        // std::endl;
        children_.erase(children_.begin() + c);
        // std::cout << " => Found and removed child from mother's list of
        // children" << std::endl; std::cout << "Afrer, nb_children = " <<
        // GetNumberOfChildren() << std::endl;
        break;
      }
    }
    if (!found) {
      Log::Warning("Person::RemoveChild()",
                   "Child to be removed not found in mother's list of "
                   "children. Age = ",
                   child->age_, " Mother:", this->GetAgentPtr(),
                   " Age mother:", this->age_,
                   " Num children:", children_.size());
    }
  }

  void SetPartner(AgentPointer<Person> partner) {
    partner_ = partner;
    // Symetric relation
    partner_->partner_ = this->GetAgentPtr<Person>();
  }

  void SeparateFromPartner() {
    if (hasPartner()) {
      // Symetric relation. Start with partner while not nullptr
      partner_->partner_ = nullptr;
      // Set partner to nullptr
      partner_ = nullptr;
    } else {
      Log::Warning("Person::SeparateFromPartner()", "Person is single");
    }
  }

  void Relocate(size_t new_location) {
    location_ = new_location;

    if (sex_ == Sex::kFemale) {
      // Children (under 15yo) migrate with their mother
      int nb_children = GetNumberOfChildren();
      // std::cout << "I am a woman with "<< nb_children << " children and I
      // migrated to location " << location_<< std::endl;
      for (int c = 0; c < nb_children; c++) {
        if (children_[c]->age_ < 15) {
          /*if (old_location != children_[c]->location_){
              Log::Warning("RandomMigration::Run()", "child and mother had
          different locations BEFORE MIGRATION. Child's age = ",
          children_[c]->age_, " ==> ", old_location, " vs. ",
          children_[c]->location_);
          }*/
          // std::cout << " ==> I am a child (" << children_[c]->age_
          // << ") and I will migrate with my mother from location " <<
          // children_[c]->location_ << " to " << location_ <<
          // std::endl;
          children_[c]->location_ = location_;
        }
      }
      // DEBUG : Check that all children migrated with Mother
      for (int c = 0; c < nb_children; c++) {
        if (children_[c]->age_ < 15) {
          if (children_[c]->location_ != location_) {
            Log::Warning("RandomMigration::Run()",
                         "DEBUG: child and mother have different locations "
                         "AFTER MIGRATION. Child's age = ",
                         children_[c]->age_);
          }
        }
      }
    } else if (hasPartner()) {
      // If a man engaged in a regular partnership relocates, his female partner
      // relocates too.
      partner_->Relocate(new_location);
    }
  }

  bool IsParentOf(AgentPointer<Person> child) {
    bool found = false;
    for (int c = 0; c < GetNumberOfChildren(); c++) {
      if (children_[c] == child) {
        found = true;
        break;
      }
    }
    return found;
  }

  bool IsChildOf(AgentPointer<Person> mother) { return mother_ == mother; }

  bool hasPartner() { return partner_ != nullptr; }
  bool IsPartnerOf(AgentPointer<Person> partner) { return partner_ == partner; }

  int GetNumberOfChildren() { return children_.size(); }

  // Restet the counter of casual partners to zero
  void ResetCasualPartners() { no_casual_partners_ = 0; }

  // Activates the protection of an agent against death.
  void LockProtection() { protected_ = true; }
  // Deactivates the protection of an agent against death.
  void UnockProtection() { protected_ = false; }
  // Returns if an agent is protected against death.
  bool IsProtected() { return protected_; }
};

}  // namespace hiv_malawi
}  // namespace bdm

#endif  // PERSON_H_
