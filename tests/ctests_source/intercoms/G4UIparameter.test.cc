//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "G4UIparameter.hh"

#include "G4UIcommandStatus.hh"

#include <gtest/gtest.h>

// Identified problems
// - Parameter can be untyped!
// - Parameter can have range and candidates, but range will always win!
//   having both in non-sensical anyway!
// - Setting candidates/range non-sensical for boolean parameters
// - Only named parameters can have ranges

TEST(G4UIparameter, Basic)
{
  // Basic integral parameter;
  {
    G4UIparameter p('i');
    // Must be same type
    ASSERT_EQ(p.GetParameterType(), 'i');

    // No default value is o.k., but should be empty string
    EXPECT_EQ(p.GetDefaultValue(), "");

    // Must get G4UIcommandStatus::fParameterUnreadable if we pass a string
    // that isn't parsable to and integer
    // NB: since these log to G4Cout/cerr, will want a test harness
    // to intercept these (e.g. dev/null them if wnat to suppress, or a checker
    // that outputs are as expected.)
    auto vals = {"12.", "12f", "-21.", "21 3", " 32", "45 "};
    for (const auto& val : vals)
    {
      EXPECT_EQ(p.CheckNewValue(val), G4UIcommandStatus::fParameterUnreadable);
    }
    // Equally, with no range, no candidates, a valid parameter must return 0
    EXPECT_EQ(p.CheckNewValue("42"), 0);
  }
}

TEST(G4UIparameter, Candidates)
{
  // Integral parameter with candidate list
  G4UIparameter p('i');
  p.SetParameterCandidates("1 2 3");

  // CheckNewValue must return 0 for all values on list
  for (const char* c : {"1", "2", "3"})
  {
    EXPECT_EQ(p.CheckNewValue(c), 0);
  }

  // CheckNewValue must return G4UIcommandStatus::fParameterOutOfCandidates
  EXPECT_EQ(p.CheckNewValue("-1"), G4UIcommandStatus::fParameterOutOfCandidates);
}

TEST(G4UIparameter, Ranges)
{
  // Integral parameter with range
  G4UIparameter p('i');
  p.SetParameterName("id");
  p.SetParameterRange("id > 0 && id < 42");

  // Good value must return 0
  EXPECT_EQ(p.CheckNewValue("23"), 0);

  // Parameter outside of range must return G4UIcommandStatus::fParameterOutOfRange
  EXPECT_EQ(p.CheckNewValue("0"), G4UIcommandStatus::fParameterOutOfRange);
  EXPECT_EQ(p.CheckNewValue("938"), G4UIcommandStatus::fParameterOutOfRange);
}
