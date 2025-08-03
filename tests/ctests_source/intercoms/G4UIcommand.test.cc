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

#include "G4UIcommand.hh"

#include "G4UIcommandStatus.hh"
#include "G4UImessenger.hh"

#include <gtest/gtest.h>

// Identified issues
// - command without parameter is to call member functions etc that don't
//   take parameters.
// - Implicit limit of 99 parameters due to range error codes.
// Not issues, but clarity
// - Command with parameters can have ranges set on both parameters and combination
//   Possible to construct these such that there is no possibly valid value, and no
//   way to check overall range for that. However, is logical to have this construct
//   e.g. Params x, y. x > 5, y > 10, x < y (though then should be possible to construct
//   overall constraint!)

TEST(G4UIcommand, Parameters)
{
  G4UImessenger mock;
  G4UIcommand c("", &mock);
  {
    auto xp = new G4UIparameter('i');
    xp->SetParameterName("x");
    xp->SetParameterRange("x > 8");
    c.SetParameter(xp);

    auto yp = new G4UIparameter('i');
    yp->SetParameterName("y");
    yp->SetParameterRange("y > 11");
    c.SetParameter(yp);
  }
  c.SetRange("x < y");

  ASSERT_EQ(c.GetParameterEntries(), 2);
  // Range satisfying both is fine
  EXPECT_EQ(c.DoIt("11 12"), 0);
  // Parameter range failure is 300 + param index
  EXPECT_EQ(c.DoIt("9 10"), 301);
  // Failure of command range is 399
  EXPECT_EQ(c.DoIt("24 20"), 399);
}
