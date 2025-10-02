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
// ----------------------------------------------------------------------
#include "G4OrderedTable.hh"

#include <gtest/gtest.h>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(G4OrderedTable, Readback)
{
  //
  // Create a G4OrderedTable object
  //
  const size_t Imax = 10;
  G4OrderedTable aTable(Imax);
  auto pl = aTable.begin();
  for (size_t I = 0; I < Imax; I++) {
    auto* aVector = new G4DataVector();
    *pl++ = aVector;
    for (size_t J = 0; J <= I; J++)
      aVector->push_back(G4double(J));
  }

  // Now access the data contained in the table
  // Store in file in ascii mode
  EXPECT_TRUE(aTable.Store("OrdTable.txt", true)) << "G4OrderedTable could not write table to file";

  // Create new table and fill from file
  G4OrderedTable txtTable;
  EXPECT_TRUE(txtTable.Retrieve("OrdTable.txt", true)) << "G4OrderedTable could not read from file";

  // Readback should succeed
  EXPECT_EQ(aTable.size(), txtTable.size()) << "Read back G4OrderedTable does not have same size";
  // Contents are pointers, so must step through each one and vector compare
  // NOTE: to be simplified with matchers etc
  for (size_t I = 0; I < Imax; I++) {
    EXPECT_EQ(*aTable[I], *txtTable[I])
      << "Read back G4OrderedTable has mismatched elements in table " << I;
  }

  // Store in file in binary mode
  aTable.Store("OrdTable.bin", false);

  // Create new table and fill from file
  G4OrderedTable binTable;
  binTable.Retrieve("OrdTable.bin", false);
  // Readback should succeed
  EXPECT_EQ(aTable.size(), binTable.size()) << "Read back G4OrderedTable does not have same size";
  // Contents are pointers, so must step through each one and vector compare
  // NOTE: to be simplified with matchers etc
  for (size_t I = 0; I < Imax; I++) {
    EXPECT_EQ(*aTable[I], *binTable[I])
      << "Read back G4OrderedTable has mismatched elements in table " << I;
  }

  // Cleanup
  aTable.clearAndDestroy();
  txtTable.clearAndDestroy();
  binTable.clearAndDestroy();
}
