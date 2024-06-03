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
/*
 * G4MoleculeTableTest.cc
 *
 *  Created on: Jun 23, 2015
 *      Author: mkaramit
 */

#include "G4MoleculeTable.hh"
#include "G4Electron_aq.hh"
#include "G4MolecularConfiguration.hh"

int main()
{
  G4MoleculeTable* moleculeTable =
      G4MoleculeTable::Instance()->GetMoleculeTable();

  //---------------------------------------------------------------------------
  // Test 1

  G4MolecularConfiguration* e_aq =
      moleculeTable->CreateConfiguration("e_aq", G4Electron_aq::Definition());

  G4MolecularConfiguration* e_aq2 = moleculeTable->GetConfiguration("e_aq");

  assert(e_aq == e_aq2);

  G4MoleculeDefinition* molDef1 =
      moleculeTable->CreateMoleculeDefinition("mol1", 3);

  G4MolecularConfiguration* mol1 = moleculeTable->CreateConfiguration("mol1",
                                                                      molDef1);

  G4MolecularConfiguration* mol1_bis = moleculeTable->GetConfiguration("mol1");

  assert(mol1 == mol1_bis);

  assert(molDef1->GetCharge() == 0);
  assert(mol1_bis->GetCharge() == 0);

  G4cout << mol1_bis->GetName() << G4endl;
  assert(mol1_bis->GetName() == "mol1^0");

  //---------------------------------------------------------------------------
  // Test 2
  G4MoleculeDefinition* molDef2 =
      moleculeTable->CreateMoleculeDefinition("mol2", 5);
  G4MolecularConfiguration* mol2 = moleculeTable->CreateConfiguration("mol2",
                                                                      molDef2,
                                                                      -1,
                                                                      9.);
  assert(mol2->GetCharge() == -1);
  assert(mol2->GetName() == "mol2^-1");

  //---------------------------------------------------------------------------
  // Test 3
  G4MoleculeDefinition* molDef3 =
      moleculeTable->CreateMoleculeDefinition("mol3", 5);
  G4MolecularConfiguration* mol3 =
      moleculeTable->CreateConfiguration("mol3_bis", molDef3, "labelMol3", -1);

  G4cout << mol3->GetCharge() << G4endl;
  assert(mol3->GetCharge() == -1);
  G4cout << mol3->GetName() << G4endl;
  assert(mol3->GetName() == "mol3^-1");
  assert(mol3->GetLabel() == "labelMol3");

  G4cout << "test OK" << G4endl;
}

