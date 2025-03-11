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
 * G4MolecularConfigurationTest.cc
 *
 *  Created on: Jun 24, 2015
 *      Author: mkaramit
 */

#include "G4MolecularConfiguration.hh"
#include "G4Electron_aq.hh"

int main()
{
  G4Electron_aq* e_aq = G4Electron_aq::Definition();

  const G4ElectronOccupancy* eOcc = e_aq->GetGroundStateElectronOccupancy();

  G4ElectronOccupancy occ = *eOcc;
  occ.AddElectron(0,1);

  G4MolecularConfiguration* conf1 = G4MolecularConfiguration::GetOrCreateMolecularConfiguration(e_aq,occ);

  assert(*(conf1->GetElectronOccupancy()) == occ);

  G4ElectronOccupancy occ2 = *eOcc;
  occ2.AddElectron(0,2); // does not make any sens, just for testing

  G4MolecularConfiguration* conf2 = G4MolecularConfiguration::GetOrCreateMolecularConfiguration(e_aq,occ2);

  assert(*(conf2->GetElectronOccupancy()) == occ2);

  G4cout << "test OK" << G4endl;
}


