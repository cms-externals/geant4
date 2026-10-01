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
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP_BERT_HP_VI
//
// Author: V. Ivanchenko (2023)
// Description: Alternative neutron HP models and processes
//
//----------------------------------------------------------------------------
//

#include "QGSP_BERT_HP_VI.hh"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysicsVI.hh"
#include "G4HadronPhysicsQGSP_BERT_HP_VI.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <iomanip>

QGSP_BERT_HP_VI::QGSP_BERT_HP_VI(G4int ver)
{
  if (ver > 0)
  {
    G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT_HP_VI" << G4endl;
    G4cout << G4endl;
  }

  defaultCutValue = 0.7 * CLHEP::mm;
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics(new G4EmStandardPhysics(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));
  RegisterPhysics(new G4RadioactiveDecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysicsVI(ver));

  // Hadron Physics
  RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP_VI(ver));

  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  RegisterPhysics(new G4IonElasticPhysics(ver));
}

void QGSP_BERT_HP_VI::SetCuts()
{
  if (verboseLevel > 1)
  {
    G4cout << "QGSP_BERT_HP_VI::SetCuts:";
  }
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types

  SetCutsWithDefault();

  // Set proton cut value to 0 for producing low energy recoil nucleus
  SetCutValue(0, "proton");
}
