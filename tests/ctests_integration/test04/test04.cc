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

#include "G4ParticleHPManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4ThermalNeutrons.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ios.hh"

#include "FTFP_BERT_HP.hh"
#include "Tst04ActionInitialization.hh"
#include "Tst04DetectorConstruction.hh"
#include "Tst04PhysicsList.hh"

int main(int argc, char** argv)
{
  // Run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();

  if (G4RunManagerFactory::GetDefault() != "Serial")
  {
    G4cout << "Warning: forcing number of threads to be 3" << G4endl;
    runManager->SetNumberOfThreads(3);
  }

  G4VModularPhysicsList* phys = new FTFP_BERT_HP();
  phys->RegisterPhysics(new G4ThermalNeutrons(1));

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst04DetectorConstruction);
  runManager->SetUserInitialization(phys);
  runManager->SetUserInitialization(new Tst04ActionInitialization);

  // Activate production of fission fragments in neutronHP
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments(true);

  if (argc == 1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  G4Random::showEngineStatus();
  delete runManager;
  return 0;
}
