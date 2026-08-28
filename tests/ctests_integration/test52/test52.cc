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

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

int main(int argc, char** argv)
{
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  PhysicsList* physics = new PhysicsList();
  runManager->SetUserInitialization(physics);

  PrimaryGenerator* source = new PrimaryGenerator();
  runManager->SetUserAction(source);

  RunAction* runAction = new RunAction();
  runManager->SetUserAction(runAction);

  EventAction* eventAction = new EventAction();
  runManager->SetUserAction(eventAction);

  SteppingAction* steppingAction = new SteppingAction();
  runManager->SetUserAction(steppingAction);

  TrackingAction* trackingAction = new TrackingAction();
  runManager->SetUserAction(trackingAction);

  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc == 2)
  {
    G4String fileName = argv[1];
    G4cout << "INFORMATION:  Commands in file " << fileName << " used to control simulation."
           << G4endl;
    UI->ApplyCommand("/control/execute " + fileName);
  }

  delete runManager;
  return 0;
}
