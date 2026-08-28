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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingVerbose.hh"

int main(int argc, char** argv)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  G4String macro;
  G4int phys_option = 2;
  for (G4int ii = 1; ii < argc; ii = ii + 2)
  {
    if (G4String(argv[ii]) == "-m")
    {
      macro = argv[ii + 1];
    }
    else if (G4String(argv[ii]) == "-p")
    {
      phys_option = G4UIcommand::ConvertToInt(argv[ii + 1]);
    }
  }
  // Set mandatory initialization classes
  runManager->SetUserInitialization(new PhysicsList(phys_option));
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new ActionInitialization());
  runManager->SetUserInitialization(new ActionInitialization());

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc == 1)  // Define UI session for interactive mode.
  {
#ifdef _WIN32
    G4UIsession* session = new G4UIterminal();
#else
    G4UIsession* session = new G4UIterminal(new G4UItcsh);
#endif
    UI->ApplyCommand("/control/execute test60.mac");
    session->SessionStart();
    delete session;
  }
  else
  {
    G4String command = "/control/execute ";
    UI->ApplyCommand(command + macro);
  }
  delete runManager;
  return 0;
}
