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
// D Chiappara, GAP Cirrone, G Petringa
//
// ** RADIOBIO  example **
// make easy the calculation of the main radiobiology quantities used in
// charged particles radiation therapy: depth dose, Fluence in depth,
// LET (Linear Energy Transfer) and RBE (Relative Biological Effectiveness)
//

#include "G4RunManagerFactory.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "ActionInitialization.hh"
#include "BeamFluence.hh"
#include "DetectorConstruction.hh"
#include "Dose.hh"
#include "LET.hh"
#include "Manager.hh"
#include "PhysicsList.hh"
#include "RBE.hh"

// #ifdef G4VIS_USE
#include "G4VisExecutive.hh"
// #endif

// #ifdef G4UI_USE
#include "G4UIExecutive.hh"
// #endif

#include <cstdlib>
#include <ctime>

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Instantiate the G4Timer object, to monitor the CPU time spent for
  // the entire execution
  G4Timer* theTimer = new G4Timer();
  // Start the benchmark
  theTimer->Start();

  // Set the Random engine
  // The following guarantees random generation also for different runs
  // in multithread
  CLHEP::RanluxEngine defaultEngine(1234567, 4);
  G4Random::setTheEngine(&defaultEngine);
  const G4long seed = static_cast<G4long>(std::time(nullptr));
  G4Random::setTheSeed(seed);

  // Create the run manager. Uses by default serial or parallel computation.
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Create Scoring manager
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
  scoringManager->SetVerboseLevel(1);

  // Creating PhysicsList
  G4VModularPhysicsList* phys = new RadioBio::PhysicsList();

  // Set mandatory initialization classes
  auto* det = new RadioBio::DetectorConstruction();
  runManager->SetUserInitialization(det);

  // Creation of Manager
  RadioBio::Manager* radioBioManager = RadioBio::Manager::CreateInstance();

  // Create and Register Radiobiological quantities
  auto* dose = new RadioBio::Dose();
  radioBioManager->Register(dose, "Dose");

  auto* beamFluence = new RadioBio::BeamFluence();
  radioBioManager->Register(beamFluence, "BeamFluence");

  auto* let = new RadioBio::LET();
  radioBioManager->Register(let, "LET");

  auto* rbe = new RadioBio::RBE();
  radioBioManager->Register(rbe, "RBE");

  // Initialisation of physics
  runManager->SetUserInitialization(phys);

  // Initialisation of the Actions
  runManager->SetUserInitialization(new RadioBio::ActionInitialization(det));

  // Initialise the Visualisation
  // #ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  // #endif

  //** Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4int commandStatus = 0;
  if (ui == nullptr) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    commandStatus = UImanager->ApplyCommand(command + fileName);
  }

  else {
    commandStatus = UImanager->ApplyCommand("/control/execute macro/visualize.mac");

    ui->SessionStart();
    delete ui;
  }

  // Save data in RadioBioManager
  if (commandStatus == 0) {
    radioBioManager->StoreAll();
  }

  // Stop the benchmark here
  theTimer->Stop();

  G4cout << "The simulation took: " << theTimer->GetRealElapsed() << " s to run (real time)"
         << G4endl;

  delete theTimer;
  delete visManager;
  delete runManager;
  delete radioBioManager;

  return commandStatus == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
