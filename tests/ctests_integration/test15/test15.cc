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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//            Underground Advanced example main program
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// main program
// --------------------------------------------------------------

#include "G4GeometryManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisExecutive.hh"

#include "Test15DetectorConstruction.hh"
#include "Test15ShellDetectorConstruction.hh"
// qel #include "Test15PhysicsList.hh"
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"

#include "Test15ActionInitialization.hh"
#include "Test15EventAction.hh"
#include "Test15StackingAction.hh"
#include "Test15SteppingAction.hh"

// Files specific for biasing and scoring
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"

#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "LBE.hh"
#include "QBBC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_AllHP.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_INCLXX.hh"
#include "QGSP_INCLXX_HP.hh"
#include "QGS_BIC.hh"
#include "Shielding.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
void PrintUsage()
{
  G4cerr << " Usage: " << G4endl;
#ifdef G4MULTITHREADED
  G4cerr << " Test15 [-m macro ] [-p physics] [-t nThreads] [-s seed]" << G4endl;
#else
  G4cerr << " Test15 [-m macro ] [-p physics] [-s seed]" << G4endl;
  G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
#endif
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  G4String macro;
  G4String physics = "QGSP_BIC_HP";
  G4String gdml_out = "";
  G4int seed = 3517801;
  G4int nThreads = 2;
  for (G4int i = 1; i < argc; i = i + 2)
  {
    if (G4String(argv[i]) == "-m") macro = argv[i + 1];
#ifdef G4MULTITHREADED
    else if (G4String(argv[i]) == "-t")
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
#endif
    else if (G4String(argv[i]) == "-p")
      physics = argv[i + 1];
    else if (G4String(argv[i]) == "-s")
      seed = G4UIcommand::ConvertToInt(argv[i + 1]);
    else
    {
      PrintUsage();
      return 1;
    }
  }

  G4Random::setTheSeed(seed);

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if (!macro.size())
  {
    ui = new G4UIExecutive(argc, argv);
  }

  G4Random::saveEngineStatus("RandomSeed.conf");

  auto runManager = G4RunManagerFactory::CreateRunManager();

  if (G4RunManagerFactory::GetDefault() != "Serial")
  {
    G4cout << "Warning: forcing number of threads to be " << nThreads << G4endl;
    runManager->SetNumberOfThreads(nThreads);
  }

  // set mandatory initialization classes
  Test15DetectorConstruction* detector = new Test15DetectorConstruction;
  runManager->SetUserInitialization(detector);

  G4String parallelName("ParallelShellWorld");
  Test15ShellDetectorConstruction* pdet = new Test15ShellDetectorConstruction(parallelName);

  detector->RegisterParallelWorld(pdet);

  G4GeometrySampler pgs(pdet->GetWorldVolume(), "neutron");

  pgs.SetParallel(true);

  G4VModularPhysicsList* physlist;
  if (physics == "FTFP_BERT")
    physlist = new FTFP_BERT;
  else if (physics == "FTFP_BERT_HP")
    physlist = new FTFP_BERT_HP;
  else if (physics == "QGSP_BERT")
    physlist = new QGSP_BERT;
  else if (physics == "QGSP_BIC_HP")
    physlist = new QGSP_BIC_HP;
  else if (physics == "QGSP_BIC_AllHP")
    physlist = new QGSP_BIC_AllHP;
  else if (physics == "QGSP_BIC")
    physlist = new QGSP_BIC;
  else if (physics == "QGS_BIC")
    physlist = new QGS_BIC;
  else if (physics == "QGSP_INCLXX")
    physlist = new QGSP_INCLXX;
  else if (physics == "QGSP_INCLXX_HP")
    physlist = new QGSP_INCLXX_HP;
  else if (physics == "QBBC")
    physlist = new QBBC;
  else if (physics == "LBE")
    physlist = new LBE;
  else if (physics == "Shielding")
    physlist = new Shielding;
  else
  {
    G4cout << G4endl << G4endl << " Physics list: " << physics << " not recognised!!! Mis-typed??? "
           << G4endl << G4endl << G4endl;
    PrintUsage();
    return 1;
  }

  G4cout << "!!!!!!! TEST15 physlist: " << physics << G4endl;

  physlist->RegisterPhysics(new G4ImportanceBiasing(&pgs, parallelName));
  physlist->RegisterPhysics(new G4ParallelWorldPhysics(parallelName));

  runManager->SetUserInitialization(physlist);

  // Set user action classes through Worker Initialization
  //
  Test15ActionInitialization* actions = new Test15ActionInitialization;
  runManager->SetUserInitialization(actions);

  runManager->Initialize();

  pdet->CreateImportanceStore();

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if (macro.size())
  {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else
  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI())
    {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

#ifdef Test15ENV_GPS_USE
  G4cout << " Using GPS and not Test15 gun " << G4endl;
#else
  G4cout << " Using the Test15 gun " << G4endl;
#endif

  // open geometry for clean biasing stores clean-up
  //
  G4GeometryManager::GetInstance()->OpenGeometry();

  pgs.ClearSampling();

  delete visManager;
  delete runManager;

  return 0;
}
