/// \file Hadr08.cc
/// \brief Main program of the hadronic/Hadr08 example
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "FTFP_INCLXX.hh"
#include "G4UImanager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "ActionInitialization.hh"
#include "G4UIExecutive.hh"
#include "G4GenericBiasingPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main( int argc, char** argv ) {

  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  #ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4);
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|              Constructing MT run manager              |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  #else
  G4RunManager * runManager = new G4RunManager;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|        Constructing sequential run manager            |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  #endif

  G4VUserDetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization( detector );

  FTFP_BERT* physicsList = new FTFP_BERT;
  //FTFP_INCLXX* physicsList = new FTFP_INCLXX;

  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics;
  biasingPhysics->Bias( "proton" );
  biasingPhysics->Bias( "neutron" );
  biasingPhysics->Bias( "pi+" );
  biasingPhysics->Bias( "pi-" );
  physicsList->RegisterPhysics( biasingPhysics );

  runManager->SetUserInitialization( physicsList );
  runManager->SetUserInitialization( new ActionInitialization );

  runManager->Initialize();

  if ( ui ) {
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
