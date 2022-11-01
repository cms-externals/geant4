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
//

#include "ExDivDetectorConstruction.hh"
#include "ExDivActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "ExDivTesterBox.hh"
#include "ExDivTesterTubs.hh"
#include "ExDivTesterCons.hh"
#include "ExDivTesterTrd.hh"
#include "ExDivTesterPara.hh"
#include "ExDivTesterPolycone.hh"
#include "ExDivTesterPolyhedra.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Select interactive mode and define UI session
  //
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);

  G4String theSolidTypeStr="box";
  G4String thePVTypeStr="division";
  G4String thePosTypeStr="normal";
  std::vector<G4String> theExtraPars;

  // First argument is the type of divisioning replica or division
  // the second argument is the type of solid
  std::vector<std::string> vsarg;
  for( G4int jj = 0; jj < argc; jj ++ )
  {
    vsarg.push_back( std::string(argv[jj] ) );
  } 
  
  G4int narg = vsarg.size();
  if( narg == 1 )
  {
    G4cout << "!!! No input division type provided. Defaulting to 'division' "
           << G4endl;
    G4cout << "!!! No positioning type provided. Defaulting to 'normal' "
           << G4endl;
    G4cout << "!!! No input solid type provided. Defaulting to 'box' "
           << G4endl;
  }
  else
  {
    if( narg == 2 )
    {
      G4cout << "!!! No positioning type provided. Defaulting to 'normal' "
             << G4endl;
      G4cout << "!!! No input solid type provided. Defaulting to 'box' "
             << G4endl;
    }
    else
    {
      if( narg == 3 )
      {
        G4cout << "!!! No input solid type provided. Defaulting to 'box' "
               << G4endl;
      }
      else
      {
        theSolidTypeStr = G4String(vsarg[3]);
      }
        thePosTypeStr = G4String(vsarg[2]);
    }
    thePVTypeStr = G4String(vsarg[1]);
  }
  if( narg > 4 )
  {
    for( G4int ii = 4; ii < narg; ii++ )
    {
      theExtraPars.push_back( vsarg[ii] );
    }
  }

  // Use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);
  
  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // UserInitialization classes (mandatory)
  ExDivDetectorConstruction* ExDivdetector =
    new ExDivDetectorConstruction(theSolidTypeStr, thePVTypeStr,
                                  thePosTypeStr, theExtraPars );
  runManager->SetUserInitialization(ExDivdetector);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // Set user action classes
  runManager->SetUserInitialization(new ExDivActionInitialization());

  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Initialise G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager 
  G4UImanager * UIman = G4UImanager::GetUIpointer();  

  // Define (G)UI terminal for interactive mode
  if(argc>=1)
  { 
    UIman->ApplyCommand("/control/execute vis.mac");    
    ui->SessionStart();
    delete ui;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UIman->ApplyCommand(command+fileName);
  }

  delete visManager;
  delete runManager;

  return 0;
}
