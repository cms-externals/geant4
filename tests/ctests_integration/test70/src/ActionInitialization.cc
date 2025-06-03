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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "ActionInitialization.hh"
#include <G4Scheduler.hh>
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "G4RunManager.hh"
#include "StackingAction.hh"
#include "G4DNAChemistryManager.hh"
#include "ITTrackingInteractivity.hh"
#include "ITSteppingAction.hh"
#include "ITTrackingAction.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4MoleculeCounter.hh"
#include "TimeStepAction.hh"
#include "G4H2O.hh"
#include "EventAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization()
  : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  // In MT mode, to be clearer, the RunAction class for the master thread might
  // be different than the one used for the workers. This RunAction will be
  // called before and after starting the workers. For more details, please
  // refer to :
  // https://twiki.cern.ch/twiki/bin/view/Geant4/Geant4MTForApplicationDevelopers
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  G4cout << "Build for = " << G4RunManager::GetRunManager()->GetRunManagerType()
         << G4endl;
  BuildMoleculeCounter();  // fUse should be localthread
  auto primGenAction = new PrimaryGeneratorAction();
  SetUserAction(primGenAction);
  //------------------------------------------------------------------
  // Set optional user action classes

  bool chemistryFlag =
    G4DNAChemistryManager::Instance()->IsChemistryActivated();

  SetUserAction(new RunAction());
  SetUserAction(new TrackingAction());
  SetUserAction(new EventAction());
  SetUserAction(new SteppingAction(primGenAction));
  SetUserAction(new StackingAction());
  if(chemistryFlag)
  {
    G4Scheduler::Instance()->SetUserAction(new TimeStepAction());
    G4Scheduler::Instance()->SetVerbose(1);
    auto itInteractivity = new ITTrackingInteractivity();
    itInteractivity->SetUserAction(new ITSteppingAction());
    itInteractivity->SetUserAction(new ITTrackingAction());
    G4Scheduler::Instance()->SetInteractivity(itInteractivity);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildMoleculeCounter() const
{
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeEvent(true);
  G4MoleculeCounterManager::Instance()->SetResetCountersBeforeRun(true);
  G4MoleculeCounterManager::Instance()->SetAccumulateCounterIntoMaster(false);

  auto counter = std::make_unique<G4MoleculeCounter>();
  counter->SetTimeComparer(G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(1 * ps));
  counter->IgnoreMolecule(G4H2O::Definition());
  G4MoleculeCounterManager::Instance()->RegisterCounter(std::move(counter));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
