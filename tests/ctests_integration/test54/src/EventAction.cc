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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "EventMessenger.hh"
#include "G4RunManager.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:fDrawFlag("none")
{
  fEventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*evt*/)
{
  fEnergyDeposit  = 0.;
  fTrackLenCharged = fTrackLenNeutral = 0.; 
  fNbStepsCharged = fNbStepsNeutral = 0;
  fTransmitFlag   = fReflectFlag    = 0;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  Run* run = static_cast<Run*>(
	      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  run->AddEnergy(fEnergyDeposit);
  run->AddTrakLenCharg(fTrackLenCharged);
  run->AddTrakLenNeutr(fTrackLenNeutral);
  
  run->CountStepsCharg(fNbStepsCharged);
  run->CountStepsNeutr(fNbStepsNeutral);
  
  run->CountTransmit (fTransmitFlag);
  run->CountReflect  (fReflectFlag);
  
  run->FillHisto(1,fEnergyDeposit);
  
  if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
      for (G4int i=0; i<n_trajectories; i++) 
	{ G4Trajectory* trj = (G4Trajectory*)
	    ((*(evt->GetTrajectoryContainer()))[i]);
	  if (fDrawFlag == "all") trj->DrawTrajectory();
	  else if ((fDrawFlag == "charged")&&(trj->GetCharge() != 0.))
	    trj->DrawTrajectory();
	  else if ((fDrawFlag == "neutral")&&(trj->GetCharge() == 0.))
	    trj->DrawTrajectory();
	  
	}
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

