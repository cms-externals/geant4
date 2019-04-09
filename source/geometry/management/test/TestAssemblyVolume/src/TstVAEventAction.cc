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
// --------------------------------------------------------------

#include "TstVAEventAction.hh"

#include "TstVAEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

TstVAEventAction::TstVAEventAction()
 : calorimeterCollID(-1), drawFlag("all"),
   printModulo(1), eventMessenger(NULL)
{
  eventMessenger = new TstVAEventActionMessenger(this);
}

TstVAEventAction::~TstVAEventAction()
{
  delete eventMessenger;
}

void TstVAEventAction::BeginOfEventAction(const G4Event*)
{
}

void TstVAEventAction::EndOfEventAction(const G4Event* evt)
{

  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");

    G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer)
    {
      n_trajectories = trajectoryContainer->entries();
    }
    
    for( int i = 0; i < n_trajectories; i++ )
    {
      (*trajectoryContainer)[i]->DrawTrajectory();
    }
  }


/*
  G4int evtNb = evt->GetEventID();
  
  // extracted from hits, compute the total energy deposit (and total charged
  // track length) in all absorbers and in all gaps

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();

  ExN03CalorHitsCollection* CHC = NULL;
  G4int n_hit = 0;
  G4double totEAbs=0, totLAbs=0, totEGap=0, totLGap=0;
    
  if (HCE) CHC = (ExN03CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
    {
     n_hit = CHC->entries();
     for (G4int i=0;i<n_hit;i++)
	{
	  totEAbs += (*CHC)[i]->GetEdepAbs(); 
	  totLAbs += (*CHC)[i]->GetTrakAbs();
	  totEGap += (*CHC)[i]->GetEdepGap(); 
	  totLGap += (*CHC)[i]->GetTrakGap();
	}
     }
   
   // print this information event by event (modulo n)  	
	  
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;	

    G4cout
       << "   Absorber: total energy: " << std::setw(7) << G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << std::setw(7) << G4BestUnit(totLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7) << G4BestUnit(totEGap,"Energy")
       << "       total track length: " << std::setw(7) << G4BestUnit(totLGap,"Length")
       << G4endl;
	  
    G4cout << "\n     " << n_hit
	   << " hits are stored in ExN03CalorHitsCollection." << G4endl;  	     
  }
  // extract the trajectories and draw them
  
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") trj->DrawTrajectory();
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory();
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                  trj->DrawTrajectory();
        }
  }             
*/
}
