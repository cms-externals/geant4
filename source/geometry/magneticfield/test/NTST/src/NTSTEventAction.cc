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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Timer.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "NTSTEventAction.hh"
#include "NTSTEventActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTEventAction::NTSTEventAction()
  : EventTime(new G4Timer()),
   MeanUserEventTime(0),  RmsUserEventTime(0),
   MeanRealEventTime(0),  RmsRealEventTime(0), 
   NumberOfEvents(0), 
   MeanVertices(0), RmsVertices(0),
   MeanTracks(0), RmsTracks(0), 
   drawFlag("all"), eventMessenger(NULL)
{
  eventMessenger = new NTSTEventActionMessenger(this);
}

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTEventAction::~NTSTEventAction()
{
  if (NumberOfEvents>0) {
    G4cout << "### Processed number of events for all runs: " 
	   << NumberOfEvents << G4endl;
    MeanUserEventTime = MeanUserEventTime / NumberOfEvents;
    MeanRealEventTime = MeanRealEventTime / NumberOfEvents;
    RmsUserEventTime = RmsUserEventTime / NumberOfEvents;
    RmsRealEventTime = RmsRealEventTime / NumberOfEvents;
    G4double ErrUserEventTime = 
      std::sqrt((RmsUserEventTime - MeanUserEventTime*MeanUserEventTime)
	   /NumberOfEvents);
    G4double ErrRealEventTime = 
      std::sqrt((RmsRealEventTime - MeanRealEventTime*MeanRealEventTime)
	   /NumberOfEvents);
    G4cout << std::setprecision(3) 
	   << "### Event user time = " << std::setw(6) << MeanUserEventTime
	   << std::setprecision(3) 
	   << " +- " << std::setw(6)
           << ErrUserEventTime << " (sec) " << G4endl;   
    G4cout << std::setprecision(3)
	   << "### Event real time = " << std::setw(6) << MeanRealEventTime
	   << std::setprecision(3) 
	   << " +- " << std::setw(6)
           << ErrRealEventTime << " (sec) " << G4endl;   
  
    MeanVertices = MeanVertices / NumberOfEvents;
    RmsVertices = RmsVertices / NumberOfEvents;
    MeanTracks = MeanTracks / NumberOfEvents;
    RmsTracks = RmsTracks / NumberOfEvents;
    G4double ErrVertices = 
      std::sqrt((RmsVertices - MeanVertices*MeanVertices)/NumberOfEvents);
    G4double ErrTracks =
      std::sqrt((RmsTracks - MeanTracks*MeanTracks)/NumberOfEvents);

    G4cout << std::setprecision(3)
	   << "### Number of Vertices = " << std::setw(6) << MeanVertices << " +- "
	   << std::setw(6) << ErrVertices
	   << " Number of Tracks = " << std::setw(6) << MeanTracks << " +- "
	   << std::setw(6) << ErrTracks << G4endl;

				
  }
  delete eventMessenger;
  delete EventTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NTSTEventAction::BeginOfEventAction(const G4Event* ) // evt)
{  
  EventTime->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NTSTEventAction::EndOfEventAction(const G4Event* evt)
{
  EventTime->Stop();
  G4cout << "### Event " << std::setprecision(3) 
	 << evt->GetEventID()+1 << " " << *EventTime << G4endl;

  // event statistics
  G4double ElapsedUserTime = EventTime->GetUserElapsed();
  G4double ElapsedRealTime = EventTime->GetRealElapsed();
  MeanUserEventTime += ElapsedUserTime;
  MeanRealEventTime += ElapsedRealTime;
  RmsUserEventTime += ElapsedUserTime*ElapsedUserTime;
  RmsRealEventTime += ElapsedRealTime*ElapsedRealTime;
  NumberOfEvents++;

  // vertex, track statistics
  G4int Vertices = evt->GetNumberOfPrimaryVertex();
  MeanVertices+=Vertices;
  RmsVertices+=Vertices*Vertices;
  G4int Tracks=0;
  for (G4int iv=0; iv<Vertices; iv++){
    Tracks+=evt->GetPrimaryVertex(iv)->GetNumberOfParticle();
  }
  MeanTracks+=Tracks;
  RmsTracks+=Tracks*Tracks;

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
    { n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories
	 << " trajectories stored in this event." << G4endl;


  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    for(G4int i=0; i<n_trajectories; i++) {
      G4VTrajectory* trj = (*(evt->GetTrajectoryContainer()))[i];
      //           if (drawFlag == "all") trj->DrawTrajectory();
      //           else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
      trj->DrawTrajectory(); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
