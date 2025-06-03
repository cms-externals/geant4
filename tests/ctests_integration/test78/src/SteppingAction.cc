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
/// \file electromagnetic/TestEmXX/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 84208 2014-10-10 14:44:50Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();

  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4bool transmit = (endPoint->GetStepStatus() <= fGeomBoundary);
  if (transmit) { run->CountProcesses(procName); }
  else {
    //count real processes and sum track length
    G4double stepLength = aStep->GetStepLength();
    run->CountProcesses(procName);
    run->SumTrack(stepLength);
  }

  //plot final state
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //scattered primary particle
  //
  //G4int id = 1;
  if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
    //G4double energy = endPoint->GetKineticEnergy();
    //analysisManager->FillH1(id,energy);

    //id = 2;
    G4ThreeVector direction = endPoint->GetMomentumDirection();
    G4double trPhi=direction.phi();
    if (trPhi<0) trPhi=trPhi+CLHEP::twopi;
//if azimuth is 0 or pi or 2pi, then this is the parallel plane
  if (trPhi<=0.02 || abs(trPhi-CLHEP::pi)<=0.02 || abs(trPhi-2*CLHEP::pi)<=0.02 )
	{
	  analysisManager->FillH1(1,acos(direction.z())*180./CLHEP::pi);
  	}
//if azimuth is pi/2 or 3pi/2, then this is the perpendicular plane
  if (abs(trPhi-1.*CLHEP::pi/2.)<=0.02 || abs(trPhi-1.5*CLHEP::pi)<=0.02)
  	{
	 analysisManager->FillH1(2,acos(direction.z())*180./CLHEP::pi);
 	}
/*
Filling the 2D histogram with theta, scattered photon polarization.
The azimuth can be selected depending on application where G4JAEAPolarizedElasticScattering is used.
Here, the plane parallel to polarization of incident photon is selected (phi=0).
Most accurate case when phi=0 exactly. But a window is selected to see good results at reasonable beamOn.
*/
if (trPhi<=0.005 || abs(trPhi-CLHEP::pi)<=0.005 || abs(trPhi-2*CLHEP::pi)<=0.005 )
    {
   	analysisManager->FillH2(1,acos(direction.z())*180./CLHEP::pi,endPoint->GetPolarization().x());
    }
  }



  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
