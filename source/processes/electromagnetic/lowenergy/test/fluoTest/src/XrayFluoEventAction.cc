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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoEventAction.hh"
#include "XrayFluoSensorHit.hh"
#include "XrayFluoEventActionMessenger.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoDataSet.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif

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
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoEventAction::XrayFluoEventAction()
  :drawFlag("all"),
   HPGeCollID(-1),
   eventMessenger(0),
   printModulo(1)
 
{
  eventMessenger = new XrayFluoEventActionMessenger(this);
  runManager = new XrayFluoRunAction();
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoEventAction::~XrayFluoEventAction()
{
   delete eventMessenger;
   eventMessenger = 0;
   delete  runManager;
   runManager = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  if (HPGeCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      HPGeCollID = SDman->GetCollectionID("HPGeCollection");
      //the pointer points to the ID number of the sensitive detector
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoEventAction::EndOfEventAction(const G4Event* evt)
{
  // extracted from hits, compute the total energy deposit (and total charged
  // track length) 
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  XrayFluoSensorHitsCollection* HPGeHC = 0;
  G4int n_hit = 0;
  G4double totEnergyDetect=0., totEnergy=0, energyD=0.;
  
  if (HCE) HPGeHC = (XrayFluoSensorHitsCollection*)(HCE->GetHC(HPGeCollID));
  if(HPGeHC)
    {
      n_hit = HPGeHC->entries();
      for (G4int i=0;i<n_hit;i++)
	{
	  totEnergy += (*HPGeHC)[i]->GetEdepTot(); 
	  
	  energyD = ResponseFunction(totEnergy);
#ifdef G4ANALYSIS_USE
	    XrayFluoAnalysisManager* analysis = XrayFluoAnalysisManager::getInstance();
	    analysis->analyseEnergyDep(energyD);
#endif
	    totEnergyDetect += energyD;
	    
	    
	}
    }
  
  // extract the trajectories and draw them
  
  if (G4VVisManager::GetConcreteInstance())
    {
      
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->size();
      
      for (G4int i=0; i<n_trajectories; i++) 
	{ G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	if (drawFlag == "all") trj->DrawTrajectory();
	else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	  trj->DrawTrajectory();
	else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
	  trj->DrawTrajectory();				   
	}
    }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



G4double XrayFluoEventAction::RandomCut(G4double energy)
  
{
  G4double efficiency = 1.;
  G4double F = 0.15;
  G4double epsilon = 2.96 * eV;
  G4double deltaE = 220 * eV;
  G4double EdepDetect = 0.;
  //const XrayFluoDataSet* dataSet = runManager->GetEfficiencySet();
     
  //G4double id = 0;
 
  //efficiency = dataSet->FindValue(energy,id); 
 
  G4double  Random = G4UniformRand(); 

    if ( Random<efficiency )
      {
	G4double sigma = std::sqrt(F*epsilon*energy+std::pow(deltaE/2355,2));

	EdepDetect = G4RandGauss::shoot(energy, sigma );

  }
    else {EdepDetect = 0.;}
    return   EdepDetect;
    
};
G4double XrayFluoEventAction::ResponseFunction(G4double energy)
{
  G4double eMin = 1* keV;
  G4double eMax = 10*keV; 
  G4double value = 0.;
  G4double efficiency = 1.;
  
  const XrayFluoDataSet* dataSet = runManager->GetEfficiencySet();
  G4double id = 0;
  
  G4double random = G4UniformRand();
 
  if (energy>=eMin && energy <=eMax)
    {
      G4double infEnergy = (G4int)(energy/keV)* keV;
      
      G4double supEnergy = ((G4int)(energy/keV) + 1)*keV;
      
      
      
      G4double infData = runManager->GetInfData(energy, random);
      
      G4double supData = runManager->GetSupData(energy,random);
      
      value = (std::log10(infData)*std::log10(supEnergy/energy) +
	       std::log10(supData)*std::log10(energy/infEnergy)) / 
	std::log10(supEnergy/infEnergy);
      value = std::pow(10,value);
    }
  else if (energy<eMin)
    { 
      G4double infEnergy = eMin;
      G4double supEnergy = eMin/keV +1*keV;
 
      G4double infData = runManager->GetInfData(eMin, random);
      G4double supData = runManager->GetSupData(eMin,random);
      value = (std::log10(infData)*std::log10(supEnergy/eMin) +
	       std::log10(supData)*std::log10(eMin/infEnergy)) / 
	std::log10(supEnergy/infEnergy);
      value = std::pow(10,value);
      value = value-eMin+ energy;


    }
 else if (energy>eMax)
    { 
      G4double infEnergy = eMax/keV - 1. *keV;
      G4double supEnergy = eMax;
 
      G4double infData = runManager->GetInfData(eMax, random);
      G4double supData = runManager->GetSupData(eMax,random);
      value = (std::log10(infData)*std::log10(supEnergy/eMax) +
	       std::log10(supData)*std::log10(eMax/infEnergy)) / 
	std::log10(supEnergy/infEnergy);
      value = std::pow(10,value);
      value = value+energy- eMax;
    }
  G4double  RandomNum = G4UniformRand(); 
  
  efficiency = dataSet->FindValue(value,id);
  if ( RandomNum>efficiency )
    {
      value = 0.;
    }
 
  return value;

}
