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

#include "DataManager.hh"
#include <iostream>
#include <fstream>

DataManager::DataManager(G4double zLow, G4double zUp) : 
    zAxisLowerBound(zLow),
    zAxisUpperBound(zUp),
    eps(0.0000000001) {
}

DataManager::~DataManager() {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();
 
  for(;iter != iter_end; iter++) {
      delete iter -> second;
  }

  iter = garbage.begin();
  iter_end = garbage.end();
 
  for(;iter != iter_end; iter++) {
      delete iter -> second;
  }
}

void DataManager::ScoreEnergyDeposit(G4double en, 
                                     G4double x, 
                                     G4double y, 
                                     G4double z,
                                     const G4String& type) {

  if(z < zAxisLowerBound || z > zAxisUpperBound) {
     return;
  }

  if(dataCollectors.begin() == dataCollectors.end()) return;

  collector::iterator iter = std::upper_bound(dataCollectors.begin(),
                                              dataCollectors.end(),
                                              z,zCompare());

  if(iter != dataCollectors.begin()) 
       (iter-1) -> second -> ScoreEnergyDeposit(en,x,y,z,type);
  
} 


void DataManager::ScoreParticleEnergy(G4double en, 
                                      G4double x, 
                                      G4double y, 
                                      G4double z,
                                      const G4String& type) {

  if(z < zAxisLowerBound || z >= zAxisUpperBound) return;

  collector::iterator iter = std::lower_bound(dataCollectors.begin(),
                                              dataCollectors.end(),
                                              z,zCompare());
  
  if(iter != dataCollectors.end()) {
     (iter-1) -> second -> ScoreParticleEnergy(en,x,y,z,type);
  }
} 


void DataManager::PrintResults() {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {
     iter -> second -> PrintResults();
  }

}


void DataManager::AddDataCollector(DataManager* comp) {
 
  G4double zLow = comp -> GetLowerBound();
  G4double zUp  = comp -> GetUpperBound();

  if(!IsContained(zLow,zUp)) {
     G4cout << "Error. Slab with lower bound z="     << zLow 
	    << " exceeds boundary of mother volume." << G4endl;
     garbage.push_back(std::make_pair(zLow,comp));
     return;
  }

  if(!OverlapsWithOtherChild(zLow,zUp)) {
     dataCollectors.push_back(std::make_pair(zLow,comp));

     G4double thickn = zUp - zLow;
     G4double zPos = zLow + 0.5 * thickn;
     G4cout << "INFORMATION. Slab with center at z="  << zPos 
	    << " and thickness " << thickn <<" added." << G4endl; 
  }
  else {
     G4cout << "Error. Slab with lower bound z="  << zLow 
	    << " overlaps with other slab." << G4endl;
     garbage.push_back(std::make_pair(zLow,comp));
     return;
  }

  std::sort(dataCollectors.begin(),dataCollectors.end(),zCompare());
}


bool DataManager::HasDataCollector(G4double z) {

  if(std::binary_search(dataCollectors.begin(),
                        dataCollectors.end(),z,zCompare())) {
     return true;
  }
 
  return false;
}


bool DataManager::OverlapsWithOtherChild(G4double zLow, G4double zUp) {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {

     if(zLow < iter -> second -> GetLowerBound() && 
        zUp > (iter -> second -> GetLowerBound() + eps)) return true;

     if(zLow < (iter -> second -> GetUpperBound() - eps) && 
        zUp > iter -> second -> GetUpperBound()) return true;

     if(iter -> second -> IsContained(zLow,zUp)) return true;
  }

  return false;
}


bool DataManager::IsContained(G4double zLow, G4double zUp){

  if(zLow >= zAxisLowerBound && zUp <= zAxisUpperBound) return true;  

  return false;
}


DataManager* DataManager::MatchingChildDataCollector(G4double zLow, G4double zUp) {

  collector::iterator iter = dataCollectors.begin();
  collector::iterator iter_end = dataCollectors.end();

  for(;iter != iter_end; iter++) {
     if(iter -> second -> IsContained(zLow,zUp)) return iter -> second;
  }

  return 0;
}
