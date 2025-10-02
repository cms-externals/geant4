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

#include <sstream>
#include <fstream>
#include <cmath>
#include "CalDataCollector.hh"

CalDataCollector::CalDataCollector(G4double zLow,G4double zUp, G4double) :
    DataManager(zLow,zUp),
    totEnergyDeposit(0),
    totEnergyDepositSumSquares(0),
    elecEnergyDeposit(0),
    photEnergyDeposit(0),
    nmbEnergyDeposits(0) {
}

CalDataCollector::~CalDataCollector() {
}

void CalDataCollector::ScoreEnergyDeposit(G4double en, 
                                          G4double x, 
                                          G4double y, 
					  G4double z,
                                          const G4String& ptype) {
  /*
  if(latEnergyDeposit == 0) {
     G4double center = 
              GetLowerBound() + (GetUpperBound() - GetLowerBound()) * 0.5;
     
     std::stringstream s;
     s << center;
     G4String histName = "z=" + s.str();
     G4String histTitle = "Lateral energy distr. in calorimeter at " + histName;

     latEnergyDeposit =
         histogramFactory -> createHistogram2D(histName,
                                               histTitle,
                                               50,-radius,radius,
                                               50,-radius,radius);
  }
*/

  if(z >= GetLowerBound() && z < GetUpperBound()) {

     if(ptype == "electron" || ptype == "e-")  elecEnergyDeposit += en;
     if(ptype == "gamma")     photEnergyDeposit += en;

     totEnergyDeposit += en;
     totEnergyDepositSumSquares += en * en;
     nmbEnergyDeposits++;
     //latEnergyDeposit -> fill(x, y, en);
     
     DataManager::ScoreEnergyDeposit(en, x, y, z, ptype);
  }
}


void CalDataCollector::PrintResults() {

  G4cout << "--------------------------------------------------------" 
	 << G4endl;
  G4cout << "CALORIMETER: " 
	 << GetLowerBound() << " to "
	 << GetUpperBound() << G4endl;
  G4cout << "--------------------------------------------------------" 
	 << G4endl;
 
  G4double totEnergyDepositErrSqu = 0.0;
  if(nmbEnergyDeposits > 0) {
    totEnergyDepositErrSqu =         
           totEnergyDepositSumSquares / (totEnergyDeposit * 
           totEnergyDeposit) - 1.0 / G4double(nmbEnergyDeposits); 
  } 
  G4double totEnergyDepositRelErr = std::sqrt(totEnergyDepositErrSqu);  

  G4double thickness = GetUpperBound() - GetLowerBound();
  G4double center = GetLowerBound() + 0.5 * thickness;

  G4cout << "Center and Thickness:  "
       << center << " " 
       << thickness 
       << "   Total Energy Deposit: "
       << totEnergyDeposit 
       << "   Rel. Error: "
       << totEnergyDepositRelErr
       << G4endl
       << "  Electron Energy Deposit: " 
       << elecEnergyDeposit
       << G4endl            
       << "  Photon Energy Deposit:   " 
       << photEnergyDeposit
       << G4endl;            

  DataManager::PrintResults();
}
