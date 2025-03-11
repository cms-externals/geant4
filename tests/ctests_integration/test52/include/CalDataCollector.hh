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
#ifndef CALDATACOLLECTOR_HH
#define CALDATACOLLECTOR_HH

#include "globals.hh"
#include "DataManager.hh"


class CalDataCollector : public DataManager {

 public:
   CalDataCollector(G4double zLow, G4double zUp, G4double radi);
   ~CalDataCollector();

   void ScoreEnergyDeposit(G4double en, 
                           G4double x, 
                           G4double y, 
                           G4double z,
                           const G4String& type="");
   void ScoreParticleEnergy(G4double, 
                            G4double, 
                            G4double, 
                            G4double,
                            const G4String&) {}
   void PrintResults();

 private:

   G4double totEnergyDeposit;
   G4double totEnergyDepositSumSquares;
   G4double elecEnergyDeposit;
   G4double photEnergyDeposit;
   G4int    nmbEnergyDeposits;
};

#endif // CALDATACOLLECTOR_HH
