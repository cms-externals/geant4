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

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"

#include "G4ThreeVector.hh"
#include "G4AnalysisManager.hh"
#include "globals.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
public:

  Run(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*, G4bool);
  ~Run();

  virtual void Merge(const G4Run*);
  void EndOfRun(std::vector<G4double>edepTrue, std::vector<G4double>rmsTrue, std::vector<G4double>limitTrue);

  void fillPerEvent(G4int,G4double,G4double);
    
  void sumEnergyFlow(G4int plane, G4double Eflow)
                                            { fEnergyFlow[plane]  += Eflow; };
  void sumLateralEleak(G4int cell, G4double Eflow)
                                            { fLateralEleak[cell] += Eflow; };

  void fillEnergyPerAbsorber(G4int histoId, G4double energy);
  void fillEdepProfilePerAbsorber(G4int histoId, G4double layerId, G4double energy);
    
private:
  
  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fPrimary;    
  HistoManager*           fHistoManager;
  G4AnalysisManager* fAnalysisManager;

  G4int fMaxAbsor;
  G4int fNEvt;
  std::vector<G4double> fSumEAbs, fSum2EAbs ; 
  std::vector<G4double> fSumLAbs, fSum2LAbs ;
    
  std::vector<G4double> fEnergyFlow;
  std::vector<G4double> fLateralEleak;
  std::vector<std::vector<G4double> > fEnergyDeposit;

  G4bool fApplyLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

