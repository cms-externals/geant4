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

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
public:
  
  HistoManager(G4int, G4int);
  ~HistoManager();

  void SetHisto (G4int,G4int,G4double,G4double,const G4String& unit="none");  
  void FillHisto(G4int id, G4double bin, G4double weight = 1.0);
  void Normalize(G4int ih, G4double norm);
  
  G4int GetMaxHisto() const { return fMaxHisto; }
    
private:
  
  void Book();
  
  G4int fNbAbsor;
  G4int fMaxAbsor;
  G4int fMaxHisto;
  std::vector<G4int> fHistoId;  

  std::vector<G4String> fLabel;
  std::vector<G4String> fTitle;
  std::vector<G4int>    fNbins;
  std::vector<G4double> fVmin ;
  std::vector<G4double> fVmax ;        
  std::vector<G4double> fUnit ;
  std::vector<G4bool>   fExist;
  
private:
  void saveAscii();         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

