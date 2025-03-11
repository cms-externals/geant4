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
#include "globals.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4AnalysisManager.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class HistoManager;
class G4ParticleDefinition;
class TestSeries;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
public:

  Run(DetectorConstruction*, HistoManager*);
  ~Run();

  virtual void Merge(const G4Run*);
  void EndOfRun(TestSeries* t1, TestSeries* t2 );

  void FillTallyEdep(G4int n, G4double e)  { fTallyEdep[n] += e; };
  void FillEdep(G4double de, G4double eni) { fEdepTot += de; fEniel += eni; };
  void FillEnIncoming(G4double in) { fEIncomingTot += in; ++fNmbein; }
  void FillEnOutgoing(G4double out) { fEOutgoingTot += out; ++fNmbeout; }
       
  G4double GetBinLength() { return fBinLength; };
  G4double GetLength()    { return fLength; };
  G4double GetOffsetX()   { return fOffsetX; }
     
  void AddProjRange (G4double x) 
  { fProjRange += x; fProjRange2 += x*x; ++nRange; };
  void AddPrimaryStep() { ++nPrimarySteps; };

  void FillHisto(G4int histoId, G4double v1, G4double v2 =1.);
    
private:
  
  DetectorConstruction*   fDetector;
  HistoManager*           fHistoManager;
  G4AnalysisManager*      fAnalysisManager;

  G4int fNEvt;

  std::vector<G4double>   fTallyEdep;   
  G4double                fBinLength;
  G4double                fOffsetX;
  G4double                fLength;

  G4double                fProjRange, fProjRange2;
  G4double                fEdepTot, fEniel;
  G4double                fEIncomingTot, fEOutgoingTot;
  G4int                   fNmbein, fNmbeout;
  G4int                   nPrimarySteps;
  G4int                   nRange;

  G4double                fEbeamCumul;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

