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

#ifndef HistoManager_h
#  define HistoManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Utility class to hold and manipulate histograms/nTuples
//
// Author:      V.Ivanchenko 30/10/03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#  include "G4AnalysisManager.hh"
#  include "G4DataVector.hh"
#  include "G4DynamicParticle.hh"
#  include "G4Track.hh"
#  include "G4VPhysicalVolume.hh"
#  include "globals.hh"

#  include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction;

class HistoManager
{
  public:

    HistoManager(DetectorConstruction*);
    ~HistoManager();

    void Book();
    // Book predefined histogramms

    // In this method histogramms are predefined
    void Add1D(G4int, const G4String&, G4int nb, G4double x1, G4double x2, G4String u1,
               G4String u2 = "none");

    // It change bins and boundaries
    void SetHisto1D(G4int, G4int, G4double, G4double, G4double);

    // Histogramms are filled
    void FillHisto(G4int, G4double, G4double w = 1.);

    // Scale histogramm
    void Scale(G4int, G4double);

  private:

    DetectorConstruction* fDetector;

    G4int fNbHisto;
    G4bool fDefaultAct;
    G4bool fVerbose;
    G4String fFileName;
    std::vector<G4int> fHistoId;
    std::vector<G4bool> fActive;
    std::vector<G4int> fNbins;
    std::vector<G4double> fXmin;
    std::vector<G4double> fXmax;
    std::vector<G4double> fUnit1;
    std::vector<G4double> fUnit2;
    std::vector<G4String> fIds;
    std::vector<G4String> fTitle;
};

#endif
