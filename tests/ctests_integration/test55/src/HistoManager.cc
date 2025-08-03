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
//---------------------------------------------------------------------------
//
// ClassName:   Histo - Generic histogram/ntuple manager class
//
//
// Author:      V.Ivanchenko 30.10.03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager(DetectorConstruction* det)
  : fDetector(det)
{
  fNbHisto     = 0;
  fDefaultAct  = true;
  fVerbose     = false;
  fFileName    = "";

  int maxHisto=5;
  fHistoId.resize(maxHisto);

  fNbins.clear();
  fXmin.clear();
  fXmax.clear();
  fUnit1.clear();
  fUnit2.clear();
  fTitle.clear();
  fActive.clear();
  fIds.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  if(fFileName==""&&analysisManager->GetFileName()!="") {
    fFileName=analysisManager->GetFileName();
  }
  if(fFileName!=""&&analysisManager->GetFileName()=="") {
    analysisManager->SetFileName(fFileName);
  }
  // Create or get analysis manager
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(false);    // enable inactivation of histograms

  // Define the histogramms 
  if(analysisManager->GetNofH1s()==0) {
    fNbHisto = 0;
    Add1D(0,"Dummy",100, 0, 100, "mm");
    Add1D(1,"Edep (MeV/mm) along absorber (mm)", 100, 0, 100, "mm", "MeV");
    Add1D(2,"Edep (MeV/mm) along absorber zoomed (mm)", 100, 0, 100, "mm", "MeV");
    Add1D(3,"Projectile range (mm)", 100, 0, 100, "mm");

    // define "1" histogram binning
    G4double length  = fDetector->GetAbsorSizeX();

    // histogram "1" is defined by the length of the target
    // zoomed histograms are defined by UI command
    SetHisto1D(1, 100, 0, length, CLHEP::mm);

    // Creating an 1-dimensional histograms in the root directory of the tree
    for(G4int i=0; i<fNbHisto; ++i) {
      fHistoId[i] = analysisManager->CreateH1(fIds[i], fTitle[i], fNbins[i], 
					      fXmin[i], fXmax[i]);
      analysisManager->SetH1Activation(fHistoId[i],fActive[i]);
    }
  }
  
  // Added to catch the SetActivation parameters set througj UI interface
  for (G4int k=0; k<fNbHisto; k++) {
    fActive[k] = analysisManager->GetH1Activation(fHistoId[k]);
  }
  // Check if a filename is set
  if(analysisManager->GetFileName()=="") return;
  
  // Activate the analysisManage ronly if a filename is set
  analysisManager->SetActivation(true);     // enable inactivation of histograms
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void HistoManager::Add1D(G4int histoId, const G4String& name, G4int nb, 
			 G4double x1, G4double x2, G4String u1, G4String u2)
{
  std::stringstream sg;
  sg<<histoId;
  const G4String id = sg.str();

  if(fVerbose) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name 
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u1<<" "<<u2 
           << G4endl;
  }
  fNbHisto++;

  double vUnit1= (u1=="none")? 1. : (G4UnitDefinition::GetValueOf(u1));
  fUnit1.push_back(vUnit1);
  double vUnit2= (u2=="none")? 1. : (G4UnitDefinition::GetValueOf(u2));
  fUnit2.push_back(vUnit2);

  x1 /= vUnit1;
  x2 /= vUnit1;
  fNbins.push_back(nb);
  fXmin.push_back(x1);
  fXmax.push_back(x2);
  fTitle.push_back(name);
  fIds.push_back(id);
  G4int active =(name=="Dummy")?false:fDefaultAct;
  fActive.push_back(active);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double unit)
{
  if(i>=0 && i<fNbHisto) {
    if(fVerbose) {
      G4cout << "Update histogram: #" << i  
             << "  " << nb << "  " << x1 << "  " << x2 << "  " 
             << G4endl;
    }
    fNbins[i] = nb;
    fUnit1[i] = unit;
    fXmin[i] = x1/fUnit1[i];
    fXmax[i] = x2/fUnit1[i];
  } else {
    G4cout << "Histo::setHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::FillHisto(G4int ih, G4double x, G4double w)
{
  if(fVerbose) {
    G4cout << "fill histogram: #" << ih << " at x= " << x 
           << "  weight= " << w
           << G4endl;   
  }
  if(!G4AnalysisManager::Instance()->IsActive()) return;
  if(!fActive[ih]) return;
  G4AnalysisManager::Instance()->FillH1( fHistoId[ih], x/fUnit1[ih], w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Scale(G4int ih, G4double fac)
{
  if(fVerbose) {
    G4cout << "Scale histogram: #" << ih << " by factor " << fac << G4endl;   
  }
  if(!G4AnalysisManager::Instance()->IsActive()) return;
  if(!fActive[ih]) return;
  G4AnalysisManager::Instance()->GetH1(fHistoId[ih])->scale(fac);
}


