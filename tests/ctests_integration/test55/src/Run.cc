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

#include "Run.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmCalculator.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "TestSeries.hh"

#include "Randomize.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, HistoManager* histoMgr)
  :fDetector(det), fHistoManager(histoMgr)
{
  fAnalysisManager = G4AnalysisManager::Instance();

  //initialisation
  fNEvt = 0;

  fTallyEdep.resize(kMaxTally, 0.);
  fBinLength = 1.*mm; 
  fOffsetX = 0.;

  //initialize projected range, tallies, Ebeam, and book histograms
  //
  nPrimarySteps = 0;
  nRange = 0;
  fProjRange = fProjRange2 = 0.;
  fEdepTot = fEniel = 0.;
  fEIncomingTot = fEOutgoingTot = 0.;
  fNmbein = fNmbeout = 0;
  fEbeamCumul = 0.;

  // define "1" histogram binning
  fLength  = fDetector->GetAbsorSizeX();
  fOffsetX = 0.5*fLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fNEvt += localRun->GetNumberOfEvent();

  nPrimarySteps  += localRun->nPrimarySteps;
  nRange += localRun-> nRange;
  fProjRange += localRun->fProjRange;
  fProjRange2 += localRun-> fProjRange2;
  fEdepTot += localRun->fEdepTot;
  fEniel += localRun->fEniel;
  fEIncomingTot += localRun->fEIncomingTot;
  fEOutgoingTot += localRun->fEOutgoingTot;
  fNmbein += localRun->fNmbein;
  fNmbeout += localRun->fNmbeout;

  for(size_t i=0; i<fTallyEdep.size(); ++i) {
    fTallyEdep[i] += localRun->fTallyEdep[i];
  }

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun(TestSeries* testRange, TestSeries* testEnergyLoss )
{
  G4int NbofEvents = fNEvt;
  if (NbofEvents == 0) return;

  //run conditions
  //  
  const G4Material* material = fDetector->GetAbsorMaterial();
  G4double density = material->GetDensity();
   
  G4String particle = fDetector->GetBeamParticle()->GetParticleName();    
  G4double energy = fDetector->GetBeamEnergy();
  G4cout << "\n The run consists of " << NbofEvents << " "<< particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
	 << G4BestUnit(fDetector->GetAbsorSizeX(),"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
	 
  //compute projected range and straggling
  //
  if(nRange > 0) {
    fProjRange /= nRange; 
    fProjRange2 /= nRange;
  }
  G4double rms = fProjRange2 - fProjRange*fProjRange;        
  rms = (rms > 0.0) ? std::sqrt(rms) : 0.;

  G4double nstep = G4double(nPrimarySteps)/G4double(NbofEvents);

  G4cout.precision(6);       
  G4cout << "\n Projected Range= "<< G4BestUnit(fProjRange,"Length")
         << "   rms= "            << G4BestUnit( rms,"Length")
         << G4endl;
  G4cout << " Mean number of primary steps = "<< nstep << G4endl;

  //compute energy deposition and NIEL
  //
  fEdepTot /= NbofEvents; 
  G4cout << " Total energy deposit= "<< G4BestUnit(fEdepTot,"Energy")
         << G4endl;
  fEniel /= NbofEvents; 
  G4cout << " NIEL energy deposit = "<< G4BestUnit(fEniel,"Energy")
         << G4endl;
     
  //print dose in tallies
  //
  G4int tallyNumber = fDetector->GetTallyNumber();
  if (tallyNumber > 0) {
    G4double Ebeam = fEbeamCumul;
    G4cout << "\n---------------------------------------------------------\n";
    G4cout << " Cumulated Doses : \tEdep      \tEdep/Ebeam \tDose" << G4endl;
    for (G4int j=0; j<tallyNumber; ++j) {
      G4double Edep = fTallyEdep[j], ratio = 100*Edep/Ebeam;
      G4double tallyMass = fDetector->GetTallyMass(j);
      G4double Dose = Edep/tallyMass;
      G4cout << " tally " << j << ": \t \t"
             << G4BestUnit(Edep,"Energy") << "\t"
	     << ratio << " % \t"
	     << G4BestUnit(Dose,"Dose")   << G4endl;
    }
    G4cout << "\n---------------------------------------------------------\n";
    G4cout << G4endl; 
  }

  // Pass computed values to tests
  testRange->EndOfRunAction(fProjRange);
  if(fNmbein > 0 && fNmbeout >0) {
    testEnergyLoss->EndOfRunAction(fEIncomingTot/fNmbein
				   -fEOutgoingTot/fNmbeout);
  }

  // save histograms
  if(fAnalysisManager) {

    if(fAnalysisManager->IsActive()) {      
      // normalize histograms
      //
      G4double fac = (mm/MeV)/(NbofEvents *  fBinLength);
      for (G4int j=0; j<3; ++j) { fHistoManager->Scale(j, fac); }
      
      // Write histogram file
      if(!fAnalysisManager->Write()) {
        G4Exception ("Histo::Save()", "hist01", FatalException, 
                     "Cannot write ROOT file.");
      }
      
      G4cout << "### Histo::Save: Histograms are saved" << G4endl;
      if(fAnalysisManager->CloseFile()) {
        G4cout << "                 File is closed" << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillHisto(G4int histoId, G4double v1, G4double v2)
{
  if(fAnalysisManager) {
    fHistoManager->FillHisto(histoId, v1, v2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

