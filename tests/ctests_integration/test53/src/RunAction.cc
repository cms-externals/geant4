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

#include "RunAction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "HistoManager.hh"
#include "EmAcceptance.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Run.hh"

#include "G4UnitsTable.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim, 
		     HistoManager*)
 :fDetector(det), fPrimary(prim)
{
  // Book predefined histograms
  fHistoManager = new HistoManager(fDetector->GetNbOfAbsor(),fDetector->GetMaxAbsor());

  // RunMessenger
  fRunMessenger = new RunActionMessenger(this);

  fMaxAbsor = fDetector->GetMaxAbsor();
  fEdepTrue.resize(fMaxAbsor);
  fRmsTrue.resize(fMaxAbsor);
  fLimitTrue.resize(fMaxAbsor);

  for (G4int k=0; k<fMaxAbsor; k++)
  { 
    fEdepTrue[k] = fRmsTrue[k] = 1.;
    fLimitTrue[k] = DBL_MAX;
  }

  fApplyLimit = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHistoManager;
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector,fPrimary,fHistoManager,fApplyLimit);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //histograms
  //        
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() )
  {
    analysisManager->OpenFile();
  } 

  //example of print dEdx tables
  //  PrintDedxTables();

  // save Rndm status
  //
  // if (!isMaster) { G4Random::showEngineStatus(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{  
  // print Run summary
  //
  if (isMaster) {
    fRun->EndOfRun(fEdepTrue, fRmsTrue, fLimitTrue);    
      
    // show Rndm status
    G4cout << "RunAction::EndOfRunAction engine status" << G4endl;
    G4Random::showEngineStatus(); 
  }
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {    
    analysisManager->Write();
    analysisManager->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetEdepAndRMS(G4int i, G4double edep, G4double rms, G4double lim)
{
  if (i>=0 && i<fMaxAbsor)
  {
    fEdepTrue [i] = edep;
    fRmsTrue  [i] = rms;
    fLimitTrue[i] = lim;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintDedxTables()
{
  //Print dE/dx tables with binning identical to the Geant3 JMATE bank.
  //The printout is readable as Geant3 ffread data cards (by the program g4mat).
  //
  const G4double tkmin=10*keV, tkmax=10*TeV;
  const G4int nbin=90;
  G4double tk[nbin];

  const G4int ncolumn = 5;

  //compute the kinetic energies
  //
  const G4double dp = std::log10(tkmax/tkmin)/nbin;
  const G4double dt = std::pow(10.,dp);
  tk[0] = tkmin;
  for (G4int i=1; i<nbin; ++i) tk[i] = tk[i-1]*dt;

  //print the kinetic energies
  //
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(3);

  G4cout << "\n kinetic energies \n ";
  for (G4int j=0; j<nbin; ++j)
  {
    G4cout << G4BestUnit(tk[j],"Energy") << "\t";
    if ((j+1)%ncolumn == 0) G4cout << "\n ";
  }
  G4cout << G4endl;

  //print the dE/dx tables
  //
  G4cout.setf(std::ios::scientific,std::ios::floatfield);

  G4ParticleDefinition*
  part = fPrimary->GetParticleGun()->GetParticleDefinition();
  
  G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;

  for (G4int iab=1;iab <= fDetector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = fDetector->GetAbsorMaterial(iab);
      G4int index = 0;
      for (size_t i=0; i<numOfCouples; i++)
      {
         couple = theCoupleTable->GetMaterialCutsCouple(i);
         if (couple->GetMaterial() == mat) {index = i; break;}
      }
      G4cout << "\nLIST";
      G4cout << "\nC \nC  dE/dx (MeV/cm) for " << part->GetParticleName()
             << " in " << mat ->GetName() << "\nC";
      G4cout << "\nKINE   (" << part->GetParticleName() << ")";
      G4cout << "\nMATE   (" << mat ->GetName() << ")";
      G4cout.precision(2);
      G4cout << "\nERAN  " << tkmin/GeV << " (ekmin)\t"
                           << tkmax/GeV << " (ekmax)\t"
			   << nbin      << " (nekbin)";
      G4double cutgam =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index];
      if (cutgam < tkmin) cutgam = tkmin;
      if (cutgam > tkmax) cutgam = tkmax;
      G4double cutele =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index];
      if (cutele < tkmin) cutele = tkmin;
      if (cutele > tkmax) cutele = tkmax;
      G4cout << "\nCUTS  " << cutgam/GeV << " (cutgam)\t"
                           << cutele/GeV << " (cutele)";

      G4cout.precision(6);
      G4cout << "\nG4VAL \n ";
      for (G4int l=0;l<nbin; ++l)
         {
           G4double dedx = G4LossTableManager::Instance()
	                                       ->GetDEDX(part,tk[l],couple);
           G4cout << dedx/(MeV/cm) << "\t";
	   if ((l+1)%ncolumn == 0) G4cout << "\n ";
         }
      G4cout << G4endl;
     }

  G4cout.precision(prec);
  G4cout.setf(mode,std::ios::floatfield);
}

