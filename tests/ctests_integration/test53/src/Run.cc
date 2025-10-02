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

#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "HistoManager.hh"
#include "EmAcceptance.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* prim, 
	 HistoManager* histoMgr, G4bool applyLimit)
  :fDetector(det), fPrimary(prim), fHistoManager(histoMgr), fApplyLimit(applyLimit)
{

  fMaxAbsor = fDetector->GetMaxAbsor();
  fNEvt = 0;

  //initialize cumulative quantities
  //
  fSumEAbs.resize(fMaxAbsor);
  fSum2EAbs.resize(fMaxAbsor);
  fSumLAbs.resize(fMaxAbsor); 
  fSum2LAbs .resize(fMaxAbsor);
  fEnergyDeposit.resize(fMaxAbsor);

  for (G4int k=0; k<fMaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k]  = fSumLAbs[k] = fSum2LAbs[k] = 0.;
    fEnergyDeposit[k].clear();  
  }
  
  //initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  fEnergyFlow.resize(nbPlanes);
  fLateralEleak.resize(nbPlanes);
  for (G4int k=0; k<nbPlanes; k++) { fEnergyFlow[k] = fLateralEleak[k] = 0.; }

  fAnalysisManager = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fNEvt += localRun->GetNumberOfEvent();

  //initialize cumulative quantities
  //
  for (G4int k=0; k<fMaxAbsor; k++) {
    fSumEAbs[k] += localRun->fSumEAbs[k];
    fSum2EAbs[k]  += localRun->fSum2EAbs[k];
    fSumLAbs[k] += localRun->fSumLAbs[k];
    fSum2LAbs[k] += localRun->fSum2LAbs[k];
    
    for(G4int j=0; j<(int)fEnergyDeposit[k].size(); j++)
      fEnergyDeposit[k][j] += localRun->fEnergyDeposit[k][j];

  }
  
  //initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor()) + 2;
  for (G4int k=0; k<nbPlanes; k++) 
    {
      fEnergyFlow[k] += localRun->fEnergyFlow[k];
      fLateralEleak[k] += localRun->fLateralEleak[k];
    }
  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun(std::vector<G4double>edepTrue, std::vector<G4double>rmsTrue, 
                       std::vector<G4double>limitTrue)
{

  // set the number of event only in case of a non MT Geant4 release 
  //       ( set in Run::Merge otherwise )
#ifndef G4MULTITHREADED
  fNEvt += this->GetNumberOfEvent();
#endif

  if (fNEvt == 0) return;

  G4cout << "Histo: End of run actions are started" << G4endl;

  G4double  norm = G4double(fNEvt);
  if(norm > 0) norm = 1./norm;
  G4double qnorm = std::sqrt(norm);

  //compute and print statistic
  //
  G4double beamEnergy = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4double sqbeam = std::sqrt(beamEnergy/GeV);

  G4double MeanEAbs,MeanEAbs2,rmsEAbs,resolution,rmsres;
  G4double MeanLAbs,MeanLAbs2,rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw(14) << "material"
         << std::setw(17) << "Edep       RMS"
	 << std::setw(33) << "sqrt(E0(GeV))*rmsE/Emean"
	 << std::setw(23) << "total tracklen \n \n";

  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++)
    {
      MeanEAbs  = fSumEAbs[k]*norm;
      MeanEAbs2 = fSum2EAbs[k]*norm;
      rmsEAbs  = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      //G4cout << "k= " << k << "  RMS= " <<  rmsEAbs 
      //     << "  applyLimit: " << applyLimit << G4endl;
      if(fApplyLimit) {
        G4int    nn    = 0;
        G4double sume  = 0.0;
        G4double sume2 = 0.0;
	// compute trancated means  
        G4double lim   = rmsEAbs * 2.5;
        for(G4int i=0; i<fNEvt; i++) {
          G4double e = (fEnergyDeposit[k])[i];
          if(std::abs(e - MeanEAbs) < lim) {
            sume  += e;
            sume2 += e*e;
            nn++;
	  }
	}
        G4double norm1 = G4double(nn);
        if(norm1 > 0.0) norm1 = 1.0/norm1;
	MeanEAbs  = sume*norm1;
	MeanEAbs2 = sume2*norm1;
	rmsEAbs  = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      }

      resolution= 100.*sqbeam*rmsEAbs/MeanEAbs;
      rmsres    = resolution*qnorm;

      // Save mean and RMS
      fSumEAbs[k] = MeanEAbs;
      fSum2EAbs[k] = rmsEAbs;

      MeanLAbs  = fSumLAbs[k]*norm;
      MeanLAbs2 = fSum2LAbs[k]*norm;
      rmsLAbs  = std::sqrt(std::abs(MeanLAbs2 - MeanLAbs*MeanLAbs));

      //print
      //
      G4cout
       << std::setw(14) << fDetector->GetAbsorMaterial(k)->GetName() << ": "
       << std::setprecision(5)
       << std::setw(6) << G4BestUnit(MeanEAbs,"Energy") << " :  "
       << std::setprecision(4)
       << std::setw(5) << G4BestUnit( rmsEAbs,"Energy")  
       << std::setw(10) << resolution  << " +- " 
       << std::setw(5) << rmsres << " %"
       << std::setprecision(3)
       << std::setw(10) << G4BestUnit(MeanLAbs,"Length")  << " +- "
       << std::setw(4) << G4BestUnit( rmsLAbs,"Length")
       << G4endl;
    }
  G4cout << "\n------------------------------------------------------------\n";

  G4cout << " Beam particle " 
	 << fPrimary->GetParticleGun()->GetParticleDefinition()->GetParticleName()
	 << "  E = " << G4BestUnit(beamEnergy,"Energy") << G4endl;
  
   //Energy flow
   //
  G4int Idmax = (fDetector->GetNbOfLayers())*(fDetector->GetNbOfAbsor());
  for (G4int Id=1; Id<=Idmax+1; Id++) {
    fHistoManager->FillHisto(2*fMaxAbsor+1, (G4double)Id, fEnergyFlow[Id]);
    fHistoManager->FillHisto(2*fMaxAbsor+2, (G4double)Id, fLateralEleak[Id]);
  }
  
  //Energy deposit from energy flow balance
  //
  std::vector<G4double> EdepTot;
  EdepTot.assign(fMaxAbsor,0.);
  
  G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  for (G4int Id=1; Id<=Idmax; Id++) {
    G4int iAbsor = Id%nbOfAbsor; if (iAbsor==0) iAbsor = nbOfAbsor;
    EdepTot [iAbsor] += (fEnergyFlow[Id] - fEnergyFlow[Id+1] - fLateralEleak[Id]);
  }
  
  G4cout << "\n Energy deposition from Energy flow balance : \n"
         << std::setw(10) << "  material \t Total Edep \n \n";
  G4cout.precision(6);
  
  for (G4int k=1; k<=nbOfAbsor; k++) {
    EdepTot [k] *= norm;
    G4cout << std::setw(10) << fDetector->GetAbsorMaterial(k)->GetName() << ":"
           << "\t " << G4BestUnit(EdepTot [k],"Energy") << "\n";
  }
  
  G4cout << "\n------------------------------------------------------------\n" 
         << G4endl;
    
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // Acceptance
  EmAcceptance acc;
  G4bool isStarted = false;
  for (G4int j=1; j<=fDetector->GetNbOfAbsor(); j++) {
    if (limitTrue[j] < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Sampling Calorimeter",fNEvt);
	isStarted = true;
      }
      MeanEAbs = fSumEAbs[j];
      rmsEAbs  = fSum2EAbs[j];
      G4String mat = fDetector->GetAbsorMaterial(j)->GetName();
      acc.EmAcceptanceGauss("Edep"+mat, fNEvt, MeanEAbs,
                             edepTrue[j], rmsTrue[j], limitTrue[j]);
      acc.EmAcceptanceGauss("Erms"+mat, fNEvt, rmsEAbs,
                             rmsTrue[j], rmsTrue[j], 2.0*limitTrue[j]);
    }
  }
  if(isStarted) acc.EndOfAcceptance();


  if(fAnalysisManager) {

    if(fAnalysisManager->IsActive()) {

      //normalize histograms
      G4int maxHisto = fHistoManager->GetMaxHisto();
      for (G4int ih = fMaxAbsor+1; ih < maxHisto; ih++) {
	fHistoManager->Normalize(ih,norm);
      }
      
      // // Write histogram file
      // if(!fAnalysisManager->Write()) {
      //   G4Exception ("Histo::Save()", "hist01", FatalException, 
      //                "Cannot write ROOT file.");
      // }

      // G4cout << "### Histo::Save: Histograms are saved" << G4endl;
      // if(fAnalysisManager->CloseFile()) {
      //   G4cout << "                 File is closed" << G4endl;
      // }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::fillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{
  //accumulate statistic with restriction
  //
  if(fApplyLimit) fEnergyDeposit[kAbs].push_back(EAbs);
  fSumEAbs[kAbs]  += EAbs;  fSum2EAbs[kAbs]  += EAbs*EAbs;
  fSumLAbs[kAbs]  += LAbs;  fSum2LAbs[kAbs]  += LAbs*LAbs;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::fillEnergyPerAbsorber(G4int histoId, G4double energy)
{
  fHistoManager->FillHisto(histoId, energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::fillEdepProfilePerAbsorber(G4int histoId, G4double layer, G4double energy)
{
  fHistoManager->FillHisto(histoId, layer, energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

