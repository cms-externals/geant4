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

#include "G4EmCalculator.hh"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* prim, HistoManager* histoMgr)
  : fDetector(det), fPrimary(prim), fHistoManager(histoMgr)
{
  fAnalysisManager = G4AnalysisManager::Instance();

  // initialisation
  fNEvt = 0;
  fEnergyDeposit = fEnergyDeposit2 = 0.;
  fTrackLenCharged = fTrackLenCharged2 = 0.;
  fTrackLenNeutral = fTrackLenNeutral2 = 0.;
  fNbStepsCharged = fNbStepsCharged2 = 0.;
  fNbStepsNeutral = fNbStepsNeutral2 = 0.;
  fMscProjecTheta = fMscProjecTheta2 = 0.;
  fMscThetaCentral = 3 * ComputeMscHighland();

  fNbGamma = fNbElect = fNbPosit = 0;

  fTransmit[0] = fTransmit[1] = fReflect[0] = fReflect[1] = 0;

  fMscEntryCentral = 0;

  fEnergyLeak[0] = fEnergyLeak[1] = fEnergyLeak2[0] = fEnergyLeak2[1] = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fNEvt += localRun->GetNumberOfEvent();

  fEnergyDeposit += localRun->fEnergyDeposit;
  fEnergyDeposit2 += localRun->fEnergyDeposit2;
  fTrackLenCharged += localRun->fTrackLenCharged;
  fTrackLenCharged2 += localRun->fTrackLenCharged2;
  fTrackLenNeutral += localRun->fTrackLenNeutral;
  fTrackLenNeutral2 += localRun->fTrackLenNeutral2;
  fNbStepsCharged += localRun->fNbStepsCharged;
  fNbStepsCharged2 += localRun->fNbStepsCharged2;
  fNbStepsNeutral += localRun->fNbStepsNeutral;
  fNbStepsNeutral2 += localRun->fNbStepsNeutral2;
  fMscProjecTheta += localRun->fMscProjecTheta;
  fMscProjecTheta2 += localRun->fMscProjecTheta2;

  fMscThetaCentral = 3 * localRun->ComputeMscHighland();

  fNbGamma += localRun->fNbGamma;
  fNbElect += localRun->fNbElect;
  fNbPosit += localRun->fNbPosit;

  fTransmit[0] += localRun->fTransmit[0];
  fTransmit[1] += localRun->fTransmit[1];
  fReflect[0] += localRun->fReflect[0];
  fReflect[1] += localRun->fReflect[1];

  fMscEntryCentral += localRun->fMscEntryCentral;

  fEnergyLeak[0] += localRun->fEnergyLeak[0];
  fEnergyLeak[1] += localRun->fEnergyLeak[1];
  fEnergyLeak2[0] += localRun->fEnergyLeak2[0];
  fEnergyLeak2[1] += localRun->fEnergyLeak2[1];

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4cout << "Histo: End of run actions are started" << G4endl;

  // set the number of event only in case of a non MT Geant4 release
  //       ( set in Run::Merge otherwise )
#ifndef G4MULTITHREADED
  fNEvt += this->GetNumberOfEvent();
#endif

  // compute mean and rms
  //
  G4int TotNbofEvents = fNEvt;
  if (TotNbofEvents == 0) return;

  G4double EnergyBalance = fEnergyDeposit + fEnergyLeak[0] + fEnergyLeak[1];
  EnergyBalance /= TotNbofEvents;

  fEnergyDeposit /= TotNbofEvents;
  fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit * fEnergyDeposit;
  if (rmsEdep > 0.)
    rmsEdep = std::sqrt(rmsEdep / TotNbofEvents);
  else
    rmsEdep = 0.;

  fTrackLenCharged /= TotNbofEvents;
  fTrackLenCharged2 /= TotNbofEvents;
  G4double rmsTLCh = fTrackLenCharged2 - fTrackLenCharged * fTrackLenCharged;
  if (rmsTLCh > 0.)
    rmsTLCh = std::sqrt(rmsTLCh / TotNbofEvents);
  else
    rmsTLCh = 0.;

  fTrackLenNeutral /= TotNbofEvents;
  fTrackLenNeutral2 /= TotNbofEvents;
  G4double rmsTLNe = fTrackLenNeutral2 - fTrackLenNeutral * fTrackLenNeutral;
  if (rmsTLNe > 0.)
    rmsTLNe = std::sqrt(rmsTLNe / TotNbofEvents);
  else
    rmsTLNe = 0.;

  fNbStepsCharged /= TotNbofEvents;
  fNbStepsCharged2 /= TotNbofEvents;
  G4double rmsStCh = fNbStepsCharged2 - fNbStepsCharged * fNbStepsCharged;
  if (rmsStCh > 0.)
    rmsStCh = std::sqrt(rmsTLCh / TotNbofEvents);
  else
    rmsStCh = 0.;

  fNbStepsNeutral /= TotNbofEvents;
  fNbStepsNeutral2 /= TotNbofEvents;
  G4double rmsStNe = fNbStepsNeutral2 - fNbStepsNeutral * fNbStepsNeutral;
  if (rmsStNe > 0.)
    rmsStNe = std::sqrt(rmsTLCh / TotNbofEvents);
  else
    rmsStNe = 0.;

  G4double Gamma = (G4double)fNbGamma / TotNbofEvents;
  G4double Elect = (G4double)fNbElect / TotNbofEvents;
  G4double Posit = (G4double)fNbPosit / TotNbofEvents;

  G4double transmit[2];
  transmit[0] = 100. * fTransmit[0] / TotNbofEvents;
  transmit[1] = 100. * fTransmit[1] / TotNbofEvents;

  G4double reflect[2];
  reflect[0] = 100. * fReflect[0] / TotNbofEvents;
  reflect[1] = 100. * fReflect[1] / TotNbofEvents;

  G4double rmsMsc = 0., tailMsc = 0.;
  if (fMscEntryCentral > 0)
  {
    fMscProjecTheta /= fMscEntryCentral;
    fMscProjecTheta2 /= fMscEntryCentral;
    rmsMsc = fMscProjecTheta2 - fMscProjecTheta * fMscProjecTheta;
    if (rmsMsc > 0.) rmsMsc = std::sqrt(rmsMsc);
    tailMsc = 100. - (100. * fMscEntryCentral) / (2 * fTransmit[1]);
  }

  fEnergyLeak[0] /= TotNbofEvents;
  fEnergyLeak2[0] /= TotNbofEvents;
  G4double rmsEl0 = fEnergyLeak2[0] - fEnergyLeak[0] * fEnergyLeak[0];
  if (rmsEl0 > 0.)
    rmsEl0 = std::sqrt(rmsEl0 / TotNbofEvents);
  else
    rmsEl0 = 0.;

  fEnergyLeak[1] /= TotNbofEvents;
  fEnergyLeak2[1] /= TotNbofEvents;
  G4double rmsEl1 = fEnergyLeak2[1] - fEnergyLeak[1] * fEnergyLeak[1];
  if (rmsEl1 > 0.)
    rmsEl1 = std::sqrt(rmsEl1 / TotNbofEvents);
  else
    rmsEl1 = 0.;

  // Stopping Power from input Table.
  //
  G4Material* material = fDetector->GetAbsorberMaterial();
  G4double length = fDetector->GetAbsorberThickness();
  G4double density = material->GetDensity();

  G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();

  G4EmCalculator emCalculator;
  G4double dEdxTable = 0., dEdxFull = 0.;
  if (particle->GetPDGCharge() != 0.)
  {
    dEdxTable = emCalculator.GetDEDX(energy, particle, material);
    dEdxFull = emCalculator.ComputeTotalDEDX(energy, particle, material);
  }
  G4double stopTable = dEdxTable / density;
  G4double stopFull = dEdxFull / density;

  // Stopping Power from simulation.
  //
  G4double meandEdx = fEnergyDeposit / length;
  G4double stopPower = meandEdx / density;

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(3);

  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(energy, "Energy") << " through " << G4BestUnit(length, "Length") << " of "
         << material->GetName() << " (density: " << G4BestUnit(density, "Volumic Mass") << ")"
         << G4endl;

  G4cout.precision(4);

  G4cout << "\n Total energy deposit in absorber per event = "
         << G4BestUnit(fEnergyDeposit, "Energy") << " +- " << G4BestUnit(rmsEdep, "Energy")
         << G4endl;

  G4cout << "\n -----> Mean dE/dx = " << meandEdx / (MeV / cm) << " MeV/cm"
         << "\t(" << stopPower / (MeV * cm2 / g) << " MeV*cm2/g)" << G4endl;

  G4cout << "\n From formulas :" << G4endl;
  G4cout << "   restricted dEdx = " << dEdxTable / (MeV / cm) << " MeV/cm"
         << "\t(" << stopTable / (MeV * cm2 / g) << " MeV*cm2/g)" << G4endl;

  G4cout << "   full dEdx       = " << dEdxFull / (MeV / cm) << " MeV/cm"
         << "\t(" << stopFull / (MeV * cm2 / g) << " MeV*cm2/g)" << G4endl;

  G4cout << "\n Leakage :  primary = " << G4BestUnit(fEnergyLeak[0], "Energy") << " +- "
         << G4BestUnit(rmsEl0, "Energy")
         << "   secondaries = " << G4BestUnit(fEnergyLeak[1], "Energy") << " +- "
         << G4BestUnit(rmsEl1, "Energy") << G4endl;

  G4cout << " Energy balance :  edep + eleak = " << G4BestUnit(EnergyBalance, "Energy") << G4endl;

  G4cout << "\n Total track length (charged) in absorber per event = "
         << G4BestUnit(fTrackLenCharged, "Length") << " +- " << G4BestUnit(rmsTLCh, "Length")
         << G4endl;

  G4cout << " Total track length (neutral) in absorber per event = "
         << G4BestUnit(fTrackLenNeutral, "Length") << " +- " << G4BestUnit(rmsTLNe, "Length")
         << G4endl;

  G4cout << "\n Number of steps (charged) in absorber per event = " << fNbStepsCharged << " +- "
         << rmsStCh << G4endl;

  G4cout << " Number of steps (neutral) in absorber per event = " << fNbStepsNeutral << " +- "
         << rmsStNe << G4endl;

  G4cout << "\n Number of secondaries per event : Gammas = " << Gamma << ";   electrons = " << Elect
         << ";   positrons = " << Posit << G4endl;

  G4cout << "\n Number of events with the primary particle transmitted = " << transmit[1] << " %"
         << G4endl;

  G4cout << " Number of events with at least  1 particle transmitted "
         << "(same charge as primary) = " << transmit[0] << " %" << G4endl;

  G4cout << "\n Number of events with the primary particle reflected = " << reflect[1] << " %"
         << G4endl;

  G4cout << " Number of events with at least  1 particle reflected "
         << "(same charge as primary) = " << reflect[0] << " %" << G4endl;

  if (fAnalysisManager)
  {
    if (fAnalysisManager->IsActive())
    {
      // compute width of the Gaussian central part of the MultipleScattering
      //
      if (fHistoManager->HistoExist(13))
      {
        G4cout << "\n MultipleScattering:"
               << "\n  rms proj angle of transmit primary particle = " << rmsMsc / mrad
               << " mrad (central part only)" << G4endl;

        G4cout << "  computed theta0 (Highland formula)          = " << ComputeMscHighland() / mrad
               << " mrad" << G4endl;

        G4cout << "  central part defined as +- " << fMscThetaCentral / mrad << " mrad; "
               << "  Tail ratio = " << tailMsc << " %" << G4endl;
      }

      G4cout.precision(prec);

      // normalize histograms
      //
      G4int ih = 1;
      G4double binWidth = fHistoManager->GetHistoBinWidth(ih);
      G4double unit = fHistoManager->GetHistoUnit(ih);
      G4double fac = unit / (TotNbofEvents * binWidth);
      fHistoManager->Scale(ih, fac);

      ih = 10;
      binWidth = fHistoManager->GetHistoBinWidth(ih);
      unit = fHistoManager->GetHistoUnit(ih);
      fac = unit / (TotNbofEvents * binWidth);
      fHistoManager->Scale(ih, fac);

      // Write histogram file
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

G4double Run::ComputeMscHighland() const
{
  // compute the width of the Gaussian central part of the MultipleScattering
  // projected angular distribution.
  // Eur. Phys. Jour. C15 (2000) page 166, formule 23.9

  G4double t =
    (fDetector->GetAbsorberThickness()) / (fDetector->GetAbsorberMaterial()->GetRadlen());
  if (t < DBL_MIN) return 0.;

  G4ParticleGun* particle = fPrimary->GetParticleGun();
  G4double T = particle->GetParticleEnergy();
  G4double M = particle->GetParticleDefinition()->GetPDGMass();
  G4double z = std::abs(particle->GetParticleDefinition()->GetPDGCharge() / eplus);

  G4double bpc = T * (T + 2 * M) / (T + M);
  G4double teta0 = 13.6 * MeV * z * std::sqrt(t) * (1. + 0.038 * std::log(t)) / bpc;
  return teta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillHisto(G4int histoId, G4double v)
{
  fHistoManager->FillHisto(histoId, v);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillHisto(G4int histoId, G4double v1, G4double v2)
{
  fHistoManager->FillHisto(histoId, v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Run::GetHistoBinWidth(G4int histoId) const
{
  return fHistoManager->GetHistoBinWidth(histoId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Run::GetHistoUnit(G4int histoId) const
{
  return fHistoManager->GetHistoUnit(histoId);
}
