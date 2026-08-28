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

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* prim, HistoManager* histoMgr)
  : fDetector(det), fPrimary(prim), fHistoManager(histoMgr)
{
  fAnalysisManager = G4AnalysisManager::Instance();

  // initialisation
  fNEvt = 0;

  fEnergyDeposit = 0.;
  fTrackLength = 0.;
  fEnergyCharged = 0.;
  fEnergyNeutral = 0.;
  fEmin[0] = fEmin[1] = DBL_MAX;
  fEmax[0] = fEmax[1] = 0.;

  fNbSteps = 0;
  fNbCharged = 0;
  fNbNeutral = 0;

  fHistoManager->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  fNEvt += localRun->GetNumberOfEvent();

  fEnergyDeposit += localRun->fEnergyDeposit;
  fTrackLength += localRun->fTrackLength;
  fEnergyCharged += localRun->fEnergyCharged;
  fEnergyNeutral += localRun->fEnergyNeutral;
  fEmin[0] = std::min(fEmin[0], localRun->fEmin[0]);
  fEmin[1] = std::min(fEmin[1], localRun->fEmin[1]);
  fEmax[0] = std::max(fEmax[0], localRun->fEmax[0]);
  fEmax[1] = std::max(fEmax[1], localRun->fEmax[1]);

  fNbSteps += localRun->fNbSteps;
  fNbCharged += localRun->fNbCharged;
  fNbNeutral += localRun->fNbNeutral;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  // set the number of event only in case of a non MT Geant4 release
  //       ( set in Run::Merge otherwise )
#ifndef G4MULTITHREADED
  fNEvt += this->GetNumberOfEvent();
#endif

  G4int nbEvents = fNEvt;
  if (nbEvents == 0) return;

  G4Material* material = fDetector->GetMaterial();
  G4double length = fDetector->GetSize();
  G4double density = material->GetDensity();

  G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double eprimary = fPrimary->GetParticleGun()->GetParticleEnergy();

  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(eprimary, "Energy") << " through " << G4BestUnit(length, "Length") << " of "
         << material->GetName() << " (density: " << G4BestUnit(density, "Volumic Mass") << ")";
  G4cout << "\n ===========================================================\n";
  G4cout << G4endl;

  G4cout.precision(5);

  // track length
  //
  G4double trackLPerEvent = fTrackLength / nbEvents;
  G4double nbStepPerEvent = G4double(fNbSteps) / nbEvents;
  G4double stepSize = fTrackLength / fNbSteps;

  G4cout << "\n trackLength= " << G4BestUnit(trackLPerEvent, "Length")
         << "\t nb of steps= " << nbStepPerEvent << "  stepSize= " << G4BestUnit(stepSize, "Length")
         << G4endl;

  // charged secondaries (ionization, direct pair production)
  //
  G4double energyPerEvent = fEnergyCharged / nbEvents;
  G4double nbPerEvent = G4double(fNbCharged) / nbEvents;
  G4double meanEkin = 0.;
  if (fNbCharged) meanEkin = fEnergyCharged / fNbCharged;

  G4cout << "\n d-rays  : eLoss/primary= " << G4BestUnit(energyPerEvent, "Energy")
         << "\t  nb of d-rays= " << nbPerEvent << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
         << "  Tmin= " << G4BestUnit(fEmin[0], "Energy")
         << "  Tmax= " << G4BestUnit(fEmax[0], "Energy") << G4endl;

  // neutral secondaries (bremsstrahlung)
  //
  energyPerEvent = fEnergyNeutral / nbEvents;
  nbPerEvent = G4double(fNbNeutral) / nbEvents;
  meanEkin = 0.;
  if (fNbNeutral) meanEkin = fEnergyNeutral / fNbNeutral;

  G4cout << "\n brems   : eLoss/primary= " << G4BestUnit(energyPerEvent, "Energy")
         << "\t  nb of gammas= " << nbPerEvent << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
         << "  Tmin= " << G4BestUnit(fEmin[1], "Energy")
         << "  Tmax= " << G4BestUnit(fEmax[1], "Energy") << G4endl;

  // Computations below only for charged particles
  if (particle->GetPDGCharge() == 0.) return;

  G4EmCalculator emCal;

  // local energy deposit
  //
  energyPerEvent = fEnergyDeposit / nbEvents;
  //
  G4double r0 = emCal.GetRangeFromRestricteDEDX(eprimary, particle, material);
  G4double r1 = r0 - trackLPerEvent;
  G4double etry = eprimary - energyPerEvent;
  G4double efinal = 0.;
  if (r1 > 0. && etry > 0.0) efinal = GetEnergyFromRestrictedRange(r1, particle, material, etry);
  G4double dEtable = eprimary - efinal;
  G4double ratio = 0.;
  if (dEtable > 0.) ratio = energyPerEvent / dEtable;

  G4cout << "\n deposit : eLoss/primary= " << G4BestUnit(energyPerEvent, "Energy")
         << "\t <dEcut > table= " << G4BestUnit(dEtable, "Energy")
         << "   ---> simul/reference= " << ratio << G4endl;

  // total energy transferred
  //
  G4double energyTotal = fEnergyDeposit + fEnergyCharged + fEnergyNeutral;
  energyPerEvent = energyTotal / nbEvents;
  //
  r0 = emCal.GetCSDARange(eprimary, particle, material);
  r1 = r0 - trackLPerEvent;
  etry = eprimary - energyPerEvent;
  efinal = 0.;
  // G4cout << "r0= " << r0 << "  r1= " << r1 << "  " << particle->GetParticleName()
  //	 << " etry= " << etry << "  " << material->GetName() << " e0= " << eprimary << G4endl;
  if (r1 > 0.0 && etry > 0.0) efinal = GetEnergyFromCSDARange(r1, particle, material, etry);
  dEtable = eprimary - efinal;
  ratio = 0.;
  if (dEtable > 0.) ratio = energyPerEvent / dEtable;

  G4cout << "\n total   : eLoss/primary= " << G4BestUnit(energyPerEvent, "Energy")
         << "\t <dEfull> table= " << G4BestUnit(dEtable, "Energy")
         << "   ---> simul/reference= " << ratio << G4endl;

  G4cout.precision(prec);

  G4int Z = 29;
  for (G4int ii = 1; ii < 4; ++ii)
  {
    G4double e = ii * MeV;

    G4double xs = emCal.ComputeShellIonisationCrossSectionPerAtom("proton", Z, fKShell, e);
    G4cout << "K-shell x-section for proton in barns " << xs / barn << "  E(MeV)= " << e
           << "  material " << material->GetName() << G4endl;
    G4double xs1 = emCal.ComputeShellIonisationCrossSectionPerAtom("proton", Z, fL1Shell, e);
    G4double xs2 = emCal.ComputeShellIonisationCrossSectionPerAtom("proton", Z, fL2Shell, e);
    G4double xs3 = emCal.ComputeShellIonisationCrossSectionPerAtom("proton", Z, fL3Shell, e);
    G4cout << "L1-shell x-section in barns " << xs1 / barn << G4endl;
    G4cout << "L2-shell x-section in barns " << xs2 / barn << G4endl;
    G4cout << "L3-shell x-section in barns " << xs3 / barn << G4endl;

    G4double xn = xs * 0.01 * mm * material->GetTotNbOfAtomsPerVolume();
    G4cout << "N K-shell vacancies in 0.01 mm absorber " << xn << G4endl;
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillHisto(G4int histoId, G4double v1, G4double v2)
{
  if (fAnalysisManager) fHistoManager->FillHisto(histoId, v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Run::GetEnergyFromRestrictedRange(G4double range, G4ParticleDefinition* particle,
                                           G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;

  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err = 1., errmax = 0.00001;
  G4int iter = 0, itermax = 10;
  while (err > errmax && iter < itermax)
  {
    iter++;
    Energy -= dE;
    r = emCal.GetRangeFromRestricteDEDX(Energy, particle, material);
    dr = r - range;
    dEdx = emCal.GetDEDX(Energy, particle, material);
    dE = dEdx * dr;
    err = std::abs(dE) / Energy;
  }
  if (iter == itermax)
  {
    G4cout << "\n  ---> warning: RunAction::GetEnergyFromRestRange() did not converge"
           << "   Etry = " << G4BestUnit(Etry, "Energy")
           << "   Energy = " << G4BestUnit(Energy, "Energy") << "   err = " << err
           << "   iter = " << iter << G4endl;
  }

  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Run::GetEnergyFromCSDARange(G4double range, G4ParticleDefinition* particle,
                                     G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;

  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err = 1., errmax = 0.00001;
  G4int iter = 0, itermax = 10;
  while (err > errmax && iter < itermax)
  {
    iter++;
    Energy -= dE;
    r = emCal.GetCSDARange(Energy, particle, material);
    dr = r - range;
    dEdx = emCal.ComputeTotalDEDX(Energy, particle, material);
    dE = dEdx * dr;
    err = std::abs(dE) / Energy;
  }
  if (iter == itermax)
  {
    G4cout << "\n  ---> warning: RunAction::GetEnergyFromCSDARange() did not converge"
           << "   Etry = " << G4BestUnit(Etry, "Energy")
           << "   Energy = " << G4BestUnit(Energy, "Energy") << "   err = " << err
           << "   iter = " << iter << G4endl;
  }

  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
