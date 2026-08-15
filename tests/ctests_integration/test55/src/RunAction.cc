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

#include "G4Alpha.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4IonTable.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PhysicsList.hh"
#include "Run.hh"
#include "RunActionMessenger.hh"
#include "StepMax.hh"
#include "TestSeries.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det) : fDetector(det)
{
  fMessenger = new RunActionMessenger(this);

  // Book predefined histograms
  fHistoManager = new HistoManager(fDetector);

  if (isMaster)
  {
    testRange = new TestSeries("range", "Length", fDetector);
    testEnergyLoss = new TestSeries("mean energy loss primary", "Energy", fDetector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHistoManager;
  delete fMessenger;
  delete testRange;
  delete testEnergyLoss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(fDetector, fHistoManager);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  fHistoManager->Book();

  // histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (analysisManager->IsActive())
  {
    analysisManager->OpenFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{
  // print Run summary
  //
  if (isMaster)
  {
    fRun->EndOfRun(testRange, testEnergyLoss);
    if (!testDone && testDNA)
    {
      TestDNAStopping();
      testDone = true;
    }
    testRange->Print();
    testEnergyLoss->Print();
    testRange->Reset();
    testEnergyLoss->Reset();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CreateRangeTest(G4double refRange, G4double relError)
{
  testRange->CreateTestForRun(refRange, relError);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CreateEnergyLossTest(G4double refLoss, G4double relError)
{
  testEnergyLoss->CreateTestForRun(refLoss, relError);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestDNAStopping()
{
  G4StateManager* g4State = G4StateManager::GetStateManager();
  auto state = g4State->GetCurrentState();
  g4State->SetNewState(G4State_Init);

  auto proton = G4Proton::Proton();
  auto alpha = G4Alpha::Alpha();
  auto c12 = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(6, 12);
  G4double aRate = alpha->GetPDGMass() / CLHEP::proton_mass_c2;
  G4double cRate = c12->GetPDGMass() / CLHEP::proton_mass_c2;
  G4DNARuddIonisationExtendedModel modp(proton, "dedxp");
  G4DNARuddIonisationExtendedModel moda(alpha, "dedxa");
  G4DNARuddIonisationExtendedModel modc(c12, "dedxc12");
  G4DataVector v;
  modp.Initialise(proton, v);
  moda.Initialise(alpha, v);
  modc.Initialise(c12, v);
  G4int nd = 50;
  G4double emin = 100 * CLHEP::eV;
  G4double emax = 10 * CLHEP::MeV;
  G4double fac = std::log(emax / emin) / nd;
  fac = std::exp(fac);
  G4double f = 1.0 / (CLHEP::cm * CLHEP::cm);
  G4double e = emin;
  auto water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  G4int pre = G4cout.precision(5);
  G4cout << "======  TestDNAStopping ============" << G4endl;
  for (G4int i = 0; i <= nd; ++i)
  {
    G4double y1 = modp.CrossSectionPerVolume(water, proton, e, 0.0, emax);
    G4double y2 = moda.CrossSectionPerVolume(water, alpha, e, 0.0, emax);
    G4double y21 = moda.CrossSectionPerVolume(water, alpha, e * aRate, 0.0, emax);
    G4double y3 = modc.CrossSectionPerVolume(water, c12, e, 0.0, emax);
    G4double y31 = modc.CrossSectionPerVolume(water, c12, e * cRate, 0.0, emax);
    G4cout << std::setw(2) << i << "." << std::setw(5) << " E(keV)=" << std::setw(7)
           << e / CLHEP::keV << " sigP(cm-2)=" << std::setw(7) << y1 * f
           << " sigHe4(cm-2)=" << std::setw(7) << y2 * f << " sigC12(cm-2)=" << std::setw(7)
           << y3 * f;
    if (y1 > 0.0 && y1 < DBL_MAX)
    {
      G4cout << " He4/P=" << std::setw(7) << y21 / (y1 * 4) << " C12/P=" << y31 / (y1 * 36);
    }
    G4cout << G4endl;
    e *= fac;
  }
  G4cout.precision(pre);
  G4cout << "====================================" << G4endl;
  g4State->SetNewState(state);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
