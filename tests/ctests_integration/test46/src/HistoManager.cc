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
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 13/07/08
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if (fManager == nullptr)
  {
    static HistoManager manager;
    fManager = &manager;
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
{
  verbose = 0;
  n_evt = -1;
  maxEnergy = 50. * CLHEP::GeV;
  nBins = 100;
  histo = new Histo();
  nmax = 3;
  factorEcal = 1.03;
  factorHcal = 100.0;
  worldZ = 100. * mm;
  sizeXY = 0.0;
  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  nHisto = 20;
  G4double emax = maxEnergy / CLHEP::GeV;
  G4cout << "### HistoManager::bookHisto() emax(GeV)= " << emax << G4endl;
  histo->Add1D("0", "e0, Normalized Evis in central crystal", nBins, 0.55, 1.05);
  histo->Add1D("1", "e9, Normalized Evis in 3x3", 2 * nBins, 0.55, 1.05);
  histo->Add1D("2", "e25, Normalized  Evis in 5x5", 2 * nBins, 0.55, 1.05);
  histo->Add1D("3", "E0/E3x3;", nBins, 0.55, 1.05);
  histo->Add1D("4", "E0/E5x5", 2 * nBins, 0.55, 1.05);
  histo->Add1D("5", "E3x3/E5x5", 2 * nBins, 0.55, 1.05);
  histo->Add1D("6", "Normalized energy ECAL", 2 * nBins, 0.5, 1.2);
  histo->Add1D("7", "Energy (GeV) Ehcal", nBins, 0., emax * 0.1);
  histo->Add1D("8", "Energy (GeV) Eehcal", nBins, 0., emax * 0.1);
  histo->Add1D("9", "Energy (GeV) Eabshcal", nBins, 0., emax);
  histo->Add1D("10", "Energy computed (GeV)", nBins, 0., 2 * emax);
  histo->Add1D("11", "Normalized reconstructed energy", nBins, 0.8, 1.2);
  histo->Add1D("12", "Normalized reconstructed energy", nBins, 0., 2.);
  histo->Add1D("13", "ECAL hits log10(edep/MeV)", 60, -4., 2.);
  histo->Add1D("14", "Time (ns) ECAL", 20, 0., 20.);
  histo->Add1D("15", "Time (ns) HCAL", 50, 0., 50.);
  histo->Add1D("16", "Gamma energy at creation log10(E/MeV)", 200, -5., 5.);
  histo->Add1D("17", "Electron energy at creation log10(E/MeV)", 200, -5., 5.);
  histo->Add1D("18", "Proton energy at creation log10(E/MeV)", 200, -5., 5.);
  histo->Add1D("19", "Neutron energy at creation log10(E/MeV)", 200, -5., 5.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookNTuple()
{
  if (histo->IsNtupleActive())
  {
    histo->AddTupleI("n_charged");
    histo->AddTupleI("n_gamma");
    histo->AddTupleF("z1interact");
    histo->AddTupleF("beam");
    histo->AddTupleF("ecal0");
    histo->AddTupleF("ecal9");
    histo->AddTupleF("ecal25");
    histo->AddTupleF("hcalZ0_0");
    histo->AddTupleF("hcalZ0_9");
    histo->AddTupleF("hcalZ0_25");
    histo->AddTupleF("hcalZ1_0");
    histo->AddTupleF("hcalZ1_9");
    histo->AddTupleF("hcalZ1_25");
    histo->AddTupleF("hcalZ2_0");
    histo->AddTupleF("hcalZ2_9");
    histo->AddTupleF("hcalZ2_25");
    histo->AddTupleF("hcalZ3_0");
    histo->AddTupleF("hcalZ3_9");
    histo->AddTupleF("hcalZ3_25");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  n_evt = 0;
  n_step = 0;
  n_lowe = 0;

  m_gamma = 0.;
  m_e = 0.;
  m_h = 0.;
  m_n = 0.;

  for (G4int i = 0; i < 6; ++i)
  {
    ecal[i] = 0;
    erms[i] = 0;
    stat[i] = 0;
  }
  Eecal = 0;
  eecal = 0;
  hcal = 0;
  ehcal = 0;
  abshcal = 0;
  edepSum = 0;
  edepSum2 = 0;
  etotSum = 0;
  etotSum2 = 0;

  histo->SetVerbose(verbose);
  histo->Book();
  bookNTuple();

  if (verbose > 0)
    G4cout << "HistoManager: Histograms are booked and run has been started" << G4endl
           << "  BeginOfRun (After histo->book)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{
  G4cout << "HistoManager: End of run actions are started" << G4endl;
  G4String nam[8] = {"1x1 (GeV)", "3x3 (GeV)", "5x5 (GeV)",   "E1/E9    ",
                     "E1/E25   ", "E9/E25   ", "Erec(GeV)  ", "Etot(GeV)  "};

  // Average values
  G4cout << "========================================================" << G4endl;
  G4double x = (G4double)n_evt;
  if (n_evt > 0)
  {
    x = 1.0 / x;
  }

  G4double xe = x * n_lowe;
  G4double xs = x * n_step;

  G4cout << "Beam particle                           " << primaryDef->GetParticleName() << G4endl;
  G4cout << "Beam Energy(MeV)                        " << primaryKineticEnergy / MeV << G4endl;
  G4cout << "Number of events                        " << n_evt << G4endl;
  G4cout << std::setprecision(4) << "Average number of MIPS (Edep < 0.8 GeV) " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps                 " << xs << G4endl;
  G4cout << "Average number of neutron steps         " << x * m_n << G4endl;
  G4cout << "Average number of hadron steps          " << x * m_h << G4endl;
  G4cout << "Average number of gamma steps           " << x * m_gamma << G4endl;
  G4cout << "Average number of e+- steps             " << x * m_e << G4endl;
  G4cout << "==============  ECAL  ====================================" << G4endl;
  for (G4int j = 0; j < 6; ++j)
  {
    G4double xx = stat[j];
    if (xx > 0.0) xx = 1.0 / xx;
    G4double e = edep[j] * xx;
    edep[j] = e;
    G4double y = erms[j] * xx - e * e;
    G4double r = 0.0;
    G4double f = 1.0;
    if (j <= 2) f = 1.0 / CLHEP::GeV;
    if (y > 0.0) r = std::sqrt(y);
    erms[j] = r;
    G4cout << "  " << nam[j] << " = " << std::setw(10) << e * f << " +- " << std::setw(10)
           << f * r * std::sqrt(xx) << "    RMS= " << f * r << G4endl;
  }
  G4cout << "==============  HCAL  ====================================" << G4endl;

  G4double sum = edepSum * x;
  G4double y = edepSum2 * x - sum * sum;
  if (y > 0.)
    y = std::sqrt(y);
  else
    y = 0.0;
  G4double r = y * std::sqrt(x);

  G4double sum1 = etotSum * x;
  G4double y1 = etotSum2 * x - sum1 * sum1;
  if (y1 > 0.)
    y1 = std::sqrt(y1);
  else
    y1 = 0.0;
  G4double r1 = y1 * std::sqrt(x);

  G4cout << std::setprecision(4) << "Visible HCAL Edep(GeV)     =   " << x * hcal / CLHEP::GeV
         << G4endl;
  G4cout << std::setprecision(4) << "Visible HCAL e-Edep(GeV)   =   " << x * ehcal / CLHEP::GeV
         << G4endl;
  G4cout << std::setprecision(4)
         << "Total HCAL Edep(GeV)       =   " << x * (ehcal + abshcal) / CLHEP::GeV << " +- "
         << r1 / CLHEP::GeV << G4endl;
  G4cout << "==========================================================" << G4endl;
  G4cout << nam[6] << " = " << std::setw(10) << sum / CLHEP::GeV << " +- " << std::setw(10)
         << r / CLHEP::GeV << "    RMS= " << y / CLHEP::GeV << G4endl;
  G4cout << nam[7] << " = " << std::setw(10) << sum1 / CLHEP::GeV << " +- " << std::setw(10)
         << r1 / CLHEP::GeV << "    RMS= " << y1 / CLHEP::GeV << G4endl;
  G4cout << "==========================================================" << G4endl;
  G4double norm = primaryKineticEnergy;
  if (primaryDef->GetBaryonNumber() == 0) norm += primaryDef->GetPDGMass();
  G4cout << "Ecal/E0=   " << edep[2] / norm << "  RMS/E0(%)= " << erms[2] * 100. / norm << G4endl;
  G4cout << "Hcal/E0=   " << std::setw(10) << x * (ehcal + abshcal) / norm << " +- " << r1 / norm
         << G4endl;
  G4cout << "Erec/E0=   " << std::setw(10) << sum / norm << " +- " << r / norm << G4endl;
  G4cout << "Etot/E0=   " << std::setw(10) << sum1 / norm << " +- " << r1 / norm << G4endl;
  G4cout << "==========================================================" << G4endl;
  G4cout << G4endl;

  // normalise histograms
  for (G4int i = 0; i < nHisto; ++i)
  {
    histo->ScaleH1(i, x);
  }

  histo->Save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  e9 = 0.0;
  e25 = 0.0;
  e19 = 0.0;
  e125 = 0.0;
  e925 = 0.0;
  Eecal = 0.0;
  Ehcal = 0.0;
  Eehcal = 0.0;
  Eabshcal = 0.0;
  for (G4int i = 0; i < 25; ++i)
  {
    E[i] = 0.0;
  }
  for (G4int j = 0; j < 4; ++j)
  {
    for (G4int i = 0; i < 3; ++i)
    {
      HCAL[j][i] = 0.0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  ++n_evt;

  // ECAL ------------------------------------------------------------

  for (G4int i = 0; i < 25; ++i)
  {
    e25 += E[i];
    //    G4cout << "### E= " << E[i] << G4endl;
    if ((6 <= i && 8 >= i) || (11 <= i && 13 >= i) || (16 <= i && 18 >= i)) e9 += E[i];
    //    G4cout << "### e9= " << e9 << G4endl;
  }

  if (e25 < 0.8 * CLHEP::GeV)
  {
    ++n_lowe;
    if (verbose > 1) G4cout << "### in the event# " << n_evt << "  E25= " << e25 << G4endl;
  }
  G4double e0 = E[12];

  if (e9 > 0.0)
  {
    // compute energies
    edep[0] += e0;
    erms[0] += e0 * e0;
    edep[1] += e9;
    erms[1] += e9 * e9;
    edep[2] += e25;
    erms[2] += e25 * e25;
    stat[0] += 1;
    stat[1] += 1;
    stat[2] += 1;
    // compute ratios
    e19 = e0 / e9;
    e125 = e0 / e25;
    e925 = e9 / e25;
    edep[3] += e19;
    erms[3] += e19 * e19;
    edep[4] += e125;
    erms[4] += e125 * e125;
    edep[5] += e925;
    erms[5] += e925 * e925;
    stat[3] += 1;
    stat[4] += 1;
    stat[5] += 1;
  }
  // HCAL ------------------------------------------------------------
  hcal += Ehcal;
  ehcal += Eehcal;
  abshcal += Eabshcal;

  // Sum of ECAl + HCAL
  G4double edep0 = e25 * factorEcal + Ehcal * factorHcal;
  edepSum += edep0;
  edepSum2 += edep0 * edep0;

  // Total sum
  G4double etot = e25 + Ehcal + Eabshcal;
  etotSum += etot;
  etotSum2 += etot * etot;

  // Fill histo
  histo->Fill(0, e0 / primaryKineticEnergy);
  histo->Fill(1, e9 / primaryKineticEnergy);
  histo->Fill(2, e25 / primaryKineticEnergy);
  histo->Fill(3, e19);
  histo->Fill(4, e125);
  histo->Fill(5, e925);
  histo->Fill(6, Eecal / primaryKineticEnergy);
  histo->Fill(7, Ehcal / CLHEP::GeV);
  histo->Fill(8, Eehcal / CLHEP::GeV);
  histo->Fill(9, (Eabshcal + Ehcal) / CLHEP::GeV);
  histo->Fill(10, edep0 / CLHEP::GeV);
  histo->Fill(11, edep0 / primaryKineticEnergy);
  histo->Fill(12, edep0 / primaryKineticEnergy);
  if (histo->IsNtupleActive())
  {
    histo->FillTupleF(2, (G4float)e0);
    histo->FillTupleF(3, (G4float)e9);
    histo->FillTupleF(4, (G4float)e25);
    G4int idx = 4;
    for (G4int j = 0; j < 4; ++j)
    {
      for (G4int i = 0; i < 3; ++i)
      {
        histo->FillTupleF(++idx, (G4float)HCAL[j][i]);
      }
    }
    histo->AddRow();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* track)
{
  const G4ParticleDefinition* pd = track->GetDefinition();
  G4double e = track->GetKineticEnergy();
  G4double w = track->GetWeight();
  //  G4cout << "### KineticEnergy= " << e << G4endl;

  // Primary track
  if (0 == track->GetParentID())
  {
    primaryDef = pd;
  }
  else
  {
    e = std::log10(e / MeV);
    if (pd == G4Gamma::Gamma())
    {
      histo->Fill(16, e, w);
    }
    else if (pd == G4Electron::Electron())
    {
      histo->Fill(17, e, w);
    }
    else if (pd == G4Proton::Proton())
    {
      histo->Fill(18, e, w);
    }
    else if (pd == G4Neutron::Neutron())
    {
      histo->Fill(19, e, w);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddEcalHit(const G4ParticleDefinition* part, G4int copyNo, G4double de,
                              G4double time)
{
  n_step++;
  E[copyNo] += de;
  //  G4cout << "### edep= " << de << "   #copyNo =" << copyNo << G4endl;
  if (part->GetPDGMass() < MeV)
  {
    Eecal += de;
  }
  //  G4cout << "### Eecal= " << Eecal << G4endl;
  histo->Fill(13, std::log10(de / MeV));
  histo->Fill(14, time / CLHEP::ns, de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddHcalHit(const G4ParticleDefinition* part, G4int copyNo, G4double de,
                              G4double time, const G4ThreeVector& pos)
{
  ++n_step;
  Ehcal += de;
  if (part->GetPDGMass() < MeV)
  {
    Eehcal += de;
  }
  histo->Fill(15, time / CLHEP::ns, de);

  G4double x = std::abs(pos.x());
  G4double y = std::abs(pos.y());

  G4int i = 2;
  if (x <= xyHcalSell && y <= xyHcalSell)
  {
    i = 0;
  }
  else if (x <= 3 * xyHcalSell && y <= 3 * xyHcalSell)
  {
    i = 1;
  }

  G4int j = 0;
  if (copyNo >= 10)
  {
    j = 3;
  }
  else if (copyNo >= 5)
  {
    j = 2;
  }
  else if (copyNo >= 2)
  {
    j = 1;
  }
  HCAL[j][i] += de;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddHcalAbsorberHit(const G4ParticleDefinition*, G4double de)
{
  ++n_step;
  Eabshcal += de;
  //  G4cout << "### edepAbs= " << de << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddStep(const G4ParticleDefinition* part)
{
  if (part == G4Neutron::Neutron())
  {
    m_n += 1.0;
  }
  else if (part == G4Gamma::Gamma())
  {
    m_gamma += 1.0;
  }
  else if (part == G4Electron::Electron())
  {
    m_e += 1.0;
  }
  else
  {
    m_h += 1.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetBeamEnergy(G4double val)
{
  if (histo->IsNtupleActive())
  {
    histo->FillTupleF(1, (G4float)val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetFirstInteraction(G4double z, G4int n_char, G4int n_neut)
{
  if (histo->IsNtupleActive())
  {
    histo->FillTupleF(0, (G4float)(z - zHcalStart));
    histo->FillTupleI(0, n_char);
    histo->FillTupleI(1, n_neut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)
{
  verbose = val;
  histo->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int HistoManager::GetVerbose() const
{
  return verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetWorldLength(G4double val)
{
  worldZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double HistoManager::GetWorldLength() const
{
  return worldZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetBeamSizeXY(G4double val)
{
  sizeXY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double HistoManager::GetBeamSizeXY() const
{
  return sizeXY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetHcalStart(G4double x)
{
  zHcalStart = x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetHcalCell(G4double x)
{
  xyHcalSell = x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Fill(G4int id, G4double x, G4double w)
{
  histo->Fill(id, x, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetNbins(G4int val)
{
  if (val > 0.0) nBins = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetFactor1(G4double val)
{
  if (val > 0.0) factorEcal = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetFactor2(G4double val)
{
  if (val > 0.0) factorHcal = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetMaxEnergy(G4double val)
{
  if (val > 0.0)
  {
    maxEnergy = val;
    primaryKineticEnergy = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
