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
// ClassName:   test31Histo
//
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31Histo.hh"

#include "G4ASTARStopping.hh"
#include "G4Alpha.hh"
#include "G4AtomicShells.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include "G4EmCorrections.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4NucleiProperties.hh"
#include "G4PSTARStopping.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Pow.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4WaterStopping.hh"

#include "Histo.hh"

#include <fstream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::fManager = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::GetPointer()
{
  if (!fManager)
  {
    fManager = new test31Histo();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::test31Histo()
{
  verbose = 0;
  nHisto = 0;
  maxEnergy = 0.0;
  nTuple = false;
  tables = true;
  extra = nullptr;
  histo = new Histo();
  histoBooked = false;
  G4NistManager* mman = G4NistManager::Instance();
  mman->FindOrBuildMaterial("G4_WATER");
  for (G4int z = 1; z < 93; ++z)
  {
    auto elm = mman->FindOrBuildElement(z, false);
    G4String nam = "G4_" + elm->GetSymbol();
    mman->FindOrBuildMaterial(nam, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::~test31Histo()
{
  delete extra;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::BeginOfHisto(G4int num)
{
  if (0 < verbose) G4cout << "test31Histo # " << num << " started " << G4endl;
  zend = 0.0;
  zend2 = 0.0;
  zEvt = 0.0;
  etot = 0.0;

  n_evt = 0;
  n_elec = 0;
  n_posit = 0;
  n_gam = 0;
  n_step = 0;

  n_charged_leak = 0;
  n_gam_leak = 0;
  n_charged_back = 0;
  n_gam_back = 0;

  n_mumu = 0;
  n_pipi = 0;

  if (nullptr == extra)
  {
    extra = new G4EnergyLossForExtrapolator(1);
  }
  bookHisto();

  if (verbose > 0)
  {
    G4cout << "test31Histo: Histograms are booked and run has been started" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::EndOfHisto()
{
  G4cout << "test31Histo: End of run actions" << G4endl;
  //  G4cout << "1.5kT(eV)= " << 1.5*k_Boltzmann*STP_Temperature/eV << G4endl;

  G4double Spin = beamParticle->GetPDGSpin();
  G4double Charge = beamParticle->GetPDGCharge();

  // Zend average
  G4cout << "====================================================================" << G4endl;
  G4cout << "Initial particle               " << beamParticle->GetParticleName()
         << "   Ekin(GeV)= " << beamEnergy / GeV << " Spin= " << Spin << " Charge= " << Charge
         << G4endl;
  if (zEvt > 0.0)
  {
    zend /= zEvt;
    zend2 /= zEvt;
    zend2 -= zend * zend;
    G4double sig = 0.0;
    if (zend2 > 0.) sig = std::sqrt(zend2);
    zend2 = sig / std::sqrt(zEvt);
    G4cout << std::setprecision(5) << "Range(mm)= " << zend / mm << "; Stragling(mm)= " << sig / mm
           << std::setprecision(2) << " +- " << zend2 / mm << "    " << zEvt << " events for range"
           << G4endl;
  }
  G4double x = (G4double)n_evt;
  if (n_evt > 0) x = 1.0 / x;
  etot *= x;
  G4double xe = x * (G4double)n_elec;
  G4double xg = x * (G4double)n_gam;
  G4double xp = x * (G4double)n_posit;
  G4double xs = x * (G4double)n_step;
  G4double xcl = x * (G4double)n_charged_leak;
  G4double xgl = x * (G4double)n_gam_leak;
  G4double xcb = x * (G4double)n_charged_back;
  G4double xgb = x * (G4double)n_gam_back;
  G4double xmu = x * (G4double)n_mumu;
  G4double xpi = x * (G4double)n_pipi;
  G4cout << "Number of events               " << n_evt << G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit         " << etot / MeV << " MeV"
         << G4endl;
  G4cout << std::setprecision(4) << "Average number of e-           " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma        " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+           " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps        " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak charged " << xcl << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak gamma   " << xgl << G4endl;
  G4cout << std::setprecision(4) << "Average number of back charged " << xcb << G4endl;
  G4cout << std::setprecision(4) << "Average number of back gamma   " << xgb << G4endl;
  G4cout << std::setprecision(4) << "Average number of mu+mu-       " << xmu << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi+pi-       " << xpi << G4endl;
  G4cout << "====================================================================" << G4endl;

  if (tables)
  {
    TableControl();
    tables = false;
  }
  MuonTest();
  //  ElectronTest();
  if (0 < nHisto)
  {
    // normalise histograms
    for (G4int i = 0; i < nHisto; ++i)
    {
      histo->ScaleH1(i, x);
    }
    histo->Save();
  }
  IonLSData();

  G4cout << "=========   End of tets31Histo  ============================" << G4endl;
  /*
  G4cout << G4endl;
  for(G4int Z=1; Z<101; ++Z) {

    G4int nn = G4AtomicShells::GetNumberOfElectrons(Z, 0);
    G4cout << nn << ",  ";
    if(Z/10*10 == Z) {G4cout << G4endl;}
  }
  G4cout << G4endl;

  G4cout << "Z=92, A=238 M(GeV)= "
   <<G4NucleiProperties::GetNuclearMass(238,92)/GeV << G4endl;
  G4cout << "Z=82, A=208 M(GeV)= "
   <<G4NucleiProperties::GetNuclearMass(208,82)/GeV << G4endl;
  G4cout << "Z=54, A=136 M(GeV)= "
   <<G4NucleiProperties::GetNuclearMass(136,54)/GeV << G4endl;
  G4cout << "Z=26, A=56 M(GeV)= "
   <<G4NucleiProperties::GetNuclearMass(56,26)/GeV << G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveEvent()
{
  // if(nTuple) histo->addRow(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& /*parname*/, G4double /*val*/)
{
  // if(nTuple) histo->fillTuple(0, parname, val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::bookHisto()
{
  zmax = (AbsorberThickness + gap) * NumberOfAbsorbers / mm;
  G4cout << "test31Histo: "
         << " AbsThick(mm)= " << AbsorberThickness / mm << " Nabs= " << NumberOfAbsorbers
         << " zmax= " << zmax << " nHisto= " << nHisto << G4endl;

  // Creating an 1-dimensional histograms in the root directory of the tree

  G4double em = maxEnergy / MeV;

  if (histoBooked && 0 < nHisto)
  {
    histo->SetHisto1D(0, NumberOfAbsorbers, 0.0, zmax / mm, mm);
    histo->SetHisto1D(1, 50, 0.0, em, MeV);
    histo->SetHisto1D(2, 36, 0.0, 180., degree);
    histo->SetHisto1D(3, 50, 0.0, em, MeV);
    histo->SetHisto1D(4, 36, 0.0, 180., degree);
    histo->SetHisto1D(5, 50, 0.0, 10., degree);
    histo->SetHisto1D(6, 100, -em, em, MeV);
  }
  else
  {
    histo->Add1D("10", "Energy deposit (MeV) in absorber (mm)", NumberOfAbsorbers, 0.0, zmax / mm,
                 mm);
    histo->Add1D("11", "Energy (MeV) of secondary electrons", 50, 0.0, em, MeV);
    histo->Add1D("12", "Theta (degrees) of delta-electrons", 36, 0.0, 180., degree);
    histo->Add1D("13", "Energy (MeV) of secondary gamma", 50, 0.0, em, MeV);
    histo->Add1D("14", "Theta (degrees) of secondary gamma", 36, 0.0, 180., degree);
    histo->Add1D("15", "Theta (degrees) of primary", 50, 0.0, 10., degree);
    histo->Add1D("16", "Delta Energy (MeV) of Reconstruction", 100, -em, em, MeV);
  }
  for (G4int i = 0; i < nHisto; ++i)
  {
    histo->Activate(i, true);
  }

  // if(nTuple){
  //   histo->addTuple( "100", "Range/Energy",
  //    "float tkin mass beta xend, yend, zend, ltpk, tend, teta, loss, dedx, back, leak, edep" );
  // }

  histoBooked = true;
  histo->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEnergy(G4double edep, G4double z)
{
  etot += edep;
  histo->Fill(0, z, edep / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEndPoint(G4double z)
{
  zend += z;
  zend2 += z * z;
  zEvt += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddDeltaElectron(const G4DynamicParticle* elec)
{
  n_elec++;
  histo->Fill(1, elec->GetKineticEnergy(), 1.0);
  histo->Fill(2, elec->GetMomentumDirection().theta(), 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddPhoton(const G4DynamicParticle* ph)
{
  n_gam++;
  histo->Fill(3, ph->GetKineticEnergy(), 1.0);
  histo->Fill(4, (ph->GetMomentumDirection()).theta(), 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleLeak(const G4Track* track)
{
  const G4DynamicParticle* dp = track->GetDynamicParticle();
  if (dp->GetDefinition() == G4Gamma::Gamma())
  {
    ++n_gam_leak;
  }
  else if (dp->GetCharge() != 0.0)
  {
    n_charged_leak++;
    if (track->GetTrackID() == 1)
    {
      G4double tet = (dp->GetMomentumDirection()).theta();
      histo->Fill(5, tet, 1.0);
      G4double e0 = track->GetVertexKineticEnergy();
      G4double e1 = dp->GetKineticEnergy();
      G4double e2 = extra->EnergyBeforeStep(e1, zmax, absMaterial, dp->GetDefinition());
      if (n_evt < 10)
        G4cout << "Extrapolation of primary " << dp->GetDefinition()->GetParticleName()
               << " E0(MeV)= " << e0 / MeV << " E1(MeV)= " << e1 / MeV
               << " Erec-E0(MeV)= " << (e2 - e0) / MeV << G4endl;
      histo->Fill(6, e2 - e0, 1.0);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleBack(const G4Track* track)
{
  if (track->GetDynamicParticle()->GetDefinition() == G4Gamma::Gamma())
  {
    n_gam_back++;
  }
  else if (track->GetDynamicParticle()->GetCharge() != 0.0)
  {
    n_charged_back++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::TableControl()
{
  /*
  G4cout << ">>>>>> Mass test 2n:  "
   << G4NucleiProperties::GetNuclearMass(2,0) - 2*neutron_mass_c2
   << " MeV" << G4endl;
  G4cout << ">>>>>> Mass test 2p:  "
   << G4NucleiProperties::GetNuclearMass(2,2) - 2*proton_mass_c2
   << " MeV" << G4endl;
  */

  G4NistManager* mman = G4NistManager::Instance();
  G4EmCorrections* emc = G4LossTableManager::Instance()->EmCorrections();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  G4double etest = 9.66 * CLHEP::MeV;
  G4double ecut = 0.15 * CLHEP::MeV;
  G4int Zt = 79;
  // mman->SetVerbose(2);
  const G4Element* elm = mman->FindOrBuildElement(Zt, false);
  // const G4Element* elm1 = mman->FindOrBuildElement(62,true);
  G4double xsb = cal.ComputeCrossSectionPerAtom(etest, "e-", "eBrem", elm, ecut);
  // const G4Material* mat79 = mman->FindOrBuildMaterial("G4_Au");
  G4cout << "##### e- bremstrahlung x-section off " << elm->GetName()
         << " for E(MeV)=" << etest / CLHEP::MeV << " Ecut(MeV)=" << ecut / CLHEP::MeV
         << " sig(b)= " << xsb / barn << G4endl;

  // parameters
  // G4double tmin = 1.*keV;
  // G4double tmax = 1.*GeV;
  // G4int    nbin = 60;
  G4String ion_name = "ionIoni";
  G4String e_name = "eIoni";
  G4String h_name = "hIoni";
  G4String mu_name = "muIoni";
  G4String proc_name = "eIoni";
  G4String part_name = "proton";

  const G4ParticleDefinition* part = cal.FindParticle(part_name);
  if (nullptr == part) return;
  // cal.PrintDEDXTable(part);
  // cal.PrintRangeTable(part);
  // cal.PrintInverseRangeTable(part);

  G4String mat_name = "G4_WATER";
  G4Material* mat = mman->FindOrBuildMaterial(mat_name);
  // mat->SetChemicalFormula("H_2O");
  G4double density = mat->GetDensity();
  G4double fact = gram / (MeV * cm2 * density);
  /*
  G4cout << "========================================" << G4endl;
  G4cout << std::setprecision(6) << mat->GetName()
   << "  amu_c2= " << amu_c2 <<  " amu= " << amu << G4endl;
  G4cout << "GetZ= " << mat->GetZ() << " GetA= " << mat->GetA()
   << " GetA()*mole/g= " <<  mat->GetA()*mole/g
   << " GetAeffective()= " <<  mat->GetAtomicMassAmu()
   << G4endl;

  G4cout << elm->GetName() << G4endl;
  G4cout << "79 GetZ= " << elm->GetZ() << " GetA= " << elm->GetA()
   << " GetA(mole/g)= " << elm->GetA()*mole/g
   << " GetA(amu)= " << elm->GetA()/(Avogadro*amu)
   << " GetN()= " <<  elm->GetN() << G4endl;
  G4cout << "Mat79 GetZ= " << mat79->GetZ() << " GetA= " << mat79->GetA()
   << " GetA= " << mat79->GetA()*mole/g
   << " GetAeff()= " << mat79->GetAtomicMassAmu() << G4endl;
  G4cout << "62 GetZ= " << elm1->GetZ() << " GetA= " << elm1->GetA()
   << " GetA= " << elm1->GetA()*mole/g
   << " GetN()= " <<  elm1->GetN() << G4endl;
  G4cout << "NIST GetZ(Au)= " << mman->GetZ("Au")
   << " GetA(Au)= " << mman->GetAtomicMassAmu("Au")
   << " GetAtomicMassAmu(79)= " << mman->GetAtomicMassAmu(79)
   << " GetA(Sm)= " << mman->GetAtomicMassAmu("Sm")
   << " GetAtomicMassAmu(Z)= " << mman->GetAtomicMassAmu(62)
   << G4endl;
  G4cout << " GetAtomicMassAmu(12)= " << mman->GetAtomicMassAmu(6)
         << " GetAtomicMassAmu(1)= " << mman->GetAtomicMassAmu(1)
   << std::setprecision(4) << G4endl;
  G4cout << "========================================" << G4endl;
  */

  const G4int ne = 57;
  G4double e0[ne] = {0.0001, 0.001,  0.0025, 0.005,  0.008,  0.01,   0.02,   0.025, 0.03,  0.04,
                     0.05,   0.08,   0.1,    0.2,    0.3,    0.4,    0.5,    0.8,   1.0,   1.1,
                     1.2,    1.3,    1.4,    1.5,    1.6,    1.7,    1.8,    1.9,   2.0,   2.1,
                     2.2,    2.3,    2.4,    2.5,    3.0,    5.0,    10.,    12.,   15.0,  20.0,
                     30.,    50.,    100.,   200.,   300.,   400.,   500.,   1000., 2000., 3000.,
                     5000.,  10000., 15000., 20000., 30000., 50000., 100000.};
  const G4int np = 10;
  G4String namep[np] = {"e-",    "mu+", "pi-",  "proton", "anti_proton",
                        "alpha", "C12", "Ar40", "Pb208",  "Rn222"};

  G4double xxx = 80 * keV;
  G4IonTable::GetIonTable()->GetIon(6, 12, 0.0);
  G4IonTable::GetIonTable()->GetIon(18, 40, 0.0);
  G4IonTable::GetIonTable()->GetIon(82, 208, 0.0);
  const G4ParticleDefinition* part3 = G4IonTable::GetIonTable()->GetIon(9, 19, 0.0);
  const G4ParticleDefinition* part4 = G4IonTable::GetIonTable()->GetIon(86, 222, 0.0);
  if (nullptr != part3)
  {
    G4cout << "***** Ion: " << part3->GetParticleName() << " mass: " << part3->GetPDGMass()
           << " Charge: " << part3->GetPDGCharge() << " E(MeV)=" << xxx / CLHEP::MeV << " off "
           << mat->GetName() << G4endl;
    G4cout << "      ElecDEDX= " << cal.ComputeElectronicDEDX(xxx, part3, mat)
           << "  NucDEDX= " << cal.ComputeNuclearDEDX(xxx, part3, mat) << G4endl;
  }
  if (part4)
  {
    G4cout << "***** Ion: " << part4->GetParticleName() << " mass: " << part4->GetPDGMass()
           << " Charge: " << part4->GetPDGCharge() << " E(MeV)=" << xxx / CLHEP::MeV << " off "
           << mat->GetName() << G4endl;
    G4cout << "      ElecDEDX= " << cal.ComputeElectronicDEDX(xxx, part4, mat)
           << "  NucDEDX= " << cal.ComputeNuclearDEDX(xxx, part4, mat) << G4endl;
  }
  G4cout << std::setprecision(5) << G4endl;
  G4int ii1 = 0;
  G4int ii2 = 7;

  G4ASTARStopping ast;
  G4PSTARStopping pst;
  G4WaterStopping wst;
  G4int idxa = ast.GetIndex(mat);
  G4int idxp = pst.GetIndex(mat);
  G4double fact1 = mm / MeV;

  for (G4int ii = ii1; ii < ii2; ++ii)
  {
    const G4ParticleDefinition* part1 = cal.FindParticle(namep[ii]);
    if (nullptr == part1) break;

    if (ii == 0)
      proc_name = e_name;
    else if (ii == 1)
      proc_name = mu_name;
    else if (ii == 2)
      proc_name = h_name;
    else if (ii == 5)
    {
      proc_name = ion_name;
    }

    G4double AA = 1.0;
    G4int ZZ = 1;
    if (ii >= 6)
    {
      ZZ = part1->GetAtomicNumber();
      AA = mman->GetAtomicMassAmu(ZZ);
    }
    else if (ii == 5)
    {
      ZZ = 2;
      AA = mman->GetAtomicMassAmu(2);
    }

    G4cout << "================================================================================"
           << G4endl;
    G4cout << "   Tables control for " << namep[ii] << " Z=" << ZZ << " A=" << AA << "  Material "
           << mat_name << G4endl;
    G4cout << "================================================================================"
           << G4endl;

    G4cout << "  N   E(MeV)  Esc(MeV)  NIST/dEdx "
           << "dEdx(MeV/mm)  NIST(MeV/mm)  dEdx(MeV*cm^2/g)" << G4endl;

    for (G4int ij = 0; ij < ne; ++ij)
    {
      G4double e = e0[ij];
      G4double e1 = e;
      if (ii >= 5) e1 *= AA;
      G4double dedx = cal.ComputeElectronicDEDX(e1, part1, mat, e1);
      G4double dedx2 = 0.0;
      if (ii == 3)
      {
        dedx2 = pst.GetElectronicDEDX(idxp, e1) * density;
      }
      else if (ii == 5)
      {
        dedx2 = ast.GetElectronicDEDX(idxa, e1) * density;
      }
      else if (ii > 5)
      {
        dedx2 = wst.GetElectronicDEDX(ZZ, e1);
      }

      G4cout << std::setw(3) << ij << "." << std::setw(8) << e / MeV << std::setw(11) << e1 / MeV
             << std::setw(10) << dedx2 / dedx << std::setw(10) << dedx << std::setw(14) << dedx2
             << std::setw(14) << dedx * fact << G4endl;
    }
  }
  G4String mat_name1 = "G4_Si";
  G4Material* mat1 = mman->FindOrBuildMaterial(mat_name1);
  fact = gram / (MeV * cm2 * mat1->GetDensity());
  for (G4int ii = 0; ii < np; ++ii)
  {
    const G4ParticleDefinition* part2 = cal.FindParticle(namep[ii]);
    if (nullptr == part2) break;
    G4double AA = 1.0;
    G4int ZZ = 1;
    if (ii >= 6)
    {
      ZZ = part2->GetAtomicNumber();
      AA = mman->GetAtomicMassAmu(ZZ);
    }
    else if (ii == 5)
    {
      ZZ = 2;
      AA = mman->GetAtomicMassAmu(2);
    }
    G4cout << "======================================================" << G4endl;
    G4cout << "   Material " << mat_name1 << "  Projectile " << part2->GetParticleName()
           << " Z=" << ZZ << " A=" << AA << G4endl;
    G4cout << "======================================================" << G4endl;
    G4cout << "  N   E(MeV) dEdx(MeV/mm) dEdx(MeV*cm^2/g) "
           << "NIEL(MeV/mm) NIEL(MeV*cm^2/g)" << G4endl;

    for (G4int ij = 0; ij < ne; ++ij)
    {
      G4double e = e0[ij];
      G4double dedx = cal.ComputeElectronicDEDX(e, part2, mat1, e);
      G4double dedxn = (ii >= 3) ? cal.ComputeNuclearDEDX(e, part2, mat1) : 0.0;
      G4cout << std::setw(3) << ij << "." << std::setw(8) << e / MeV << std::setw(12)
             << dedx * fact1 << std::setw(13) << dedx * fact << std::setw(12) << dedxn * fact1
             << std::setw(14) << dedxn * fact << G4endl;
    }
  }
  G4Material* mat2 = mman->FindOrBuildMaterial("G4_Cr");
  auto part31 = cal.FindParticle(namep[3]);
  auto part41 = cal.FindParticle(namep[4]);
  G4cout << "======================================================" << G4endl;
  G4cout << "   Material " << mat2->GetName() << "  Projectiles p, pbar " << G4endl;
  G4cout << "======================================================" << G4endl;
  G4cout << "  N   E(MeV)  p dEdx(MeV/mm)   pbar dEdx(MeV/mm) " << G4endl;
  G4cout << "======================================================" << G4endl;

  for (G4int ij = 0; ij < ne; ++ij)
  {
    G4double e = e0[ij];
    G4double dedxp = cal.ComputeElectronicDEDX(e, part31, mat2, e) * fact1;
    G4double dedxb = cal.ComputeElectronicDEDX(e, part41, mat2, e) * fact1;
    G4cout << std::setw(3) << ij << "." << std::setw(8) << e / MeV << std::setw(12) << dedxp
           << std::setw(16) << dedxb << G4endl;
  }

  //    G4bool icorr = true;
  G4bool icorr = false;
  if (icorr)
  {
    G4cout << "================================================================" << G4endl;
    G4cout << "             Ionisation Corrections" << G4endl;
    G4cout << "================================================================" << G4endl;

    const G4int nmm = 7;
    const G4String nmat[nmm] = {"G4_H", "G4_C", "G4_Al", "G4_Cu", "G4_Ag", "G4_Au", "G4_Pb"};
    const G4int kkk = 25;
    G4double ek[kkk] = {0.3,  1.0,   1.25,  1.5,    1.75,    2.0,      2.5,      3.0, 3.5,
                        4.0,  4.5,   5.0,   6.5,    8.0,     12.5,     20.0,     30., 100.,
                        300., 1000., 3000., 10000., 100000., 1000000., 10000000.};

    G4double mass = part->GetPDGMass();

    G4double aL, L0, L1, L2, KS, LS, S, del, mk, dedx, fac, fs(0);

    for (G4int i = 0; i < nmm; ++i)
    {
      mat = mman->FindOrBuildMaterial(nmat[i]);
      fac =
        2.0 * twopi_mc2_rcl2 * (mat->GetElectronDensity()) * gram / (MeV * cm2 * mat->GetDensity());
      G4cout << "   New Material  " << mat->GetName() << G4endl;
      for (G4int j = 0; j < kkk; ++j)
      {
        G4double e = ek[j] * MeV;
        G4double tau = e / mass;
        G4double gamma = 1.0 + tau;
        G4double beta2 = tau * (tau + 2.0) / (gamma * gamma);
        L0 = emc->Bethe(part, mat, e);
        // Spin = emc->SpinCorrection(part,mat,e);
        KS = emc->KShellCorrection(part, mat, e);
        LS = emc->LShellCorrection(part, mat, e);
        S = emc->ShellCorrection(part, mat, e);
        // S0 = emc->ShellCorrectionSTD(part,mat,e);
        L1 = emc->BarkasCorrection(part, mat, e);
        L2 = emc->BlochCorrection(part, mat, e);
        del = -0.5 * emc->DensityCorrection(part, mat, e);
        mk = 0.5 * emc->MottCorrection(part, mat, e);
        aL = L0 + L1 + L2 + fs + del - S;
        dedx = aL * fac / beta2;

        G4cout << j + 1 << ". " << ek[j] << " MeV "
               << " L0= "
               << L0
               //		 << " Spin= " << Spin
               << " KSh= " << KS << " LSh= " << LS << " Sh= "
               << S
               //		 << " Sh0= " << S0
               << " L1= " << L1 << " L2= " << L2 << " dn= " << del << " mott= "
               << mk
               //	    	 << " fs= " << fs
               << " L= " << aL << " dedx= "
               << dedx
               //     << " dedx0= " << dedx0
               << G4endl;
      }
      G4cout << "==============================================================" << G4endl;
    }
  }

  G4cout << "==============  End of Table Control ===============================" << G4endl;

  G4bool ish = false;
  if (ish)
  {
    const G4ParticleDefinition* proton = cal.FindParticle("proton");
    G4BraggModel bragg;
    G4BetheBlochModel bethe;
    G4DataVector empty;
    bragg.Initialise(proton, empty);
    bethe.Initialise(proton, empty);
    const G4Element* elm1;
    const G4Material* ma;

    G4double e = 2.0 * MeV;
    G4double tau = e / proton_mass_c2;
    G4double gam = tau + 1.0;
    G4double bg2 = tau * (tau + 2.0);
    G4double beta2 = bg2 / (gam * gam);
    G4double eta = beta2 / (fine_structure_const * fine_structure_const);
    G4double fact2 = 0.5 * beta2 * eta * eta / twopi_mc2_rcl2;
    G4double dedx0, dedx1;

    for (G4int z = 1; z < 93; ++z)
    {
      elm1 = mman->FindOrBuildElement(z, false);
      G4String nam = "G4_" + elm1->GetSymbol();
      ma = mman->FindOrBuildMaterial(nam, false);
      //      G4cout << "Elm " << elm << "  mat " << ma << "   " << nam << G4endl;
      dedx0 = bragg.ComputeDEDXPerVolume(ma, proton, e, GeV);
      dedx1 = bethe.ComputeDEDXPerVolume(ma, proton, e, GeV);
      G4cout << " " << (dedx1 - dedx0) * fact2 / (ma->GetElectronDensity()) << ",";
      if (z / 10 * 10 == z) G4cout << G4endl;
    }
    G4cout << G4endl;
    G4double fe = 8.0 * MeV;
    G4double ftau = fe / proton_mass_c2;
    G4double fgam = ftau + 1.0;
    G4double fbg2 = ftau * (ftau + 2.0);
    G4double fbeta2 = fbg2 / (fgam * fgam);
    G4double ffact = 0.5 * fbeta2 / twopi_mc2_rcl2;
    G4double fdedx0, fdedx1, s0, s1;
    G4cout << G4endl;
    for (G4int z = 1; z < 93; z++)
    {
      elm1 = mman->FindOrBuildElement(z, false);
      G4String nam = "G4_" + elm1->GetSymbol();
      ma = mman->FindOrBuildMaterial(nam, false);
      fdedx0 = bragg.ComputeDEDXPerVolume(ma, proton, e, GeV) * ffact / ma->GetElectronDensity();
      fdedx1 = bethe.ComputeDEDXPerVolume(ma, proton, e, GeV) * ffact / ma->GetElectronDensity();
      s0 = fdedx1 - fdedx0;
      s1 = emc->ShellCorrectionSTD(proton, ma, fe);

      G4cout << " " << (s1 * std::log(tau) - s0 * std::log(ftau)) / (ftau - tau) << ",";
      // << " s0= " << s0 << " s1= " << s1 << G4endl;;
      if (z / 10 * 10 == z) G4cout << G4endl;
    }
    G4cout << G4endl;
  }

  G4bool ihist = false;
  if (ihist)
  {
    G4cout << "=================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "=================================================================" << G4endl;

    G4String nmk[7] = {"G4_Be", "G4_Al", "G4_Si", "G4_Ge", "G4_Fe", "G4_Ag", "G4_Au"};
    const G4String partc[2] = {"proton", "alpha"};
    const G4String proc[2] = {"hIoni", "ionIoni"};

    G4double e, se, sn, st, mce, mcn, mct;
    char line[200];

    for (G4int ii = 0; ii < 2; ++ii)
    {
      const G4ParticleDefinition* part2 = cal.FindParticle(partc[ii]);

      for (G4int i = 0; i < 7; ++i)
      {
        mat = mman->FindOrBuildMaterial(nmk[i]);
        G4cout << "  Particle  " << partc[ii] << " in  Material  " << mat->GetName() << G4endl;
        G4double fact3 = gram / (MeV * cm2 * mat->GetDensity());
        std::ifstream* fin = new std::ifstream();
        std::string fname = "stopping/" + partc[ii] + "_" + nmk[i];
        std::string fnamef = fname + ".txt";
        fin->open(fnamef.c_str());
        if (!fin->is_open())
        {
          G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
          exit(1);
        }
        // G4int i1 = histo->addCloud1D(fname);
        for (G4int j = 0; j < 8; ++j)
        {
          fin->getline(line, 200);
        }
        G4int i2 = 0;
        do
        {
          i2++;
          (*fin) >> e >> se >> sn >> st;
          e *= MeV;
          mce = fact3 * cal.ComputeDEDX(e, part2, proc[ii], mat);
          mcn = fact3 * cal.ComputeNuclearDEDX(e, part2, mat);
          mct = mce + mcn;
          G4double diff = 100. * (mct / st - 1.0);
          /*
                G4cout << e/MeV
                       << " NIST:  dedx= " << se
                       << " nuc= " << sn
                       << " tot= " << st
                       << " G4:  dedx= " << mce
                       << " nuc= " << mcn
                       << " tot= " << mct
                       << " diff= " << diff << " %"
                       << G4endl;
          */
          if (ii == 0 && i == 0)
            G4cout << " " << e;
          else
            G4cout << " " << -diff;
          // histo->Fill(i1,e,diff);

        } while (std::fabs(e - 1000.) > MeV);
        G4cout << G4endl;
        G4cout << "========= n= " << i2
               << " ===========================================================" << G4endl;
        fin->close();
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::MuonTest()
{
  G4NistManager* mman = G4NistManager::Instance();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  const G4ParticleDefinition* part = cal.FindParticle("mu+");

  G4bool imu = true;
  if (imu)
  {
    G4cout << "====================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "====================================================================" << G4endl;

    G4String nmk[4] = {"G4_WATER", "G4_Al", "G4_Fe", "G4_He"};

    G4double energy1[25] = {0.0001, 0.0002, 0.0003, 0.0005, 0.0007, 0.001, 0.002, 0.003, 0.005,
                            0.007,  0.01,   0.02,   0.03,   0.05,   0.07,  0.1,   0.2,   0.3,
                            0.5,    0.7,    1.,     2.,     3.,     5.,    7.};

    G4double energy[43] = {
      10.,       14.,       20.,       30.,      40.,       80.,       100.,      140.,
      200.,      300.,      400.,      800.,     1000.,     1400.,     2000.,     3000.,
      4000.,     8000.,     10000.,    14000.,   20000.,    30000.,    40000.,    80000.,
      100000.,   140000.,   200000.,   300000.,  400000.,   800000.,   1000000.,  1400000.,
      2000000.,  3000000.,  4000000.,  8000000., 10000000., 14000000., 20000000., 30000000.,
      40000000., 80000000., 100000000.};
    G4double dedx[4][43] = {
      {7.965,  6.213,  4.852,  3.764,  3.214,  2.413,  2.270,   2.116,   2.026,   1.992,  1.999,
       2.075,  2.109,  2.166,  2.229,  2.300,  2.351,  2.470,   2.507,   2.564,   2.625,  2.701,
       2.760,  2.942,  3.020,  3.166,  3.376,  3.709,  4.039,   5.352,   6.014,   7.330,  9.327,
       12.654, 16.018, 29.603, 36.462, 50.192, 70.964, 105.740, 140.768, 282.138, 353.358},
      {6.188,  4.849, 3.802,  2.961, 2.533,  1.908,   1.798,   1.688,   1.630,   1.616,  1.630,
       1.711,  1.745, 1.799,  1.858, 1.925,  1.971,   2.082,   2.117,   2.172,   2.233,  2.312,
       2.377,  2.594, 2.694,  2.885, 3.167,  3.628,   4.093,   5.964,   6.914,   8.805,  11.628,
       16.474, 21.32, 40.865, 50.72, 70.429, 100.206, 149.938, 199.946, 401.196, 502.351},
      {5.494,  4.321,  3.399,  2.654,  2.274,   1.717,   1.616,   1.516,   1.463,   1.453,  1.467,
       1.548,  1.582,  1.637,  1.697,  1.767,   1.816,   1.936,   1.975,   2.039,   2.113,  2.214,
       2.303,  2.623,  2.777,  3.082,  3.543,   4.304,   5.079,   8.221,   9.820,   13.013, 17.877,
       25.974, 34.162, 67.147, 83.761, 116.947, 167.024, 250.537, 334.408, 671.133, 840.063},
      {7.709, 5.998, 4.673,  3.616, 3.083,  2.305,  2.165,  2.026,  1.954,   1.939,  1.961,
       2.082, 2.134, 2.219,  2.315, 2.429,  2.511,  2.712,  2.776,  2.874,   2.972,  3.062,
       3.117, 3.250, 3.299,  3.386, 3.501,  3.679,  3.848,  4.506,  4.832,   5.481,  6.461,
       8.097, 9.753, 16.472, 19.88, 26.734, 37.147, 54.704, 72.476, 144.871, 181.606}};
    G4double dedxn[4][43] = {
      {0.,    0.,    0.,    0.,    0.,    0.,     0.,    0.,     0.,     0.,    0.,
       0.,    0.,    0.001, 0.001, 0.001, 0.002,  0.004, 0.005,  0.007,  0.009, 0.013,
       0.018, 0.034, 0.042, 0.059, 0.084, 0.125,  0.167, 0.337,  0.423,  0.601, 0.870,
       1.332, 1.803, 3.763, 4.773, 6.854, 10.051, 15.6,  21.296, 45.199, 57.59},
      {0.,    0.,    0.,    0.,    0.,    0.,    0.,     0.,     0.,     0.,    0.,
       0.,    0.,    0.001, 0.001, 0.001, 0.002, 0.004,  0.005,  0.006,  0.009, 0.013,
       0.017, 0.033, 0.040, 0.056, 0.080, 0.120, 0.160,  0.323,  0.405,  0.575, 0.832,
       1.274, 1.723, 3.59,  4.551, 6.528, 9.562, 14.821, 20.213, 42.793, 54.48},
      {0.,    0.,    0.,    0.,    0.,    0.,    0.,     0.,     0.,     0.,    0.,
       0.,    0.,    0.001, 0.001, 0.001, 0.002, 0.003,  0.004,  0.006,  0.008, 0.012,
       0.016, 0.031, 0.038, 0.054, 0.076, 0.114, 0.152,  0.307,  0.386,  0.547, 0.791,
       1.211, 1.637, 3.406, 4.316, 6.184, 9.05,  14.009, 19.089, 40.329, 51.31},
      {0.,    0.,    0.,    0.,    0.,    0.,    0.,     0.,    0.,     0.,    0.,
       0.,    0.001, 0.001, 0.001, 0.002, 0.002, 0.004,  0.005, 0.007,  0.010, 0.014,
       0.019, 0.036, 0.045, 0.062, 0.088, 0.132, 0.175,  0.354, 0.445,  0.632, 0.916,
       1.404, 1.902, 3.979, 5.05,  7.26,  10.66, 16.578, 22.66, 48.256, 61.55}};

    G4cout << "###  Energy (MeV) ### n= 43" << G4endl;
    G4int i, j;
    for (i = 0; i < 43; ++i)
    {
      G4cout << energy[i] << " ";
    }
    G4cout << G4endl;
    G4cout << G4endl;

    for (i = 0; i < 4; ++i)
    {
      const G4Material* mat = mman->FindOrBuildMaterial(nmk[i]);
      G4double fact4 = gram / (MeV * cm2 * mat->GetDensity());
      G4cout << "###  Material ### " << mat->GetName() << " Data" << G4endl;
      for (j = 0; j < 43; ++j)
      {
        dedx[i][j] -= dedxn[i][j];
        dedxn[i][j] = fact4 * (cal.ComputeDEDX(energy[j], part, "muIoni", mat));
        dedxn[i][j] += fact4 * (cal.ComputeDEDX(energy[j], part, "muBrems", mat));
        dedxn[i][j] += fact4 * (cal.ComputeDEDX(energy[j], part, "muPairProd", mat));
        G4cout << dedx[i][j] << " ";
      }
      G4cout << G4endl;
      G4cout << "### Geant4 " << G4endl;
      for (j = 0; j < 43; ++j)
      {
        G4cout << dedxn[i][j] << " ";
      }
      G4cout << G4endl;
      G4cout << "### 1 - Geant4/Data (%)" << G4endl;
      for (j = 0; j < 43; ++j)
      {
        G4cout << (1.0 - dedxn[i][j] / dedx[i][j]) * 100 << " ";
      }
      G4cout << G4endl;
      G4cout << "### Geant4 low energy" << G4endl;
      for (j = 0; j < 25; ++j)
      {
        G4cout << fact4 * (cal.ComputeDEDX(energy1[j], part, "muIoni", mat)) << " ";
      }
      G4cout << G4endl;
    }
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::ElectronTest()
{
  G4NistManager* mman = G4NistManager::Instance();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  // const G4ParticleDefinition* part = cal.FindParticle("e-");
  const G4ParticleDefinition* part = cal.FindParticle("e+");

  G4bool imu = true;
  if (imu)
  {
    G4cout << "====================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "====================================================================" << G4endl;

    G4String nmk[3] = {"G4_WATER", "G4_Si", "G4_W"};
    const G4int nmax = 30;

    G4double e = MeV;

    G4cout << "###  Electron test for Energy (MeV) = " << e << G4endl;
    G4int i, j;
    G4double cs;
    G4double cut[nmax];

    for (j = 0; j < nmax; ++j)
    {
      cut[j] = std::pow(10.0, -3.0 + 0.1 * G4double(j));
      G4cout << cut[j] << " ";
    }
    G4cout << G4endl;
    for (i = 0; i < 3; ++i)
    {
      const G4Material* mat = mman->FindOrBuildMaterial(nmk[i]);
      //      G4double fact = gram/(barn*cm3*mat->GetDensity());
      G4double fact = 1.0;
      G4cout << "###  Material ### " << mat->GetName() << "   Cross Sections: " << G4endl;
      for (j = 0; j < 30; ++j)
      {
        cs = fact * (cal.ComputeCrossSectionPerVolume(e, part, "eIoni", mat, cut[j]));
        G4cout << cs << " ";
      }
      G4cout << G4endl;
    }
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::CountProcess(const G4String& name)
{
  //  G4cout << "### " << name << G4endl;
  if (name == "AnnihiToMuPair")
    ++n_mumu;
  else if (name == "ee2hadr")
    ++n_pipi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::IonLSData()
{
  G4cout << "### Lindhard-Sorensen data " << G4endl;
  const G4double deltaLvsGammaZ1[25][2] = {
    {0.01436024, 0.00263658},    {0.02393100, 0.00430234},    {0.03987898, 0.00388275},
    {0.06645603, 0.00438721},    {0.11075155, 0.00814783},    {0.18456661, 0.01034251},
    {0.30755487, 0.00814138},    {0.51253753, 0.01030063},    {0.85412787, 0.01170064},
    {1.42336359, 0.01254833},    {2.37191935, 0.01224110},    {3.95260172, 0.01181356},
    {6.58667436, 0.01138603},    {10.97586072, 0.00956405},   {18.29031040, 0.00905230},
    {30.47895289, 0.00808340},   {50.78685427, 0.00351957},   {84.62338056, -0.00262961},
    {140.96650648, -0.02350989}, {234.63049810, -0.09073405}, {362.30586618, -0.20764796},
    {490.19352099, -0.34578343}, {610.80514427, -0.48198184}, {736.40234574, -0.61997669},
    {891.36657041, -0.76914050}};
  const G4double deltaLvsGammaZ10[30][2] = {
    {0.01202207, -0.19854843},   {0.02099562, -0.10489606},   {0.03782180, -0.04106323},
    {0.06119283, -0.00364427},   {0.10202443, 0.02455497},    {0.17008341, 0.04679504},
    {0.27950914, 0.06870928},    {0.44431066, 0.08325708},    {0.79232428, 0.09783465},
    {1.39350789, 0.10733428},    {2.14593637, 0.11410545},    {3.68714589, 0.11690894},
    {6.77755927, 0.11831730},    {11.87526209, 0.11566406},   {18.09865075, 0.11171429},
    {32.07701302, 0.10292498},   {55.94968126, 0.06616626},   {88.00460405, -0.00129972},
    {120.76998331, -0.09534659}, {166.90277002, -0.25372356}, {211.37558468, -0.40267391},
    {235.70203465, -0.47316627}, {281.20760851, -0.62281767}, {349.72829144, -0.80641512},
    {424.82121389, -0.98519941}, {481.67794615, -1.10489968}, {564.03352570, -1.25820702},
    {657.79581552, -1.41145176}, {766.09357076, -1.55842034}, {862.85817722, -1.67966419}};
  const G4double deltaLvsGammaZ18[29][2] = {
    {0.01368051, -0.45942297},   {0.02194677, -0.30305129},   {0.03523195, -0.17703268},
    {0.06383301, -0.05847129},   {0.10542129, 0.01356253},    {0.16617781, 0.06283975},
    {0.28469867, 0.11249538},    {0.47813759, 0.15056925},    {0.80899004, 0.17918734},
    {1.35684705, 0.19833635},    {2.23038168, 0.20957980},    {3.80405513, 0.21624204},
    {6.40220186, 0.21847218},    {10.62228614, 0.21651724},   {18.09393997, 0.20807300},
    {30.79993967, 0.18592100},   {54.00295889, 0.12553152},   {79.58737186, 0.03627608},
    {113.70053738, -0.12092116}, {148.66143472, -0.28060710}, {184.17154674, -0.44533159},
    {222.65574495, -0.59201451}, {272.71832992, -0.77410848}, {344.08630260, -0.99006350},
    {428.23216159, -1.20446930}, {508.47772213, -1.38216281}, {602.98208369, -1.54839338},
    {718.39739422, -1.72602869}, {832.24156046, -1.86939536}};
  const G4double deltaLvsGammaZ36[32][2] = {
    {0.01317123, -1.01728294},   {0.01871739, -0.84641558},   {0.02609706, -0.68307609},
    {0.03708605, -0.51747388},   {0.05342089, -0.35515102},   {0.07828659, -0.19893214},
    {0.12999426, -0.01770073},   {0.17880707, 0.07510724},    {0.31565813, 0.21578290},
    {0.58986626, 0.32570520},    {1.10108133, 0.39856874},    {2.16513655, 0.44211026},
    {3.58935883, 0.45897830},    {6.09099454, 0.46535519},    {10.22788938, 0.45976556},
    {17.13991851, 0.44155645},   {28.90121362, 0.39469247},   {46.53755751, 0.29981448},
    {67.01025933, 0.16255136},   {87.38162673, 0.00875576},   {105.67704531, -0.11764334},
    {128.49783242, -0.28071123}, {152.88880120, -0.43145426}, {178.59796732, -0.57611267},
    {207.69630639, -0.71029649}, {243.53276268, -0.86098010}, {285.86419173, -1.02264871},
    {332.82891953, -1.17165546}, {379.90978773, -1.31707932}, {430.16875099, -1.43965337},
    {496.33689703, -1.58970170}, {572.61603408, -1.72702584}};
  const G4double deltaLvsGammaZ54[37][2] = {
    {0.01276025, -1.41505623},   {0.01728549, -1.24674624},   {0.02344324, -1.07236209},
    {0.02858124, -0.96322100},   {0.03912607, -0.79022042},   {0.05046641, -0.64777877},
    {0.06414268, -0.50890215},   {0.08483619, -0.35569446},   {0.11481326, -0.19009809},
    {0.16289135, -0.00657670},   {0.22311738, 0.13152180},    {0.30210259, 0.27841392},
    {0.44004810, 0.39999659},    {0.61976323, 0.50376402},    {0.94408685, 0.59164981},
    {1.67037489, 0.67563035},    {2.90467510, 0.71924390},    {5.05499305, 0.73547666},
    {8.63202269, 0.72932379},    {14.71959008, 0.70026652},   {25.91934362, 0.61908457},
    {38.45595993, 0.50573543},   {60.49972575, 0.27941988},   {81.97433171, 0.04688405},
    {97.67836112, -0.10468711},  {117.24947728, -0.26045982}, {140.73897166, -0.43328414},
    {162.53425060, -0.58464879}, {192.65452766, -0.74807022}, {223.05738178, -0.90047556},
    {258.74289704, -1.03810962}, {301.28672417, -1.20764386}, {347.59370087, -1.35140361},
    {398.70997540, -1.49620496}, {463.60490108, -1.64987506}, {535.55668182, -1.79760312},
    {616.83830284, -1.93613689}};
  const G4double deltaLvsGammaZ66[40][2] = {
    {0.01278948, -1.60033866},   {0.01696744, -1.44707852},   {0.02171071, -1.30563913},
    {0.02695706, -1.17696670},   {0.03382574, -1.03300579},   {0.04417117, -0.88039869},
    {0.05364450, -0.75698838},   {0.07318593, -0.55259761},   {0.10057273, -0.36149869},
    {0.12636166, -0.21271331},   {0.15292897, -0.09533867},   {0.19362111, 0.03858734},
    {0.25535362, 0.19540667},    {0.33987529, 0.35602813},    {0.46820920, 0.50424383},
    {0.63416166, 0.62213350},    {1.03333257, 0.76210888},    {1.76089900, 0.86478668},
    {2.93708994, 0.91623537},    {4.89597547, 0.93383047},    {8.15764123, 0.92590445},
    {13.58421624, 0.88482970},   {22.71621649, 0.79212245},   {33.11638874, 0.66896213},
    {46.20183130, 0.48792338},   {57.92187457, 0.33817138},   {69.94044899, 0.17879205},
    {86.55051725, 0.00833958},   {99.35801832, -0.12608069},  {116.80654186, -0.27866177},
    {137.22395051, -0.44187629}, {158.41458371, -0.58512484}, {184.54009805, -0.74063532},
    {219.30422497, -0.92715276}, {247.95206761, -1.05205101}, {285.83290353, -1.20209066},
    {326.42688088, -1.34627894}, {369.56340312, -1.47179257}, {423.53799302, -1.62005403},
    {493.60707268, -1.78017235}};
  const G4double deltaLvsGammaZ79[44][2] = {
    {0.01268307, -1.78313045},   {0.01740260, -1.60578389},   {0.02219116, -1.46471990},
    {0.02780693, -1.32669752},   {0.03158530, -1.24217882},   {0.04010501, -1.08712997},
    {0.05258857, -0.90780164},   {0.06314499, -0.78460243},   {0.08089141, -0.61415961},
    {0.10304992, -0.45174288},   {0.12618317, -0.29975919},   {0.15467069, -0.14472039},
    {0.19629153, 0.01668428},    {0.24509496, 0.15178213},    {0.30799735, 0.31116338},
    {0.37296559, 0.44502142},    {0.46577787, 0.57229678},    {0.64515101, 0.74822672},
    {0.92427957, 0.88603601},    {1.39779152, 1.01603743},    {2.50189979, 1.12021265},
    {4.20038281, 1.15989655},    {7.01805876, 1.15532517},    {11.70177258, 1.10296244},
    {18.76064170, 0.99766019},   {28.00237663, 0.84963148},   {37.00547378, 0.69697499},
    {46.35800095, 0.54791280},   {54.72655174, 0.40466052},   {64.79025672, 0.25026659},
    {79.41569629, 0.06889349},   {93.93078003, -0.08265986},  {109.86892676, -0.23545368},
    {124.97940662, -0.36789983}, {143.39388393, -0.51717617}, {166.06494807, -0.66532605},
    {196.15070865, -0.83942240}, {232.12898874, -1.02470616}, {272.36500491, -1.21088628},
    {311.34753127, -1.35055706}, {356.10832067, -1.49364874}, {414.41212072, -1.65147816},
    {470.02002347, -1.78614884}, {537.21274879, -1.91854572}};
  const G4double deltaLvsGammaZ92[47][2] = {
    {0.01257207, -1.93718925},   {0.01672941, -1.77917182},   {0.02372074, -1.56777495},
    {0.03213168, -1.37139905},   {0.04750962, -1.10892169},   {0.06028701, -0.93719455},
    {0.07349109, -0.80036056},   {0.09113865, -0.64697183},   {0.11230138, -0.48771053},
    {0.13619570, -0.34099724},   {0.16389200, -0.18691223},   {0.20175431, -0.01321294},
    {0.23748156, 0.11494242},    {0.28184184, 0.25443931},    {0.33868154, 0.39103330},
    {0.39562472, 0.51762970},    {0.50100885, 0.68100788},    {0.60606278, 0.80478269},
    {0.83546750, 0.99198662},    {1.14505824, 1.13497122},    {1.79211375, 1.27660323},
    {2.97566787, 1.36248726},    {4.93497518, 1.38108748},    {8.21829270, 1.34347154},
    {13.55932346, 1.24575615},   {20.53379599, 1.10743160},   {27.61135733, 0.94666942},
    {34.04776974, 0.80830552},   {41.70403225, 0.66333841},   {49.78754043, 0.50048655},
    {57.80657230, 0.36381615},   {69.11721857, 0.19831657},   {84.96128702, 0.01534661},
    {97.72328290, -0.10995011},  {114.08996019, -0.25808770}, {133.57986879, -0.42360971},
    {156.39785083, -0.57034020}, {183.61019708, -0.73909299}, {206.05064588, -0.84776905},
    {238.38416721, -1.00494558}, {288.50208323, -1.21819213}, {330.12705018, -1.35797844},
    {379.41268288, -1.51019318}, {443.60060783, -1.66243610}, {509.19591284, -1.78976451},
    {579.31510336, -1.92750315}, {617.45512259, -1.99783924}};
  const G4double deltaLvsGammaZ109[48][2] = {
    {0.01769748, -1.91890296},   {0.02254740, -1.78135199},   {0.02855958, -1.63624759},
    {0.04101947, -1.39497274},   {0.05292755, -1.21866773},   {0.06962793, -1.01105164},
    {0.09010080, -0.81411925},   {0.11078015, -0.65010733},   {0.13107173, -0.50447083},
    {0.15850521, -0.33805424},   {0.19012462, -0.17104479},   {0.22625240, -0.01068787},
    {0.26154433, 0.12631697},    {0.31103732, 0.29291130},    {0.35748406, 0.41821815},
    {0.40717130, 0.53571422},    {0.49052508, 0.67629924},    {0.56510305, 0.79670822},
    {0.73262997, 0.99841072},    {1.00906350, 1.20459504},    {1.30518854, 1.34548199},
    {1.99302363, 1.49775889},    {3.44296079, 1.59397579},    {5.87444180, 1.59312643},
    {10.22882791, 1.50603219},   {14.32399265, 1.40888328},   {20.16579053, 1.25099176},
    {25.04298752, 1.11531198},   {31.51395886, 0.96192129},   {38.47728245, 0.79629636},
    {45.43204397, 0.64551404},   {52.54860121, 0.49877992},   {58.87411258, 0.36438190},
    {67.46869346, 0.22757635},   {82.25395110, 0.02844372},   {97.33796272, -0.12901102},
    {117.46154819, -0.33402683}, {138.87573153, -0.51534248}, {159.59856501, -0.67348002},
    {180.91445587, -0.81434060}, {215.31289549, -1.01644066}, {259.19871817, -1.19178342},
    {296.78537436, -1.34208086}, {340.94200720, -1.49385375}, {393.95922992, -1.64935359},
    {448.00363561, -1.78315964}, {525.77066748, -1.94113109}, {563.52771793, -2.00098092}};
  //
  std::size_t i;
  G4PhysicsFreeVector* v[9];
  v[0] = new G4PhysicsFreeVector(25, false);
  for (i = 0; i < 25; ++i)
  {
    v[0]->PutValues(i, deltaLvsGammaZ1[i][0], deltaLvsGammaZ1[i][1]);
  }
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double fact = 1. / g4calc->Z23(10);
  v[1] = new G4PhysicsFreeVector(30, true);
  for (i = 0; i < 30; ++i)
  {
    v[1]->PutValues(i, deltaLvsGammaZ10[i][0], deltaLvsGammaZ10[i][1] * fact);
  }
  v[1]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(18);
  v[2] = new G4PhysicsFreeVector(29, true);
  for (i = 0; i < 29; ++i)
  {
    v[2]->PutValues(i, deltaLvsGammaZ18[i][0], deltaLvsGammaZ18[i][1] * fact);
  }
  v[2]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(36);
  v[3] = new G4PhysicsFreeVector(32, true);
  for (i = 0; i < 32; ++i)
  {
    v[3]->PutValues(i, deltaLvsGammaZ36[i][0], deltaLvsGammaZ36[i][1] * fact);
  }
  v[3]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(54);
  v[4] = new G4PhysicsFreeVector(37, true);
  for (i = 0; i < 37; ++i)
  {
    v[4]->PutValues(i, deltaLvsGammaZ54[i][0], deltaLvsGammaZ54[i][1] * fact);
  }
  v[4]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(66);
  v[5] = new G4PhysicsFreeVector(40, true);
  for (i = 0; i < 40; ++i)
  {
    v[5]->PutValues(i, deltaLvsGammaZ66[i][0], deltaLvsGammaZ66[i][1] * fact);
  }
  v[5]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(79);
  v[6] = new G4PhysicsFreeVector(44, true);
  for (i = 0; i < 44; ++i)
  {
    v[6]->PutValues(i, deltaLvsGammaZ79[i][0], deltaLvsGammaZ79[i][1] * fact);
  }
  v[6]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(92);
  v[7] = new G4PhysicsFreeVector(47, true);
  for (i = 0; i < 47; ++i)
  {
    v[7]->PutValues(i, deltaLvsGammaZ92[i][0], deltaLvsGammaZ92[i][1] * fact);
  }
  v[7]->FillSecondDerivatives();

  fact = 1. / g4calc->Z23(109);
  v[8] = new G4PhysicsFreeVector(48, true);
  for (i = 0; i < 48; ++i)
  {
    v[8]->PutValues(i, deltaLvsGammaZ109[i][0], deltaLvsGammaZ109[i][1] * fact);
  }
  v[8]->FillSecondDerivatives();

  G4cout << "  {" << std::setprecision(8) << G4endl;
  G4double xmin = 0.02;
  G4double xmax = std::pow(10., 2.5);
  fact = std::exp(std::log(xmax / xmin) / 40.);
  G4double x = 0.0;
  for (G4int j = 0; j < 9; ++j)
  {
    x = xmin;
    for (i = 0; i <= 40; ++i)
    {
      G4double y = v[j]->Value(x);
      if (0 == i)
      {
        G4cout << "  {" << y << "," << G4endl;
      }
      else if (40 == i)
      {
        G4cout << std::setw(12) << y << "}," << G4endl;
      }
      else
      {
        G4cout << std::setw(12) << y << ", ";
        if (i / 5 * 5 == i)
        {
          G4cout << " // " << i - 5 << "-" << i << G4endl;
          G4cout << "  ";
        }
      }
      x *= fact;
    }
  }
  G4cout << "\n   xmin= " << xmin << "  xmax= " << xmax << "  x= " << x << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
