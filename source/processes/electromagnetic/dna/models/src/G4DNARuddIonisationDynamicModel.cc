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
// Modified by Z. Francis, S. Incerti to handle HZE 
// && inverse rudd function sampling 26-10-2010
//
// Rewitten by V.Ivanchenko 21.05.2023
//

#include "G4EmCorrections.hh"
#include "G4DNARuddIonisationDynamicModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNACrossSectionDataSet* G4DNARuddIonisationDynamicModel::xsdata = nullptr;
const std::vector<G4double>* G4DNARuddIonisationDynamicModel::fpWaterDensity = nullptr;

namespace
{
  const G4double scaleFactor = CLHEP::m*CLHEP::m;
  const G4double tolerance = 1*CLHEP::eV;
  const G4double Ry = 13.6*CLHEP::eV;

  // Following values provided by M. Dingfelder (priv. comm)
  const G4double Bj[5] = {12.60*CLHEP::eV, 14.70*CLHEP::eV, 18.40*CLHEP::eV,
                          32.20*CLHEP::eV, 539*CLHEP::eV};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationDynamicModel::G4DNARuddIonisationDynamicModel(const G4ParticleDefinition*,
                                                                 const G4String& nam)
  : G4VEmModel(nam)
{
  fEmCorrections = G4LossTableManager::Instance()->EmCorrections();
  fGpow = G4Pow::GetInstance();
  fLowestEnergy = 100*CLHEP::eV;

  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  // Define default angular generator
  SetAngularDistribution(new G4DNARuddAngle());

  if (nullptr == xsdata) {
    isFirst = true;
    LoadData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationDynamicModel::~G4DNARuddIonisationDynamicModel()
{  
  if (isFirst) { delete xsdata; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::LoadData()
{
  // initialisation of static data once
  G4String filename = "dna/sigma_ionisation_p_rudd";
  xsdata = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, CLHEP::eV, scaleFactor);
  xsdata->LoadData(filename);

  // to avoid possible threading problem fill this vector only once
  auto water = G4NistManager::Instance()->FindMaterial("G4_WATER");
  fpWaterDensity =
    G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(water);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::Initialise(const G4ParticleDefinition* p,
						 const G4DataVector&)
{
  if (p != fParticle) { SetParticle(p); }

  // particle change object may be externally set
  if (nullptr == fParticleChangeForGamma) {
    fParticleChangeForGamma = GetParticleChangeForGamma();
  }

  // the same definition of generic ion as in G4VEmProcess class
  if (p->GetParticleType() == "nucleus" && p->GetParticleSubType() == "generic") {
    G4String pname = p->GetParticleName();
    if (pname != "deuteron" && pname != "triton" &&
	pname != "He3" && pname != "alpha" && pname != "alpha+" &&
	pname != "helium" && pname != "hydrogen") {
      isIon = true;
    }
  }

  // initialisation once in each thread
  if (!isInitialised) {
    isInitialised = true;

    // defined stationary mode
    statCode = G4EmParameters::Instance()->DNAStationary();

    // initialise atomic de-excitation
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

    if (verbose > 0) {
      G4cout << "### G4DNARuddIonisationDynamicModel::Initialise(..) "
	     << fParticle->GetParticleName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::SetParticle(const G4ParticleDefinition* p)
{
  fParticle = p;
  fMass = p->GetPDGMass();
  fMassRate = CLHEP::proton_mass_c2/fMass; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationDynamicModel::StartTracking(G4Track* track)
{
  fDynParticle = track->GetDynamicParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4DNARuddIonisationDynamicModel::CrossSectionPerVolume(const G4Material* material,
						       const G4ParticleDefinition* part,
						       G4double kinE,
						       G4double, G4double)
{
  // check if model is applicable for given material
  G4double density = (material->GetIndex() < fpWaterDensity->size())
    ? (*fpWaterDensity)[material->GetIndex()] : 0.0;
  if (0.0 == density) { return 0.0; }

  const G4double q = fDynParticle->GetCharge()*inveplus;
  if (0.0 == q) { return 0.0; }
  
  // check on kinetic energy (not scaled energy) to stop low-energy ion
  if (kinE < fLowestEnergy) { return DBL_MAX; }

  // ion may be different
  if (fParticle != part) { SetParticle(part); }

  // cross section for scaled energy
  const G4double e = kinE*fMassRate;
  G4double sigma = (e > fLowestEnergy) ? xsdata->FindValue(e)
    : xsdata->FindValue(fLowestEnergy) * e / fLowestEnergy;

  sigma *= q * q * density;

  if (verbose > 1) {
    G4cout << "G4DNARuddIonisationDynamicModel for " << part->GetParticleName() 
           << " Ekin(keV)=" << kinE/CLHEP::keV 
           << " sigma(cm^2)=" << sigma/CLHEP::cm2 << G4endl;
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNARuddIonisationDynamicModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                    const G4MaterialCutsCouple* couple,
                                                    const G4DynamicParticle* dpart,
                                                    G4double, G4double)
{
  const G4ParticleDefinition* pd = dpart->GetDefinition();
  if (fParticle != pd) { SetParticle(pd); }

  // stop ion with energy below low energy limit
  G4double kinE = dpart->GetKineticEnergy();
  // ion shoud be stopped - check on kinetic energy and not scaled energy
  if (kinE <= fLowestEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopButAlive);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE);
    return;
  }

  const G4double eScaled = kinE*fMassRate;
  G4int shell = SelectShell(eScaled);
  G4double bindingEnergy = (useDNAWaterStructure)
    ? waterStructure.IonisationEnergy(shell) : Bj[shell];

  //Si: additional protection if tcs interpolation method is modified
  if (kinE < bindingEnergy) { return; }
  
  G4double esec = SampleElectronEnergy(eScaled, shell);
  G4double esum = 0.0;

  // sample deexcitation
  // here we assume that H2O electronic levels are the same as Oxygen.
  // this can be considered true with a rough 10% error in energy on K-shell,
  G4int Z = 8;	
  G4ThreeVector deltaDir = 
    GetAngularDistribution()->SampleDirectionForShell(dpart, esec, Z, shell, couple->GetMaterial());

  // SI: only atomic deexcitation from K shell is considered
  if(fAtomDeexcitation != nullptr && shell == 4) {
    auto as = G4AtomicShellEnumerator(0);
    auto ashell = fAtomDeexcitation->GetAtomicShell(Z, as);
    fAtomDeexcitation->GenerateParticles(fvect, ashell, Z, 0, 0);

    // compute energy sum from de-excitation
    for (auto const & ptr : *fvect) {
      esum += ptr->GetKineticEnergy();
    }
  }
  // check energy balance
  // remaining excitation energy of water molecule
  G4double exc = bindingEnergy - esum;

  // remaining projectile energy
  G4double scatteredEnergy = kinE - bindingEnergy - esec;
  if(scatteredEnergy < -tolerance || exc < -tolerance) {
    G4cout << "G4DNARuddIonisationDynamicModel::SampleSecondaries: "
           << "negative final E(keV)=" << scatteredEnergy/CLHEP::keV << " Ein(keV)="
           << kinE/CLHEP::keV << "  " << pd->GetParticleName()
           << " Edelta(keV)=" << esec/CLHEP::keV << " MeV, Exc(keV)=" << exc/CLHEP::keV
	   << G4endl;
  }

  // projectile
  if (!statCode) {
    fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(exc);
  } else {
    fParticleChangeForGamma->SetProposedKineticEnergy(kinE);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(kinE - scatteredEnergy);
  }

  // delta-electron
  auto  dp = new G4DynamicParticle(G4Electron::Electron(), deltaDir, esec);
  fvect->push_back(dp);

  // create radical
  const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule, shell,
							 theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationDynamicModel::SelectShell(const G4double e)
{
  G4double sum = 0.0;
  G4double xs;
  for (G4int i=0; i<5; ++i) {
    auto ptr = xsdata->GetComponent(i);
    xs = (e > fLowestEnergy) ? ptr->FindValue(e)
      : ptr->FindValue(fLowestEnergy)*e/fLowestEnergy;
    sum += xs;
    fTemp[i] = sum;
  }
  sum *= G4UniformRand();
  for (G4int i=0; i<5; ++i) {
    if (sum <= fTemp[i]) { return i; }
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNARuddIonisationDynamicModel::MaxEnergy(const G4double kine, const G4int shell)
{
  // kinematic limit
  G4double tau = kine/CLHEP::proton_mass_c2;
  G4double gam = 1.0 + tau;
  G4double emax = 2.0*CLHEP::electron_mass_c2*tau*(tau + 2.0);

  // Initialisation of sampling
  G4double A1, B1, C1, D1, E1, A2, B2, C2, D2;
  if (shell == 4) {
    //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
    A1 = 1.25;
    B1 = 0.5;
    C1 = 1.00;
    D1 = 1.00;
    E1 = 3.00;
    A2 = 1.10;
    B2 = 1.30;
    C2 = 1.00;
    D2 = 0.00;
    alphaConst = 0.66;
  } else {
    //Data For Liquid Water from Dingfelder (Protons in Water)
    A1 = 1.02;
    B1 = 82.0;
    C1 = 0.45;
    D1 = -0.80;
    E1 = 0.38;
    A2 = 1.07;
    // Value provided by M. Dingfelder (priv. comm)
    B2 = 11.6;
    C2 = 0.60;
    D2 = 0.04;
    alphaConst = 0.64;
  }
  bEnergy = Bj[shell];
  G4double v2 = 0.25*emax/(bEnergy*gam*gam);
  v = std::sqrt(v2);
  u = Ry/bEnergy;
  wc = 4.*v2 - 2.*v - 0.25*u;

  G4double L1 = (C1 * fGpow->powA(v, D1)) / (1. + E1 * fGpow->powA(v, (D1 + 4.)));
  G4double L2 = C2 * fGpow->powA(v, D2);
  G4double H1 = (A1 * G4Log(1. + v2)) / (v2 + (B1 / v2));
  G4double H2 = (A2 / v2) + (B2 / (v2 * v2));

  F1 = L1 + H1;
  F2 = (L2 * H2) / (L2 + H2);
  return emax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNARuddIonisationDynamicModel::SampleElectronEnergy(const G4double kine,
						      const G4int shell)
{
  // sampling is performed for proton projectile
  G4double emax = MaxEnergy(kine, shell);
  // compute cumulative probability function
  G4double step = 1*CLHEP::eV;
  auto nn = (G4int)(emax/step);
  nn = std::min(std::max(nn, 10), 100);
  step = emax/(G4double)nn;

  // find max probability
  G4double pmax = ProbabilityFunction(kine, 0.0, shell);
  //G4cout << "## E(keV)=" << kine/keV << " emax=" << emax/keV
  //       << " pmax(0)=" << pmax << " shell=" << shell << " nn=" << nn << G4endl;

  G4double e0 = 0.0; // energy with max probability
  // 2 areas after point with max probability
  G4double e1 = emax;
  G4double e2 = emax;
  G4double p1 = 0.0;
  G4double p2 = 0.0;
  const G4double f = 0.25;

  // find max probability
  G4double e = 0.0;
  G4double p = 0.0;
  for (G4int i=0; i<nn; ++i) {
    e += step;
    p = ProbabilityFunction(kine, e, shell);
    if (p > pmax) {
      pmax = p;
      e0 = e;
    } else {
      break;
    }
  }
  // increase step to be more effective
  step *= 2.0;
  // 2-nd area
  for (G4int i=0; i<nn; ++i) {
    e += step;
    if (std::abs(e - emax) < step) {
      e1 = emax;
      break;
    }
    p = ProbabilityFunction(kine, e, shell);
    if (p < f*pmax) {
      p1 = p;
      e1 = e;
      break;
    }
  }
  // 3-d area
  if (e < emax) {
    for (G4int i=0; i<nn; ++i) {
      e += step;
      if (std::abs(e - emax) < step) {
        e2 = emax;
	break;
      }
      p = ProbabilityFunction(kine, e, shell);
      if (p < f*p1) {
	p2 = p;
	e2 = e;
        break;
      }
    }
  }
  pmax *= 1.05;
  // regression method with 3 regions
  G4double s0 = pmax*e1;
  G4double s1 = s0 + p1 * (e2 - e1);
  G4double s2 = s1 + p2 * (emax - e2);
  s0 = (s0 == s1) ? 1.0 : s0 / s2;
  s1 = (s1 == s2) ? 1.0 : s1 / s2;

  //G4cout << "pmax=" << pmax << " e1(keV)=" << e1/keV << " p1=" << p1 << " e2(keV)=" << e2/keV
  //	 << " p2=" << p2 << " s0=" << s0 << " s1=" << s1 << " s2=" << s2 << G4endl;

  // sampling
  G4int count = 0;
  G4double ymax, y, deltae;
  for (G4int i = 0; i<100000; ++i) {
    G4double q = G4UniformRand();
    if (q <= s0) {
      ymax = pmax;
      deltae = e1 * q / s0;
    } else if (q <= s1) {
      ymax = p1;
      deltae = e1 + (e2 - e1) * (q - s0) / (s1 - s0);
    } else {
      ymax = p2;
      deltae = e2 + (emax - e2) * (q - s1) / (1.0 - s1);
    }
    y = ProbabilityFunction(kine, deltae, shell);
    //G4cout << "    " << i << ".  deltae=" << deltae/CLHEP::keV 
    //       << " y=" << y << " ymax=" << ymax << G4endl; 
    if (y > ymax && count < 5) {
      ++count;
      G4cout << "G4DNARuddIonisationDynamicModel::SampleElectronEnergy warning: "
	     << fParticle->GetParticleName() << " Escaled(keV)=" << kine/CLHEP::keV
	     << " Edelta(keV)=" << deltae/CLHEP::keV 
	     << " y=" << y << " ymax=" << ymax << " n=" << i << G4endl; 
    }
    if (ymax * G4UniformRand() <= y) {
      return deltae;
    }
  }
  deltae = std::min(e0 + step, 0.5*emax);
  return deltae;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::ProbabilityFunction(const G4double kine,
                                                              const G4double deltae,
                                                              const G4int shell)
{
  // Shells ids are 0 1 2 3 4 (4 is k shell)
  // !!Attention, "energyTransfer" here is the energy transfered to the electron which means
  //             that the secondary kinetic energy is w = energyTransfer - bindingEnergy
  //
  //   ds            S                F1(nu) + w * F2(nu)
  //  ---- = G(k) * ----     -------------------------------------------
  //   dw            Bj       (1+w)^3 * [1 + exp{alpha * (w - wc) / nu}]
  //
  // w is the secondary electron kinetic Energy in eV
  //
  // All the other parameters can be found in Rudd's Papers
  //
  // M.Eugene Rudd, 1988, User-Friendly model for the energy distribution of
  // electrons from protons or electron collisions. Nucl. Tracks Rad. Meas.Vol 16 N0 2/3 pp 219-218
  //
  G4double w = deltae/bEnergy;
  G4double x = alphaConst*(w - wc)/v;
  G4double y = (x > -15.) ? 1.0 + G4Exp(x) : 1.0;

  G4double res = CorrectionFactor(kine, shell) * (F1 + w*F2) /
    (fGpow->powN((1. + w)/u, 3) * y);

  return std::max(res, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationDynamicModel::ComputeProbabilityFunction(
         const G4ParticleDefinition* p, G4double e, G4double deltae, G4int shell)
{
  if (fParticle != p) { SetParticle(p); }
  MaxEnergy(e, shell);
  return ProbabilityFunction(e, deltae, shell);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4DNARuddIonisationDynamicModel::CorrectionFactor(G4double kine, G4int shell) 
{
  // ZF Shortened
  G4double res = 1.0;
  if (shell < 4) {
    const G4double ln10 = fGpow->logZ(10);
    G4double x = 2.0*((G4Log(kine/CLHEP::eV)/ln10) - 4.2);
    // The following values are provided by M. Dingfelder (priv. comm)
    res = 0.6/(1.0 + G4Exp(x)) + 0.9;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
