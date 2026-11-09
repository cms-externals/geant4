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
// ClassName:    G4EnergyLossForExtrapolator
//
// Description:  This class provide calculation of energy loss, fluctuation,
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko
//
// Modification:
// 08-04-05 Rename Propogator -> Extrapolator (V.Ivanchenko)
// 16-03-06 Add muon tables and fix bug in units (V.Ivanchenko)
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 12-05-06 SEt linLossLimit=0.001 (VI)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossForExtrapolator.hh"

#include "G4AutoLock.hh"
#include "G4Electron.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Positron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace
{
G4Mutex extrMutex = G4MUTEX_INITIALIZER;
}

G4TablesForExtrapolator* G4EnergyLossForExtrapolator::tables = nullptr;

static std::once_flag applyOnce;

G4EnergyLossForExtrapolator::G4EnergyLossForExtrapolator(G4int verb)
  : maxEnergyTransfer(DBL_MAX), verbose(verb)
{
  emin = 1. * CLHEP::MeV;
  emax = 100. * CLHEP::TeV;
  if (nullptr == tables)
  {
    G4AutoLock l(&extrMutex);
    if (nullptr == tables)
    {
      tables = new G4TablesForExtrapolator(verbose, nbins, emin, emax);
    }
    l.unlock();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator::~G4EnergyLossForExtrapolator()
{
  if (isInitializer)
  {
    delete tables;
    tables = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyAfterStep(G4double kinEnergy, G4double stepLength,
                                                      const G4Material* mat,
                                                      const G4ParticleDefinition* part)
{
  G4double kinEnergyFinal = kinEnergy;
  if (SetupKinematics(part, mat, kinEnergy))
  {
    G4double step = TrueStepLength(kinEnergy, stepLength, mat, part);
    G4double r = ComputeRange(kinEnergy, part, mat);
    if (r <= step)
    {
      kinEnergyFinal = 0.0;
    }
    else if (step < linLossLimit * r)
    {
      kinEnergyFinal -= step * ComputeDEDX(kinEnergy, part, mat);
    }
    else
    {
      G4double r1 = r - step;
      kinEnergyFinal = ComputeEnergy(r1, part, mat);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyBeforeStep(G4double kinEnergy, G4double stepLength,
                                                       const G4Material* mat,
                                                       const G4ParticleDefinition* part)
{
  // G4cout << "G4EnergyLossForExtrapolator::EnergyBeforeStep" << G4endl;
  G4double kinEnergyFinal = kinEnergy;

  if (SetupKinematics(part, mat, kinEnergy))
  {
    G4double step = TrueStepLength(kinEnergy, stepLength, mat, part);
    G4double r = ComputeRange(kinEnergy, part, mat);

    if (step < linLossLimit * r)
    {
      kinEnergyFinal += step * ComputeDEDX(kinEnergy, part, mat);
    }
    else
    {
      G4double r1 = r + step;
      kinEnergyFinal = ComputeEnergy(r1, part, mat);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::TrueStepLength(G4double kinEnergy, G4double stepLength,
                                                     const G4Material* mat,
                                                     const G4ParticleDefinition* part)
{
  G4double res = stepLength;
  // G4cout << "## G4EnergyLossForExtrapolator::TrueStepLength L= " << res
  //          <<  "  " << part->GetParticleName() << G4endl;
  if (SetupKinematics(part, mat, kinEnergy))
  {
    if (part == electron || part == positron)
    {
      const G4double x =
        stepLength * ComputeValue(kinEnergy, GetPhysicsTable(fMscElectron), FindIndex(mat));
      // G4cout << " x= " << x << G4endl;
      if (x < 0.2)
      {
        res *= (1.0 + 0.5 * x + x * x / 3.0);
      }
      else if (x < 0.9999)
      {
        res = -G4Log(1.0 - x) * stepLength / x;
      }
      else
      {
        res = ComputeRange(kinEnergy, part, mat);
      }
    }
    else
    {
      res = ComputeTrueStep(mat, part, kinEnergy, stepLength);
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EnergyLossForExtrapolator::SetupKinematics(const G4ParticleDefinition* part,
                                                    const G4Material* mat, G4double kinEnergy)
{
  if (mat->GetNumberOfMaterials() != nmat)
  {
    Initialisation();
  }
  if (nullptr == part || nullptr == mat || kinEnergy < CLHEP::keV)
  {
    return false;
  }

  G4bool redo = false;
  if (part != currentParticle)
  {
    currentParticle = part;
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    charge2 = q * q;
    redo = true;
  }
  if (mat != currentMaterial)
  {
    std::size_t i = mat->GetIndex();
    if (i >= nmat)
    {
      G4cout << "### G4EnergyLossForExtrapolator WARNING: material index i= " << i
             << " above number of materials " << nmat << G4endl;
      return false;
    }
    else
    {
      currentMaterial = mat;
      electronDensity = mat->GetElectronDensity();
      radLength = mat->GetRadlen();
    }
  }
  if (redo || kinEnergy != kineticEnergy)
  {
    kineticEnergy = kinEnergy;
    G4double mass = part->GetPDGMass();
    G4double tau = kinEnergy / mass;

    gam = tau + 1.0;
    bg2 = tau * (tau + 2.0);
    beta2 = bg2 / (gam * gam);
    tmax = kinEnergy;
    if (part == electron)
      tmax *= 0.5;
    else if (part != positron)
    {
      G4double r = CLHEP::electron_mass_c2 / mass;
      tmax = 2.0 * bg2 * CLHEP::electron_mass_c2 / (1.0 + 2.0 * gam * r + r * r);
    }
    tmax = std::min(tmax, maxEnergyTransfer);
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EnergyLossForExtrapolator::FindParticle(const G4String& name)
{
  currentParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
  if (nullptr == currentParticle)
  {
    G4cout << "### G4EnergyLossForExtrapolator WARNING: "
           << "FindParticle() fails to find " << name << G4endl;
  }
  return currentParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4EnergyLossForExtrapolator::FindIndex(const G4Material* mat)
{
  auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t nCouples = theCoupleTable->GetTableSize();
  for (std::size_t i = 0; i < nCouples; ++i)
  {
    if (mat == theCoupleTable->GetMaterialCutsCouple((G4int)i)->GetMaterial())
    {
      return i;
    }
  }
  // the material is not in the geometry
  G4ExceptionDescription ed;
  ed << " FindIndex() for " << mat->GetName()
     << " fail, because there is no such material in the geometry";

  G4Exception("G4EnergyLossForExtrapolator::FindIndex", "em004", FatalException, ed, "");
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeDEDX(G4double ekin, const G4ParticleDefinition* part,
                                                  const G4Material* mat)
{
  if (mat->GetNumberOfMaterials() != nmat)
  {
    Initialisation();
  }
  G4double x = 0.0;
  auto idx = FindIndex(mat);
  if (part == electron)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxElectron), idx);
  }
  else if (part == positron)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxPositron), idx);
  }
  else if (part == muonPlus || part == muonMinus)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxMuon), idx);
  }
  else
  {
    G4double e = ekin * CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    x = ComputeValue(e, GetPhysicsTable(fDedxProton), idx) * q * q;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeRange(G4double ekin, const G4ParticleDefinition* part,
                                                   const G4Material* mat)
{
  if (mat->GetNumberOfMaterials() != nmat)
  {
    Initialisation();
  }
  G4double x = 0.0;
  auto idx = FindIndex(mat);
  if (part == electron)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fRangeElectron), idx);
  }
  else if (part == positron)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fRangePositron), idx);
  }
  else if (part == muonPlus || part == muonMinus)
  {
    x = ComputeValue(ekin, GetPhysicsTable(fRangeMuon), idx);
  }
  else
  {
    G4double massratio = CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double e = ekin * massratio;
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    x = ComputeValue(e, GetPhysicsTable(fRangeProton), idx) / (q * q * massratio);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeEnergy(G4double range,
                                                    const G4ParticleDefinition* part,
                                                    const G4Material* mat)
{
  if (mat->GetNumberOfMaterials() != nmat)
  {
    Initialisation();
  }
  G4double x = 0.0;
  auto idx = FindIndex(mat);
  if (part == electron)
  {
    x = ComputeValue(range, GetPhysicsTable(fInvRangeElectron), idx);
  }
  else if (part == positron)
  {
    x = ComputeValue(range, GetPhysicsTable(fInvRangePositron), idx);
  }
  else if (part == muonPlus || part == muonMinus)
  {
    x = ComputeValue(range, GetPhysicsTable(fInvRangeMuon), idx);
  }
  else
  {
    G4double massratio = CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    G4double r = range * massratio * q * q;
    x = ComputeValue(r, GetPhysicsTable(fInvRangeProton), idx) / massratio;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyDispersion(G4double kinEnergy, G4double stepLength,
                                                       const G4Material* mat,
                                                       const G4ParticleDefinition* part)
{
  G4double sig2 = 0.0;
  if (SetupKinematics(part, mat, kinEnergy))
  {
    G4double step = ComputeTrueStep(mat, part, kinEnergy, stepLength);
    sig2 = (1.0 / beta2 - 0.5) * CLHEP::twopi_mc2_rcl2 * tmax * step * electronDensity * charge2;
  }
  return sig2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::AverageScatteringAngle(G4double kinEnergy,
                                                             G4double stepLength,
                                                             const G4Material* mat,
                                                             const G4ParticleDefinition* part)
{
  G4double theta = 0.0;
  if (SetupKinematics(part, mat, kinEnergy))
  {
    G4double t = stepLength / radLength;
    G4double y = std::max(0.001, t);
    theta = 19.23 * CLHEP::MeV * std::sqrt(charge2 * t) * (1.0 + 0.038 * G4Log(y))
            / (beta2 * gam * part->GetPDGMass());
  }
  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::Initialisation()
{
  // initialisation for the 1st run
  std::call_once(applyOnce, [this]() { isInitializer = true; });

  if (verbose > 0)
  {
    G4cout << "### G4EnergyLossForExtrapolator::Initialisation tables= " << tables
           << "  initializer:" << isInitializer << "  ThreadID=" << G4Threading::G4GetThreadId()
           << G4endl;
  }
  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton = G4Proton::Proton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus = G4MuonMinus::MuonMinus();

  if (G4Material::GetNumberOfMaterials() != nmat)
  {
    // set lock only to create tables
    G4AutoLock l(&extrMutex);
    if (G4Material::GetNumberOfMaterials() != nmat)
    {
      tables->Initialisation();
      if (verbose > 0)
      {
        G4cout << "### G4EnergyLossForExtrapolator::BuildTables for " << nmat
               << " materials Nbins= " << nbins << " Emin(MeV)= " << emin << "  Emax(MeV)= " << emax
               << G4endl;
      }
    }
    l.unlock();
  }
  nmat = G4Material::GetNumberOfMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
