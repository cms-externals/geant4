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
// V.Ivanchenko 2013/10/19
// step limiter and killer for e+,e- and other charged particles
//
#include "ElectronLimiter.hh"

#include "G4DummyModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TransportationProcessType.hh"

ElectronLimiter::ElectronLimiter(const G4ParticleDefinition* part)
  : G4VEmProcess("eLimiter", fGeneral),
    particle(part),
    regionEcal(nullptr),
    regionHcal(nullptr),
    limitEcal(DBL_MAX),
    factEcal(1.0),
    rmsEcal(0.0),
    limitHcal(DBL_MAX),
    factHcal(1.0),
    rmsHcal(0.0),
    insideEcal(false)
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(STEP_LIMITER));
}

ElectronLimiter::~ElectronLimiter() {}

void ElectronLimiter::InitialiseProcess(const G4ParticleDefinition*)
{
  AddEmModel(0, new G4DummyModel());
  regionEcal = G4RegionStore::GetInstance()->GetRegion("EcalRegion");
  regionHcal = G4RegionStore::GetInstance()->GetRegion("HcalRegion");
}

G4double ElectronLimiter::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double,
                                                               G4ForceCondition* cond)
{
  *cond = NotForced;
  G4double limit = DBL_MAX;

  G4double kinEnergy = aTrack.GetKineticEnergy();
  const G4Region* reg = aTrack.GetVolume()->GetLogicalVolume()->GetRegion();
  if (reg == regionEcal && kinEnergy < limitEcal)
  {
    insideEcal = true;
    limit = 0.0;
  }
  else if (reg == regionHcal && kinEnergy < limitHcal)
  {
    insideEcal = false;
    limit = 0.0;
  }
  return limit;
}

G4VParticleChange* ElectronLimiter::PostStepDoIt(const G4Track& track, const G4Step&)
{
  fParticleChange.Initialize(track);
  G4double ekin = track.GetKineticEnergy();
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  if (insideEcal)
  {
    ekin *= factEcal * G4RandGauss::shoot(rndmEngine, 1.0, rmsEcal);
  }
  else
  {
    ekin *= factHcal * G4RandGauss::shoot(rndmEngine, 1.0, rmsHcal);
  }
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  fParticleChange.ProposeLocalEnergyDeposit(std::max(ekin, 0.0));
  fParticleChange.SetProposedKineticEnergy(0.0);
  return &fParticleChange;
}

G4bool ElectronLimiter::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

void ElectronLimiter::StartTracking(G4Track*) {}
