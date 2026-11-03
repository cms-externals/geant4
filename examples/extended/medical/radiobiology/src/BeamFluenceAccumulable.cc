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
/// \file radiobiology/src/BeamFluenceAccumulable.cc
/// \brief Implementation of the RadioBio::BeamFluenceAccumulable class

#include "BeamFluenceAccumulable.hh"

#include "G4ParticleDefinition.hh"

#include "Hit.hh"
#include "Manager.hh"
#include "VoxelizedSensitiveDetector.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamFluenceAccumulable::BeamFluenceAccumulable() : VRadiobiologicalAccumulable("BeamFluence") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluenceAccumulable::Merge(const G4VAccumulable& rhs)
{
  if (GetVerboseLevel() > 1) G4cout << "BeamFluenceAccumulable::Merge()" << G4endl;

  const auto& other = dynamic_cast<const BeamFluenceAccumulable&>(rhs);

  fTotalParticleCount += other.GetTotalParticleCount();
  fPrimaryParticleCount += other.GetPrimaryParticleCount();
  fSecondaryParticleCount += other.GetSecondaryParticleCount();

  // Different worker threads may observe different particle species.
  for (const auto& entry : other.GetParticleCountBySpecies()) {
    const G4String& particleName = entry.first;
    const count_array_type& otherCounts = entry.second;

    auto found = fParticleCountBySpecies.find(particleName);
    if (found == fParticleCountBySpecies.end()) {
      fParticleCountBySpecies.emplace(particleName, otherCounts);
    }
    else {
      found->second += otherCounts;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluenceAccumulable::Reset()
{
  if (GetVerboseLevel() > 0) G4cout << "BeamFluenceAccumulable::Reset()" << G4endl;
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluenceAccumulable::Accumulate(Hit* hit)
{
  if (!hit->IsVoxelEntry()) return;

  auto* voxelizedDetector = VoxelizedSensitiveDetector::GetInstance();
  const G4int voxel =
    voxelizedDetector->GetThisVoxelNumber(hit->GetXindex(), hit->GetYindex(), hit->GetZindex());

  if (voxel < 0 || voxel >= voxelizedDetector->GetTotalVoxelNumber()) {
    G4Exception("BeamFluenceAccumulable::Accumulate", "InvalidVoxelIndex", JustWarning,
                "Ignoring a particle entry with an invalid voxel index.");
    return;
  }

  ++fTotalParticleCount[voxel];

  if (hit->GetParentID() == 0) {
    ++fPrimaryParticleCount[voxel];
  }
  else {
    ++fSecondaryParticleCount[voxel];
  }

  const G4ParticleDefinition* particle = hit->GetPartType();
  const G4String particleName = particle ? particle->GetParticleName() : G4String("unknown");

  auto found = fParticleCountBySpecies.find(particleName);
  if (found == fParticleCountBySpecies.end()) {
    const auto nVoxels = voxelizedDetector->GetTotalVoxelNumber();
    found =
      fParticleCountBySpecies.emplace(particleName, count_array_type(G4long{0}, nVoxels)).first;
  }

  ++found->second[voxel];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int BeamFluenceAccumulable::GetVerboseLevel() const
{
  return Manager::GetInstance()->GetQuantity("BeamFluence")->GetVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluenceAccumulable::Initialize()
{
  const auto nVoxels = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();

  fTotalParticleCount = count_array_type(G4long{0}, nVoxels);
  fPrimaryParticleCount = count_array_type(G4long{0}, nVoxels);
  fSecondaryParticleCount = count_array_type(G4long{0}, nVoxels);
  fParticleCountBySpecies.clear();
}

}  // namespace RadioBio
