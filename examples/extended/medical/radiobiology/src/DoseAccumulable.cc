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
/// \file radiobiology/src/DoseAccumulable.cc
/// \brief Implementation of the RadioBio::DoseAccumulable class

#include "DoseAccumulable.hh"

#include "G4ParticleDefinition.hh"

#include "Hit.hh"
#include "Manager.hh"
#include "VoxelizedSensitiveDetector.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DoseAccumulable::DoseAccumulable() : VRadiobiologicalAccumulable("Dose") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Merge(const G4VAccumulable& rhs)
{
  if (GetVerboseLevel() > 1) G4cout << "DoseAccumulable::Merge()" << G4endl;

  const auto& other = dynamic_cast<const DoseAccumulable&>(rhs);

  fEnDep += other.GetEnDeposit();
  fPrimaryEnDep += other.GetPrimaryEnDeposit();
  fSecondaryEnDep += other.GetSecondaryEnDeposit();

  // Different worker threads may observe different secondary species.
  for (const auto& entry : other.GetSecondaryEnDepositByParticle()) {
    const G4String& particleName = entry.first;
    const array_type& otherEnergy = entry.second;

    auto found = fSecondaryEnDepByParticle.find(particleName);
    if (found == fSecondaryEnDepByParticle.end()) {
      fSecondaryEnDepByParticle.emplace(particleName, otherEnergy);
    }
    else {
      found->second += otherEnergy;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Reset()
{
  if (GetVerboseLevel() > 0) G4cout << "DoseAccumulable::Reset()" << G4endl;
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Accumulate(Hit* hit)
{
  if (hit->GetDeltaE() <= 0.) return;

  if (GetVerboseLevel() > 1) G4cout << "DoseAccumulable::Accumulate()" << G4endl;

  auto* voxelizedDetector = VoxelizedSensitiveDetector::GetInstance();
  const G4int voxel =
    voxelizedDetector->GetThisVoxelNumber(hit->GetXindex(), hit->GetYindex(), hit->GetZindex());
  if (voxel < 0 || voxel >= voxelizedDetector->GetTotalVoxelNumber()) {
    G4Exception("DoseAccumulable::Accumulate", "InvalidVoxelIndex", JustWarning,
                "Ignoring an energy deposit with an invalid voxel index.");
    return;
  }

  const G4double depositedEnergy = hit->GetDeltaE();

  // The original total dose already contains energy deposited by every track.
  fEnDep[voxel] += depositedEnergy;

  // In Geant4, primary tracks have ParentID == 0. Every descendant has ParentID > 0.
  if (hit->GetParentID() == 0) {
    fPrimaryEnDep[voxel] += depositedEnergy;
    return;
  }

  fSecondaryEnDep[voxel] += depositedEnergy;

  const G4ParticleDefinition* particle = hit->GetPartType();
  const G4String particleName = particle ? particle->GetParticleName() : G4String("unknown");

  auto found = fSecondaryEnDepByParticle.find(particleName);
  if (found == fSecondaryEnDepByParticle.end()) {
    const auto nVoxels = voxelizedDetector->GetTotalVoxelNumber();
    found = fSecondaryEnDepByParticle.emplace(particleName, array_type(0.0, nVoxels)).first;
  }
  found->second[voxel] += depositedEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DoseAccumulable::GetVerboseLevel() const
{
  return Manager::GetInstance()->GetQuantity("Dose")->GetVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseAccumulable::Initialize()
{
  if (GetVerboseLevel() > 0) G4cout << "DoseAccumulable::Initialize(): " << G4endl;

  const auto nVoxels = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();
  fEnDep = array_type(0.0, nVoxels);
  fPrimaryEnDep = array_type(0.0, nVoxels);
  fSecondaryEnDep = array_type(0.0, nVoxels);
  fSecondaryEnDepByParticle.clear();
}

}  // namespace RadioBio
