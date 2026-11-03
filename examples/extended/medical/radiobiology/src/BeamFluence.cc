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
/// \file radiobiology/src/BeamFluence.cc
/// \brief Implementation of the RadioBio::BeamFluence class

#include "BeamFluence.hh"

#include "BeamFluenceAccumulable.hh"
#include "BeamFluenceMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

#include <cctype>
#include <fstream>
#include <iomanip>

namespace
{
constexpr std::streamsize columnWidth = 30;

G4String MakeColumnSafe(G4String name)
{
  for (auto& character : name) {
    const unsigned char c = static_cast<unsigned char>(character);
    if (!(std::isalnum(c) || character == '+' || character == '-')) character = '_';
  }
  return name;
}
}  // namespace

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamFluence::BeamFluence()
{
  fPath = "RadioBio_BeamFluence.out";
  fVerboseLevel = 0;
  fMessenger = new BeamFluenceMessenger(this);
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamFluence::~BeamFluence()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::Initialize()
{
  const auto nVoxels = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();

  fTotalParticleCount = count_array_type(G4long{0}, nVoxels);
  fPrimaryParticleCount = count_array_type(G4long{0}, nVoxels);
  fSecondaryParticleCount = count_array_type(G4long{0}, nVoxels);
  fParticleCountBySpecies.clear();

  fInitialized = true;
  fCalculated = false;
  fSaved = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::Compute()
{
  if (!fCalculationEnabled) {
    return;
  }

  // The accumulated quantities are particle-entry counts; no conversion is required.
  fCalculated = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::Store()
{
  if (!fCalculationEnabled) {
    return;
  }

  if (fSaved) {
    G4Exception("BeamFluence::Store", "BeamFluenceOverwrite", JustWarning,
                "Overwriting beam-fluence file. For multiple runs, change filename.");
  }

  Compute();

  std::ofstream output(fPath);
  if (!output.is_open()) {
    G4Exception("BeamFluence::Store", "BeamFluenceFileNotOpen", JustWarning,
                ("Cannot open beam-fluence output file: " + fPath).c_str());
    return;
  }

  output << std::left << std::setw(12) << "x_index" << std::setw(12) << "y_index" << std::setw(12)
         << "z_index" << std::setw(columnWidth) << "N_enter_total" << std::setw(columnWidth)
         << "N_enter_prim" << std::setw(columnWidth) << "N_enter_sec";

  for (const auto& entry : fParticleCountBySpecies) {
    const G4String header = "N_enter_" + MakeColumnSafe(entry.first);
    output << std::setw(columnWidth) << header;
  }
  output << G4endl;

  auto* voxelizedDetector = VoxelizedSensitiveDetector::GetInstance();

  for (G4int i = 0; i < voxelizedDetector->GetVoxelNumberAlongX(); ++i) {
    for (G4int j = 0; j < voxelizedDetector->GetVoxelNumberAlongY(); ++j) {
      for (G4int k = 0; k < voxelizedDetector->GetVoxelNumberAlongZ(); ++k) {
        const G4int voxel = voxelizedDetector->GetThisVoxelNumber(i, j, k);

        output << std::left << std::setw(12) << i << std::setw(12) << j << std::setw(12) << k
               << std::setw(columnWidth) << fTotalParticleCount[voxel] << std::setw(columnWidth)
               << fPrimaryParticleCount[voxel] << std::setw(columnWidth)
               << fSecondaryParticleCount[voxel];

        for (const auto& entry : fParticleCountBySpecies) {
          output << std::setw(columnWidth) << entry.second[voxel];
        }
        output << G4endl;
      }
    }
  }

  if (fVerboseLevel > 0) {
    G4cout << "BeamFluence: voxel-resolved particle-entry counts written to " << fPath << G4endl;
  }

  fSaved = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::AddParticleCountBySpecies(const particle_map_type& countsBySpecies)
{
  for (const auto& entry : countsBySpecies) {
    auto found = fParticleCountBySpecies.find(entry.first);
    if (found == fParticleCountBySpecies.end()) {
      fParticleCountBySpecies.emplace(entry.first, entry.second);
    }
    else {
      found->second += entry.second;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::AddFromAccumulable(G4VAccumulable* genericAccumulable)
{
  auto* accumulable = dynamic_cast<BeamFluenceAccumulable*>(genericAccumulable);
  if (accumulable == nullptr) {
    G4Exception("BeamFluence::AddFromAccumulable", "WrongAccumulable", FatalException,
                "BeamFluence received an accumulable of the wrong type.");
    return;
  }

  fTotalParticleCount += accumulable->GetTotalParticleCount();
  fPrimaryParticleCount += accumulable->GetPrimaryParticleCount();
  fSecondaryParticleCount += accumulable->GetSecondaryParticleCount();
  AddParticleCountBySpecies(accumulable->GetParticleCountBySpecies());

  fCalculated = false;
  fSaved = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::Reset()
{
  fTotalParticleCount = G4long{0};
  fPrimaryParticleCount = G4long{0};
  fSecondaryParticleCount = G4long{0};
  fParticleCountBySpecies.clear();
  fCalculated = false;
  fSaved = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluence::PrintParameters()
{
  G4cout << "************************************************" << G4endl
         << "****** Parameters of class BeamFluence ********" << G4endl
         << "************************************************" << G4endl
         << "Counts each track once in every voxel it visits." << G4endl
         << "Results are stored for every (x,y,z) voxel and" << G4endl
         << "separated dynamically by particle species." << G4endl;
  PrintVirtualParameters();
}

}  // namespace RadioBio
