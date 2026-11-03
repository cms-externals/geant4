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
/// \file radiobiology/src/Dose.cc
/// \brief Implementation of the RadioBio::Dose class

#include "Dose.hh"

#include "G4SystemOfUnits.hh"

#include "DoseAccumulable.hh"
#include "DoseMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

#include <cctype>
#include <fstream>
#include <iomanip>

namespace
{
constexpr std::streamsize columnWidth = 28;

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

Dose::Dose()
{
  fPath = "RadioBio_Dose.out";
  fMessenger = new DoseMessenger(this);
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Dose::~Dose()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Initialize()
{
  if (fVerboseLevel > 0) G4cout << "Dose::Initialize() called" << G4endl;

  const auto nVoxels = VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber();

  fEnDep.resize(nVoxels);
  fPrimaryEnDep.resize(nVoxels);
  fSecondaryEnDep.resize(nVoxels);
  fDose.resize(nVoxels);
  fPrimaryDose.resize(nVoxels);
  fSecondaryDose.resize(nVoxels);

  fSecondaryEnDepByParticle.clear();
  fSecondaryDoseByParticle.clear();
  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Compute()
{
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "Dose::Compute() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fCalculated) return;

  if (fVerboseLevel > 0) G4cout << "Dose::Compute() called" << G4endl;

  auto* voxSensDet = VoxelizedSensitiveDetector::GetInstance();
  if (voxSensDet == nullptr) {
    G4Exception("Dose::Compute", "VoxNotInit", FatalException,
                "Calculating dose without voxelized geometry pointer!");
  }

  const G4double voxelMass = voxSensDet->GetVoxelMass();

  for (G4int voxel = 0; voxel < voxSensDet->GetTotalVoxelNumber(); ++voxel) {
    fDose[voxel] = fEnDep[voxel] / voxelMass;
    fPrimaryDose[voxel] = fPrimaryEnDep[voxel] / voxelMass;
    fSecondaryDose[voxel] = fSecondaryEnDep[voxel] / voxelMass;
  }

  fSecondaryDoseByParticle.clear();
  for (const auto& entry : fSecondaryEnDepByParticle) {
    array_type particleDose(0.0, entry.second.size());
    particleDose = entry.second / voxelMass;
    fSecondaryDoseByParticle.emplace(entry.first, particleDose);
  }

  fCalculated = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Store()
{
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "Dose::Store() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fSaved) {
    G4Exception("Dose::Store", "DoseOverwrite", JustWarning,
                "Overwriting Dose file. For multiple runs, change filename.");
  }

  Compute();

  std::ofstream ofs(fPath);
  if (!ofs.is_open()) {
    G4Exception("Dose::Store", "DoseFileNotOpen", JustWarning,
                ("Cannot open dose output file: " + fPath).c_str());
    return;
  }

  ofs << std::left << std::setw(12) << "x_index" << std::setw(12) << "y_index" << std::setw(12)
      << "z_index" << std::setw(columnWidth) << "Dose_total(Gy)" << std::setw(columnWidth)
      << "Dose_prim(Gy)" << std::setw(columnWidth) << "Dose_sec(Gy)";

  for (const auto& entry : fSecondaryDoseByParticle) {
    const G4String header = "Dose_sec_" + MakeColumnSafe(entry.first) + "(Gy)";
    ofs << std::setw(columnWidth) << header;
  }
  ofs << G4endl;

  auto* voxSensDet = VoxelizedSensitiveDetector::GetInstance();

  ofs << std::scientific << std::setprecision(10);
  for (G4int i = 0; i < voxSensDet->GetVoxelNumberAlongX(); ++i) {
    for (G4int j = 0; j < voxSensDet->GetVoxelNumberAlongY(); ++j) {
      for (G4int k = 0; k < voxSensDet->GetVoxelNumberAlongZ(); ++k) {
        const G4int voxel = voxSensDet->GetThisVoxelNumber(i, j, k);

        ofs << std::left << std::setw(12) << i << std::setw(12) << j << std::setw(12) << k
            << std::setw(columnWidth) << fDose[voxel] / gray << std::setw(columnWidth)
            << fPrimaryDose[voxel] / gray << std::setw(columnWidth) << fSecondaryDose[voxel] / gray;

        for (const auto& entry : fSecondaryDoseByParticle) {
          ofs << std::setw(columnWidth) << entry.second[voxel] / gray;
        }
        ofs << G4endl;
      }
    }
  }

  if (fVerboseLevel > 0) G4cout << "Dose: Dose written to " << fPath << G4endl;
  fSaved = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::AddSecondaryEnergyDepositByParticle(const particle_map_type& depByParticle)
{
  for (const auto& entry : depByParticle) {
    auto found = fSecondaryEnDepByParticle.find(entry.first);
    if (found == fSecondaryEnDepByParticle.end()) {
      fSecondaryEnDepByParticle.emplace(entry.first, entry.second);
    }
    else {
      found->second += entry.second;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::AddFromAccumulable(G4VAccumulable* genericAccumulable)
{
  auto* accumulable = dynamic_cast<DoseAccumulable*>(genericAccumulable);
  if (accumulable == nullptr) {
    G4Exception("Dose::AddFromAccumulable", "WrongAccumulable", FatalException,
                "Dose received an accumulable of the wrong type.");
    return;
  }

  AddEnergyDeposit(accumulable->GetEnDeposit());
  AddPrimaryEnergyDeposit(accumulable->GetPrimaryEnDeposit());
  AddSecondaryEnergyDeposit(accumulable->GetSecondaryEnDeposit());
  AddSecondaryEnergyDepositByParticle(accumulable->GetSecondaryEnDepositByParticle());
  fCalculated = false;
  fSaved = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::Reset()
{
  if (fVerboseLevel > 1) G4cout << "Dose::Reset(): ";

  fEnDep = 0.0;
  fPrimaryEnDep = 0.0;
  fSecondaryEnDep = 0.0;
  fDose = 0.0;
  fPrimaryDose = 0.0;
  fSecondaryDose = 0.0;
  fSecondaryEnDepByParticle.clear();
  fSecondaryDoseByParticle.clear();
  fCalculated = false;
  fSaved = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dose::PrintParameters()
{
  G4cout << "*******************************************" << G4endl
         << "****** Parameters of the class Dose *******" << G4endl
         << "*******************************************" << G4endl;
  PrintVirtualParameters();
}

}  // namespace RadioBio
