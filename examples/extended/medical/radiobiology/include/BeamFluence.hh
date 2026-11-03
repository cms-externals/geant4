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
/// \file BeamFluence.hh
/// \brief Definition of the RadioBio::BeamFluence class

#ifndef RADIOBIOLOGYBEAMFLUENCE_H
#define RADIOBIOLOGYBEAMFLUENCE_H

#include "VRadiobiologicalQuantity.hh"

#include <map>
#include <valarray>

namespace RadioBio
{

class BeamFluenceMessenger;

/// Stores voxel-resolved unique track-entry counts.
///
/// A track is counted once in each voxel it visits, including the voxel in
/// which a secondary track is created. The output contains total, primary,
/// and secondary counts, plus one dynamically generated column for each
/// particle species observed during the run. These counts are not normalized
/// by area and therefore are not a differential fluence estimator.
class BeamFluence : public VRadiobiologicalQuantity
{
  public:
    using count_array_type = std::valarray<G4long>;
    using particle_map_type = std::map<G4String, count_array_type>;

    BeamFluence();
    ~BeamFluence() override;

    void AddFromAccumulable(G4VAccumulable*) override;
    void Initialize() override;
    void Compute() override;
    void Reset() override;
    void Store() override;
    void PrintParameters() override;

    const count_array_type& GetTotalParticleCount() const { return fTotalParticleCount; }
    const count_array_type& GetPrimaryParticleCount() const { return fPrimaryParticleCount; }
    const count_array_type& GetSecondaryParticleCount() const { return fSecondaryParticleCount; }
    const particle_map_type& GetParticleCountBySpecies() const { return fParticleCountBySpecies; }

  private:
    void AddParticleCountBySpecies(const particle_map_type& countsBySpecies);

    count_array_type fTotalParticleCount{};
    count_array_type fPrimaryParticleCount{};
    count_array_type fSecondaryParticleCount{};
    particle_map_type fParticleCountBySpecies{};
};

}  // namespace RadioBio

#endif  // RADIOBIOLOGYBEAMFLUENCE_H
