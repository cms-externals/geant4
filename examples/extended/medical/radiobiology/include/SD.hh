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
/// \file SD.hh
/// \brief Definition of the RadioBio::SD class

#ifndef RADIOBIOLOGYSD_H
#define RADIOBIOLOGYSD_H

#include "G4VSensitiveDetector.hh"

#include "Hit.hh"

#include <set>
#include <tuple>
#include <vector>

class G4Step;
class G4HCofThisEvent;

namespace RadioBio
{

/// Creates hits and dispatches them to enabled accumulables.
class SD : public G4VSensitiveDetector
{
  public:
    SD(const G4String& name, const G4String& hitsCollectionName);
    ~SD() override = default;

    // Methods from base class
    void Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void EndOfEvent(G4HCofThisEvent* hitCollection) override;

  private:
    RadioBioHitsCollection* fHitsCollection = nullptr;

    // Track IDs are unique only within one event. This set is therefore
    // cleared in Initialize() and prevents multiple counting of the same
    // track in the same voxel after backscattering or re-entry.
    // Key: (trackID, xIndex, yIndex, zIndex).
    std::set<std::tuple<G4int, G4int, G4int, G4int>> fCountedVoxelEntries{};
};

}  // namespace RadioBio

#endif  // RADIOBIOLOGYSD_H
