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
// G4ParameterisationChecker implementation
//
// Author: John Apostolakis (CERN), July 2026
// --------------------------------------------------------------------

#include "G4ParameterisationChecker.hh"

#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

#include <set>

// ********************************************************************
// CheckAndMark
// ********************************************************************
//
G4int G4ParameterisationChecker::CheckAndMark(G4bool verbose)
{
  // Marking a logical volume mutable creates shared/master-side data
  // structures; this must be done once, on the master thread, before the
  // per-thread copies are made. Do nothing on workers.
  if (!G4Threading::IsMasterThread())
  {
    return 0;
  }

  G4int nFixed = 0;

  for (auto* pv : *G4PhysicalVolumeStore::GetInstance())
  {
    if (pv == nullptr || !pv->IsParameterised())
    {
      continue;
    }
    G4VPVParameterisation* param = pv->GetParameterisation();
    if (param == nullptr)
    {
      continue;
    }
    const G4int nCopies = pv->GetMultiplicity();

    // Snapshot the physical volume's mutable state: the exploratory
    // ComputeTransformation() calls below reassign its translation/rotation/
    // name/logical-volume, and we must leave the volume exactly as we found
    // it so that later navigation recomputes from a clean starting point.
    G4LogicalVolume* savedLV = pv->GetLogicalVolume();
    const G4ThreeVector savedTranslation = pv->GetTranslation();
    G4RotationMatrix* savedRotation = pv->GetRotation();
    const G4String savedName = pv->GetName();

    // Discover every distinct logical volume this parameterisation binds.
    // A parameterisation may reassign the logical volume from any of its
    // three per-copy hooks (ComputeTransformation(), ComputeSolid(),
    // ComputeMaterial()), so all three must be exercised - probing only
    // ComputeTransformation() would silently miss reassignments made from
    // the other two.
    std::set<G4LogicalVolume*> boundVolumes;
    auto captureBoundVolume = [&pv, &boundVolumes]() {
      G4LogicalVolume* lv = pv->GetLogicalVolume();
      if (lv != nullptr)
      {
        boundVolumes.insert(lv);
      }
    };
    const G4bool isNested = param->IsNested();
    for (G4int copyNo = 0; copyNo < nCopies; ++copyNo)
    {
      param->ComputeTransformation(copyNo, pv);
      captureBoundVolume();

      // Discard the returned solid: the checker only needs to discover
      // which logical volume this hook binds, not to install a solid.
      (void)param->ComputeSolid(copyNo, pv);
      captureBoundVolume();

      // ComputeMaterial()'s contract requires it to cope with a valid
      // parent touchable when the parameterisation is nested (used by the
      // navigator's SetupHierarchy()). No such touchable exists here at
      // CloseGeometry() time, so nested parameterisations are not probed
      // via this hook - a documented gap, not exercised for correctness.
      if (!isNested)
      {
        // Discard the returned material for the same reason as above.
        (void)param->ComputeMaterial(copyNo, pv, nullptr);
        captureBoundVolume();
      }
    }

    // Restore the physical volume to its pre-check state.
    pv->SetLogicalVolume(savedLV);
    pv->SetTranslation(savedTranslation);
    pv->SetRotation(savedRotation);
    pv->SetName(savedName);

    // Mark any bound logical volume that is not yet mutable.
    for (auto* lv : boundVolumes)
    {
      if (!lv->IsMaterialMutable())
      {
        lv->MarkMaterialMutable();
        ++nFixed;
        if (verbose)
        {
          G4cout << "G4ParameterisationChecker: marked logical volume '"
                 << lv->GetName() << "' as mutable - reassigned by parameterised"
                 << " physical volume '" << savedName << "'." << G4endl;
        }
      }
    }
  }

  // Always report the number of volumes repaired.
  G4cout << "G4ParameterisationChecker: " << nFixed
         << " logical volume(s) marked mutable due to parameterisation reassignment."
         << G4endl;

  return nFixed;
}
