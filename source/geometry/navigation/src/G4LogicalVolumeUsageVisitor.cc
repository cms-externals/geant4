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
// G4LogicalVolumeUsageVisitor class implementation
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------

#include "G4LogicalVolumeUsageVisitor.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "geomdefs.hh"

namespace
{
  // Null-safe names for the diagnostic printout: a logical volume may have no
  // solid or no material assigned (e.g. a partially-built or mother-only volume).
  G4String SolidName(G4LogicalVolume* lv)
  {
    const G4VSolid* solid = lv->GetSolid();
    return (solid != nullptr) ? solid->GetName() : G4String("<none>");
  }

  G4String MaterialName(G4LogicalVolume* lv)
  {
    const G4Material* material = lv->GetMaterial();
    return (material != nullptr) ? material->GetName() : G4String("<none>");
  }
}

// ---------------------------------------------------------------------------
// Internal helpers - index the per-type count vectors by LV instance ID
// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::IncrementCount(std::vector<std::uint32_t>& counts,
                                                 G4LogicalVolume* lv)
{
  if (lv == nullptr) return;
  G4int instanceID = lv->GetInstanceID();
  if (instanceID < 0) return;

  // The logical-volume instance ID indexes the count vector, which is grown
  // on demand to be large enough to hold this ID.
  if (instanceID >= (G4int)counts.size())
  {
    counts.resize(instanceID + 1, 0);
  }
  counts[instanceID]++;
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetCount(const std::vector<std::uint32_t>& counts,
                                                    G4LogicalVolume* lv) const
{
  if (lv == nullptr) return 0;
  G4int instanceID = lv->GetInstanceID();
  if (instanceID < 0 || instanceID >= (G4int)counts.size())
  {
    return 0;
  }
  return counts[instanceID];
}

// ---------------------------------------------------------------------------
// G4VVolumeVisitor Callbacks
// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::VisitPlacement(G4PVPlacement* pv, G4int /*depth*/)
{
  if (pv != nullptr)
  {
    IncrementCount(fPlacementCounts, pv->GetLogicalVolume());
  }
}

// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::VisitReplica(G4PVReplica* pv, G4int /*depth*/)
{
  if (pv != nullptr)
  {
    IncrementCount(fReplicaCounts, pv->GetLogicalVolume());
  }
}

// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::VisitParameterised(G4PVReplica* pv, G4int /*depth*/)
{
  // 'pv' may actually be a G4PVDivision or G4ReplicatedSlice; G4PVReplica
  // is their real common base with G4PVParameterised - see G4VVolumeVisitor.hh.
  if (pv != nullptr)
  {
    IncrementCount(fParameterisedCounts, pv->GetLogicalVolume());
  }
}

// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::VisitUnknown(G4VPhysicalVolume* pv, G4int /*depth*/)
{
  // Physical volumes of an unknown or external type (G4VExternalPhysicalVolume).
  if (pv != nullptr)
  {
    IncrementCount(fExternalCounts, pv->GetLogicalVolume());
  }
}

// ---------------------------------------------------------------------------
// Query Methods
// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::IsUsedAsPlacement(G4LogicalVolume* lv) const
{
  return GetCount(fPlacementCounts, lv) > 0;
}

// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::IsUsedAsReplica(G4LogicalVolume* lv) const
{
  return GetCount(fReplicaCounts, lv) > 0;
}

// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::IsUsedAsParameterised(G4LogicalVolume* lv) const
{
  return GetCount(fParameterisedCounts, lv) > 0;
}

// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::IsUsedAsExternal(G4LogicalVolume* lv) const
{
  return GetCount(fExternalCounts, lv) > 0;
}

// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::IsUsedAsReplicaOrParameterised(G4LogicalVolume* lv) const
{
  return IsUsedAsReplica(lv) || IsUsedAsParameterised(lv);
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetPlacementCount(G4LogicalVolume* lv) const
{
  return GetCount(fPlacementCounts, lv);
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetReplicaCount(G4LogicalVolume* lv) const
{
  return GetCount(fReplicaCounts, lv);
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetParameterisedCount(G4LogicalVolume* lv) const
{
  return GetCount(fParameterisedCounts, lv);
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetExternalCount(G4LogicalVolume* lv) const
{
  return GetCount(fExternalCounts, lv);
}

// ---------------------------------------------------------------------------

std::uint32_t G4LogicalVolumeUsageVisitor::GetReplicaOrParameterisedCount(G4LogicalVolume* lv) const
{
  return GetReplicaCount(lv) + GetParameterisedCount(lv);
}

// ---------------------------------------------------------------------------

G4bool G4LogicalVolumeUsageVisitor::ContainsOnlyPlacements(G4LogicalVolume* lv) const
{
  if (lv == nullptr) return false;
  G4int nDaughters = (G4int)lv->GetNoDaughters();
  if (nDaughters == 0) return false;

  // Check if all daughters are normal placements (not replicas/parameterisations)
  for (G4int i = 0; i < nDaughters; ++i)
  {
    G4VPhysicalVolume* daughter = lv->GetDaughter(i);
    if (daughter == nullptr || daughter->VolumeType() != kNormal)
    {
      return false;
    }
  }
  return true;
}

// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::Clear()
{
  fPlacementCounts.clear();
  fReplicaCounts.clear();
  fParameterisedCounts.clear();
  fExternalCounts.clear();
}

// ---------------------------------------------------------------------------

void G4LogicalVolumeUsageVisitor::PrintSummary() const
{
  // Query all logical volumes in the store and group them by usage type. A
  // logical volume used in more than one physical-volume type is collected in
  // 'multiple' rather than in any single-type list.
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  if (store == nullptr) return;

  std::vector<G4LogicalVolume*> placements;
  std::vector<G4LogicalVolume*> replicas;
  std::vector<G4LogicalVolume*> parameterised;
  std::vector<G4LogicalVolume*> externals;
  std::vector<G4LogicalVolume*> multiple;

  // Totals per type, counting every use (a 'multiple' volume adds to each
  // type it participates in).
  std::size_t nPlacement = 0, nReplica = 0, nParameterised = 0, nExternal = 0;

  for (G4LogicalVolume* lv : *store)
  {
    const G4bool isPlacement     = IsUsedAsPlacement(lv);
    const G4bool isReplica       = IsUsedAsReplica(lv);
    const G4bool isParameterised = IsUsedAsParameterised(lv);
    const G4bool isExternal      = IsUsedAsExternal(lv);

    if (isPlacement)     ++nPlacement;
    if (isReplica)       ++nReplica;
    if (isParameterised) ++nParameterised;
    if (isExternal)      ++nExternal;

    const G4int nTypes = (isPlacement ? 1 : 0) + (isReplica ? 1 : 0)
                       + (isParameterised ? 1 : 0) + (isExternal ? 1 : 0);

    if (nTypes == 0)          { continue; }            // not seen in the walked tree(s)
    else if (nTypes > 1)      { multiple.push_back(lv); }
    else if (isPlacement)     { placements.push_back(lv); }
    else if (isReplica)       { replicas.push_back(lv); }
    else if (isParameterised) { parameterised.push_back(lv); }
    else                      { externals.push_back(lv); }
  }

  auto printLine = [](G4LogicalVolume* lv, std::uint32_t count)
  {
    G4cout << "  - " << lv->GetName() << " [Solid: " << SolidName(lv)
           << ", Material: " << MaterialName(lv) << "] (Count: " << count << ")";
  };

  G4cout << "\n========== G4LogicalVolumeUsageVisitor Summary ==========" << G4endl;
  G4cout << "Total logical volumes in store: " << store->size() << G4endl;
  G4cout << "Used as Placements               : " << nPlacement << G4endl;
  G4cout << "Used as Replicas                 : " << nReplica << G4endl;
  G4cout << "Used as Parameterised            : " << nParameterised << G4endl;
  G4cout << "Used as External/Unknown         : " << nExternal << G4endl;
  G4cout << "Used in multiple types           : " << multiple.size() << G4endl;

  G4cout << "\nLogical Volumes used as Placements:" << G4endl;
  for (auto* lv : placements)
  {
    printLine(lv, GetPlacementCount(lv));
    if (ContainsOnlyPlacements(lv)) { G4cout << " (Mother containing only placements)"; }
    G4cout << G4endl;
  }

  G4cout << "\nLogical Volumes used as Replicas:" << G4endl;
  for (auto* lv : replicas)
  {
    printLine(lv, GetReplicaCount(lv));
    G4cout << G4endl;
  }

  G4cout << "\nLogical Volumes used as Parameterised:" << G4endl;
  for (auto* lv : parameterised)
  {
    printLine(lv, GetParameterisedCount(lv));
    G4cout << G4endl;
  }

  G4cout << "\nLogical Volumes used as External/Unknown:" << G4endl;
  for (auto* lv : externals)
  {
    printLine(lv, GetExternalCount(lv));
    G4cout << G4endl;
  }

  G4cout << "\nLogical Volumes used in multiple types:" << G4endl;
  for (auto* lv : multiple)
  {
    G4cout << "  - " << lv->GetName() << " [Solid: " << SolidName(lv)
           << ", Material: " << MaterialName(lv) << "]"
           << " placements="    << GetPlacementCount(lv)
           << " replicas="      << GetReplicaCount(lv)
           << " parameterised=" << GetParameterisedCount(lv)
           << " external="      << GetExternalCount(lv) << G4endl;
  }

  if (!multiple.empty())
  {
    G4cout << "\nWARNING: The following logical volumes are used in multiple"
           << " physical-volume types:" << G4endl;
    for (auto* lv : multiple)
    {
      G4cout << "  * " << lv->GetName() << G4endl;
    }
  }
  else
  {
    G4cout << "\nNo logical volumes are used in multiple physical-volume types." << G4endl;
  }
  G4cout << "=========================================================\n" << G4endl;
}
