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
// G4LogicalVolumeUsageVisitor
//
// Class description:
//
// G4VVolumeVisitor that gathers usage statistics for G4LogicalVolume
// instances across the geometry hierarchy, counting separately every type of
// physical volume in which each logical volume is used (placement, replica,
// parameterised, and external/unknown), indexed by logical volume instance ID.
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------
#ifndef G4LOGICALVOLUMEUSAGEVISITOR_HH
#define G4LOGICALVOLUMEUSAGEVISITOR_HH

#include "G4VVolumeVisitor.hh"
#include "globals.hh"
#include <vector>
#include <cstdint>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4PVReplica;

/**
 * @brief G4LogicalVolumeUsageVisitor checks usage of logical volumes in the geometry tree.
 *
 * Checks if each G4LogicalVolume is used as a placement (G4PVPlacement) or as a
 * replica/parameterised volume (G4PVReplica), gathering statistics across the
 * geometry hierarchy. It optimizes this lookup by indexing counts via the logical
 * volume's instance ID.
 *
 * This class walks the tree of a geometry hierarchy once for each subtree
 * (not once per touchable) and counts the volumes as it encounters them.
 *
 * Notes:
 *  - It is expected that a particular PhysicalVolume instance will always refer to
 *    the same G4LogicalVolume instance throughout the lifetime of an application -
 *    certainly during a 'Run'.
 *  - If the simulation uses 'User' methods or the methods of G4PVParameterisation to
 *    change the logical volume, this will not be evident from this count.
 *    G4PVParameterised volume's parameterisations which do this are 'unorthodox'
 *    and require separate treatment (see G4ParameterisationChecker).
 *
 * @ingroup geometry_navigation
 */

class G4LogicalVolumeUsageVisitor : public G4VVolumeVisitor
{
  public:

    /// Default constructor.
    G4LogicalVolumeUsageVisitor() = default;

    /// Destructor.
    ~G4LogicalVolumeUsageVisitor() override = default;

    // -- G4VVolumeVisitor Callbacks -------------------------------------------

    /**
     * @brief Callback invoked when visiting a placed physical volume (G4PVPlacement).
     * @param pv Pointer to the physical volume.
     * @param depth Current depth in the geometry tree.
     */
    void VisitPlacement(G4PVPlacement* pv, G4int depth) override;

    /**
     * @brief Callback invoked when visiting a replica physical volume (G4PVReplica).
     * @param pv Pointer to the replica physical volume.
     * @param depth Current depth in the geometry tree.
     */
    void VisitReplica(G4PVReplica* pv, G4int depth) override;

    /**
     * @brief Callback invoked when visiting a parameterised, divided, or replicated slice volume.
     * @param pv Pointer to the physical volume (represented as G4PVReplica).
     * @param depth Current depth in the geometry tree.
     */
    void VisitParameterised(G4PVReplica* pv, G4int depth) override;

    /**
     * @brief Callback invoked when visiting a physical volume of an unknown or external type.
     * @param pv Pointer to the physical volume.
     * @param depth Current depth in the geometry tree.
     */
    void VisitUnknown(G4VPhysicalVolume* pv, G4int depth) override;

    // -- Queries --------------------------------------------------------------

    /**
     * @brief Check if a logical volume is used as a placement (kNormal) anywhere.
     * @param lv Pointer to the logical volume to query.
     * @return True if the logical volume is used as a placement.
     */
    G4bool IsUsedAsPlacement(G4LogicalVolume* lv) const;

    /**
     * @brief Check if a logical volume is used as a plain replica (kReplica).
     * @param lv Pointer to the logical volume to query.
     * @return True if used as a replica.
     */
    G4bool IsUsedAsReplica(G4LogicalVolume* lv) const;

    /**
     * @brief Check if a logical volume is used as a parameterised volume
     *        (kParameterised, including G4PVDivision and G4ReplicatedSlice).
     * @param lv Pointer to the logical volume to query.
     * @return True if used as a parameterised volume.
     */
    G4bool IsUsedAsParameterised(G4LogicalVolume* lv) const;

    /**
     * @brief Check if a logical volume is used as an external/unknown volume (kExternal).
     * @param lv Pointer to the logical volume to query.
     * @return True if used as an external or unknown physical-volume type.
     */
    G4bool IsUsedAsExternal(G4LogicalVolume* lv) const;

    /**
     * @brief Check if a logical volume is used as a replica or parameterised volume.
     * @param lv Pointer to the logical volume to query.
     * @return True if used as either a replica or a parameterised volume.
     */
    G4bool IsUsedAsReplicaOrParameterised(G4LogicalVolume* lv) const;

    /**
     * @brief Check if a logical volume contains only placement daughter volumes.
     * @param lv Pointer to the logical volume to query.
     * @return True if all daughter volumes are placements (type kNormal).
     */
    G4bool ContainsOnlyPlacements(G4LogicalVolume* lv) const;

    /**
     * @brief Get the count of placement (kNormal) instances referencing a logical volume.
     * @param lv Pointer to the logical volume.
     * @return Number of placement instances.
     */
    std::uint32_t GetPlacementCount(G4LogicalVolume* lv) const;

    /**
     * @brief Get the count of replica (kReplica) instances referencing a logical volume.
     * @param lv Pointer to the logical volume.
     * @return Number of replica instances.
     */
    std::uint32_t GetReplicaCount(G4LogicalVolume* lv) const;

    /**
     * @brief Get the count of parameterised instances referencing a logical volume.
     * @param lv Pointer to the logical volume.
     * @return Number of parameterised instances.
     */
    std::uint32_t GetParameterisedCount(G4LogicalVolume* lv) const;

    /**
     * @brief Get the count of external/unknown instances referencing a logical volume.
     * @param lv Pointer to the logical volume.
     * @return Number of external/unknown instances.
     */
    std::uint32_t GetExternalCount(G4LogicalVolume* lv) const;

    /**
     * @brief Get the combined count of replica and parameterised instances.
     * @param lv Pointer to the logical volume.
     * @return Number of replica plus parameterised instances.
     */
    std::uint32_t GetReplicaOrParameterisedCount(G4LogicalVolume* lv) const;

    /**
     * @brief Clear all collected usage statistics.
     */
    void Clear();

    /**
     * @brief Print a summary report of logical volume usage to G4cout.
     */
    void PrintSummary() const;

  private:

    /**
     * @brief Increment (in the given per-type count vector) the usage count
     *        for a logical volume, growing the vector on demand.
     * @param counts The per-type count vector, indexed by LV instance ID.
     * @param lv Pointer to the logical volume.
     */
    void IncrementCount(std::vector<std::uint32_t>& counts, G4LogicalVolume* lv);

    /**
     * @brief Read (from the given per-type count vector) the usage count for a
     *        logical volume; returns 0 if it was never recorded.
     * @param counts The per-type count vector, indexed by LV instance ID.
     * @param lv Pointer to the logical volume.
     */
    std::uint32_t GetCount(const std::vector<std::uint32_t>& counts,
                           G4LogicalVolume* lv) const;

  private:

    // One count vector per physical-volume type (EVolume), indexed by the
    // logical volume's instance ID.
    std::vector<std::uint32_t> fPlacementCounts;      // kNormal
    std::vector<std::uint32_t> fReplicaCounts;        // kReplica
    std::vector<std::uint32_t> fParameterisedCounts;  // kParameterised
    std::vector<std::uint32_t> fExternalCounts;       // kExternal / unknown
};

#endif
