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
// G4GeometryMemoryAdvisor
//
// Class description:
//
// Walks a constructed geometry (via G4VolumeTreeWalker) and reports the
// mutable/shared split of logical volumes (per G4LogicalVolume::
// MarkMaterialMutable()), the per-thread memory this saves, and an
// inventory of physical-volume types with their per-thread G4PVData cost.
// Sizes are reported both as platform-independent "reference" values and
// as "this platform" sizeof() values; see the .cc file for details and
// caveats (e.g. sizeof() excludes any heap memory a solid owns internally).
//
// Usage:
//   G4GeometryMemoryAdvisor advisor;
//   advisor.Analyse(worldPhysicalVolume);
//   advisor.PrintReport(G4cout);
//   advisor.SpaceOfPhysicalVolumes(G4cout);
//
// Author: John Apostolakis (CERN), July 2026
// --------------------------------------------------------------------
#ifndef G4GEOMETRYMEMORYADVISOR_HH
#define G4GEOMETRYMEMORYADVISOR_HH

#include "G4VVolumeVisitor.hh"
#include "G4Types.hh"

#include <iosfwd>
#include <set>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4PVReplica;

/**
 * @brief G4GeometryMemoryAdvisor analyses a constructed geometry and
 * estimates the per-thread memory saved automatically by the runtime
 * mutable/shared logical-volume split (G4LogicalVolume::
 * MarkMaterialMutable()).
 * @ingroup geometry_navigation
 */
class G4GeometryMemoryAdvisor : public G4VVolumeVisitor
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4GeometryMemoryAdvisor() = default;
    ~G4GeometryMemoryAdvisor() override = default;

    /**
     * Walks the subtree rooted at @p world and gathers statistics.
     * Resets any statistics from a previous call first.
     */
    void Analyse(G4VPhysicalVolume* world);

    /**
     * G4VVolumeVisitor callbacks, used internally by Analyse().
     */
    void VisitPlacement(G4PVPlacement* pv, G4int depth) override;
    void VisitReplica(G4PVReplica* pv, G4int depth) override;
    void VisitParameterised(G4PVReplica* pv, G4int depth) override;
    void VisitUnknown(G4VPhysicalVolume* pv, G4int depth) override;
    void EnterLogicalVolume(G4LogicalVolume* lv, G4int depth) override;

    /**
     * Queries, valid after Analyse().
     */
    std::size_t GetNumberOfLogicalVolumes() const { return fLogicalVolumesSeen.size(); }

    /**
     * Distinct logical volumes used as the daughter of a
     * G4PVParameterised/G4PVDivision/G4ReplicatedSlice - flagged by
     * G4LogicalVolume::MarkMaterialMutable() for a per-thread copy.
     */
    std::size_t GetNumberOfMutableLogicalVolumes() const { return fMutableLogicalVolumesSeen.size(); }

    /**
     * Distinct logical volumes not used inside a parameterised/division
     * placement - these share their solid/material/mass/cuts-couple
     * across all threads.
     */
    std::size_t GetNumberOfSharedLogicalVolumes() const
    {
      return GetNumberOfLogicalVolumes() - GetNumberOfMutableLogicalVolumes();
    }

    std::size_t GetNumberOfPlacements() const { return fPlacementCount; }
    std::size_t GetNumberOfReplicaInstances() const { return fReplicaInstanceCount; }
    std::size_t GetNumberOfParameterisedOrDivisionInstances() const { return fParamInstanceCount; }
    std::size_t GetNumberOfUnknownInstances() const { return fUnknownCount; }

    /**
     * Total number of physical-volume instances of any kind found.
     */
    std::size_t GetTotalPhysicalVolumes() const;

    /**
     * Writes a human-readable report to @p os: geometry composition,
     * mutable vs shared logical-volume counts, and the estimated per-thread
     * memory saving from the runtime mutable/shared split.
     *  @param[in] nThreads Worker-thread count; if <= 0, uses
     *             G4Threading::G4GetNumberOfCores().
     */
    void PrintReport(std::ostream& os, G4int nThreads = 0) const;

    /**
     * Writes an inventory of the physical-volume types found in the global
     * G4PhysicalVolumeStore to @p os: count, master-thread instance size,
     * and additional per-thread G4PVData cost, per type.
     *  @param[in] nThreads Worker-thread count; if <= 0, uses
     *             G4Threading::G4GetNumberOfCores().
     */
    void SpaceOfPhysicalVolumes(std::ostream& os, G4int nThreads = 0) const;

  private:

    /**
     * Helper for PrintReport(): the logical-volume-side footprint table
     * (reference vs. this-platform sizes) and per-thread split overhead.
     */
    void ReportLogicalVolumeFootprint(std::ostream& os, G4int nThreads) const;

    std::set<const G4LogicalVolume*> fLogicalVolumesSeen;
    std::set<const G4LogicalVolume*> fMutableLogicalVolumesSeen;
    std::size_t fPlacementCount = 0;
    std::size_t fReplicaInstanceCount = 0;
    std::size_t fParamInstanceCount = 0;
    std::size_t fUnknownCount = 0;
};

#endif
