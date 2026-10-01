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
// G4ParameterisationCheckVisitor
//
// Class description:
//
// Concrete G4VTouchableVisitor that exercises every parameterisation found
// during a G4TouchableWalker traversal with a real, live parent touchable -
// including parameterisations nested inside a replica (level 2) or inside
// two levels of replica (level 3). This is the complement of
// G4ParameterisationChecker: that class iterates the flat physical-volume
// store and cannot supply a parent touchable, so it must skip
// ComputeMaterial() for nested parameterisations (IsNested() == true). Here
// the walker builds a genuine G4TouchableHistory of the actual ancestor
// chain before every ComputeMaterial() call, so nested parameterisations are
// triggered exactly as the navigator triggers them at tracking time.
//
// Detects:
//  - a parameterisation that reassigns its physical volume's logical volume
//    (via ComputeTransformation(), ComputeSolid() or ComputeMaterial()) to a
//    logical volume not marked mutable - the case most likely to be missed,
//    since the walker's own G4LogicalVolume::SetSolid()/UpdateMaterial() are
//    silently applied to the *old* logical volume when this happens;
//  - materials returned by a nested ComputeMaterial() that are absent from
//    its declared G4VVolumeMaterialScanner list, which G4Region::ScanVolumeTree
//    relies on exclusively to build production-cut couples;
//  - a nested parameterisation whose material never appears to depend on its
//    parent touchable, suggesting the nesting is not genuinely exercised;
//  - null materials and other basic per-instance sanity conditions.
//
// The walk is O(total touchables) and allocates a G4TouchableHistory per
// node, so by default it samples rather than exhaustively visiting: see
// SetMaxInstancesPerVolume(), SetMaxParentSamples(), SetMaxTotalVisits() and
// SetPruneUnparameterisedSubtrees(). Set any cap to -1 for exhaustive
// coverage of that dimension.
//
// CheckGeometry() is the entry point; it also owns the required snapshot and
// restoration of physical- and logical-volume state around the walk, since
// (unlike G4ParameterisationChecker) G4TouchableWalker does not restore what
// it changes, and a non-mutable logical volume's state is overwritten with
// no diagnostic when this runs on the master thread (see
// G4LogicalVolume::CheckGeometryConfiguration()).
//
// Author: John Apostolakis (CERN), August 2026
// --------------------------------------------------------------------
#ifndef G4PARAMETERISATIONCHECKVISITOR_HH
#define G4PARAMETERISATIONCHECKVISITOR_HH

#include "G4Types.hh"
#include "G4VTouchableVisitor.hh"

#include <iosfwd>
#include <map>
#include <set>

class G4Material;
class G4VPhysicalVolume;

/**
 * @brief Concrete G4VTouchableVisitor that drives a G4TouchableWalker
 * traversal to exercise every parameterisation - including those nested
 * inside one or two levels of replica - with a real, live parent touchable.
 * @ingroup geometry_navigation
 */

class G4ParameterisationCheckVisitor : public G4VTouchableVisitor
{
  public:

    /** Accumulated findings for one parameterised physical volume. */
    struct VolumeReport
    {
      /** The parameterised physical volume this report is about. */
      G4VPhysicalVolume* fPV = nullptr;

      /** True if this volume's parameterisation is a G4VNestedParameterisation. */
      G4bool fIsNested = false;

      /** True if reported via VisitRegularStructure() rather than expanded. */
      G4bool fIsRegularStructure = false;

      /** Number of copies GetMultiplicity() declares for this volume. */
      G4long fInstancesDeclared = 0;

      /** Number of instances actually visited (may be less, if capped). */
      G4long fInstancesVisited = 0;

      /** Shallowest and deepest history depth at which an instance was seen. */
      G4int fMinDepth = 0;
      G4int fMaxDepth = 0;

      /** Every logical volume observed bound to this physical volume. */
      std::set<G4LogicalVolume*> fBoundLVs;

      /** Bound logical volumes that differ from the one declared at construction. */
      std::set<G4LogicalVolume*> fReboundLVs;

      /** Every non-null material observed on the bound logical volume. */
      std::set<G4Material*> fMaterialsSeen;

      /** Materials observed that are absent from the nested parameterisation's
       *  declared G4VVolumeMaterialScanner list (and from GetBaseMaterial()). */
      std::set<G4Material*> fUndeclaredMaterials;

      /** Instances with a null bound material (excluding rebound instances,
       *  where the reading is indeterminate - see .cc for rationale). */
      G4long fNullMaterialCount = 0;

      /** Nested-parameterisation instances visited at history depth 0, where
       *  there is no parent replica to key the material on. */
      G4long fShallowNestedCount = 0;

      /** Number of distinct (parent physical volume, parent replica number)
       *  contexts this instance's material computation was exercised under. */
      std::size_t fDistinctParentContexts = 0;

      /** True if, for at least one replica number, the material differed
       *  across distinct parent contexts - positive evidence that a nested
       *  parameterisation's ComputeMaterial() genuinely used its parent
       *  touchable rather than ignoring it. Always false for non-nested. */
      G4bool fMaterialVariesWithParent = false;

      /** Set if any cost-control cap curtailed coverage of this volume. */
      G4bool fTruncated = false;
    };

    explicit G4ParameterisationCheckVisitor(G4bool verbose = true);
    ~G4ParameterisationCheckVisitor() override = default;

    /**
     * Snapshots physical- and logical-volume state, walks the subtree rooted
     * at @p startVol exercising every parameterisation encountered, restores
     * the snapshot, then repairs and/or reports. Master thread only (returns
     * 0 immediately on a worker thread).
     *  @param[in] startVol Root of the subtree to check; nullptr => world volume.
     *  @param[in] repair   If true and the geometry is still open,
     *             G4LogicalVolume::MarkMaterialMutable() is called on every
     *             newly-discovered non-mutable logical volume. Forced to
     *             false if the geometry is already closed (repair is then
     *             illegal); in that case a single G4Exception is raised
     *             naming every offender and the remedy.
     *  @param[in] verbose  Print the per-volume report (only for volumes with
     *             a nested parameterisation, or with at least one problem).
     *  @returns The number of problems found; 0 means clean.
     */
    static G4int CheckGeometry(G4VPhysicalVolume* startVol = nullptr, G4bool repair = true,
                               G4bool verbose = true);

    // -- Cost control (pass -1 to any setter to disable that cap) ----------

    /** Stop recursing into a parameterised volume's daughters once this many
     *  of its instances have been visited. Does not stop ComputeMaterial()
     *  from being exercised for every declared copy. Default: 64. */
    void SetMaxInstancesPerVolume(G4long n) { fMaxInstancesPerVolume = n; }

    /** For a replicated/parameterised mother, recurse into daughters only for
     *  the first N copies plus the last copy. Default: 3. */
    void SetMaxParentSamples(G4int n) { fMaxParentSamples = n; }

    /** Hard stop on the total number of logical-volume entries walked.
     *  Default: 100000. */
    void SetMaxTotalVisits(G4long n) { fMaxTotalVisits = n; }

    /** If true (default), a logical volume whose subtree contains no
     *  parameterised physical volume is not descended into. */
    void SetPruneUnparameterisedSubtrees(G4bool b) { fPruneUnparameterisedSubtrees = b; }

    // -- G4VTouchableVisitor overrides --------------------------------------

    void VisitParameterised(const G4TouchableHistoryHandle& touchable, G4PVParameterised* pv,
                            G4int replicaNo, G4int depth) override;

    void VisitRegularStructure(const G4TouchableHistoryHandle& touchable, G4PVParameterised* pv,
                               G4int nInstances, G4int depth) override;

    G4bool DescendInto(const G4TouchableHistoryHandle& touchable, G4LogicalVolume* lv,
                       G4int depth) override;

    // -- Setup and results ---------------------------------------------------

    /**
     * Records, for each parameterised physical volume, the logical volume it
     * was declared with before the walk began. Used to detect rebinding.
     * Called by CheckGeometry(); exposed for direct use of the visitor.
     */
    void SetDeclaredLogicalVolumes(std::map<G4VPhysicalVolume*, G4LogicalVolume*> declared);

    /** Returns the accumulated per-volume reports, keyed by physical volume. */
    const std::map<G4VPhysicalVolume*, VolumeReport>& GetReports() const { return fReports; }

    /** Returns the number of volumes for which a problem was found. */
    G4int CountProblems() const;

    /** True if any cost-control cap curtailed the walk. */
    G4bool IsTruncated() const { return fGlobalTruncated; }

    /** Prints the report: a statistics block for every volume with a nested
     *  parameterisation, and a summary line for every other volume that has
     *  at least one problem. */
    void PrintReport(std::ostream& os) const;

  private:

    G4bool HasProblem(const VolumeReport& rep) const;
    G4bool SubtreeHasParameterised(G4LogicalVolume* lv);
    VolumeReport& GetOrCreateReport(G4VPhysicalVolume* pv);

    G4bool fVerbose = true;

    G4long fMaxInstancesPerVolume = 64;
    G4int fMaxParentSamples = 3;
    G4long fMaxTotalVisits = 100000;
    G4bool fPruneUnparameterisedSubtrees = true;

    G4long fTotalVisits = 0;
    G4bool fGlobalTruncated = false;

    std::map<G4VPhysicalVolume*, G4LogicalVolume*> fDeclaredLV;
    std::map<G4VPhysicalVolume*, VolumeReport> fReports;

    /** Per parameterised PV: replica number -> materials seen across all
     *  distinct parent contexts encountered for that replica number. Used to
     *  compute VolumeReport::fMaterialVariesWithParent. */
    std::map<G4VPhysicalVolume*, std::map<G4int, std::set<G4Material*>>> fMaterialByReplica;

    /** Per parameterised PV: (parent physical volume, parent replica number)
     *  contexts seen, used to compute VolumeReport::fDistinctParentContexts. */
    std::map<G4VPhysicalVolume*, std::set<std::pair<G4VPhysicalVolume*, G4int>>> fParentContexts;

    /** Per nested parameterised PV: its declared material list, built once
     *  from G4VVolumeMaterialScanner on first use. */
    std::map<G4VPhysicalVolume*, std::set<G4Material*>> fDeclaredMaterials;

    /** Memoises whether a logical volume's subtree contains a parameterised
     *  physical volume, for SetPruneUnparameterisedSubtrees(). */
    std::map<G4LogicalVolume*, G4bool> fSubtreeHasParamCache;
};

#endif
