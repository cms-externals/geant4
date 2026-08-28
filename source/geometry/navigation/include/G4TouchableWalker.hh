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
// G4TouchableWalker
//
// Class description:
//
// Performs a depth-first traversal of a Geant4 physical-volume subtree,
// visiting each physical instance (touchable) of the tree, and
// dispatching to a G4VTouchableVisitor at each node (touchable).
//
// Usage:
//   G4TouchableCollector     touchableCollector;
//   G4TouchableWalker walker(touchableCollector);
//   walker.Walk(startPhysicalVolume);
//
// Traversal order for each node:
//   1. VisitPlacement / VisitReplica / VisitParameterised  (physical vol)
//   2. EnterLogicalVolume                                   (logical vol)
//   3. Recurse into daughters (depth-first)
//   4. LeaveLogicalVolume
//
// In most uses we expect the Walk to start at the world volume.  This is
//   the default behaviour if no volume is provided.
//
// It is possible to walk a subtree anchored by a different 'start' volume
// - The simpler case is that the starting volume is a placement volume.
//
// - It can be a replica or parameterised volume which is the direct
//   daughter of a placement volume.
//   In this case an artificial PVPlacement volume is constructed
//   to hold it and this is registered in G4PhysicalVolumeStore.
//   The walker destructor then removes it from the store.
//   This creation and destruction of PVs affect any operations with/on
//   the store.
//
// Author: John Apostolakis (CERN), May 2026
// --------------------------------------------------------------------
#ifndef G4TOUCHABLEWALKER_HH
#define G4TOUCHABLEWALKER_HH

#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchableVisitor.hh"

#include "geomdefs.hh"

class G4PVPlacement;
class G4PVReplica;
class G4PVParameterised;

/**
 * @brief Performs a depth-first traversal of a Geant4 physical-volume
 * subtree, visiting every touchable instance and dispatching to a
 * G4VTouchableVisitor at each node. Replicas and parameterised volumes are
 * expanded: every copy is visited individually. If a volume has daughters,
 * every daughter of every copy is visited recursively.
 * @ingroup geometry_navigation
 */

class G4TouchableWalker
{
  public:

    /**
     * Constructor. Registers @p visitor as the callback object for this walk.
     *  @param[in] visitor  The visitor to invoke at each node.
     */
    explicit G4TouchableWalker(G4VTouchableVisitor& visitor);

    /**
     * Destructor. Releases any artificial placement volume created during
     * a walk that started from a replica or parameterised volume.
     */
    virtual ~G4TouchableWalker();

    /**
     * Performs a depth-first traversal of the subtree rooted at @p startVol,
     * invoking visitor callbacks at every node.
     *  @param[in] startVol  Root of the subtree to walk. If nullptr the world
     *             volume returned by G4TransportationManager is used. For a
     *             non-world kReplica or kParameterised start volume an artificial
     *             placement (whose logical volume is the replica's mother) is
     *             created internally so that all replica instances are visited.
     *  @note Geometry must be closed before calling Walk().
     */
    void Walk(G4VPhysicalVolume* startVol = nullptr);

    /**
     * Returns a reference-counted G4TouchableHistoryHandle encapsulating a
     * snapshot of the current navigation history. Intended to be called by
     * the visitor inside a traversal callback.
     *  @returns Handle wrapping a newly allocated G4TouchableHistory.
     */
    G4TouchableHistoryHandle GetTouchableHistoryHandle() const
    {
      return G4TouchableHistoryHandle(new G4TouchableHistory(fNavHistory));
    }

    /**
     * Called when a physical volume of an unrecognised type (e.g. kExternal)
     * is encountered. The default implementation emits a diagnostic when
     * verbose mode is on. Subclasses may override this to handle external types.
     *  @param[in] pv     The physical volume of unknown type.
     *  @param[in] depth  Depth of @p pv in the traversal tree (root = 0).
     */
    virtual void DescendUnknown(G4VPhysicalVolume* pv, G4int depth);

  protected:

    /**
     * Hook called after all daughters of a logical volume have been visited.
     * Currently empty; provided for potential subclass use.
     *  @param[in] lv     The logical volume being left.
     *  @param[in] depth  Depth of the owning physical volume in the tree.
     */
    virtual void LeaveLogicalVolume(G4LogicalVolume* lv, G4int depth);

    /**
     * Returns a string of spaces proportional to @p depth, for use in
     * formatted diagnostic output.
     *  @param[in] depth  Indentation level.
     *  @returns  Indentation string (2 * @p depth spaces).
     */
    G4String Indent(G4int depth) const;

  private:

    /**
     * Calls EnterLogicalVolume on the visitor, recurses into all daughters of
     * @p lv (dispatching to DescendPlacement, DescendAllReplicas, or
     * DescendAllParameterised as appropriate), then calls LeaveLogicalVolume.
     *  @param[in] lv     The logical volume to enter.
     *  @param[in] depth  Depth of the owning physical volume in the tree.
     */
    void DescendLogicalVolume(G4LogicalVolume* lv, G4int depth);

    /**
     * Pops one level from the navigation history via BackLevel(). Called after
     * the full subtree rooted at a non-root physical volume has been visited,
     * to restore the history to the state before that volume was entered.
     *  @param[in] pv     The volume being departed (used only for diagnostics).
     *  @param[in] depth  Depth of @p pv (used only for diagnostics).
     */
    void DepartPhysicalVolume(G4VPhysicalVolume* pv, G4int depth);

    /**
     * Returns the world volume from G4TransportationManager. If @p fatalOnFailure
     * is true a fatal G4Exception is thrown when the world is not found;
     * otherwise a warning is issued and nullptr is returned.
     *  @param[in] fatalOnFailure  Treat a missing world volume as a fatal error.
     *  @returns  Pointer to the world G4VPhysicalVolume, or nullptr.
     */
    G4VPhysicalVolume* GetWorldVolume(G4bool fatalOnFailure = false);

    /**
     * Pushes @p pv onto the navigation history, invokes VisitPlacement on the
     * visitor, descends into the logical volume, then pops @p pv from the history.
     *  @param[in] pv     The placement volume to descend into.
     *  @param[in] depth  Depth of @p pv in the traversal tree.
     */
    void DescendPlacement(G4PVPlacement* pv, G4int depth);

    /**
     * Pushes replica copy @p replicaNo of @p pv onto the navigation history,
     * invokes VisitReplica on the visitor, descends into the logical volume,
     * then pops the level. Visits a single replica instance.
     *  @param[in] pv        The replica volume.
     *  @param[in] replicaNo Copy number of the specific instance to visit.
     *  @param[in] depth     Depth of @p pv in the traversal tree.
     */
    void DescendReplica(G4PVReplica* pv, G4int replicaNo, G4int depth);

    /**
     * Computes the solid, dimensions, transformation, and material for
     * parameterised copy @p replicaNo of @p pv, pushes it onto the navigation
     * history, invokes VisitParameterised on the visitor, descends into the
     * logical volume, then pops the level. Visits a single parameterised instance.
     *  @param[in] pv        The parameterised volume.
     *  @param[in] replicaNo Copy number of the specific instance to visit.
     *  @param[in] depth     Depth of @p pv in the traversal tree.
     */
    void DescendParameterised(G4PVParameterised* pv, G4int replicaNo, G4int depth);

    /**
     * Iterates over all replica copies of @p pv by calling DescendReplica for
     * each copy number in [0, nReplicas). Called from DescendLogicalVolume when
     * the sole daughter of a logical volume is a replica.
     *  @param[in] pv     The replica volume whose instances are to be visited.
     *  @param[in] depth  Depth of @p pv in the traversal tree.
     */
    void DescendAllReplicas(G4PVReplica* pv, G4int depth);

    /**
     * Iterates over all parameterised copies of @p pv by calling
     * DescendParameterised for each copy number in [0, nReplicas). Called from
     * DescendLogicalVolume when the sole daughter of a logical volume is a
     * parameterised volume.
     *  @param[in] pv     The parameterised volume whose instances are to be visited.
     *  @param[in] depth  Depth of @p pv in the traversal tree.
     */
    void DescendAllParameterised(G4PVParameterised* pv, G4int depth);

    /**
     * Handles the root physical volume of the walk. Dispatches to
     * DescendPlacement, DescendReplica, or DescendParameterised as appropriate.
     * For replica and parameterised roots, Walk() will have already primed the
     * navigation history with an artificial placement and @p pv will be that
     * placement.
     *  @param[in] pv     The root physical volume of the walk.
     *  @param[in] depth  Depth to assign to the root (typically 0).
     */
    void DescendFirstPhysical(G4VPhysicalVolume* pv, G4int depth);

  private:

    /** Enables verbose diagnostic output when true. */
    G4bool fVerbose = false;

    /** The visitor whose callbacks are invoked at each node during traversal. */
    G4VTouchableVisitor& fVisitor;

    /** Navigation history maintained and updated throughout the traversal. */
    G4NavigationHistory fNavHistory;

    /** Converts an EAxis enumeration value to a human-readable string. */
    static G4String AxisName(EAxis axis);

    /**
     * Artificial placement created when Walk() is called with a replica or
     * parameterised start volume. Its logical volume is the start volume's
     * mother, providing a kNormal root from which the normal traversal logic
     * can proceed. Owned by this walker; deleted in the destructor.
     */
    std::unique_ptr<G4PVPlacement> fArtificialWorldPV;
};

#endif
