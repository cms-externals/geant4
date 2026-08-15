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
// G4TouchableCollector
//
// Class description:
//
// Concrete G4VTouchableVisitor that generates a list of
// G4TouchableHistoryHandle objects, one for each instance of a selected
// logical or physical volume found during a G4TouchableWalker traversal.
//
// Author: John Apostolakis (CERN), May-June 2026
// --------------------------------------------------------------------
#ifndef G4TOUCHABLECOLLECTOR_HH
#define G4TOUCHABLECOLLECTOR_HH

#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4VTouchableVisitor.hh"

#include "geomdefs.hh"

class G4TouchableWalker;

/**
 * @brief G4TouchableCollector is a concrete G4VTouchableVisitor that records
 * a G4TouchableHistoryHandle for every instance of a specified target volume
 * found during a G4TouchableWalker traversal.
 * The target may be either a G4LogicalVolume (matched in EnterLogicalVolume)
 * or a G4VPhysicalVolume (matched in VisitPlacement / VisitReplica /
 * VisitParameterised). All instances, including every replica copy and every
 * parameterised copy, are recorded as separate entries.
 *
 * Usage:
 * @code
 *   G4TouchableCollector collector(targetLogicalVolume);
 *   G4TouchableWalker    walker(collector);
 *   walker.Walk();  // walks from the world volume
 *   const auto& list = collector.GetListOfMatchingTouchables();
 * @endcode
 *
 * @note Walk() may be called repeatedly with different start volumes; results
 *       accumulate. Call ClearList() before re-walking to start fresh.
 * @note Touchable handles are only valid for navigation when the walk was
 *       started from the world volume.
 * @ingroup geometry_navigation
 */

class G4TouchableCollector : public G4VTouchableVisitor
{
  public:

    /**
     * Constructor. Collects a touchable handle for every instance of the
     * logical volume @p lv encountered during the walk.
     *  @param[in] lv  Target logical volume whose instances are to be collected.
     */
    explicit G4TouchableCollector(const G4LogicalVolume* lv);

    /**
     * Constructor. Collects a touchable handle for every occurrence of the
     * physical volume @p targetPV encountered during the walk (including each
     * replica or parameterised copy).
     *  @param[in] targetPV  Target physical volume to search for.
     */
    explicit G4TouchableCollector(const G4VPhysicalVolume* targetPV);

    ~G4TouchableCollector() override = default;

    /**
     * Called by the walker for each placement physical volume visited.
     * Records a touchable handle if @p pv matches the target physical volume.
     *  @param[in] touchable Navigation history at the moment of the call.
     *  @param[in] pv     The placement volume being visited.
     *  @param[in] depth  Depth of @p pv in the traversal tree.
     */
    void VisitPlacement(const G4TouchableHistoryHandle& touchable, G4PVPlacement* pv,
                        G4int depth) override;

    /**
     * Called by the walker for each replica copy visited.
     * Records a touchable handle if @p pv matches the target physical volume.
     *  @param[in] touchable Navigation history at the moment of the call.
     *  @param[in] pv        The replica volume being visited.
     *  @param[in] replicaNo Copy number of the current replica instance.
     *  @param[in] depth     Depth of @p pv in the traversal tree.
     */
    void VisitReplica(const G4TouchableHistoryHandle& touchable, G4PVReplica* pv, G4int replicaNo,
                      G4int depth) override;

    /**
     * Called by the walker for each parameterised copy visited.
     * Records a touchable handle if @p pv matches the target physical volume.
     *  @param[in] touchable Navigation history at the moment of the call.
     *  @param[in] pv        The parameterised volume being visited.
     *  @param[in] replicaNo Copy number of the current parameterised instance.
     *  @param[in] depth     Depth of @p pv in the traversal tree.
     */
    void VisitParameterised(const G4TouchableHistoryHandle& touchable, G4PVParameterised* pv,
                            G4int replicaNo, G4int depth) override;

    /**
     * Called by the walker for physical volumes of unrecognised type.
     * Records a touchable handle if @p pv matches the target physical volume.
     *  @param[in] touchable Navigation history at the moment of the call.
     *  @param[in] pv     The physical volume of unknown type.
     *  @param[in] depth  Depth of @p pv in the traversal tree.
     */
    void VisitUnknown(const G4TouchableHistoryHandle& touchable, G4VPhysicalVolume* pv,
                      G4int depth) override;

    /**
     * Called by the walker on entering a logical volume.
     * Records a touchable handle if @p lv matches the target logical volume.
     *  @param[in] touchable Navigation history at the moment of the call.
     *  @param[in] lv     The logical volume being entered.
     *  @param[in] depth  Depth of the owning physical volume in the tree.
     */
    void EnterLogicalVolume(const G4TouchableHistoryHandle& touchable, G4LogicalVolume* lv,
                            G4int depth) override;

    /**
     * Returns the list of touchable handles collected so far, one entry per
     * matched volume instance.
     *  @returns Const reference to the internal vector of handles.
     */
    const std::vector<G4TouchableHistoryHandle>& GetListOfMatchingTouchables() const
    {
      return fListOfTouchableHandles;
    }

    /**
     * Clears the collected list of touchable handles. Call this before
     * re-walking to avoid accumulating results from previous walks.
     */
    void ClearList() { fListOfTouchableHandles.clear(); }

  protected:

    /**
     * Returns a string of spaces proportional to @p depth, for use in
     * formatted diagnostic output.
     *  @param[in] depth  Indentation level.
     *  @returns  Indentation string (2 * @p depth spaces).
     */
    G4String Indent(G4int depth) const;

  private:

    /** Target logical volume; a touchable is recorded when this LV is entered. */
    const G4LogicalVolume* fTargetLogVol = nullptr;

    /** Target physical volume; a touchable is recorded at each matching Visit call. */
    const G4VPhysicalVolume* fTargetPhysVol = nullptr;

    /** Enables verbose diagnostic output when true. */
    G4bool fVerbose = false;

    /** Accumulated list of touchable handles - the output of the collection. */
    std::vector<G4TouchableHistoryHandle> fListOfTouchableHandles;
};

#endif
