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
// G4VTouchableVisitor
//
// Class description:
//
// Abstract base class for visitors that traverse the tree of Geant4
// touchables. Concrete subclasses override only the methods they need;
// all methods have empty default implementations.
//
// The traversal is driven by G4TouchableWalker. Depth is measured in
// physical-volume levels: the root passed to Walker::Walk() has depth 0,
// its daughters depth 1, and so on.
//
// For each physical volume the walker calls exactly one of:
//   VisitPlacement / VisitEachReplica / VisitEachParameterised
// before descending. For the associated logical volume it calls
//   EnterLogicalVolume before and LeaveLogicalVolume after the descent.
//
// Author: John Apostolakis (CERN), May-June 2026
// --------------------------------------------------------------------
#ifndef G4VTOUCHABLEVISITOR_HH
#define G4VTOUCHABLEVISITOR_HH

#include "G4TouchableHistoryHandle.hh"
#include "G4Types.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4PVReplica;
class G4PVParameterised;

/**
 * @brief G4VTouchableVisitor is an abstract base class for visitors that
 * traverses the true Geant4 tree of touchable volumes. Concrete subclasses
 * override only the methods they need.
 * For each physical volume the walker calls either:
 *   VisitPlacement / VisitReplica / VisitParameterised
 * once for each instance of that volume type (replicated volumes get called
 * multiple times, equal to the number of replicas or parameterised volumes)
 * before descending.
 * For its logical volume it calls first EnterLogicalVolume and later
 *   LeaveLogicalVolume (after the descent).
 * @ingroup geometry_navigation
 */

class G4VTouchableVisitor
{
  public:

    /*
     * Default constructor and virtual destructor.
     */
    G4VTouchableVisitor() = default;
    virtual ~G4VTouchableVisitor() = default;

    // -- Per touchable volume dispatch ------------------------------------
    // Called once per physical volume before descending into its children.
    // Replicas and parameterisations represent GetMultiplicity() instances
    // that share one logical volume; the walker descends that LV once.
    // The touchable handle encapsulates the full navigation history at the
    // moment of the call; it remains valid beyond the scope of the callback.

    /*
     * Called once for each placement in mother logical volume
     *   before descending into its children.
     */
    virtual void VisitPlacement(const G4TouchableHistoryHandle& /*touchable*/,
                                G4PVPlacement* /*pv*/, G4int /*depth*/)
    {}

    /*
     * Called once for each replica id in mother logical volume
     *   before descending into its children.
     */
    virtual void VisitReplica(const G4TouchableHistoryHandle& /*touchable*/, G4PVReplica* /*pv*/,
                              G4int /*replicaNo*/, G4int /*depth*/)
    {}
    /*
     * Called once for each parameterisation's instance in mother logical
     *   volume before descending into its children.
     */
    virtual void VisitParameterised(const G4TouchableHistoryHandle& /*touchable*/,
                                    G4PVParameterised* /*pv*/, G4int /*replicaNo*/, G4int /*depth*/)
    {}

    /*
     * Called once for (mother) volume types not covered above (kExternal).
     */
    virtual void VisitUnknown(const G4TouchableHistoryHandle& /*touchable*/,
                              G4VPhysicalVolume* /*pv*/, G4int /*depth*/)
    {}

    // -- Logical volume bracketing ---------------------------------------
    // EnterLogicalVolume is called at the same depth as its owning PV,
    // immediately before visiting the LV's daughters.
    // LeaveLogicalVolume is called after all daughters have been visited.

    /*
     * Called once for each logical volume at the same depth as its owning PV,
     * immediately before visiting the LV's daughters.
     */
    virtual void EnterLogicalVolume(const G4TouchableHistoryHandle& /*touchable*/, G4LogicalVolume*,
                                    G4int /* depth */)
    {}

    /*
     * Called once for each logical volume after all daughters have been visited.
     */
    virtual void LeaveLogicalVolume(const G4TouchableHistoryHandle& /*touchable*/, G4LogicalVolume*,
                                    G4int /* depth */)
    {}
};

#endif
