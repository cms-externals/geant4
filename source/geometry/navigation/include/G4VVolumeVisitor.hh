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
// G4VVolumeVisitor
//
// Class description:
//
// Abstract base class for visitors that traverse the Geant4 physical
// volume tree. Concrete subclasses override only the methods they need;
// all methods have empty default implementations.
//
// The traversal is driven by G4VolumeTreeWalker. Depth is measured in
// physical-volume levels: the root passed to Walker::Walk() has depth 0,
// its daughters depth 1, and so on.
//
// For each physical volume the walker calls exactly one type of method
// either VisitPlacement / VisitReplica / VisitParameterised
// for each volume it contains ( each placement, each replica number,
// each instance of theparameterised volume ) before descending.
//  For the associated logical volume it calls
//    - EnterLogicalVolume before and
//    -LeaveLogicalVolume after the descent.
//
// Author: John Apostolakis (CERN), May 2026
// --------------------------------------------------------------------
#ifndef G4VVOLUMEVISITOR_HH
#define G4VVOLUMEVISITOR_HH

#include "G4Types.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVPlacement;
class G4PVReplica;
class G4PVParameterised;

/**
 * @brief G4VVolumeVisitor is an abstract base class for visitors that
 * traverses the Geant4 physical volume tree. Concrete subclasses override
 * only the methods they need.
 * For each physical volume the walker calls either:
 *   VisitPlacement / VisitReplica / VisitParameterised
 * once for each instance of that volume type (replicated volumes get called
 * only one time - independent of the number of replicas or parameterised volumes)
 * before descending.
 * For its logical volume it calls first EnterLogicalVolume and later
 *   LeaveLogicalVolume (after the descent).
 * @ingroup geometry_navigation
 */

class G4VVolumeVisitor
{
  public:

    /*
     * Default constructor and virtual destructor.
     */
    G4VVolumeVisitor() = default;
    virtual ~G4VVolumeVisitor() = default;

    // -- Physical volume dispatch ----------------------------------------
    // Called once per physical volume before descending into its children.
    // Replicas and parameterisations represent GetMultiplicity() instances
    // that share one logical volume; the walker descends that LV once.

    /*
     * Called once for each placement in mother logical volume
     *   before descending into its children.
     */
    virtual void VisitPlacement(G4PVPlacement*, G4int /*depth*/) {}

    /*
     * Called once for the whole replica before descending into its children.
     */
    virtual void VisitReplica(G4PVReplica*, G4int /*depth*/) {}

    /*
     * Called once for the whole parameterisation, before descending into its children.
     */
    virtual void VisitParameterised(G4PVParameterised*, G4int /*depth*/) {}

    /*
     * Called once for (physical) volume types not covered above (kExternal).
     */
    virtual void VisitUnknown(G4VPhysicalVolume*, G4int /*depth*/) {}

    /*
     * Called once for each placement in mother logical
     *   volume after all children have been visited.
     */
    virtual void DepartPlacement(G4PVPlacement*, G4int /*depth*/) {}

    /*
     * Called once for the whole replica/parameterisation in mother logical
     *   volume after all children have been visited.
     */
    virtual void DepartReplica(G4PVReplica*, G4int /*depth*/) {}
    virtual void DepartParameterised(G4PVParameterised*, G4int /*depth*/) {}

    /*
     * Called once for an external-defined type of physical volumes in mother
     *   logical volume after all children have been visited.
     */
    virtual void DepartUnknown(G4VPhysicalVolume*, G4int /*depth*/) {}

    // -- Logical volume bracketing ---------------------------------------

    /*
     * Called once for each logical volume at the same depth as its owning PV,
     * immediately before visiting the LV's daughters.
     */
    virtual void EnterLogicalVolume(G4LogicalVolume*, G4int /* depth */) {}

    /*
     * Called once for each logical volume after all daughters have been visited.
     */
    virtual void LeaveLogicalVolume(G4LogicalVolume*, G4int /* depth */) {}
};

#endif
