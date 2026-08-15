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
// G4VolumeTreeWalker
//
// Class description:
//
// Performs a depth-first traversal of a Geant4 physical-volume subtree,
// dispatching to a G4VVolumeVisitor at each node.
//
// Usage:
//   MyVisitor visitor;
//   G4VolumeTreeWalker walker(visitor);
//   walker.Walk(worldPhysVol);
//
// Traversal order for each node:
//   1. VisitPlacement / VisitReplica / VisitParameterised  (physical vol)
//   2. EnterLogicalVolume                                   (logical vol)
//   3. Recurse into daughters (depth-first)
//   4. LeaveLogicalVolume
//
// For replicated and parameterised volumes GetMultiplicity() > 1 copies
// share a single logical volume; the walker descends that logical volume
// exactly once, visiting the common daughter structure.
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------
#ifndef G4VOLUMETREEWALKER_HH
#define G4VOLUMETREEWALKER_HH

#include "G4Types.hh"

class G4VVolumeVisitor;
class G4VPhysicalVolume;
class G4LogicalVolume;

class G4PVPlacement;
class G4PVParameterised;
class G4PVReplica;

/*
 * @brief G4VolumeTreeWalker performs a depth-first traversal of a Geant4
 * physical-volume subtree, dispatching to a G4VVolumeVisitor at each node.
 * @ingroup geometry_navigation
 */

class G4VolumeTreeWalker
{
  public:

    /*
     * Constructor and default destructor.
     */
    explicit G4VolumeTreeWalker(G4VVolumeVisitor& visitor);
    ~G4VolumeTreeWalker() = default;

    /*
     * Walks the subtree rooted at startVol (inclusive).
     *  @note  startVol != nullptr.
     */
    void Walk(G4VPhysicalVolume* startVol);

    /*
     * Sets the verbosity level.
     */
    inline G4bool SetVerbose(G4bool val)
    {
      G4bool oldVal = val;
      fVerbose = val;
      return oldVal;
    }

  private:

    void DescendPhysical(G4VPhysicalVolume* pv, G4int depth);
    void DescendLogical(G4LogicalVolume* lv, G4int depth);
    void DescendPlacement(G4PVPlacement* pv, G4int depth);
    void DescendAllParameterised(G4PVParameterised* pv, G4int depth);
    void DescendAllReplicas(G4PVReplica* pv, G4int depth);

  private:

    G4bool fVerbose = false;
    G4VVolumeVisitor& fVisitor;
};

#endif
