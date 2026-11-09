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
// G4LVMutableData
//
// Class description:
//
// Per-thread storage for the solid, material, mass and cuts-couple of a
// G4LogicalVolume, allocated only for the (typically small) subset of
// logical volumes that are actually used as the daughter of a
// G4PVParameterised, G4PVDivision or G4ReplicatedSlice (i.e. whose
// VolumeType() is kParameterised). This decision is made entirely at
// runtime, per logical-volume instance, requiring no build-time
// configuration - see G4LogicalVolume::MarkMaterialMutable().
//
// A logical volume that is never used inside a parameterised/division
// volume instead keeps these fields shared (its own plain members,
// identical to all threads). See also G4LVData.hh, which retains only the
// fields that are always per-thread (sensitive detector, field manager).
//
// Author: John Apostolakis (CERN), 20.07.2026 - Initial version
// ------------------------------------------------------------------------
#ifndef G4LVMUTABLEDATA_HH
#define G4LVMUTABLEDATA_HH

#include "G4GeomSplitter.hh"
#include "G4Types.hh"

class G4VSolid;
class G4Material;
class G4MaterialCutsCouple;

/**
 * @brief G4LVMutableData encapsulates the fields of G4LogicalVolume that
 * must be per-thread only for the subset of logical volumes used inside a
 * parameterised/division-type placement.
 * @ingroup geometry_management
 */
class G4LVMutableData
{
  public:

    G4LVMutableData() = default;

    void initialize()
    {
      fSolid = nullptr;
      fMaterial = nullptr;
      fMass = 0.0;
      fCutsCouple = nullptr;
    }

  public:

    G4VSolid* fSolid = nullptr;
    // Pointer to solid.
    G4Material* fMaterial = nullptr;
    // Pointer to material at this node.
    G4double fMass = 0.0;
    // Mass of the logical volume tree. Plain (not atomic): each thread
    // only ever reads/writes its own slot.
    G4MaterialCutsCouple* fCutsCouple = nullptr;
    // Pointer (possibly nullptr) to associated production cuts.
};

/** G4LVMutableManager encapsulates the methods used by both the master
    thread and worker threads to allocate memory space for the fields
    encapsulated by the class G4LVMutableData. */
using G4LVMutableManager = G4GeomSplitter<G4LVMutableData>;

#endif
