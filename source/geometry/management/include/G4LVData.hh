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
// G4LVData
//
// Class description:
//
// Encapsulates the fields associated to G4LogicalVolume that are always
// per-thread under Geant4's split-class multithreading model (see
// G4GeomSplitter): the sensitive detector and field manager.
//
// The solid, material, mass and cuts-couple are NOT here: they live either
// as ordinary (shared) members of G4LogicalVolume itself, or in a per-thread
// G4LVMutableData slot, decided at runtime per logical-volume instance (see
// G4LVMutableData.hh and G4LogicalVolume::MarkMaterialMutable()).

// Author: Xin Dong (Northeastern Univ.), 01.25.2009
//          - Initial version (as part of G4LogicalVolume)
// ------------------------------------------------------------------------
#ifndef G4LVDATA_HH
#define G4LVDATA_HH

#include "G4GeomSplitter.hh"
#include "G4Types.hh"

class G4VSensitiveDetector;
class G4FieldManager;

/**
 * @brief G4LVData encapsulates the fields associated to the class
 * G4LogicalVolume that are always per-thread.
 * @ingroup geometry_management
 */
class G4LVData
{
  public:

    G4LVData() = default;

    void initialize()
    {
      fSensitiveDetector = nullptr;
      fFieldManager = nullptr;
    }

  public:

    G4VSensitiveDetector* fSensitiveDetector = nullptr;
    // Pointer to sensitive detector.
    G4FieldManager* fFieldManager = nullptr;
    // Pointer (possibly nullptr) to (magnetic or other) field manager object.
};

/** G4LVManager encapsulates the methods used by both the master thread and
    worker threads to allocate memory space for the fields encapsulated by
    the class G4LVData. */
using G4LVManager = G4GeomSplitter<G4LVData>;

#endif
