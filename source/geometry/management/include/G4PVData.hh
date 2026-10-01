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
// G4PVData
//
// Class description:
//
// Encapsulates the fields associated to G4VPhysicalVolume that are not
// read-only - they will change during simulation and must have a per-thread
// state under Geant4's split-class multithreading model (see
// G4GeomSplitter).

// Author: Xin Dong (Northeastern Univ.), 01.25.2009 
//         - Initial version (as part of G4VPhysicalVolume)
// --------------------------------------------------------------------
#ifndef G4PVDATA_HH
#define G4PVDATA_HH

#include "G4GeomSplitter.hh"
#include "G4RotationMatrix.hh"
#include "G4Types.hh"

/**
 * @brief G4PVData encapsulates the fields associated to G4VPhysicalVolume
 * that are not read-only - they will change during simulation and must have
 * a per-thread state.
 * @ingroup geometry_management
 */
class G4PVData
{
  public:

    G4PVData() = default;

    void initialize()
    {
      frot = nullptr;
      tx = 0.;
      ty = 0.;
      tz = 0.;
    }

    G4RotationMatrix* frot = nullptr;
    G4double tx = 0., ty = 0., tz = 0.;
};

/** Type defined for use of G4PVData objects. */
using G4PVManager = G4GeomSplitter<G4PVData>;

#endif
