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
// G4ParameterisationChecker
//
// Class description:
//
// Optional geometry-close consistency pass for an 'unorthodox' use of
// parameterisation: a G4VPVParameterisation whose ComputeTransformation()
// (or related hooks) reassigns the physical volume's logical volume per copy
// number, via G4VPhysicalVolume::SetLogicalVolume(), rather than varying a
// single logical volume.
//
// This checker discovers those extra logical volumes by
// Exercises each parameterisation over all of its copies - calling
// ComputeTransformation(), ComputeSolid() and ComputeMaterial() 
// and marks any that are not yet mutable. It must run on the master thread,
// while the geometry is still open and before voxelisation.
//
// ComputeMaterial() is not probed for *nested* parameterisations
//
// Author: John Apostolakis (CERN), July 2026
// --------------------------------------------------------------------
#ifndef G4PARAMETERISATIONCHECKER_HH
#define G4PARAMETERISATIONCHECKER_HH

#include "G4Types.hh"

/**
 * @brief Detects (and repairs) parameterisations that reassign the logical
 * volume of their physical volume to volumes not marked mutable, and reports.
 * @ingroup geometry_management
 */
class G4ParameterisationChecker
{
  public:

    /**
     * Walks all parameterised physical volumes, exercising each
     * parameterisation's ComputeTransformation(), ComputeSolid() and
     * ComputeMaterial() (the latter skipped for nested parameterisations)
     * over its copies to discover every logical volume it binds, and calls
     * G4LogicalVolume::MarkMaterialMutable() on any that is not yet mutable.
     * Must be called on the master thread while the geometry is still open.
     *  @param[in] verbose If true (default), prints one line per repaired
     *             logical volume, naming it and the parameterised physical
     *             volume that reassigns it. The count of volumes marked is
     *             reported regardless of this flag.
     *  @returns the number of logical volumes newly marked mutable.
     */
    static G4int CheckAndMark(G4bool verbose = true);
};

#endif
