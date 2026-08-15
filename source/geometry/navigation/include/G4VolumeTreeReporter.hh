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
// G4VolumeTreeReporter
//
// Class description:
//
// Concrete G4VVolumeVisitor that prints an indented summary of the
// physical-volume tree to G4cout.  Example output (2 spaces per level):
//
//   [PV] World  copy=0  type=Placement
//     [LV] World_LV  solid=G4Box  material=G4_AIR
//       [PV] Detector  copy=0  type=Placement
//         [LV] Detector_LV  solid=G4Tubs  material=G4_Si
//       [PV] Layers  copy=-1  type=Replica  n=10  axis=kZAxis
//         [LV] Layer_LV  solid=G4Box  material=G4_Si
//
// Usage:
//   G4VolumeTreeReporter printer;
//   G4VolumeTreeWalker walker(printer);
//   walker.Walk(worldPV);
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------
#ifndef G4VOLUMETREEREPORTER_HH
#define G4VOLUMETREEREPORTER_HH

#include "G4VVolumeVisitor.hh"

#include "geomdefs.hh"

/**
 * @brief G4VolumeTreeReporter is concrete G4VVolumeVisitor that prints
 * an indented summary of the physical-volume tree to G4cout.
 * @ingroup geometry_navigation
 */

class G4VolumeTreeReporter : public G4VVolumeVisitor
{
  public:

    /*
     * Constructor. Sets the number of spaces per depth level.
     *  @param[in] indentSize  Number of spaces per depth level (default 2).
     */
    explicit G4VolumeTreeReporter(G4int indentSize = 2);

    /*
     * Default virtual destructor.
     */
    ~G4VolumeTreeReporter() override = default;

    /*
     * Visit methods.
     */
    void VisitPlacement(G4PVPlacement* pv, G4int depth) override;
    void VisitReplica(G4PVReplica* pv, G4int depth) override;
    void VisitParameterised(G4PVParameterised* pv, G4int depth) override;
    void VisitUnknown(G4VPhysicalVolume* pv, G4int depth) override;

    void EnterLogicalVolume(G4LogicalVolume* lv, G4int depth) override;

  private:

    /** Returns a string of fIndentSize*depth spaces. */
    G4String Indent(G4int depth) const;

    /** Formats the replication axis name for output. */
    static G4String AxisName(EAxis axis);

    /** Number of spaces per depth level. */
    G4int fIndentSize;
};

#endif
