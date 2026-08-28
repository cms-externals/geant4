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
// G4VolumeTreeReporter class Implementation
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------

#include "G4VolumeTreeReporter.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ios.hh"

#include "geomdefs.hh"

// ---------------------------------------------------------------------------

G4VolumeTreeReporter::G4VolumeTreeReporter(G4int indentSize) : fIndentSize(indentSize) {}

// ---------------------------------------------------------------------------

void G4VolumeTreeReporter::VisitPlacement(G4PVPlacement* pv, G4int depth)
{
  G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
         << "  type=Placement" << G4endl;
}

// ---------------------------------------------------------------------------

void G4VolumeTreeReporter::VisitReplica(G4PVReplica* pv, G4int depth)
{
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

  G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
         << "  type=Replica"
         << "  n=" << nReplicas << "  axis=" << AxisName(axis) << "  width=" << width << G4endl;
}

// ---------------------------------------------------------------------------

void G4VolumeTreeReporter::VisitParameterised(G4PVParameterised* pv, G4int depth)
{
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

  G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
         << "  type=Parameterised"
         << "  n=" << nReplicas << "  axis=" << AxisName(axis) << G4endl;
}

// ---------------------------------------------------------------------------

void G4VolumeTreeReporter::VisitUnknown(G4VPhysicalVolume* pv, G4int depth)
{
  G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
         << "  type=Unknown(kExternal?)" << G4endl;
}

// ---------------------------------------------------------------------------

void G4VolumeTreeReporter::EnterLogicalVolume(G4LogicalVolume* lv, G4int depth)
{
  const G4String solidName = lv->GetSolid() ? lv->GetSolid()->GetName() : G4String("(none)");
  const G4String materialName =
    lv->GetMaterial() ? lv->GetMaterial()->GetName() : G4String("(none)");

  G4cout << Indent(depth) << "  [LV] " << lv->GetName() << "  solid=" << solidName
         << "  material=" << materialName << G4endl;
}

// ---------------------------------------------------------------------------

G4String G4VolumeTreeReporter::Indent(G4int depth) const
{
  return G4String(static_cast<std::size_t>(fIndentSize * depth), ' ');
}

// ---------------------------------------------------------------------------

G4String G4VolumeTreeReporter::AxisName(EAxis axis)
{
  switch (axis)
  {
    case kXAxis:
      return "kXAxis";
    case kYAxis:
      return "kYAxis";
    case kZAxis:
      return "kZAxis";
    case kRho:
      return "kRho";
    case kRadial3D:
      return "kRadial3D";
    case kPhi:
      return "kPhi";
    default:
      return "kUndefined";
  }
}
