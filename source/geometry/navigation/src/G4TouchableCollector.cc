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
// G4TouchableCollector class Implementation
//
// Author: John Apostolakis (CERN), May-June 2026
// --------------------------------------------------------------------

#include "G4TouchableCollector.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NavigationHistory.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ios.hh"

#include "geomdefs.hh"

// ---------------------------------------------------------------------------

G4TouchableCollector::G4TouchableCollector(const G4LogicalVolume* targetLogVol)
  : G4VTouchableVisitor(), fTargetLogVol(targetLogVol), fTargetPhysVol(nullptr)
{
  if (fVerbose)
  {
    G4cout << "* G4TouchableCollector: targetLV= " << targetLogVol
           << "  and its name= " << targetLogVol->GetName() << G4endl;
  }
}

// ---------------------------------------------------------------------------

G4TouchableCollector::G4TouchableCollector(const G4VPhysicalVolume* targetPhysVol)
  : G4VTouchableVisitor(), fTargetLogVol(nullptr), fTargetPhysVol(targetPhysVol)
{
  if (fVerbose)
  {
    G4cout << "* G4TouchableCollector: targetPhysical= " << targetPhysVol
           << "  and its name= " << targetPhysVol->GetName() << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableCollector::VisitPlacement(const G4TouchableHistoryHandle& touchable,
                                          G4PVPlacement* pv, G4int depth)
{
  if (pv == fTargetPhysVol)
  {
    fListOfTouchableHandles.push_back(touchable);
    if (fVerbose)
    {
      G4cout << " **> Added new entry (Placement) for Target Phys-Vol in List of Touchables."
             << G4endl;
    }
  }

  if (fVerbose)
  {
    G4VPhysicalVolume* topPV = touchable->GetVolume();
    G4cout << "  TouchableCollector:VisitPlacement level= " << touchable->GetHistoryDepth()
           << " info:  ptr= " << topPV << " = " << topPV->GetName()
           << " copyNo= " << touchable->GetCopyNumber(0) << " -- at depth = " << depth << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableCollector::VisitReplica(const G4TouchableHistoryHandle& touchable, G4PVReplica* pv,
                                        G4int replicaNo, G4int depth)
{
  if (pv == fTargetPhysVol)
  {
    fListOfTouchableHandles.push_back(touchable);
    if (fVerbose)
    {
      G4cout << " **> Added new entry (Replica) for Target Phys-Vol in List of Touchables."
             << G4endl;
    }
  }

  if (fVerbose)
  {
    EAxis axis;
    G4int nReplicas;
    G4double width, offset;
    G4bool consuming;
    pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

    G4cout << "TouchCollect: "  // << Indent(depth)
           << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo() << "  type=Replica"
           << "  repNum= " << replicaNo << " (FORCED) "
           << "  ( numInstances=" << nReplicas << " ) "
           << "  axis=" << axis  // << AxisName(axis)
           << "  width=" << width << " -- at depth= " << depth << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableCollector::VisitParameterised(const G4TouchableHistoryHandle& touchable,
                                              G4PVParameterised* pv, G4int replicaNo, G4int depth)
{
  if (pv == fTargetPhysVol)
  {
    fListOfTouchableHandles.push_back(touchable);
    if (fVerbose)
    {
      G4cout << " **> Added new entry (Parameterised) for Target Phys-Vol in List of Touchables."
             << G4endl;
    }
  }

  if (fVerbose)
  {
    G4cout << "TouchCollect: "
           << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo() << "  type=Parameterised"
           << "  n=" << replicaNo << " -- at depth = " << depth << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableCollector::VisitUnknown(const G4TouchableHistoryHandle& touchable,
                                        G4VPhysicalVolume* pv, G4int depth)
{
  if (pv == fTargetPhysVol)
  {
    fListOfTouchableHandles.push_back(touchable);
    if (fVerbose)
    {
      G4cout << " **> Added new entry (Unknown) for Target Phys-Vol in List of Touchables."
             << G4endl;
    }
  }
  if (fVerbose)
  {
    G4cout << "TouchCollect: "
           << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
           << "  type=Unknown(kExternal?)"
           << " at depth= " << depth << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableCollector::EnterLogicalVolume(const G4TouchableHistoryHandle& touchable,
                                              G4LogicalVolume* lv, G4int depth)
{
  if (lv == fTargetLogVol)
  {
    fListOfTouchableHandles.push_back(touchable);
    if (fVerbose)
    {
      G4cout << " **> Added new entry for LV in List of Touchables." << G4endl;
    }
  }

  if (fVerbose)
  {
    const G4String solidName = lv->GetSolid() ? lv->GetSolid()->GetName() : G4String("(none)");
    const G4String materialName =
      lv->GetMaterial() ? lv->GetMaterial()->GetName() : G4String("(none)");

    G4cout << "TC::ELV - "  // << Indent(depth)
           << "  [LV] " << lv->GetName() << "  solid=" << solidName << "  material=" << materialName
           << "  depth= " << depth << G4endl;
  }
}

// ---------------------------------------------------------------------------

G4String G4TouchableCollector::Indent(G4int depth) const
{
  const G4int indentSize = 2;
  return G4String(static_cast<std::size_t>(indentSize * depth), ' ');
}
