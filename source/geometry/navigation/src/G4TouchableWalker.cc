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
// G4TouchableWalker class Implementation
//
// Author: John Apostolakis (CERN), May-June 2026
// --------------------------------------------------------------------

#include "G4TouchableWalker.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NavigationHistory.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ios.hh"

#include "geomdefs.hh"

#include <cassert>
#include <iomanip>
#include <memory>

// Design Note:  On the use of the G4TouchableWalker class and thread safety
//               -----------------------------------------------------------
//
// The Walker is expected to be used outside the run/event/tracking loop,
// or in the BeginOfRun or EndOfRun action.  This is because it has
// the potential to change the state of the geometry objects in memory.
//
// If the geometry contains only placements and replicas, the Walker will not
// change the state of the geometry. So it is thread-safe by definition.
//
// For replica volumes, the walker will set the Replica Number.
//
// For parameterised volumes, the walker sets the state of the parameterised
// volume's attributes (solid pointer, material pointer, position/rotation.)
// However this state is anyway per-thread -- as it must be in order to allow
// for correct navigation in a multi-threaded application.
//
// In other words, the classes which have these attributes are split classes:
// the changeable part of the state, that which depends on the instance/replica
// number of the parameterisation, has a separate per-thread sub-object to hold
// it. In this way parameterised volume prmA's logical volume volA can point to
// solid1 and material1 for replicaNo=1 in thread 1 and to solid2 and material2
// for replicaNo=2 in thread 2.

// ---------------------------------------------------------------------------

G4TouchableWalker::G4TouchableWalker(G4VTouchableVisitor& visitor)
  : fVisitor(visitor), fNavHistory()
{}

// ---------------------------------------------------------------------------

G4TouchableWalker::~G4TouchableWalker() {}

// ---------------------------------------------------------------------------

void G4TouchableWalker::Walk(G4VPhysicalVolume* startVol)
{
  G4VPhysicalVolume* worldVol = GetWorldVolume();
  if (startVol == nullptr)
  {
    if (worldVol == nullptr)
    {
      G4Exception("G4TouchableWalker::Walk", "GeomNav001", FatalException,
                  "Both requested start volume and World Volume are nullptr");
    }
    assert(worldVol != nullptr);

    startVol = worldVol;
  }
  G4int assumedDepth = 0;

  if (startVol == worldVol)
  {
    fNavHistory.SetFirstEntry(startVol);
    DescendFirstPhysical(startVol, 0);
  }
  else
  {
    G4LogicalVolume* pMotherLogical = startVol->GetMotherLogical();

    G4ExceptionDescription ed;
    ed << "Walking from volume: '" << startVol->GetName() << "'"
       << " which is NOT the world volume."
       << " Volume type: " << startVol->VolumeType();
    G4Exception("G4TouchableWalker::Walk", "GeomNav002", JustWarning, ed);

    // Let's examine the volume type
    switch (startVol->VolumeType())
    {
      case kNormal:
        fNavHistory.SetFirstEntry(startVol);
        DescendFirstPhysical(startVol, assumedDepth);
        break;

      case kParameterised:
        fArtificialWorldPV =  // new G4PVPlacement(...
          std::make_unique<G4PVPlacement>(nullptr, G4ThreeVector(), pMotherLogical, "", nullptr,
                                          false, 0);
        fNavHistory.SetFirstEntry(fArtificialWorldPV.get());
        DescendFirstPhysical(fArtificialWorldPV.get(), assumedDepth);
        break;

      case kReplica:
        fArtificialWorldPV = std::make_unique<G4PVPlacement>(nullptr, G4ThreeVector(),
                                                             pMotherLogical, "", nullptr, false, 0);
        fNavHistory.SetFirstEntry(fArtificialWorldPV.get());
        DescendFirstPhysical(fArtificialWorldPV.get(), assumedDepth);
        break;

      default:  // kExternal or any future type
        auto touchable = GetTouchableHistoryHandle();
        fVisitor.VisitUnknown(touchable, startVol, assumedDepth);
        break;
    }
  }
}

// ---------------------------------------------------------------------------

G4VPhysicalVolume* G4TouchableWalker::GetWorldVolume(G4bool fatalOnFailure)
{
  auto transportMgr = G4TransportationManager::GetTransportationManager();

  G4VPhysicalVolume* worldPV = transportMgr->GetNavigatorForTracking()->GetWorldVolume();
  if (worldPV == nullptr)
  {
    if (fatalOnFailure)
    {
      G4Exception("G4TouchableWalker::GetWorldVolume", "G4Geom-Run007", FatalException,
                  "Cannot find World Volume - was mandatory");
    }
    else
    {
      G4Exception("G4TouchableWalker::GetWorldVolume", "G4Geom-Run008", JustWarning,
                  "Cannot find World Volume");
      return nullptr;
    }
  }

  if (fVerbose && worldPV != nullptr)
  {
    G4cout << "* G4TouchableWalker: found worldPV= " << worldPV
           << "  and its name= " << worldPV->GetName() << G4endl;
  }

  return worldPV;
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendFirstPhysical(G4VPhysicalVolume* pv, G4int depth)
{
  // The root volume has already been registered in fNavHistory by SetFirstEntry()
  // in Walk().  We must therefore visit it and descend into its daughters WITHOUT
  // calling NewLevel() (which would push a duplicate history entry) and WITHOUT
  // calling DepartPhysicalVolume() (which pops a NewLevel that was never pushed).

  if (fVerbose)
  {
    G4cout << "Walker::DescendFirstPhysical depth= " << depth << " - entering." << G4endl;
  }

  auto touchable = GetTouchableHistoryHandle();
  switch (pv->VolumeType())
  {
    case kNormal:
      fVisitor.VisitPlacement(touchable, static_cast<G4PVPlacement*>(pv), depth);
      DescendLogicalVolume(pv->GetLogicalVolume(), depth);
      break;

    case kParameterised:
      // Walk() wraps a parameterised root in an artificial kNormal placement,
      // so this branch is reached only for a bare parameterised start volume.
      fVisitor.VisitParameterised(touchable, static_cast<G4PVParameterised*>(pv), 0, depth);
      DescendLogicalVolume(pv->GetLogicalVolume(), depth);
      break;

    case kReplica:
      fVisitor.VisitReplica(touchable, static_cast<G4PVReplica*>(pv), 0, depth);
      DescendLogicalVolume(pv->GetLogicalVolume(), depth);
      break;

    default:
      DescendUnknown(pv, depth);
      break;
  }

  if (fVerbose)
  {
    G4cout << "Walker::DescendFirstPhysical depth= " << depth << " - EXITing." << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DepartPhysicalVolume(G4VPhysicalVolume*, G4int depth)
{
  fNavHistory.BackLevel();

  auto navLevel = fNavHistory.GetDepth();
  if (fVerbose)
  {
    G4cout << " DepartPhysicalVolumeCalled - popped to level " << navLevel << " vs depth= " << depth
           << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendPlacement(G4PVPlacement* pv, G4int depth)
{
  fNavHistory.NewLevel(pv, kNormal, pv->GetCopyNo());

  auto touchable = GetTouchableHistoryHandle();
  fVisitor.VisitPlacement(touchable, pv, depth);

  if (fVerbose)
  {
    G4VPhysicalVolume* topPV = fNavHistory.GetTopVolume();
    G4cout << "-G4TouchableWalker::DescendPlacement: level= " << fNavHistory.GetDepth()
           << " info:  ptr= " << topPV << " = " << topPV->GetName()
           << " copyNo= " << fNavHistory.GetReplicaNo(0) << G4endl;
  }

  DescendLogicalVolume(pv->GetLogicalVolume(), depth);
  DepartPhysicalVolume(pv, depth);
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendReplica(G4PVReplica* pv, G4int replicaNo, G4int depth)
{
  fNavHistory.NewLevel(pv, kReplica, replicaNo);
  auto touchable = GetTouchableHistoryHandle();
  fVisitor.VisitReplica(touchable, pv, replicaNo, depth);
  DescendLogicalVolume(pv->GetLogicalVolume(), depth);
  DepartPhysicalVolume(pv, depth);
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendParameterised(G4PVParameterised* pv, G4int replicaNo, G4int depth)
{
  G4TouchableHistory parentTouchable(fNavHistory);
  G4VPVParameterisation* pParam = pv->GetParameterisation();
  G4VSolid* pSolid = pParam->ComputeSolid(replicaNo, pv);

  if (pSolid == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Parameterisation failed to compute solid for replica " << replicaNo << " of volume '"
       << pv->GetName() << "'";
    G4Exception("G4TouchableWalker::DescendParameterised", "GeomNav0006", FatalException, ed);
  }

  pSolid->ComputeDimensions(pParam, replicaNo, pv);
  pParam->ComputeTransformation(replicaNo, pv);

  fNavHistory.NewLevel(pv, kParameterised, replicaNo);
  pv->SetCopyNo(replicaNo);

  G4LogicalVolume* pLogical = pv->GetLogicalVolume();
  pLogical->SetSolid(pSolid);

  G4Material* mat = pParam->ComputeMaterial(replicaNo, pv, &parentTouchable);
  pLogical->UpdateMaterial(mat);

  auto touchable = GetTouchableHistoryHandle();
  fVisitor.VisitParameterised(touchable, pv, replicaNo, depth);
  DescendLogicalVolume(pv->GetLogicalVolume(), depth);

  DepartPhysicalVolume(pv, depth);
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendAllReplicas(G4PVReplica* pv, G4int depth)
{
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

  for (G4int replicaNo = 0; replicaNo < nReplicas; ++replicaNo)
  {
    if (fVerbose)
    {
      G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
             << "  type=Replica"
             << "  replicaNo= " << replicaNo << "  ( numInstances=" << nReplicas << " ) "
             << "  axis=" << AxisName(axis) << "  width=" << width
             << " -- Method: TW:DescendAllReplicas" << G4endl;
    }
    DescendReplica(pv, replicaNo, depth);
  }
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendAllParameterised(G4PVParameterised* pv, G4int depth)
{
  EAxis axis = kUndefined;
  G4int nReplicas = -1;
  G4double width = 0., offset = 0.;
  G4bool consuming = false;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

  if (pv->GetRegularStructureId() == 0)
  {
    G4TouchableHistory parentTouchable(fNavHistory);
    G4VPVParameterisation* pParam = pv->GetParameterisation();

    for (G4int replicaNo = 0; replicaNo < nReplicas; ++replicaNo)
    {
      G4VSolid* pSolid = pParam->ComputeSolid(replicaNo, pv);
      pSolid->ComputeDimensions(pParam, replicaNo, pv);
      pParam->ComputeTransformation(replicaNo, pv);

      fNavHistory.NewLevel(pv, kParameterised, replicaNo);
      pv->SetCopyNo(replicaNo);
      //
      // Set the correct solid and material in Logical Volume
      //
      G4LogicalVolume* pLogical = pv->GetLogicalVolume();
      pLogical->SetSolid(pSolid);
      G4Material* mat = pParam->ComputeMaterial(replicaNo, pv, &parentTouchable);
      pLogical->UpdateMaterial(mat);

      assert(pv->GetCopyNo() == replicaNo);
      assert(pLogical->GetSolid() == pSolid);
      assert(pLogical->GetMaterial() == mat);

      auto touchable = GetTouchableHistoryHandle();
      fVisitor.VisitParameterised(touchable, pv, replicaNo, depth);

      DescendLogicalVolume(pv->GetLogicalVolume(), depth);
      DepartPhysicalVolume(pv, depth);

      if (fVerbose)
      {
        G4cout << Indent(depth) << "[PV-param] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
               << "  type=Parameterised"
               << "  replicaNo=" << replicaNo << "  n=" << nReplicas << "  axis=" << AxisName(axis)
               << "  -- Method: VisitAllParameterised " << G4endl;
      }
    }
  }
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendUnknown(G4VPhysicalVolume* pv, G4int depth)
{
  if (fVerbose)
  {
    G4cout << Indent(depth) << "[PV] " << pv->GetName() << "  copy=" << pv->GetCopyNo()
           << "  type=Unknown(kExternal?)" << G4endl;
  }
}

// ---------------------------------------------------------------------------

G4String G4TouchableWalker::Indent(G4int depth) const
{
  const G4int indentSize = 2;
  return G4String(static_cast<std::size_t>(indentSize * depth), ' ');
}

// ---------------------------------------------------------------------------

G4String G4TouchableWalker::AxisName(EAxis axis)
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

// ---------------------------------------------------------------------------

void G4TouchableWalker::DescendLogicalVolume(G4LogicalVolume* lv, G4int depth)
{
  if (lv == nullptr)
  {
    G4Exception("G4TouchableWalker::DescendLogicalVolume", "Geom-InvalidLogical", FatalException,
                "Entering null logical volume.");
    return;
  }

  if (fVerbose)
  {
    const G4String solidName = lv->GetSolid() ? lv->GetSolid()->GetName() : G4String("(none)");
    const G4String matName = lv->GetMaterial() ? lv->GetMaterial()->GetName() : G4String("(none)");

    G4cout << "TW::DLV - " << Indent(depth) << "  [LV] " << std::setw(15 - depth) << lv->GetName()
           << "  solid= " << std::setw(15) << solidName << "  material= " << matName << G4endl;
  }

  auto touchable = GetTouchableHistoryHandle();
  fVisitor.EnterLogicalVolume(touchable, lv, depth);

  const std::size_t nDaughters = lv->GetNoDaughters();

  if (nDaughters > 0)
  {
    G4VPhysicalVolume* firstDaughter = lv->GetDaughter(0);
    auto daughterType = firstDaughter->VolumeType();

    if (daughterType == kExternal)
    {
      DescendUnknown(firstDaughter, depth + 1);
    }
    else
    {
      if (nDaughters == 1 && daughterType != kNormal)
      {
        if (daughterType == kReplica)
        {
          DescendAllReplicas(static_cast<G4PVReplica*>(firstDaughter), depth + 1);
        }
        else
        {
          assert(daughterType == kParameterised);
          DescendAllParameterised(static_cast<G4PVParameterised*>(firstDaughter), depth + 1);
        }
      }
      else
      {
        for (std::size_t i = 0; i < nDaughters; ++i)
        {
          G4VPhysicalVolume* daughter = lv->GetDaughter(i);
          if (daughter->VolumeType() != kNormal)
          {
            G4ExceptionDescription ed;
            ed << "Daughter " << i << " ('" << daughter->GetName() << "')"
               << " of logical volume '" << lv->GetName() << "'"
               << " is not a placement (kNormal); found type=" << daughter->VolumeType() << ".\n"
               << "In a multi-daughter logical volume all daughters must be placements.";
            G4Exception("G4TouchableWalker::DescendLogicalVolume", "GeomNav0004", FatalException,
                        ed);
          }
          DescendPlacement(static_cast<G4PVPlacement*>(daughter), depth + 1);
        }
      }
    }
  }
  if (fVerbose)
  {
    G4cout << "Walker::DescendLogical depth= " << depth << " - exiting." << G4endl;
  }

  fVisitor.LeaveLogicalVolume(touchable, lv, depth);
  LeaveLogicalVolume(lv, depth);
}

// ---------------------------------------------------------------------------

void G4TouchableWalker::LeaveLogicalVolume(G4LogicalVolume*, G4int /*depth*/) {}
