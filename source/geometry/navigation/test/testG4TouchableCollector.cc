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
//
//
//
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

#include "ApproxEqual.hh"
#include <assert.h>

// Global defs
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableCollector.hh"
#include "G4TouchableWalker.hh"
#include "G4TransportationManager.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

// Sample Paramterisation
class G4LinScale : public G4VPVParameterisation
{
    virtual void ComputeTransformation(const G4int n, G4VPhysicalVolume* pRep) const
    {
      pRep->SetTranslation(G4ThreeVector(0, (n - 1) * 15, 0));
    }

    virtual void ComputeDimensions(G4Box& pBox, const G4int n, const G4VPhysicalVolume*) const
    {
      pBox.SetXHalfLength(10);
      pBox.SetYHalfLength(5 + n);
      pBox.SetZHalfLength(5 + n);
    }

    virtual void ComputeDimensions(G4Tubs&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Trd&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Trap&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Hype&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Orb&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Sphere&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Torus&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Para&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Polycone&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Polyhedra&, const G4int, const G4VPhysicalVolume*) const {}
};
G4LinScale myParam;

G4PVPlacement* worldPhys = nullptr;
G4PVPlacement* offMXYPhys = nullptr;
G4PVPlacement* offYPhys = nullptr;
G4PVPlacement* offMXMYPhys = nullptr;

G4PVParameterised* paramPhys = nullptr;
G4PVPlacement* offYZPhys = nullptr;
G4PVPlacement* offMYZPhys = nullptr;

G4LogicalVolume* worldLog = nullptr;
G4LogicalVolume* slabLog = nullptr;
G4LogicalVolume* boxLog = nullptr;
G4LogicalVolume* variLog = nullptr;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{
  G4Box* myBigBox = new G4Box("cuboid", 25, 25, 20);
  G4Box* mySlab = new G4Box("slab", 10, 25, 10);
  G4Box* myBox = new G4Box("cube", 10, 10, 10);

  G4Box* myVariableBox = new G4Box("Variable Box", 10, 5, 5);

  worldLog = new G4LogicalVolume(myBigBox, 0, "World", 0, 0, 0);
  // Logical with no material,field,
  // sensitive detector or user limits

  worldPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "World", worldLog, 0, false, 0);
  // Note: no mother pointer set

  slabLog = new G4LogicalVolume(mySlab, 0, "Crystal Slab", 0, 0, 0);
  boxLog = new G4LogicalVolume(myBox, 0, "Crystal Box", 0, 0, 0);

  offMXYPhys =
    new G4PVPlacement(0, G4ThreeVector(-15, 15, -10), "Target 1", boxLog, worldPhys, false, 0);

  offMXMYPhys =
    new G4PVPlacement(0, G4ThreeVector(-15, -15, -10), "Target 2", boxLog, worldPhys, false, 0);

  //  Use Rotation to check exit normal -- Box is unchanged !
  // G4RotationMatrix *rot1=new G4RotationMatrix();
  // rot1->rotateY(pi);

  offYPhys = new G4PVPlacement(0, G4ThreeVector(15, 0, -10),  // (rot1,
                               "Target 3", slabLog, worldPhys, false, 3);  // copy# = 3

  //  Use Rotation to check exit normal -- Box is unchanged !
  G4RotationMatrix* rot2 = new G4RotationMatrix();
  rot2->rotateZ(0.5 * pi);
  offYZPhys =
    new G4PVPlacement(rot2, G4ThreeVector(0, 15, 10), "Target 4", boxLog, worldPhys, false, 4);

  offMYZPhys =
    new G4PVPlacement(0, G4ThreeVector(0, -15, 10), "Target 5", boxLog, worldPhys, false, 5);

  variLog = new G4LogicalVolume(myVariableBox, 0, "Variable Blocks", 0, 0, 0);

  paramPhys = new G4PVParameterised("Vari' Blocks", variLog, offYPhys, kYAxis, 3, &myParam);
  return worldPhys;
}

G4bool gVerbose = true;

// ------------------------------------------------------------------------------------------
// checkCollectorForLogicalVolumes
//
G4bool checkCollectorForLogicalVolumes(G4LogicalVolume* targetLV, G4int expectedCount)
{
  assert(targetLV != nullptr);
  G4cout << G4endl;
  G4cout << "checkCollectorForLogicalVolumes: called  with targetLV = " << targetLV
         << " name = " << targetLV->GetName() << G4endl;
  G4TouchableCollector listBuilder(targetLV);
  G4TouchableWalker lbWalker(listBuilder);
  lbWalker.Walk(worldPhys);
  // G4cout << " WALKing ... Ended." << G4endl;
  auto listOfTouchables = listBuilder.GetListOfMatchingTouchables();

  G4int count = listOfTouchables.size();
  if (count != expectedCount)
  {
    G4cout << "Expected " << expectedCount << " touchables, got " << count << G4endl;
    return false;
  }

  G4cout << "Touchables collected - with matching LV are: " << G4endl;
  for (auto item : listOfTouchables)
  {
    G4int depth = item->GetHistoryDepth();
    G4VPhysicalVolume* pv = item->GetVolume();
    G4VSolid* solid = item->GetSolid();
    G4int replicaNo = item->GetReplicaNumber();
    G4LogicalVolume* lv = pv->GetLogicalVolume();

    G4cout << "Touchable: : "  // << Indent(depth)
           << "[PV] "
           << pv->GetName()
           //           << "  type= " << pv->VolumeType()
           << "  repNum= " << replicaNo << "  pv-copy=" << pv->GetCopyNo()
           << " -- at depth= " << depth << " (lv= " << lv
           << "  name= " << (lv != nullptr ? lv->GetName() : G4String("-Not-Available-")) << " )"
           << " solid = " << solid << " name= " << solid->GetName() << G4endl;

    if (lv != targetLV)
    {
      G4Exception("RunAction", "User-code", FatalException,
                  "Mismatch in logical volume - does not agree with requested Target");
    }
  }
  if (gVerbose)
  {
    G4cout << "check Collector for Logical Volumes - SUCCESS " << G4endl << G4endl;
  }
  return true;
}

// ------------------------------------------------------------------------------------------
//
// Test TouchableCollector for logical volumes
//
G4bool testG4TouchableCollector1(G4VPhysicalVolume* pTopNode)
{
  MyNavigator myNav;
  G4VPhysicalVolume* located;
  myNav.SetWorldVolume(pTopNode);

  G4cout << "-----------------------------------------------------------------" << G4endl;
  G4cout << "-Starting  testG4TouchableCollector1" << G4endl;

  assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity, 0, 0), 0, false));
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(0, 0, 0), 0, false);
  assert(located->GetName() == "World");

  G4TransportationManager::GetTransportationManager()->SetNavigatorForTracking(&myNav);

  assert(checkCollectorForLogicalVolumes(worldLog, 1));
  assert(checkCollectorForLogicalVolumes(boxLog, 4));
  assert(checkCollectorForLogicalVolumes(variLog, 3));

  assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity, 0, 0)));

  // Check relative search that causes backup one level and then search down:
  // Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 0, -10), 0, false);
  assert(located->GetName() == "Vari' Blocks");

  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(0, -15, 20));
  assert(located->GetName() == "Target 5");
  assert(ApproxEqual(myNav.CurrentLocalCoordinate(), G4ThreeVector(0, 0, 10)));
  // Check that outside point causes stack to unwind
  assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity, 0, 0)));

  // Check parameterised volumes

  // Replication 0
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -15, -10));
  assert(located->GetName() == "Vari' Blocks");
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -15, -16));
  assert(located->GetName() == "Target 3");

  // Replication 1
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 0, -10));
  assert(located->GetName() == "Vari' Blocks");
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 0, -17));
  assert(located->GetName() == "Target 3");

  // Replication 2
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 15, -10));
  assert(located->GetName() == "Vari' Blocks");
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 15, -18));
  assert(located->GetName() == "Target 3");

  G4cout << "-Ended    testG4TouchableCollector1" << G4endl;

  return true;
}

// ------------------------------------------------------------------------------------
// checkCollectorForPhysicalVolumes
//
G4bool checkCollectorForPhysicalVolumes(G4VPhysicalVolume* targetPV, G4int expectedCount)
{
  assert(targetPV != nullptr);
  G4cout << "checkCollectorForPhysicalVolumes: called  with targetPV = " << targetPV << " name = '"
         << targetPV->GetName() << "'" << G4endl;
  G4TouchableCollector listBuilder(targetPV);
  G4TouchableWalker lbWalker(listBuilder);
  lbWalker.Walk(worldPhys);
  auto listOfTouchables = listBuilder.GetListOfMatchingTouchables();

  G4int count = listOfTouchables.size();
  if (count != expectedCount)
  {
    G4cout << "checkCollectorForPhysicalVolumes: Expected " << expectedCount << " touchables, "
           << " but got " << count << G4endl;
    return false;
  }

  G4cout << "Touchables collected - with matching LV are: " << G4endl;
  for (auto item : listOfTouchables)
  {
    G4int depth = item->GetHistoryDepth();
    G4VPhysicalVolume* pv = item->GetVolume();
    G4VSolid* solid = item->GetSolid();
    G4int replicaNo = item->GetReplicaNumber();
    G4LogicalVolume* lv = pv->GetLogicalVolume();

    G4cout << "Touchable: : "  // << Indent(depth)
           << "[PV] " << pv->GetName() << "  type= " << pv->VolumeType() << "  copy/replica-No= "
           << replicaNo
           //<< "  pv-copy=" << pv->GetCopyNo()
           << "  name= " << (pv != nullptr ? pv->GetName() : G4String("-Not-Available-")) << " )"
           << " -- at depth= " << depth << " (lv= " << lv
           << "  name= " << (lv != nullptr ? lv->GetName() : G4String("-Not-Available-")) << " )"
           << " solid = " << solid << " name= " << solid->GetName() << G4endl;

    if (pv != targetPV)
    {
      G4Exception("checkCollectorForPhysicalVolumes", "Test-code", FatalException,
                  "Mismatch in physical volume - does not agree with requested Target");
    }
  }
  return true;
}

//
// Test Stepping
//
G4bool testG4TouchableCollector2(G4VPhysicalVolume* pTopNode)
{
  MyNavigator myNav;
  G4VPhysicalVolume* located;
  G4double Step, physStep, safety;
  G4ThreeVector xHat(1, 0, 0), yHat(0, 1, 0), zHat(0, 0, 1);
  G4ThreeVector mxHat(-1, 0, 0), myHat(0, -1, 0), mzHat(0, 0, -1);

  G4ThreeVector position, direction;  // memory for use in checks, move
  G4ThreeVector globalNormal;
  G4bool validNormal;

  myNav.SetWorldVolume(pTopNode);

  G4cout << "-----------------------------------------------------------------" << G4endl;
  G4cout << "-Starting testG4TouchableCollector2" << G4endl;

  G4TransportationManager::GetTransportationManager()->SetNavigatorForTracking(&myNav);

  // assert(checkCollectorForPhysicalVolumes(worldPhys, 1));
  assert(checkCollectorForPhysicalVolumes(offYZPhys, 1));
  assert(checkCollectorForPhysicalVolumes(offMYZPhys, 1));

  assert(checkCollectorForPhysicalVolumes(paramPhys, 3));

  //
  // Test location & Step computation
  //
  located = myNav.LocateGlobalPointAndSetup(position = G4ThreeVector(0, 0, -10));
  assert(located->GetName() == "World");
  physStep = kInfinity;

  Step = myNav.ComputeStep(position, direction = mxHat, physStep, safety);
  // Intersect 'Target 3' at ( 25, 0, -10)
  assert(ApproxEqual(Step, 25));
  //    assert(ApproxEqual(safety,5));
  assert(safety >= 0);

  // G4cout << " Checking Exit Normal at " << position+Step*direction << G4endl;
  globalNormal = myNav.GetGlobalExitNormal(position + Step * direction, &validNormal);
  assert(globalNormal == mxHat);
  assert(validNormal);

  // 2nd locate / compute - step
  located = myNav.LocateGlobalPointAndSetup(position = G4ThreeVector(0, 0, -10));
  assert(located->GetName() == "World");
  physStep = kInfinity;

  Step = myNav.ComputeStep(position, direction = xHat, physStep, safety);
  assert(ApproxEqual(Step, 5));
  //    assert(ApproxEqual(safety,5));
  assert(safety >= 0);
  myNav.SetGeometricallyLimitedStep();

  globalNormal = myNav.GetGlobalExitNormal(position + Step * direction, &validNormal);
  // G4cout << " TestNav1> Exit Normal (global) = " << globalNormal << G4endl;
  assert(globalNormal == direction);
  assert(validNormal);

  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(5, 0, -10), 0, true);
  assert(located->GetName() == "Vari' Blocks");

  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(0, 0, -10));
  assert(located->GetName() == "World");
  physStep = kInfinity;
  Step = myNav.ComputeStep(G4ThreeVector(0, 0, -10), direction = zHat, physStep, safety);
  assert(ApproxEqual(Step, 30));
  //    assert(ApproxEqual(safety,5));
  assert(safety >= 0);

  globalNormal = myNav.GetGlobalExitNormal(position + Step * direction, &validNormal);
  // G4cout << " TestNav1> Exit Normal (global) = " << globalNormal << G4endl;
  assert(globalNormal == direction);
  assert(validNormal);

  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(0, 0, -10));
  assert(located->GetName() == "World");
  physStep = kInfinity;
  Step = myNav.ComputeStep(G4ThreeVector(0, 0, -10), mzHat, physStep, safety);
  assert(ApproxEqual(Step, 10));
  //    assert(ApproxEqual(safety,5));
  assert(safety >= 0);

  //
  // Test stepping through common boundaries
  //
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7, 7, -20));
  assert(located->GetName() == "Target 1");
  physStep = kInfinity;
  Step = myNav.ComputeStep(G4ThreeVector(-7, 7, -20), zHat, physStep, safety);
  assert(ApproxEqual(Step, 20));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7, 7, 0));
  assert(located->GetName() == "Target 4");
  Step = myNav.ComputeStep(G4ThreeVector(-7, 7, 0), zHat, physStep, safety);
  assert(ApproxEqual(Step, 20));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7, 7, 20));
  assert(!located);

  //
  // Test mother limited Step
  //
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25, 0, 10));
  assert(located->GetName() == "World");
  physStep = kInfinity;
  Step = myNav.ComputeStep(G4ThreeVector(-25, 0, 10), xHat, physStep, safety);
  assert(ApproxEqual(Step, 50));
  assert(ApproxEqual(safety, 0));

  //
  // Test stepping through parameterised volumes
  //
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -25, -10), 0, false);
  assert(located->GetName() == "Target 3");
  physStep = kInfinity;
  Step = myNav.ComputeStep(G4ThreeVector(15, -25, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 5));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -20, -10));
  assert(located->GetName() == "Vari' Blocks");
  Step = myNav.ComputeStep(G4ThreeVector(15, -20, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 10));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -10, -10));
  assert(located->GetName() == "Target 3");
  Step = myNav.ComputeStep(G4ThreeVector(15, -10, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 4));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, -6, -10));
  assert(located->GetName() == "Vari' Blocks");
  Step = myNav.ComputeStep(G4ThreeVector(15, -6, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 12));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 6, -10));
  assert(located->GetName() == "Target 3");
  Step = myNav.ComputeStep(G4ThreeVector(15, 6, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 2));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 8, -10));
  assert(located->GetName() == "Vari' Blocks");
  Step = myNav.ComputeStep(G4ThreeVector(15, 8, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 14));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 22, -10));
  assert(located->GetName() == "Target 3");
  Step = myNav.ComputeStep(G4ThreeVector(15, 22, -10), yHat, physStep, safety);
  assert(ApproxEqual(Step, 3));
  assert(ApproxEqual(safety, 0));
  myNav.SetGeometricallyLimitedStep();
  located = myNav.LocateGlobalPointAndSetup(G4ThreeVector(15, 25, -10));
  assert(!located);

  G4cout << "-Ended    testG4TouchableCollector2" << G4endl;
  G4cout << "-----------------------------------------------------------------" << G4endl << G4endl;

  return true;
}

int main()
{
  G4VPhysicalVolume* myTopNode;
  myTopNode = BuildGeometry();  // Build the geometry
  G4GeometryManager::GetInstance()->CloseGeometry(false);
  testG4TouchableCollector1(myTopNode);
  testG4TouchableCollector2(myTopNode);
  // Repeat tests but with full voxels
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4GeometryManager::GetInstance()->CloseGeometry(true);
  testG4TouchableCollector1(myTopNode);
  testG4TouchableCollector2(myTopNode);

  G4GeometryManager::GetInstance()->OpenGeometry();
  return 0;
}
