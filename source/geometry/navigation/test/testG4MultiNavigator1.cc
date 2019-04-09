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

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4MultiNavigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "MyMultiNavigator.hh"

// Sample Paramterisation
class G4LinScale : public G4VPVParameterisation
{
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,(n-1)*15,0));
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int n,
				 const G4VPhysicalVolume*) const
  {
    pBox.SetXHalfLength(10);
    pBox.SetYHalfLength(5+n);
    pBox.SetZHalfLength(5+n);
  }

  virtual void ComputeDimensions(G4Tubs &,
				 const G4int ,
                                 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trd &, 
				 const G4int,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Cons &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trap &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Hype &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Orb &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Sphere &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Torus &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Para &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polycone &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polyhedra &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
};
G4LinScale myParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("cuboid",25,25,20);
    G4Box *myBox=new G4Box("cube",10,10,10);
    G4Box *mySlab= new G4Box("slab",10,25,10);

    G4Box *myVariableBox=new G4Box("Variable Box",10,5,5);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume *boxLog=new G4LogicalVolume(myBox,0,
						"Crystal Box",0,0,0);
    G4LogicalVolume *slabLog=new G4LogicalVolume(mySlab,0,
						 "Crystal Slab",0,0,0);
//  G4PVPlacement *offMXYPhys=
                              new G4PVPlacement(0,G4ThreeVector(-15,15,-10),
						"Target 1",boxLog,
						worldPhys,false,0);
//  G4PVPlacement *offMXMYPhys=
                               new G4PVPlacement(0,G4ThreeVector(-15,-15,-10),
						 "Target 2",boxLog,
						 worldPhys,false,0);

    G4PVPlacement *offYPhys=new G4PVPlacement(0,G4ThreeVector(15,0,-10),
					      "Target 3",slabLog,
					      worldPhys,false,0);

//  G4PVPlacement *offYZPhys=
                             new G4PVPlacement(0,G4ThreeVector(0,15,10),
					       "Target 4",boxLog,
					       worldPhys,false,0);
//  G4PVPlacement *offMYZPhys=
                              new G4PVPlacement(0,G4ThreeVector(0,-15,10),
						"Target 5",boxLog,
						worldPhys,false,0);


    G4LogicalVolume *variLog=new G4LogicalVolume(myVariableBox,0,
					         "Variable Blocks",0,0,0);
//  G4PVParameterised *paramPhys=
                                 new G4PVParameterised("Vari' Blocks",
						       variLog,
						       offYPhys,
						       kYAxis,
						       3,
						       &myParam);
    return worldPhys;
}

MyMultiNavigator& CreateMyMultiNavigator(G4VPhysicalVolume *pTopNode)
{
    MyMultiNavigator* pMultiNav= new MyMultiNavigator(); 
    pMultiNav->SetWorldVolume(pTopNode);
    pMultiNav->PrepareNavigators();     

    // G4ThreeVector zero(0.,0.,0.); 
    // pMultiNav->PrepareNewTrack( zero, zero );     

    return *pMultiNav; 
}

//
// Test LocateGlobalPointAndSetup
//
G4bool testG4MultiNavigator1(G4VPhysicalVolume *pTopNode)
{
  // MyMultiNavigator*  pNav= CreateMyMultiNavigator();
    MyMultiNavigator& myNav= CreateMyMultiNavigator(pTopNode);
    G4VPhysicalVolume *located;

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check relative search that causes backup one level and then search down:
// Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10),0,false);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-15,20));
    assert(located->GetName()=="Target 5");
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),G4ThreeVector(0,0,10)));
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check parameterised volumes

// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-16));
    assert(located->GetName()=="Target 3");

// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-17));
    assert(located->GetName()=="Target 3");

// Replication 2
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-18));
    assert(located->GetName()=="Target 3");

    delete &myNav; 

    return true;
}


//
// Test Stepping
//
G4bool testG4MultiNavigator2(G4VPhysicalVolume *pTopNode)
{
    // MyMultiNavigator& myNav= CreateMyMultiNavigator(pTopNode);
    MyMultiNavigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    myNav.SetWorldVolume(pTopNode);
    myNav.PrepareNavigators();     
  
//
// Test location & Step computation
//  
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mxHat,physStep,safety);
    assert(ApproxEqual(Step,25));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),xHat,physStep,safety);
    assert(ApproxEqual(Step,5));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(5,0,-10),0,true);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),zHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mzHat,physStep,safety);
    assert(ApproxEqual(Step,10));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);


//
// Test stepping through common boundaries
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,-20));
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,-20),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,0));
    assert(located->GetName()=="Target 4");
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,0),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,20));
    assert(!located);

//
// Test mother limited Step
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25,0,10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-25,0,10),xHat,physStep,safety);
    assert(ApproxEqual(Step,50));
    assert(ApproxEqual(safety,0));

//
// Test stepping through parameterised volumes
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-25,-10),0,false);
    assert(located->GetName()=="Target 3");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(15,-25,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-20,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-20,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,10));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-10,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,-10,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-6,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,12));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,6,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,8,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,8,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,14));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,22,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,22,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,3));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,25,-10));
    assert(!located);

    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4MultiNavigator1(myTopNode);
    testG4MultiNavigator2(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4MultiNavigator1(myTopNode);
    testG4MultiNavigator2(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}
