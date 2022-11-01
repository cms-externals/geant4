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
//   Locate & Step within simple rotated boxlike geometry, both
//   with and without voxels.

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
// #include "G4PVParameterised.hh"
// #include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("cube",50,50,50);
    G4Box *myBox=new G4Box("cuboid",5,10,15);
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
    
    G4RotationMatrix *rot1=new G4RotationMatrix();
    rot1->rotateZ(pi*0.5);
//  G4PVPlacement *boxPhys1=
                            new G4PVPlacement(rot1,G4ThreeVector(0,0,0),
					      "Target 1",boxLog,
					      worldPhys,false,0);

    G4RotationMatrix *rot2=new G4RotationMatrix();
    rot2->rotateX(pi*0.5);
//  G4PVPlacement *boxPhys2=
                            new G4PVPlacement(rot2,G4ThreeVector(-30,0,0),
					      "Target 2",boxLog,
					      worldPhys,false,0);
    G4RotationMatrix *rot3=new G4RotationMatrix();
    rot3->rotateY(pi*0.5);
//  G4PVPlacement *boxPhys3=
                            new G4PVPlacement(rot3,G4ThreeVector(30,0,0),
					      "Target 3",boxLog,
					      worldPhys,false,0);
    return worldPhys;
}


//
// Test LocateGlobalPointAndSetup
//
G4bool testG4NavigatorLocate(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="Target 1");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(10,0,0),0,false);
    assert(located->GetName()=="Target 1");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,10,0),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,5,0),0,false);
    assert(located->GetName()=="Target 1");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-30,0,0),0,false);
    assert(located->GetName()=="Target 2");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-30,10,15),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-35,15,10),0,false);
    assert(located->GetName()=="Target 2");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25,-15,-10),0,false);
    assert(located->GetName()=="Target 2");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,0,0),0,false);
    assert(located->GetName()=="Target 3");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,15,5),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,10,5),0,false);
    assert(located->GetName()=="Target 3");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(45,-10,-5),0,false);
    assert(located->GetName()=="Target 3");

    return true;
}

// Test ComputeStep
//
G4bool testG4NavigatorSteps(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector pos,dir,origin,xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    myNav.SetWorldVolume(pTopNode);

    pos=origin;
    dir=xHat;
    located=myNav.LocateGlobalPointAndSetup(pos,0,false);
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    Step=myNav.ComputeStep(origin,dir,physStep,safety);
    assert(ApproxEqual(Step,10));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,30));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    pos+=Step*dir;

    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located==0);
    return true;
}

//
//
//
G4double
CheckNextStep_And_TakeIt( G4Navigator &Navigator, 
                          const G4ThreeVector &position,
                          const G4ThreeVector &direction,
                          G4double proposedStep,
  			  G4double &resNewSafety,           //  Output !!
 			  G4double expectedStep)
{
    // G4double Step1= 0.0, Step2= 0.0; 
    G4double Step= 0.0;
    G4double safety1=0.0, safety2=0.0; 
    // Check new ComputeStep method
    Step=Navigator.CheckNextStep(position,direction,proposedStep,safety1);
    assert(ApproxEqual(Step,expectedStep));

    Step=Navigator.ComputeStep(position,direction,proposedStep,safety2);
    assert(ApproxEqual(Step,expectedStep));

    assert(ApproxEqual(safety1,safety2)); 

    resNewSafety= safety2; 
    return Step; 
}

//

//
// Test CheckNextSteps
//
G4bool testG4NavigatorCheckNextSteps(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector pos,dir,origin,xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    myNav.SetWorldVolume(pTopNode);

    G4double expectedStep=0.0;

    pos=origin;
    dir=xHat;
    located=myNav.LocateGlobalPointAndSetup(pos,0,false);
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    // Step=myNav.ComputeStep(origin,dir,physStep,safety);
    // assert(ApproxEqual(Step,10));
    Step= CheckNextStep_And_TakeIt(myNav,pos,dir,physStep,safety,expectedStep=10.0);
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    // Step=myNav.ComputeStep(pos,dir,physStep,safety);
    // assert(ApproxEqual(Step,5));
    Step= CheckNextStep_And_TakeIt(myNav,pos,dir,physStep,safety,expectedStep=5.0);
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="Target 3");

    // Step=myNav.ComputeStep(pos,dir,physStep,safety);
    // assert(ApproxEqual(Step,30));
    Step= CheckNextStep_And_TakeIt(myNav,pos,dir,physStep,safety,expectedStep=30.0);
    pos+=Step*dir;

    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");

    // Step=myNav.ComputeStep(pos,dir,physStep,safety);
    // assert(ApproxEqual(Step,5));
    Step= CheckNextStep_And_TakeIt(myNav,pos,dir,physStep,safety,expectedStep=5.0);
    // ** New: Check new ComputeStep method
    // Step=myNav.CheckNextStep(pos,dir,physStep,safety);
    // assert(ApproxEqual(Step,5));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located==0);
    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4NavigatorLocate(myTopNode);
    testG4NavigatorSteps(myTopNode);
    testG4NavigatorCheckNextSteps(myTopNode);

// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4NavigatorLocate(myTopNode);
    testG4NavigatorSteps(myTopNode);
    testG4NavigatorCheckNextSteps(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();

    G4PhysicalVolumeStore *pvs=G4PhysicalVolumeStore::GetInstance();
    for (G4int i=pvs->size()-1;i>=0;i--)
      {
	// Delete any rotation matrices
	delete (*pvs)[i]->GetRotation();
      }
    return 0;
}
