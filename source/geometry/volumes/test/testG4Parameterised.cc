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
//   Test the Navigation in geometry with parameterised volumes (which 
//    include rotations as well as translations).
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.
//   Started from testG4Navigator1.cc

#include <assert.h>
#include "G4ios.hh"
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Sample Parameterisation
class MoveNRotate : public G4VPVParameterisation
{
 public:
  MoveNRotate(G4double twistAngle)
  { 
    fTwistAngle= twistAngle;
    fRotationVec= new G4RotationMatrix();
  }

  virtual ~MoveNRotate() { delete fRotationVec; }

  G4double GetTwistAngle() { return fTwistAngle; }
  void     SetTwistAngle(G4double newAngle ) { fTwistAngle= newAngle; }

private:
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,n*100,0));
    *fRotationVec = G4RotationMatrix();             // Unit matrix
    fRotationVec->rotateZ( n * fTwistAngle );
    pRep->SetRotation( fRotationVec );
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int,
				 const G4VPhysicalVolume*) const
  {
    pBox.SetXHalfLength(10);
    pBox.SetYHalfLength(10);
    pBox.SetZHalfLength(10);
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
 private:
    G4RotationMatrix *fRotationVec;
    G4double fTwistAngle;
};

G4double    angle1= 15.0*pi/180.;
MoveNRotate myParam(angle1);


// Build simple geometry:
// 4 small cubes (G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    // The world volume
    //
    G4Box *myBigBox= new G4Box ("Big Cube", 500, 500, 500);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set


    // A set of boxes
    G4Box *myBox=new G4Box("cube",10,10,10);
    G4LogicalVolume *boxLog=new G4LogicalVolume(myBox,0,
						"Rotating Box",0,0,0);

//  G4PVParameterised *paramP=
                              new G4PVParameterised("Rotating Blocks",
						    boxLog,
						    worldPhys, //OR worldLog,
						    kYAxis,
						    3,
						    &myParam);
    // Copies 0, 1 & 2 will exist    

    return worldPhys;
}

//
// Test LocateGlobalPointAndSetup
//
G4bool testG4Navigator1(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0, false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(100,100,100),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// 
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-5,-5),0,false);
    assert(located->GetName()=="Rotating Blocks");
    assert(located->GetCopyNo()== 0);
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),G4ThreeVector(0,-5,-5)));
    G4cout << " Local coords = " << myNav.CurrentLocalCoordinate() << G4endl;

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,100,5));
    assert(located->GetName()=="Rotating Blocks");
    assert(located->GetCopyNo()== 1);
    G4cout << " Local coords = " << myNav.CurrentLocalCoordinate() << G4endl;
//    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),
//                       G4ThreeVector(0,0,10)));
    
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check parameterised volumes

// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,5,5));
    assert(located->GetName()=="Rotating Blocks");
    assert(located->GetCopyNo()== 0);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,15,15));
    assert(located->GetName()=="World");

// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,105,5));
    assert(located->GetName()=="Rotating Blocks");
    assert(located->GetCopyNo()== 1);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-17));
    assert(located->GetName()=="World");

// Replication 2
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,205,5));
    assert(located->GetName()=="Rotating Blocks");
    assert(located->GetCopyNo()== 2);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-18));
    assert(located->GetName()=="World");

    return true;
}


//
// Test Stepping
//
G4bool testG4Navigator2(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    myNav.SetWorldVolume(pTopNode);
  
//
// Test location & Step computation
//  
    G4ThreeVector  StartPoint(-50,0,-5);
    located=myNav.LocateGlobalPointAndSetup( StartPoint ); 
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( StartPoint, mxHat,physStep,safety);  // -x dir
    assert(ApproxEqual(Step,450));
    // assert(ApproxEqual(safety,40));
    // assert(safety>=0);

    StartPoint= G4ThreeVector(-15,0,-5);
    located=myNav.LocateGlobalPointAndSetup( StartPoint ); 
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( StartPoint,xHat,physStep,safety); // +x dir
    assert(ApproxEqual(Step,5));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);
    myNav.SetGeometricallyLimitedStep();
    G4ThreeVector EndPoint = StartPoint + Step * xHat;
    located=myNav.LocateGlobalPointAndSetup(EndPoint,0,true);
    assert(located->GetName()=="Rotating Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-40));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-40),zHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0, 40));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,40),mzHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);


//
// Test moving through series of volumes
//
    StartPoint= G4ThreeVector(0,-20,0);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-20,0));
    assert(located->GetName()=="World");
    
    // Replication 0 block
    //
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,-20,0),yHat,physStep,safety);
    assert(ApproxEqual(Step,10));
    EndPoint= StartPoint + Step * yHat;   //  Should be  0, -10, 0
    assert(ApproxEqual( 0, (EndPoint-G4ThreeVector(0,-10,0)).mag()) );
    // assert(ApproxEqual(safety,0));
    
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint) ;
    assert(located->GetName()=="Rotating Blocks");
    Step=myNav.ComputeStep(EndPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * yHat;   //  Should be  0, +10, 0
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Replication 1 block
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step,90.-10./std::cos(angle1)));
    EndPoint= StartPoint + Step * yHat;   //  Should near  0, 90, 0
    assert(safety<=Step);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint) ;
    assert(located->GetName()=="Rotating Blocks");
    
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step,20./std::cos(angle1)));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * yHat;   //  Should be near 0, 110, 0
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Replication 2 block
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step,100.-10.*(1./std::cos(angle1)+1./std::cos(2.*angle1))));
    EndPoint= StartPoint + Step * yHat;   //  Should near  0, 190, 0
    assert(safety<=Step);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint);
    assert(located->GetName()=="Rotating Blocks");
    
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step,20./std::cos(2.*angle1)));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * yHat;   //  Should be near 0, 110, 0
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Edge of the world
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,yHat,physStep,safety);
    assert(ApproxEqual(Step, 300. - 10./std::cos(2.*angle1) ));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * yHat;   //  Should be near 0, 110, 0
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(!located);


    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}




