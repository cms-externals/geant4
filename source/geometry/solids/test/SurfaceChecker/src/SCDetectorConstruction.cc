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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "SCDetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "SCDetectorMessenger.hh"
#include "SCMagneticField.hh"
#include "SCTrackerSD.hh"

#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"

#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tet.hh"

#include "G4Polycone.hh"

#include "G4TwistedTubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ReflectedSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
SCDetectorConstruction::SCDetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 logicTracker(0),physiTracker(0), 
 fpMagField(0), fWorldLength(0.),  fTrackerpDz(0.)
{
  fpMagField = new SCMagneticField();
  detectorMessenger = new SCDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
SCDetectorConstruction::~SCDetectorConstruction()
{
  delete fpMagField;
  delete detectorMessenger;             
}

////////////////////////////////////////////////////////////////

void SCDetectorConstruction::SwitchDetector()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(physiWorld);
}


G4VPhysicalVolume*
SCDetectorConstruction::SelectDetector( const G4String& val )
{


  G4double a, z;

  G4double density;
  G4int nel;

  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  // Print all the materials defined.
  //
  //  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4Box* b1 = new G4Box ( "b1", 100*cm, 50*cm, 50*cm );
  G4Box* b2 = new G4Box ( "b2", 50*cm, 100*cm, 50*cm );

  G4ThreeVector b1Xb2(50*cm,50*cm,50*cm);	// Offset for boolean tests

  fval = val ;

  if (val == "Sphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);
  }

  else if (val == "HalfSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 0*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 12*cm ;
      fTrackerR2 = 14*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 360*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "HalfHollowSphere")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 12*cm ;
      fTrackerR2 = 14*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 180*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 180*deg ;
      aVolume = new G4Sphere ("aHalfHollowSphere",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q1Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ1Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q2Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ2Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q3Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ3Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q4Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 0*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ4Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q5Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 0*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
else if (val == "Q6Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 90*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ5Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }  
  else if (val == "Q7Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 180*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ7Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Q8Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 270*deg ;
      fPhiSegment = 90*deg ;
      fTheta = 90*deg ;
      fThetaSegment = 90*deg ;
      aVolume = new G4Sphere ("aQ8Shell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }


  else if (val == "Shell")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 100*deg ;
      aVolume = new G4Sphere ("aShell",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if (val == "Ring")
  {
    // only solid sphere is supportet for now
      fTrackerR1 = 10*cm ;
      fTrackerR2 = 12*cm ;
      fPhi = 30*deg ;
      fPhiSegment = 120*deg ;
      fTheta = 10*deg ;
      fThetaSegment = 40*deg ;
      aVolume = new G4Sphere ("aRing",  fTrackerR1, fTrackerR2 ,fPhi, fPhiSegment, fTheta, fThetaSegment);

  }
  else if ( val == "Ellipsoid" ) {

    fSemiAxisX = 1*cm ;
    fSemiAxisY = 10*cm ;
    fSemiAxisZ = 100*cm ;

    aVolume = new  G4Ellipsoid("aEllipsoid",
			       fSemiAxisX,
			       fSemiAxisY,
			       fSemiAxisZ,
			       -fSemiAxisZ,fSemiAxisZ);


  }
  else if (val == "Orb")
  {

    fTrackerR = 11*cm ;
    aVolume = new G4Orb ( "aOrb", fTrackerR );

  }
  else if (val == "Box") 
  {         

    fTrackerpDx1 = 10*cm ;
    fTrackerpDy1 = 10*cm ;
    fTrackerpDz  = 10*cm ;

    aVolume = new G4Box ( "aBox", fTrackerpDx1, fTrackerpDy1, fTrackerpDz );
  }
  else if (val == "Cons")
  {        

    fTrackerpDz = 80*cm ;
    fTrackerR1 = 11*cm ;
    fTrackerR2 = 16*cm ;

    aVolume = new G4Cons ( "aCons", 0.8*fTrackerR1 , fTrackerR1, 0.8*fTrackerR2, fTrackerR2, fTrackerpDz, 10*deg, 120.*deg  ); 

  }
  else if (val == "manyCons")
  {        
    aVolume = new G4Cons ( "aCone", 2*cm, 6*cm, 8*cm, 14*cm,
                           10*cm, 10*deg, 300*deg ); 
  //  10*cm, 10*deg, 300*deg ); 
			   //  0., pi);

 
  }
  else if (val == "Tubs")
  {

    // only solid Tube is supported.
    fTrackerpDz = 80*cm ;
    fTrackerR = 11*cm ;

    aVolume = new G4Tubs ( "aTube",0.8*fTrackerR,fTrackerR,fTrackerpDz,40.,100*deg) ;

  }
  else if (val == "Hype")
  {
    aVolume = new G4Hype ("aHype", 7*cm, 10*cm, 40*deg, 40*deg, 40*cm );
  }
  else if (val == "Torus")
  {

    fPhi = 40*deg ;
    fPhiSegment = 100*deg ;

    fTrackerR1 = 5*cm ;
    fTrackerR2 = 6*cm ;
    fTrackerR  = 20*cm ;

    aVolume = new G4Torus("aTorus", fTrackerR1, fTrackerR2 ,fTrackerR, fPhi, fPhiSegment) ;

  }
  else if (val == "Para")
  {
    aVolume = new G4Para ("aPara", 8*cm, 10*cm, 12*cm, 30*deg, 45*deg, 60*deg);
  }
  else if (val == "Paraboloid")
  {
    aVolume = new G4Paraboloid ("aParaboloid", 8*cm, 1*cm, 12*cm);
  }
  else if (val == "Trd")
  {
    aVolume = new G4Trd ("aTrd", 80*cm, 100*cm, 70*cm, 90*cm, 100*cm);
  }
  else if (val == "b1Ub2") 
  {         
    aVolume = new G4UnionSolid("b1Ub2",b1,b2,0,b1Xb2);
  }
  else if (val == "b1Ib2") 
  {         
    aVolume = new G4IntersectionSolid("b1Ib2",b1,b2,0,b1Xb2);
  }
  else if (val == "b1Sb2") 
  {         
    aVolume = new G4SubtractionSolid("b1Sb2",b1,b2,0,b1Xb2);
  }
  else if (val == "b1Ib1") 
  {         
    aVolume = new G4IntersectionSolid("b1Ib1",b1,b1,0,b1Xb2);
  }
  else if (val == "b1Ub1") 
  {         
    aVolume = new G4UnionSolid("b1Ub1",b1,b1,0,b1Xb2);
  }
  else if (val == "b1Sb1") 
  {         
    aVolume = new G4SubtractionSolid("b1Sb1",b1,b1,0,b1Xb2);
  }
  else if ( val == "TwistedTubs" )
  {

    fTwistAngle = 20*deg ;
    fTrackerpDz = 80*cm ;
    fTrackerR1  = 5*cm ;
    fTrackerR2  = 10*cm ;
    fPhi        = 50*deg ;

    aVolume = new G4TwistedTubs("aTwistedTubs", fTwistAngle, fTrackerR1 , fTrackerR2, fTrackerpDz, fPhi ) ;
  }

  else if (val == "TwistedBox")
  {


    fTwistAngle = 20*deg ;
    fTrackerpDx1 = 11*cm ;
    fTrackerpDy1 = 8*cm ;
    fTrackerpDz  = 80*cm ;

    aVolume = new G4TwistedBox("aTwistedBox",fTwistAngle,fTrackerpDx1,fTrackerpDy1,fTrackerpDz) ;
  }
  else if (val == "TwistedTrd")
  {

    fTrackerpDx1 = 5*cm ;
    fTrackerpDx2 = 10*cm ;
    fTrackerpDy1 = 8*cm ;
    fTrackerpDy2 = 15*cm ;
    fTrackerpDz  = 80*cm ;
    fTwistAngle = 20*deg ;

    aVolume = new G4TwistedTrd("aTwistedTrd",fTrackerpDx1,fTrackerpDx2,fTrackerpDy1,fTrackerpDy2,fTrackerpDz,fTwistAngle);

  }
  else if ( val == "TwistedTrap") 
  {
    fTwistAngle = 60*deg ; 
    fTrackerpDz = 80*cm;
    fTheta = 10*deg ;
    fPhi  =  40*deg ;
    fTrackerpDy1 = 16*cm ;
    fTrackerpDx1 = 24*cm ;
    fTrackerpDx2 = 14*cm ;
    fTrackerpDy2 = 8*cm ;
    fTrackerpDx3 = 16*cm ;
    fTrackerpDx4 = 11*cm ;
    fAlph = 50*deg    ;

    aVolume = new G4TwistedTrap("aTwistedTrap",
				fTwistAngle,         // twist angle
				fTrackerpDz,         // half z length
				fTheta,              // direction between end planes
				fPhi,                // defined by polar and azimutal angles.
				fTrackerpDy1,        // half y length at -pDz
				fTrackerpDx1,        // half x length at -pDz,-pDy
				fTrackerpDx2,        // half x length at -pDz,+pDy
				fTrackerpDy2,        // half y length at +pDz
				fTrackerpDx3,        // half x length at +pDz,-pDy
				fTrackerpDx4,        // half x length at +pDz,+pDy
				fAlph                // tilt angle at +pDz
				) ;
  }
  else if ( val == "Tet" ) 
  {

      G4ThreeVector pzero(0,0,0);
      G4ThreeVector pnt1(10.*cm,0.*cm,0.*cm),pnt2(5.0*cm,10.*cm,0.*cm), pnt3(5.*cm,5.*cm,10.*cm);
      G4bool  goodTet;
      G4Tet   t1( "aTet", pzero, pnt1, pnt2, pnt3, &goodTet);
  }
  else if ( val == "Trap") 
  {
    fTrackerpDz = 80*cm;
    fTheta = 10*deg ;
    fPhi  =  40*deg ;
    fTrackerpDy1 = 16*cm ;
    fTrackerpDx1 = 24*cm ;
    fTrackerpDx2 = 14*cm ;
    fTrackerpDy2 = 8*cm ;
    fTrackerpDx3 = 16*cm ;
    fTrackerpDx4 = 11*cm ;
    fAlph = 50*deg    ;

    aVolume = new G4Trap("aTrap",
				fTrackerpDz,         // half z length
				fTheta,              // direction between end planes
				fPhi,                // defined by polar and azimutal angles.
				fTrackerpDy1,        // half y length at -pDz
				fTrackerpDx1,        // half x length at -pDz,-pDy
				fTrackerpDx2,        // half x length at -pDz,+pDy
				fAlph,                // tilt angle at +pDz
				fTrackerpDy2,        // half y length at +pDz
				fTrackerpDx3,        // half x length at +pDz,-pDy
				fTrackerpDx4,        // half x length at +pDz,+pDy
				fAlph                // tilt angle at +pDz
				) ;
  }
  else if ( val == "EllipticalCone" ) 
  {
    aVolume = new G4EllipticalCone("aEllipticalCone",
                        0.5*mm,       // xSemiAxis
                        1*mm,       // ySemiAxis
                        40*mm,      // zheight
                        25*mm) ;    // zTopCut

  }
  else if ( val == "EllipticalTube" ) 
  {
    aVolume = new G4EllipticalTube("aEllipticalTube" ,
				   2*cm,   // xSemiAxis
				   5*cm,   // ySemiAxis
				   35*cm) ;  // zheight

  }
  else
  {
    G4Exception("Sc01DetectorConstruction::SelectDetector()", "Sc01DC001",
		FatalException, "Invalid shape!");
  }

  fWorldLength= 10*m ;
   
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
  solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
				 


  G4LogicalVolume* aVolume_log = new G4LogicalVolume(aVolume, Air, "aVolume_L", 0,0,0);

  // G4VPhysicalVolume * aVolume_phys1 =
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, 0*cm),
                        aVolume_log, 
			val, 
			logicWorld, 
			false,
			0);



//--------- Visualization attributes -------------------------------


// the world is transparent
  G4VisAttributes* WorldAtt = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
  WorldAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(WorldAtt);  

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(0.0,.0,1.0,0.6));
  BoxVisAtt->SetVisibility(true);
  aVolume_log->SetVisAttributes(BoxVisAtt);
  
//--------- example of User Limits -------------------------------

  
  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* SCDetectorConstruction::Construct()
{


  //-------------------Hall ----------------------------------
  
  return SelectDetector ("Sphere");  // default is Sphere

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void SCDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
