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
// DetectorConstruction program
// --------------------------------------------------------------

#include "Test15DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// #include "G4UserLimits.hh"

#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
Test15DetectorConstruction::Test15DetectorConstruction() : olap_test(false)
{
  // create commands for interactive definition of time cuts:


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
Test15DetectorConstruction::~Test15DetectorConstruction() 
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::DefineMaterials() 
{

  G4double density,      // density
    a,                   // atomic mass
    z;                   // atomic number
  G4String name,         // name
    symbol;              // symbol
  G4int ncomponents,     // n components
    iz,                  // number of protons
    in;                  // number of nuceons
  G4double abundance,    // abundance
    temperature,         // temperature
    pressure;            // pressure

  // making vacuum
  G4Material* vacuum = new G4Material 
    (name="Vacuum", z=1., a=1.*g/mole, density=1.e-20*g/cm3,
     kStateGas, temperature=0.1*kelvin, pressure=1.e-20*bar);


  // air
  G4Element* N = new G4Element
    (name="Nitrogen",symbol="N" , z= 7., a=14.00674*g/mole);
  G4Element* O  = new G4Element
    (name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
  G4Material* Air = new G4Material
    ("AIR", 1.2929*kg/m3, 2, kStateGas, 300.00*kelvin, 1.0*atmosphere);
  Air->AddElement(N, 0.8);
  Air->AddElement(O , 0.2);

  // copper
  G4Element* Cu = new G4Element
    (name="Copper"  ,symbol="Cu" , z= 29., a=63.55*g/mole);  
  G4Material* metalCu = new G4Material
    (name="MetalCopper", density=8.960*g/cm3, ncomponents=1);
  metalCu->AddElement(Cu, 1);

  // lead
  G4Element* Pb = new G4Element
    (name="Lead",symbol="Pb" , z= 82., a=207.2*g/mole);
  G4Material* metalPb = new G4Material
    (name="MetalLead", density=11.340*g/cm3, ncomponents=1);
  metalPb->AddElement(Pb, 1);


  // lead
  G4Isotope* Lead204 = new G4Isotope(name="Lead204", iz=82, in=204, a=204.0*g/mole);
  G4Isotope* Lead206 = new G4Isotope(name="Lead206", iz=82, in=206, a=206.0*g/mole);
  G4Isotope* Lead207 = new G4Isotope(name="Lead207", iz=82, in=207, a=207.0*g/mole);
  G4Isotope* Lead208 = new G4Isotope(name="Lead208", iz=82, in=208, a=208.0*g/mole);

  G4Element* Pb204 = new G4Element
    (name="Pb204",symbol="Pb204" , ncomponents=1);
  Pb204->AddIsotope(Lead204, abundance=1);
  G4Element* Pb206 = new G4Element
    (name="Pb206",symbol="Pb206" , ncomponents=1);
  Pb206->AddIsotope(Lead206, abundance=1);
  G4Element* Pb207 = new G4Element
    (name="Pb207",symbol="Pb207" , ncomponents=1);
  Pb207->AddIsotope(Lead207, abundance=1);
  G4Element* Pb208 = new G4Element
    (name="Pb208",symbol="Pb208" , ncomponents=1);
  Pb208->AddIsotope(Lead208, abundance=1);
  G4Element* PbNat = new G4Element
    (name="PbNat",symbol="PbNat" , ncomponents=4);
  PbNat->AddIsotope(Lead204, abundance=0.014);
  PbNat->AddIsotope(Lead206, abundance=0.241);
  PbNat->AddIsotope(Lead207, abundance=0.221);
  PbNat->AddIsotope(Lead208, abundance=0.524);

//  PbNat->AddIsotope(Lead204, abundance=1.4);
//  PbNat->AddIsotope(Lead206, abundance=24.1);
//  PbNat->AddIsotope(Lead207, abundance=22.1);
//  PbNat->AddIsotope(Lead208, abundance=52.4);

  G4Material* metalPb204 = new G4Material
    (name="MetalLead204", density=11.340*g/cm3, ncomponents=1);
  metalPb204->AddElement(Pb204, 1);
  G4Material* metalPb206 = new G4Material
    (name="MetalLead206", density=11.340*g/cm3, ncomponents=1);
  metalPb206->AddElement(Pb206, 1);
  G4Material* metalPb207 = new G4Material
    (name="MetalLead207", density=11.340*g/cm3, ncomponents=1);
  metalPb207->AddElement(Pb207, 1);
  G4Material* metalPb208 = new G4Material
    (name="MetalLead208", density=11.340*g/cm3, ncomponents=1);
  metalPb208->AddElement(Pb208, 1);
  G4Material* metalPbNat = new G4Material
    (name="MetalLeadNat", density=11.340*g/cm3, ncomponents=1);
  metalPbNat->AddElement(PbNat, 1);

  G4Isotope* He3 = new G4Isotope
    (name="Helium3", iz= 2, in=3, a=3.0*g/mole);
  G4Element* He = new G4Element
    (name="Helium", "He", ncomponents=1);
  He->AddIsotope(He3, abundance=1);
  G4Material* helium3 = new G4Material
    (name="helium3", density= 0.1785*kg/m3, ncomponents=1);
  helium3->AddElement(He, 1);

  G4Element* Si = new G4Element
    (name="Silicon",symbol="Si" , z= 14., a=28.09*g/mole);
  G4Material* silicon = new G4Material
    (name="silicon", density=3.000*g/cm3, ncomponents=1);
  silicon->AddElement(Si, 1);

  world_mat = vacuum;
  lab_mat = Air;
  lead_mat = metalPbNat;
  sample_mat = metalCu;
  sampleSphere_mat = metalPbNat;
  sampleTube_mat = metalPbNat;

  //xx    sampleSphere_mat = helium3;
  //xx       lead_mat = helium3;
  //       lead_mat = metalCu;
  //    sampleTube_mat = helium3;
  //    sampleTube_mat = metalPb;
  //    sampleTube_mat = silicon;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* Test15DetectorConstruction::Construct() {

  DefineMaterials();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;
  

  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;



  // Universe - room wall - CONCRETE ****************************************

  //NB: measured INSIDE of lab, therefore have to add twice wall thickness
  G4double wallThick   = 24.*cm;
  G4double worldWidth  = 5.0*m + 2.*wallThick; // "x"
  G4double worldLength = 5.0*m + 2.*wallThick; // "y"
  G4double worldHeight = 5.0*m + 2.*wallThick; // "z"

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, world_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "world_phys", world_log, NULL, false,0);

  //  G4VisAttributes* world_vat= new G4VisAttributes(white);
  world_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //world_vat->SetVisibility(true);
  //world_vat->SetVisibility(false);
  //world_log->SetVisAttributes(world_vat);


  // Lab Space - AIR ********************************************************

  G4double labWidth  = worldWidth  - 2.*wallThick; //X
  G4double labLength = worldLength - 2.*wallThick; //Y
  G4double labHeight = worldHeight - 2.*wallThick; //Z

  G4Box* lab_box = new G4Box
     ("lab_box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight );
  lab_log  = new G4LogicalVolume(lab_box, lab_mat, "lab_log");
  lab_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lab_phys", 
			       lab_log, world_phys, false,0, olap_test);

  G4VisAttributes* lab_vat= new G4VisAttributes(white);
  //  lab_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  lab_vat->SetVisibility(true);
  lab_vat->SetVisibility(false);
  lab_log->SetVisAttributes(lab_vat);

  // Now start with detector assembly:

  // Lead Block

  G4double leadWidth  = 3.3*m; //X
  G4double leadHeight = 3.3*m; //Y
  G4double leadLength = 3.0*m; //Z

  G4Box* lead_box = new G4Box
     ("lead_box", 0.5*leadWidth, 0.5*leadHeight, 0.5*leadLength );
  lead_log  = new G4LogicalVolume(lead_box, lab_mat, "lead_log");
  lead_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lead_phys", 
     lead_log, lab_phys, false,0, olap_test);

  G4VisAttributes* lead_vat= new G4VisAttributes(grey);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  lead_vat->SetVisibility(true);
  lead_log->SetVisAttributes(lead_vat);

//--------------------------------------------------------------
//--------------------------------------------------------------


  // Lead Blocks
  // Type A
  blockWidthA  = 30.0*cm; //X
  blockHeightA = 30.0*cm; //Y
  blockLengthA = 60.0*cm; //Z
  blockWidthB  = blockWidthA; //X
  blockHeightB = blockHeightA; //Y
  blockLengthB = blockLengthA; //Z
  blockWidthC  = 15.0*cm; //X
  blockHeightC = 30.0*cm; //Y
  blockLengthC = 60.0*cm; //Z


  G4Box* blockA_box = new G4Box
    ("blockA_box", 0.5*blockWidthA-1.0*nanometer, 0.5*blockHeightA, 0.5*blockLengthA );
  //     ("blockA_box", 0.5*blockWidthA, 0.5*blockHeightA, 0.5*blockLengthA );

  blockA_log  = new G4LogicalVolume(blockA_box, lead_mat, "blockA_log");


  G4VisAttributes* blockA_vat= new G4VisAttributes(blue);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  blockA_vat->SetVisibility(true);
  //xx blockA_vat->SetForceSolid(true);
  blockA_log->SetVisAttributes(blockA_vat);

//--------------------------------------------------------------
  // Type B

  G4double blockB_holeRadius = 32.0*mm;
  G4double blockB_beamholeDiameter = 77.2*mm;
//   G4RotationMatrix rotblockBhole;
//   rotblockBhole.rotateX(90.0*deg);
  // tubs is from 0 to 360 degrees to remove problem of exact interface between block and cylinder hole
  G4Tubs* blockB_hole = new G4Tubs("blockB_hole",0.*cm, blockB_holeRadius,0.5*blockLengthA+1.0*mm,0.*deg,360.*deg);
  G4Tubs* blockB_beamhole = new G4Tubs("blockB_beamhole",0.*cm, 0.5*blockB_beamholeDiameter,0.5*blockLengthA+1.0*mm,0.*deg,360.*deg);
  //  G4SubtractionSolid* blockB_box = new G4SubtractionSolid("blockB_box",blockA_box,blockB_hole,G4Transform3D(rotblockBhole, G4ThreeVector(0., 0.5*blockWidthA, 0.)));
  G4SubtractionSolid* blockB_box = new G4SubtractionSolid("blockB_box",blockA_box,blockB_hole,G4Transform3D(G4RotationMatrix(), G4ThreeVector(-0.5*blockWidthA, 0., 0.)));
  G4SubtractionSolid* blockB_beambox = new G4SubtractionSolid("blockB_beambox",blockB_box,blockB_beamhole,G4Transform3D(G4RotationMatrix(), G4ThreeVector(0.0, 0., 0.)));

  blockB_log  = new G4LogicalVolume(blockB_box, lead_mat, "blockB_log");
  beamblockB_log  = new G4LogicalVolume(blockB_beambox, lead_mat, "beamblockB_log");


  G4VisAttributes* blockB_vat= new G4VisAttributes(yellow);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  blockB_vat->SetVisibility(true);
  //xx blockB_vat->SetForceSolid(true);
  blockB_log->SetVisAttributes(blockB_vat);

  G4VisAttributes* beamblockB_vat= new G4VisAttributes(brown);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  beamblockB_vat->SetVisibility(true);
  //xx beamblockB_vat->SetForceSolid(true);
  beamblockB_log->SetVisAttributes(beamblockB_vat);

//--------------------------------------------------------------
  // Type C

  G4Box* blockC_box = new G4Box
     ("blockC_box", 0.5*blockWidthC, 0.5*blockHeightC, 0.5*blockLengthC );

  blockC_log  = new G4LogicalVolume(blockC_box, lead_mat, "blockC_log");


  G4VisAttributes* blockC_vat= new G4VisAttributes(green);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  blockC_vat->SetVisibility(true);
  //xx blockC_vat->SetForceSolid(true);
  blockC_log->SetVisAttributes(blockC_vat);

//--------------------------------------------------------------
// ......................................................................
  // Sample Block

  G4double sampleTubeRadius  = blockB_holeRadius;
  G4double sampleTubeLength = 15.0*cm;
  //xtest  G4double sampleTubeLength = 2.0*blockB_holeRadius;
  G4double sampleSphereRadius  = blockB_holeRadius;
  G4double sampleSphereXpos = -0.45*m;
  G4double sampleSphereYpos = 0.0*m;
  G4double sampleSphereZpos = -7.5*cm;
  G4double sampleTubeXpos = -0.45*m;
  G4double sampleTubeYpos = 0.0*m;
  G4double sampleTubeZpos = 7.5*cm;
  G4double sampleSphereXpos2 = 0.0*m;
  G4double sampleSphereYpos2 = 0.0*m;
  G4double sampleSphereZpos2 = -14.4*cm;

  G4Tubs* sampleTube = new G4Tubs("sampleTube",0.*cm, sampleTubeRadius,0.5*sampleTubeLength,0.*deg,360.*deg);
  sampleTube_log  = new G4LogicalVolume(sampleTube, sampleTube_mat, "sampleTube_log");

  G4Sphere* sampleSphere = new G4Sphere("sampleSphere",0.*cm,sampleSphereRadius,0.*deg, 360.*deg, 0.*deg, 180.*deg);

  sampleSphere_log  = new G4LogicalVolume(sampleSphere, sampleSphere_mat, "sampleSphere_log");
  sampleSphere_log2  = new G4LogicalVolume(sampleSphere, sampleSphere_mat, "sampleSphere_log2");

  sample_phys = new G4PVPlacement(0, G4ThreeVector(sampleSphereXpos, sampleSphereYpos, sampleSphereZpos), "sample_phys", sampleSphere_log, lead_phys, false,0, olap_test);

//   sample_phys2 = new G4PVPlacement(0, G4ThreeVector(sampleSphereXpos2, sampleSphereYpos2, sampleSphereZpos2), "sample_phys2", sampleSphere_log2, blockB_phys, false,0);
// have to use logical mother as it's sitting inside multiple placements...........!
  sample_phys2 = new G4PVPlacement(0, G4ThreeVector(sampleSphereXpos2, sampleSphereYpos2, sampleSphereZpos2), sampleSphere_log2, "sample_phys2", blockB_log, false,0, olap_test);

  sampleTube_phys = new G4PVPlacement(0, G4ThreeVector(sampleTubeXpos,sampleTubeYpos,sampleTubeZpos), "sampleTube_phys", sampleTube_log, lead_phys, false,0, olap_test);

  G4VisAttributes* sampleSphere_vat= new G4VisAttributes(magenta);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  sampleSphere_vat->SetVisibility(true);
  sampleSphere_vat->SetForceSolid(true);
  sampleSphere_log->SetVisAttributes(sampleSphere_vat);
  sampleTube_log->SetVisAttributes(sampleSphere_vat);

//--------------------------------------------------------------

  CreateLayer1();
  CreateLayer2();
  CreateLayer3();
  CreateLayer4();
  CreateLayer5();
  CreateLayer6();
  CreateLayer7();
  CreateLayer8();
  CreateLayer9();
  CreateLayer10();
  CreateLayer11();

//--------------------------------------------------------------
  // Sample Block

  G4double sampleWidth  = 10.0*cm; //X
  G4double sampleLength = 4.0*cm; //Y
  G4double sampleHeight = 3.0*cm; //Z

  G4Box* sample_box = new G4Box
     ("sample_box", 0.5*sampleWidth, 0.5*sampleLength, 0.5*sampleHeight );
  sample_log  = new G4LogicalVolume(sample_box, sample_mat, "sample_log");
//   sample_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "sample_phys", 
//      sample_log, lead_phys, false,0);

  G4VisAttributes* sample_vat= new G4VisAttributes(red);
  //  lead_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  sample_vat->SetVisibility(true);
  sample_log->SetVisAttributes(sample_vat);

  // ......................................................................
  // attach user limits ...................................................

  
  // G4cout << G4endl << "User Limits: " << G4endl 
  // 	 << "\t theMaxTimeCuts:     " << G4BestUnit(theMaxTimeCuts,"Time")  
  // 	 << G4endl
  // 	 << "\t theMaxStepSize:     " << G4BestUnit(theMaxStepSize,"Length")
  // 	 << G4endl
  // 	 << "\t theMinEKine:        " << G4BestUnit(theMinEkine,"Energy")   
  // 	 << G4endl;

  // if (theUserLimitsForDetector != 0) delete theUserLimitsForDetector;

  // theUserLimitsForDetector = new G4UserLimits(theDetectorStepSize,
  // 					      DBL_MAX, // Track Max
  // 					      theMaxTimeCuts,
  // 					      theMinEkine);

  //     world_log->SetUserLimits(theUserLimitsForDetector);
  //       lab_log->SetUserLimits(theUserLimitsForDetector);
  //       lead_log->SetUserLimits(theUserLimitsForDetector);
  //       sample_log->SetUserLimits(theUserLimitsForDetector);


  return world_phys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer1() 
{

  //  olap_test = true;
  // olap_test = false;
  replica_idx_A = 0;
  replica_idx_B = 0;
  replica_idx_C = 0;
  replica_sample = 0;
  // layer 1 column 1,2 left
  G4double blockA_x = blockWidthB + blockWidthC + 1.5*blockWidthA;
  G4double blockA_y = -2.0*blockHeightB - 3.0*blockHeightA;
  G4double blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  //	  CheckOverlaps();	  
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x -= blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }


  // layer 1 column 1,2 right
  blockA_x = -(blockWidthB + blockWidthC + 1.5*blockWidthA);
  blockA_y = -2.0*blockHeightB - 3.0*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x += blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }


  // layer 1 column 4,5 central
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = -2.0*blockHeightB - 3.0*blockHeightA;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 4;
	    //xxx  	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys4", 
	    //xxx  						sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


  // layer 1 column 3 left,right
  G4double blockC_x = blockWidthB + 0.5*blockWidthC;
  G4double blockC_y = -2.0*blockHeightB - 3.0*blockHeightA;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockWidthB + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthB;
    }


}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer2() 
{

  // layer 2 column 1,2 left
  G4double blockA_x = 0.5*blockWidthA + 1.5*blockLengthA;
  G4double blockA_y = -1.5*blockHeightB - 2.5*blockHeightA;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x -= blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }


  // layer 2 column 1,2 right
  blockA_x = -1.0*(0.5*blockWidthA + 1.5*blockLengthA);
  blockA_y = -1.5*blockHeightB - 2.5*blockHeightA;
  blockA_z = -4.5*blockWidthA;
  // rotation matrix is additive so don't rotate again
  // rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
 	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 	  				  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x += blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }

  // layer 2 central column
  blockA_x = 0.0*m;
  blockA_y = -1.5*blockHeightB - 2.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int i=0; i<5; i++)
    {	  
      blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
				      blockA_log, lead_phys, false,replica_idx_A,olap_test);
      replica_idx_A++;
      blockA_z += blockLengthA;
    }

}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer3() 
{

  // layer 3 column 1 left,right
  G4double blockA_x = 3.5*blockWidthA + 0.5*blockLengthA;
  G4double blockA_y = -1.5*blockHeightB - 1.5*blockHeightA;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(3.5*blockWidthA + 0.5*blockLengthA);
      blockA_z = -4.5*blockWidthA;
    }

  // layer 3 column 1,2,3 left
  blockA_x = 3.0*blockWidthA;
  blockA_y = -1.5*blockHeightB - 1.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<3; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
 	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x -= blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }

  // layer 3 column 1,2,3 right
  blockA_x = -3.0*blockWidthA;
  blockA_y = -1.5*blockHeightB - 1.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<3; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
 	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x += blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }

  // layer 3 central column
  blockA_x = 0.0*m;
  blockA_y = -1.5*blockHeightB - 1.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int i=0; i<5; i++)
    {	  
      blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
				      blockA_log, lead_phys, false,replica_idx_A,olap_test);
      replica_idx_A++;
      blockA_z += blockLengthA;
    }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer4() 
{

  // layer 4 column 1 left,right
  G4double blockC_x = blockWidthB + 2.0*blockLengthA + 0.5*blockWidthC;
  G4double blockC_y = -1.0*blockHeightB - 1.0*blockHeightA;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockWidthB + 2.0*blockLengthA + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthC;
    }


  // layer 4 column 2,3 left
  G4double blockA_x = blockWidthB + 1.5*blockLengthA;
  G4double blockA_y = -1.0*blockHeightB - 1.0*blockHeightA;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x -= blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }


  // layer 4 column 2,3 right
  blockA_x = -1.0*(blockWidthB + 1.5*blockLengthA);
  blockA_y = -1.0*blockHeightB - 1.0*blockHeightA;
  blockA_z = -4.5*blockWidthA;
  // don't need to rotate again!
//   G4RotationMatrix rotblockA;
//   rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x += blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }


  // layer 4 column 4,5 central
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = -1.0*blockHeightB - 1.0*blockHeightA;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 5;
	    //xxx  	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys5", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer5() 
{

  // layer 5 column 1 left,right
  G4double blockA_x = 1.5*blockWidthA + 1.5*blockLengthA;
  G4double blockA_y = -0.5*blockHeightB - 0.5*blockHeightA;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(1.5*blockWidthA + 1.5*blockLengthA);
      blockA_z = -4.5*blockWidthA;
    }

  // layer 5 column 2 left,right
  blockA_x = blockWidthA + blockLengthA;
  blockA_y = -0.5*blockHeightB - 0.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x = -1.0*(blockWidthA + blockLengthA);
      blockA_z = -2.0*blockLengthA;
    }
     
  // layer 5 column 3 left,right
  blockA_x = 0.5*blockWidthA + 0.5*blockLengthA;
  blockA_y = -0.5*blockHeightB - 0.5*blockHeightA;
  blockA_z = -4.5*blockWidthA;
// don't need to rotate matrix again!
//   G4RotationMatrix rotblockA;
//   rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(0.5*blockWidthA + 0.5*blockLengthA);
      blockA_z = -4.5*blockWidthA;
    }

  // layer 5 column 4 central
  blockA_x = 0.0*m;
  blockA_y = -0.5*blockHeightB - 0.5*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int i=0; i<5; i++)
    {	  
      blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
				      blockA_log, lead_phys, false,replica_idx_A,olap_test);
      replica_idx_A++;
      blockA_z += blockLengthA;
    }
     
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer6() 
{

  // layer 6 column 1 left,right
  G4double blockA_x = 3.5*blockWidthB + 1.5*blockWidthA;
  G4double blockA_y = 0.0*m;
  G4double blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x = -1.0*(3.5*blockWidthB + 1.5*blockWidthA);
      blockA_z = -2.0*blockLengthA;
    }


  // layer 6 column 2,3 left
  G4double blockB_x = 3.0*blockWidthB + blockWidthA;
  G4double blockB_y = 0.0*m;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 1;
	    //xxx  sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys1", 
	    //xxx				sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }

  // layer 6 column 2,3 right
  blockB_x = -4.0*blockWidthB;
  blockB_y = 0.0*m;
  blockB_z = -2.0*blockLengthA;
  // extra 180 degreee rotation due to the previous rotation above
  rotblockB.rotateZ(180.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 12;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x-0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys12", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x += blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }

  // layer 6 column 4,5 right
  blockB_x = -2.0*blockWidthB;
  blockB_y = 0.0*m;
  blockB_z = -2.0*blockLengthA;
  // extra 180 degreee rotation not needed as following from 2,3 right above
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 10;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x-0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys10", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x += blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }

  // layer 6 column central
  blockB_x = 0.0*m;
  blockB_y = 0.0*m;
  blockB_z = -2.0*blockLengthA;
  // extra 180 degreee rotation not needed as following 5, right above
  rotblockB.rotateZ(0.0*deg);

  beamblockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "beamblockB_phys", 
				      beamblockB_log, lead_phys, false,0,olap_test);
  blockB_z += blockLengthB;
  beamblockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "beamblockB_phys", 
				      beamblockB_log, lead_phys, false,1,olap_test);
  replica_sample = 3;
  //xxx  sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys3", 
  //xxx				      sampleTube_log, lead_phys, false, replica_sample, olap_test);
  //  replica_sample++;
  blockB_z += blockLengthB;
  for (G4int i=2; i<5; i++)
    {	  
      blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
				      blockB_log, lead_phys, false,replica_idx_B,olap_test);
      replica_idx_B++;
      blockB_z += blockLengthB;
      
      G4cout << " BLOCKB: " << i << " x: " << blockB_x << " y: " << blockB_y << " z: " << blockB_z << " replica: " << replica_idx_B << G4endl;


    }


  // layer 6 column 4 left
  blockA_x = 1.5*blockWidthB + 0.5*blockWidthA;
  blockA_y = 0.0*m;
  blockA_z = -2.0*blockLengthA;
  for (G4int i=0; i<5; i++)
    {	  
      blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
				      blockA_log, lead_phys, false,replica_idx_A,olap_test);
      replica_idx_A++;
      blockA_z += blockLengthA;
    }

  // layer 6 column 5 left
  blockB_x = 1.0*blockWidthB;
  blockB_y = 0.0*m;
  blockB_z = -2.0*blockLengthA;
  // extra 180 degreee rotation not needed as following from 2,3 right above
  rotblockB.rotateZ(180.0*deg);
  for (G4int i=0; i<5; i++)
    {	  
      blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
				      blockB_log, lead_phys, false,replica_idx_B,olap_test);
      replica_idx_B++;
      blockB_z += blockLengthB;
    }


}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer7() 
{

  // layer 7 column 1 left,right
  G4double blockA_x = blockLengthA + 2.0*blockWidthB + blockWidthC;
  G4double blockA_y = 0.5*blockHeightA + 0.5*blockHeightB;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(blockLengthA + 2.0*blockWidthB + blockWidthC);
      blockA_z = -4.5*blockWidthA;
    }


  // layer 7 column 2 left,right
  G4double blockC_x = 0.5*blockLengthA + 2.0*blockWidthB + 0.5*blockWidthC;
  G4double blockC_y = 0.5*blockHeightA + 0.5*blockHeightB;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(0.5*blockLengthA + 2.0*blockWidthB + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthC;
    }


  // layer 7 column 3,4 left
  G4double blockB_x = 0.5*blockLengthA + 1.5*blockWidthB;
  G4double blockB_y =  0.5*blockHeightA + 0.5*blockHeightB;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 2;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys2", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


  // layer 7 column 3,4 right
  blockB_x = -1.0*(0.5*blockLengthA + 1.5*blockWidthB);
  blockB_y =  0.5*blockHeightA + 0.5*blockHeightB;
  blockB_z = -2.0*blockLengthA;
  rotblockB.rotateZ(180.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 11;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x-0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys11", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x += blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


  // layer 7 central
  blockA_x = 0.0*m;
  blockA_y = 0.5*blockHeightA + 0.5*blockHeightB;
  blockA_z = -4.5*blockWidthA;
  // block already rotated from column 1 above
  rotblockA.rotateY(0.0*deg);
  for (G4int i=0; i<10; i++)
    {	  
      blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
				      blockA_log, lead_phys, false,replica_idx_A,olap_test);
      replica_idx_A++;
      blockA_z += blockWidthA;
    }


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer8() 
{

  // layer 8 column 1 left,right
  G4double blockC_x = blockWidthB + 2.0*blockLengthA + 0.5*blockWidthC;
  G4double blockC_y = 1.0*blockHeightB + 1.0*blockHeightA;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockWidthB + 2.0*blockLengthA + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthC;
    }


  // layer 8 column 2,3 left
  G4double blockA_x = blockWidthB + 1.5*blockLengthA;
  G4double blockA_y = 1.0*blockHeightB + 1.0*blockHeightA;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x -= blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }


  // layer 8 column 2,3 right
  blockA_x = -1.0*(blockWidthB + 1.5*blockLengthA);
  blockA_y = 1.0*blockHeightB + 1.0*blockHeightA;
  blockA_z = -4.5*blockWidthA;
  // don't need to rotate again!
//   G4RotationMatrix rotblockA;
//   rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x += blockLengthA;
      blockA_z = -4.5*blockWidthA;
    }


  // layer 8 column 4,5 central
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = 1.0*blockHeightB + 1.0*blockHeightA;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 6;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys6", 
	    //xxx					sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer9() 
{

  // layer 9 column 1 left,right
  G4double blockA_x = blockLengthA + 2.0*blockWidthB + blockWidthC;
  G4double blockA_y = blockHeightA + 2.0*blockHeightB;
  G4double blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(blockLengthA + 2.0*blockWidthB + blockWidthC);
      blockA_z = -4.5*blockWidthA;
    }

  // layer 9 column 2 left,right
  blockA_x = 1.5*blockWidthA + blockWidthB + blockWidthC;
  blockA_y = blockHeightA + 2.0*blockHeightB;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x = -1.0*(1.5*blockWidthA + blockWidthB + blockWidthC);
      blockA_z = -2.0*blockLengthA;
    }

  // layer 9 column 3 left,right
  G4double blockC_x = blockWidthA + blockWidthB + 0.5*blockWidthC;
  G4double blockC_y = blockHeightA + 2.0*blockHeightB;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockWidthA + blockWidthB + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthC;
    }


  // layer 9 column 4 left,right
  blockA_x = 0.5*blockWidthA + blockWidthB;
  blockA_y = blockHeightA + 2.0*blockHeightB;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x = -1.0*(0.5*blockWidthA + blockWidthB);
      blockA_z = -2.0*blockLengthA;
    }


  // layer 9 column central left,right
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = blockHeightA + 2.0*blockHeightB;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 7;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phys7", 
	    //xxx				sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }



}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer10() 
{

  // layer 10 column 1 left,right
  G4double blockA_x = 0.5*blockWidthA + blockLengthA + blockWidthB + blockWidthC;
  G4double blockA_y = blockHeightA + 3.0*blockHeightB;
  G4double blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x = -1.0*(0.5*blockWidthA + blockLengthA + blockWidthB + blockWidthC);
      blockA_z = -2.0*blockLengthA;
    }


  // layer 10 column 2 left,right
  G4double blockC_x = blockLengthA + blockWidthB + 0.5*blockWidthC;
  G4double blockC_y = blockHeightA + 3.0*blockHeightB;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockLengthA + blockWidthB + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthC;
    }


  // layer 10 column 3 left,right
  blockA_x = 0.5*blockLengthA + blockWidthB;
  blockA_y = blockHeightA + 3.0*blockHeightB;
  blockA_z = -4.5*blockWidthA;
  G4RotationMatrix rotblockA;
  rotblockA.rotateY(90.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<10; i++)
	{	  
	  blockA_phys = new G4PVPlacement(G4Transform3D(rotblockA, G4ThreeVector(blockA_x,blockA_y,blockA_z)), "blockA_phys", 
 					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockWidthA;
	}
      blockA_x = -1.0*(0.5*blockLengthA + blockWidthB);
      blockA_z = -4.5*blockWidthA;
    }

  // layer 10 column central left,right
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = blockHeightA + 3.0*blockHeightB;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 8;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phy8", 
	    //xxx				sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Test15DetectorConstruction::CreateLayer11() 
{

  // layer 11 column 1,2 left
  G4double blockA_x = blockWidthB + blockWidthC + 1.5*blockWidthA;
  G4double blockA_y = 2.0*blockHeightB + 3.0*blockHeightA;
  G4double blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x -= blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }


  // layer 11 column 1,2 right
  blockA_x = -(blockWidthB + blockWidthC + 1.5*blockWidthA);
  blockA_y = 2.0*blockHeightB + 3.0*blockHeightA;
  blockA_z = -2.0*blockLengthA;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockA_phys = new G4PVPlacement(0, G4ThreeVector(blockA_x,blockA_y,blockA_z), "blockA_phys", 
					  blockA_log, lead_phys, false,replica_idx_A,olap_test);
	  replica_idx_A++;
	  blockA_z += blockLengthA;
	}
      blockA_x += blockWidthA;
      blockA_z = -2.0*blockLengthA;
    }


  // layer 11 column 4,5 central
  G4double blockB_x = 0.5*blockWidthB;
  G4double blockB_y = 2.0*blockHeightB + 3.0*blockHeightA;
  G4double blockB_z = -2.0*blockLengthA;
  G4RotationMatrix rotblockB;
  rotblockB.rotateZ(0.0*deg);
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockB_phys = new G4PVPlacement(G4Transform3D(rotblockB, G4ThreeVector(blockB_x,blockB_y,blockB_z)), "blockB_phys", 
					  blockB_log, lead_phys, false,replica_idx_B,olap_test);
	  replica_idx_B++;
	  if(j == 1) {
	    replica_sample = 9;
	    //xxx   	    sampleTube_phys = new G4PVPlacement(0, G4ThreeVector((blockB_x+0.5*blockWidthB),blockB_y,blockB_z), "sampleTube_phy9", 
	    //xxx				sampleTube_log, lead_phys, false, replica_sample, olap_test);
	    //	    replica_sample++;
	  }
	  blockB_z += blockLengthB;
	}
      rotblockB.rotateZ(180.0*deg);
      blockB_x -= blockWidthB;
      blockB_z = -2.0*blockLengthB;
    }


  // layer 11 column 3 left,right
  G4double blockC_x = blockWidthB + 0.5*blockWidthC;
  G4double blockC_y = 2.0*blockHeightB + 3.0*blockHeightA;
  G4double blockC_z = -2.0*blockLengthC;
  for (G4int j=0; j<2; j++)
    {
      for (G4int i=0; i<5; i++)
	{	  
	  blockC_phys = new G4PVPlacement(0, G4ThreeVector(blockC_x,blockC_y,blockC_z), "blockC_phys", 
					  blockC_log, lead_phys, false,replica_idx_C,olap_test);
	  replica_idx_C++;
	  blockC_z += blockLengthC;
	}
      blockC_x = -1.0*(blockWidthB + 0.5*blockWidthC);
      blockC_z = -2.0*blockLengthB;
    }


  G4cout << G4endl
	 << G4endl
	 << G4endl
	 << " Got to end of lead block construction " << G4endl
	 << G4endl
	 << G4endl
	 << G4endl
	 << G4endl
	 << " and number of blocks is: " << G4endl
	 << " Block Type A: " << replica_idx_A << G4endl
	 << " Block Type B: " << replica_idx_B << G4endl
	 << " Block Type C: " << replica_idx_C << G4endl
	 << G4endl
	 << G4endl
	 << G4endl
	 << G4endl;

}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
// void Test15DetectorConstruction::SetEnergyCut(G4double val)
// {
//   // set minimum charged particle energy cut - NB: for Xenon Detector
//   theMinEkine = val;
//   if (theUserLimitsForDetector != 0) 
//     {
//       theUserLimitsForDetector->SetUserMinEkine(val);
//       G4cout << "Changing Detector energy cut to: " << G4BestUnit(val,"Energy")
// 	     << G4endl;
//     }
// }  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
// void Test15DetectorConstruction::SetTimeCut(G4double val)
// {
//   // set detector time cut:
//   theMaxTimeCuts = val;
//   if (theUserLimitsForDetector != 0) 
//     {
//       theUserLimitsForDetector->SetUserMaxTime(val);
//       G4cout << " Changing Detector Time cut to: " << G4BestUnit(val,"Time")
// 	     << G4endl;
//     }
// }  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//void Test15DetectorConstruction::UpdateGeometry()
//{
//  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VPhysicalVolume *Test15DetectorConstruction::GetWorldVolume() {
   return world_phys;
}


G4VPhysicalVolume &Test15DetectorConstruction::GetWorldVolumeAddress() const{
  return *world_phys;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Test15DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
