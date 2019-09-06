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

#include "ExampleDetectorConstruction.hh"

ExampleDetectorConstruction::ExampleDetectorConstruction()
{
  expHall_x = 1111600.*cm;
  expHall_y = 1111600.*cm;
  expHall_z = 1111600.*cm;

  calBox_x = 50.*cm;
  calBox_y = 50.*cm;
  calBox_z = 50.*cm;
  rotAngle = 30.*deg;
  calPos = 200.*cm;

  trackerRadius = 50.*cm;
  trackerHight = 100.*cm;
  trackerPos = -200.*cm;
}

ExampleDetectorConstruction::~ExampleDetectorConstruction()
{;}

G4VPhysicalVolume* ExampleDetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall

  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);

  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);

  G4ThreeVector coords(0.0,0.0,0.0);

  G4Box *vol1_solid = new G4Box("v1_sol",10.*m,10.*m,10.*m);

  G4LogicalVolume *vol1_log = new G4LogicalVolume(vol1_solid,Lead,"v1_log",0,0,0);

  G4VPhysicalVolume *vol1_phys = new G4PVPlacement(0,coords,"v1_sol",vol1_log,
						   experimentalHall_phys,false,0);


  G4Box *vol2_solid = new G4Box("v2_sol",10.*mm,10.*mm,10.*mm);

  G4LogicalVolume *vol2_log = new G4LogicalVolume(vol2_solid,Lead,"v2_log",0,0,0);

  G4VPhysicalVolume *vol2_phys = new G4PVPlacement(0,coords,"v2_sol",vol2_log,
						   vol1_phys,false,0);

  
  return experimentalHall_phys;
}
















