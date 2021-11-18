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
// Makoto Asai - based on Long Baseline Neutrino Observatory experiment.
// Embryo detector with rotated volumes

#include "BuildGeom_Example2.hh"

#include <cmath>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"
#include "G4GeometryManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"


G4PVPlacement* calUnitPhys; 
G4PVPlacement* calRowPhys;
G4PVPlacement* calCellPhys;

G4VPhysicalVolume* BuildGeom_Example2()
{

//--------- Material definition ---------

  G4double a, iz, density;
  G4String name, symbol;
  G4int nel;

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", iz=14., a);
  a = 207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a);

  //a = 26.98*g/mole;
  //density = 2.7*g/cm3;
  //G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  //a = 55.85*g/mole;
  //density = 7.87*g/cm3;
  //G4Material* Fe = new G4Material(name="Iron", z=26., a, density);
  //a = 207.19*g/mole;
  //density = 11.35*g/cm3;
  //G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  density = 5.2*g/cm3;
  G4Material* LeadGlass = new G4Material(name="LeadGlass", density, nel=3);
  LeadGlass->AddElement(elO, .199);
  LeadGlass->AddElement(elSi, .127);
  LeadGlass->AddElement(elPb, .674);
  density = 1.0*g/cm3;
  G4Material* Water = new G4Material(name="water",density,nel=2);
  Water->AddElement(elH, (G4int)2);
  Water->AddElement(elO, (G4int)1);


// Experimental hall (world volume)

  G4Box *myWorldBox= new G4Box("WBox",750*cm,1000*cm,1200*cm);
  G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,Air,
						    "WLog", 0, 0, 0);
  G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);

  G4double waterCherenkovZpos = 100.*cm;


// LeadGlass unit

  G4Tubs* calUnitTubs 
    = new G4Tubs("CUTubs",170*cm,240*cm,30.5*cm,-22.5*deg,45.0*deg);
  G4LogicalVolume* calUnitLog
    = new G4LogicalVolume(calUnitTubs,Air,"CUlog",0,0,0);

  G4RotationMatrix* calUnitRotR = new G4RotationMatrix();
  calUnitRotR->rotateX(-90.0*deg);
  calUnitRotR->rotateY(-112.5*deg);
  G4RotationMatrix* calUnitRotL = new G4RotationMatrix();
  calUnitRotL->rotateX(-90.0*deg);
  calUnitRotL->rotateY(-67.5*deg);

  G4int calUnitCopyNo = 0;
  //G4PVPlacement *calUnitPhys; 
  //for(int iCalUnit=2; iCalUnit>=-2; iCalUnit--)
  for(int iCalUnit=2; iCalUnit>=2; iCalUnit--)
  {
    G4double yCalRow = iCalUnit*61.0*cm;
    //calUnitPhys = new G4PVPlacement(calUnitRotL,
    //           G4ThreeVector(0.*cm,yCalRow,waterCherenkovZpos),
    //           "calUnitPhys",calUnitLog,myWorldPhys,false,calUnitCopyNo++);
    calUnitPhys = new G4PVPlacement(calUnitRotR,
               G4ThreeVector(0.*cm,yCalRow,waterCherenkovZpos),
               "calUnitPhys",calUnitLog,myWorldPhys,false,calUnitCopyNo++);
  }


// LeadGlass row

  G4Tubs* calRowTubs 
    = new G4Tubs("CRTubs",170*cm,212.9*cm,6.1*cm,-22.5*deg,45.0*deg);
  G4LogicalVolume* calRowLog
    = new G4LogicalVolume(calRowTubs,Air,"CRlog",0,0,0);

  G4int calRowCopyNo = 0;
  //G4PVPlacement *calRowPhys;
  //for(int iCalRow=-2; iCalRow<=2; iCalRow++)
  for(int iCalRow=-2; iCalRow<=-2; iCalRow++)
  {
    G4double zCalRow = iCalRow*12.2*cm;
    calRowPhys = new G4PVPlacement(0,
               G4ThreeVector(0.*cm,0.*cm,zCalRow),
               "calRowPhys",calRowLog,calUnitPhys,false,calRowCopyNo++);
  }


// LeadGlass cell

  G4Trd* calCellTrd  
    = new G4Trd("CCtrd",5.65*cm,6.75*cm,6.1*cm,6.1*cm,17.0*cm);
  G4LogicalVolume* calCellLog
    = new G4LogicalVolume(calCellTrd,LeadGlass,"CCtrd",0,0,0);

  G4double calCellOpeningAngle = 3.702*deg;
  G4double calCellZpos = 191.636*cm;
  //G4PVPlacement* calCellPhys;
  for(int iCalCell=0; iCalCell<=11; iCalCell++)
  {
    G4RotationMatrix* calCellRot = new G4RotationMatrix();
    calCellRot->rotateX(90.0*deg);
    G4double cellRotAngle = (5.5-iCalCell)*calCellOpeningAngle;
    calCellRot->rotateZ(cellRotAngle+90.0*deg);
/*
    G4cout << calCellRot << G4endl;
    G4cout << calCellRot->xx() << " "
         << calCellRot->xy() << " "
         << calCellRot->xz() << " " << G4endl;
    G4cout << calCellRot->yx() << " "
         << calCellRot->yy() << " "
         << calCellRot->yz() << " " << G4endl;
    G4cout << calCellRot->zx() << " "
         << calCellRot->zy() << " "
         << calCellRot->zz() << " " << G4endl;
*/
    G4double calCellX = calCellZpos * std::cos(cellRotAngle/rad);
    G4double calCellY = calCellZpos * std::sin(cellRotAngle/rad);
    calCellPhys = new G4PVPlacement(calCellRot,
	       G4ThreeVector(calCellX,calCellY,0.*cm),
	       "calRowPhys",calCellLog,calRowPhys,false,iCalCell);
/*
    G4RotationMatrix* tmpRot = calCellPhys->GetRotation();
    G4cout << "saved : " << tmpRot << G4endl;
    G4cout << tmpRot->xx() << " "
         << tmpRot->xy() << " "
         << tmpRot->xz() << " " << G4endl;
    G4cout << tmpRot->yx() << " "
         << tmpRot->yy() << " "
         << tmpRot->yz() << " " << G4endl;
    G4cout << tmpRot->zx() << " "
         << tmpRot->zy() << " "
         << tmpRot->zz() << " " << G4endl;
*/
  }

  return myWorldPhys;
}

