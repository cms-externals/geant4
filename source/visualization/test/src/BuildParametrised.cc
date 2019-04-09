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
// Based on test code from Hans-Peter Wellisch

#include "BuildParametrised.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4PVParameterised.hh"
#include "ParametrisedBox.hh"

G4VPhysicalVolume* BuildParametrised ()
{
  G4Material* theAir=new G4Material("N2-air", 14.0, 28.0 * g / mole,
				     0.0001 *g/cm3);
  G4Material * theCopper =
    new G4Material("Quaterhard-Cu",25.0, 63.0 * g / mole,
		   6.8*g/cm3);

  G4RotationMatrix* theNull = new G4RotationMatrix ();
  G4ThreeVector* theCenter = new G4ThreeVector(0,0,0);

  // The Pit.
  G4Box* aBox = new G4Box("aBox",10*m, 10*m, 10*m);
  G4LogicalVolume*  theLogicalPit =
    new G4LogicalVolume(aBox, theAir,
			"thePit", 0, 0, 0);
  theLogicalPit->SetVisAttributes (G4VisAttributes::GetInvisible());
  G4PVPlacement* thePit =
    new G4PVPlacement(    theNull,
			  *theCenter,
			  "thePit",
			  theLogicalPit,    
			  0, 0, 0);

  // The Mother of Parametrised.
  G4VSolid * anotherBox = new G4Box("aBox",5*m, 5*m, 5*m);
  G4LogicalVolume*  theLogicalMother1 =
    new G4LogicalVolume(anotherBox, theAir,
			"aTest1", 0, 0,
			0);
  G4PVPlacement* theMother1 =
    new G4PVPlacement(    0, G4ThreeVector (-4 * m, 0 , 0),
			  "theMother1",
			  theLogicalMother1,    
			  thePit, 0, 0);

  // Parametrized volume.
  G4Box * theBox = new G4Box("aBox",1,1,1);
  G4LogicalVolume*  theLogicalBox1 =
    new G4LogicalVolume(theBox, theCopper,
			"aTest1", 0, 0,
			0);
  G4VisAttributes * logicalPAttributes1;
  logicalPAttributes1 = new G4VisAttributes(G4Colour(1.,0.,0.));
  theLogicalBox1->SetVisAttributes(logicalPAttributes1);
  ParametrisedBox* thePar1 = new ParametrisedBox(theCopper,theAir);
  new G4PVParameterised("theParametrizedtest1",
			theLogicalBox1,
			theMother1,
			kYAxis,
			2,
			thePar1);
  //  A second parameterised volume for test purposes...
  G4LogicalVolume*  theLogicalMother2 =
    new G4LogicalVolume(anotherBox, theAir,
			"aTest2", 0, 0,
			0);
  G4PVPlacement* theMother2 =
    new G4PVPlacement(    0, G4ThreeVector (0, -4 * m, 0),
			  "theMother2",
			  theLogicalMother2,    
			  thePit, 0, 0);
  G4LogicalVolume*  theLogicalBox2 =
    new G4LogicalVolume(theBox, theCopper,
			"aTest2", 0, 0,
			0);
  G4VisAttributes * logicalPAttributes2;
  logicalPAttributes2 = new G4VisAttributes(G4Colour(0.,1.,0.));
  theLogicalBox2->SetVisAttributes(logicalPAttributes2);
  ParametrisedBox* thePar2 = new ParametrisedBox(theAir,theCopper);
  new G4PVParameterised("theParametrizedtest2",
			theLogicalBox2,
			theMother2,
			kYAxis,
			5,
			thePar2);

  // The Mother of Replica.
  G4VSolid * box3 = new G4Box("aBox", 1 * m, 2 * m, 3 * m);
  G4LogicalVolume*  theLogicalRMother =
    new G4LogicalVolume(box3, theAir,
			"aTest", 0, 0,
			0);
  G4PVPlacement* theRMother =
    new G4PVPlacement(    0, G4ThreeVector (4 * m, 0, 0),
			  "theRMother",
			  theLogicalRMother,    
			  thePit, 0, 0);

  // Replicated volume.
  G4VSolid * RBox = new G4Box("aBox", 0.2 * m, 2 * m, 3 * m);
  G4LogicalVolume*  theLogicalRBox =
    new G4LogicalVolume(RBox, theCopper,
			"aTest", 0, 0,
			0);
  G4VisAttributes * logicalRAttributes;
  logicalRAttributes = new G4VisAttributes(G4Colour(0.,1.,0.));
  theLogicalRBox->SetVisAttributes(logicalRAttributes);
  new G4PVReplica("theReplicatest",
		  theLogicalRBox,
		  theRMother,
		  kXAxis,
		  5, 0.4 * m, 0);

  // The Mother of Tubs Replica.
  G4VSolid * tubs = new G4Tubs ("aTubsBox", 0.5 * m, 1 * m, 2 * m,
				90. * deg, 180. * deg);
  G4LogicalVolume*  theLogicalTMother =
    new G4LogicalVolume(tubs, theAir,
			"aTest", 0, 0,
			0);
  G4PVPlacement* theTMother =
    new G4PVPlacement(    0, G4ThreeVector (4 * m, 4 * m, 0),
			  "theTMother",
			  theLogicalTMother,    
			  thePit, 0, 0);

  // Replicated volume.
  G4VSolid * tubsPart = new G4Tubs ("aTubsPart", 0.5 * m, 1 * m, 2 * m,
				    0., 36. * deg);
  G4LogicalVolume*  theLogicalTubsPart =
    new G4LogicalVolume(tubsPart, theCopper,
			"aTest", 0, 0,
			0);
  G4VisAttributes * logicalTAttributes;
  logicalTAttributes = new G4VisAttributes(G4Colour(0.,0.,1.));
  theLogicalTubsPart->SetVisAttributes(logicalTAttributes);
  new G4PVReplica("theTReplicatest",
		  theLogicalTubsPart,
		  theTMother,
		  kPhi,
		  5, 36. * deg, 90. * deg);

  /*************** Bloggs's box to workaround bug in G4SmartVoxelHeader.

  G4Box * thebloggsBox = new G4Box("aBox",5*m, 4*m, 2*m);
  G4LogicalVolume*  theLogicalbloggs =
    new G4LogicalVolume(thebloggsBox, theCopper,
			"aTest", 0, 0,
			0);
  G4VisAttributes * bloggsAttributes = 
    new G4VisAttributes(G4Colour(0.,1.,0.));
  theLogicalbloggs->SetVisAttributes(bloggsAttributes);
  G4PVPlacement* bloggs =
    new G4PVPlacement(    theNull,
			  G4ThreeVector (10*m, -30*m, -40*m),
			  "bloggs",
			  theLogicalbloggs,
			  thePitPosition,    
			  0, 0);
			  ***************/

  return thePit;
}
