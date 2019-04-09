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
// class ExDivTesterPolyhedra Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterPolyhedra.hh"
#include "G4Polyhedra.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterPolyhedra::
ExDivTesterPolyhedra( PVType& pvtype, PlaceType& postype,
                      std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, postype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kRho );
  theAxis.push_back( kPhi );
  theAxis.push_back( kZAxis );
}

//--------------------------------------------------------------------------
void ExDivTesterPolyhedra::GenerateScanPoints()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerDiv = 2;
  numberOfPoints = theNDiv * nPointsPerDiv;
  // For division along X
  G4ThreeVector centre(0.,0.,-2*theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any Z, any Y
    G4ThreeVector pR( 0., theWorldLengthXY/100., theWorldLengthXY/100. );
    G4double X = -theWorldLengthXY + (ii+0.001) * 2*theWorldLengthXY/numberOfPoints;
    pR.setX( X );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Y
  centre = G4ThreeVector(0.,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any X, any Z
    G4ThreeVector pR( theWorldLengthXY/100., 0., theWorldLengthXY/100. );
    G4double Y = -theWorldLengthXY + (ii+0.001) * 2*theWorldLengthXY/numberOfPoints;
    pR.setY( Y );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Z
  centre = G4ThreeVector(0.,0.,2*theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any X, any Y
    G4ThreeVector pR( theWorldLengthXY/100., 0., theWorldLengthXY/100. );
    G4double Z = -theWorldLengthXY + (ii+0.001) * 2*theWorldLengthXY/numberOfPoints;
    pR.setZ( Z );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}

//--------------------------------------------------------------------------
void ExDivTesterPolyhedra::BuildParentSolids()
{
  G4int numSides = 3;
  G4int numZPlanes = 4;
  G4double* zPlane1 = new G4double[numZPlanes];
            zPlane1[0]=-theWorldLengthXY;
	    zPlane1[1]=-0.25*theWorldLengthXY;
	    zPlane1[2]= 0.5*theWorldLengthXY;
	    zPlane1[3]= theWorldLengthXY;
  G4double* rInner1 = new G4double[numZPlanes];
            rInner1[0]=0./2.;
	    rInner1[1]=0.1*theWorldLengthXY/2.;
	    rInner1[2]=0.3*theWorldLengthXY/2.;
	    rInner1[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter1  = new G4double[numZPlanes];
            rOuter1[0]=0.2*theWorldLengthXY/2.;
	    rOuter1[1]=0.4*theWorldLengthXY/2.;
	    rOuter1[2]=0.6*theWorldLengthXY/2.;
	    rOuter1[3]=0.9*theWorldLengthXY/2.;
  G4double* zPlane2 = new G4double[numZPlanes];
            zPlane2[0]=-theWorldLengthXY;
	    zPlane2[1]=-0.25*theWorldLengthXY;
	    zPlane2[2]= 0.5*theWorldLengthXY;
	    zPlane2[3]= theWorldLengthXY;
  G4double* rInner2 = new G4double[numZPlanes];
            rInner2[0]=0./2.;
	    rInner2[1]=0.1*theWorldLengthXY/2.;
	    rInner2[2]=0.3*theWorldLengthXY/2.;
	    rInner2[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter2  = new G4double[numZPlanes];
            rOuter2[0]=0.2*theWorldLengthXY/2.;
	    rOuter2[1]=0.4*theWorldLengthXY/2.;
	    rOuter2[2]=0.6*theWorldLengthXY/2.;
	    rOuter2[3]=0.9*theWorldLengthXY/2.;
  G4double* zPlane3 = new G4double[numZPlanes];
            zPlane3[0]=-theWorldLengthXY;
	    zPlane3[1]=-0.25*theWorldLengthXY;
	    zPlane3[2]= 0.5*theWorldLengthXY;
	    zPlane3[3]= theWorldLengthXY;
  G4double* rInner3 = new G4double[numZPlanes];
            rInner3[0]=0./2.;
	    rInner3[1]=0.1*theWorldLengthXY/2.;
	    rInner3[2]=0.2*theWorldLengthXY/2.;
	    rInner3[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter3  = new G4double[numZPlanes];
            rOuter3[0]=0.2*theWorldLengthXY/2.;
	    rOuter3[1]=0.4*theWorldLengthXY/2.;
	    rOuter3[2]=0.6*theWorldLengthXY/2.;
	    rOuter3[3]=0.9*theWorldLengthXY/2.;
  theParentSolids.push_back( new G4Polyhedra("parent_1", theStartPhi, theDeltaPhi,
                             numSides, numZPlanes, zPlane1, rInner1, rOuter1 ) );
  theParentSolids.push_back( new G4Polyhedra("parent_2", theStartPhi, theDeltaPhi,
                             numSides, numZPlanes, zPlane2, rInner2, rOuter2 ) );
  theParentSolids.push_back( new G4Polyhedra("parent_3", theStartPhi, theDeltaPhi,
			     numSides, numZPlanes, zPlane3, rInner3, rOuter3 ) );
}

//--------------------------------------------------------------------------
void ExDivTesterPolyhedra::BuildChildrenSolids()
{
  G4int numSides = 3;
  G4int numZPlanes = 4;
  G4double* zPlane1 = new G4double[numZPlanes];
            zPlane1[0]=-theWorldLengthXY;
	    zPlane1[1]=-0.25*theWorldLengthXY;
	    zPlane1[2]= 0.5*theWorldLengthXY;
	    zPlane1[3]= theWorldLengthXY;
  G4double* rInner1 = new G4double[numZPlanes];
            rInner1[0]=0./2.;
	    rInner1[1]=0.1*theWorldLengthXY/2.;
	    rInner1[2]=0.2*theWorldLengthXY/2.;
	    rInner1[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter1  = new G4double[numZPlanes];
            rOuter1[0]=0.2*theWorldLengthXY/2.;
	    rOuter1[1]=0.4*theWorldLengthXY/2.;
	    rOuter1[2]=0.6*theWorldLengthXY/2.;
	    rOuter1[3]=0.9*theWorldLengthXY/2.;
  G4double* zPlane2 = new G4double[numZPlanes];
            zPlane2[0]=-theWorldLengthXY;
	    zPlane2[1]=-0.25*theWorldLengthXY;
	    zPlane2[2]= 0.5*theWorldLengthXY;
	    zPlane2[3]= theWorldLengthXY;
  G4double* rInner2 = new G4double[numZPlanes];
            rInner2[0]=0./2.;
	    rInner2[1]=0.1*theWorldLengthXY/2.;
	    rInner2[2]=0.2*theWorldLengthXY/2.;
	    rInner2[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter2  = new G4double[numZPlanes];
            rOuter2[0]=0.2*theWorldLengthXY/2.;
	    rOuter2[1]=0.4*theWorldLengthXY/2.;
	    rOuter2[2]=0.6*theWorldLengthXY/2.;
	    rOuter2[3]=0.9*theWorldLengthXY/2.;
  G4double* zPlane3 = new G4double[numZPlanes];
            zPlane3[0]=-theWorldLengthXY;
	    zPlane3[1]=-0.25*theWorldLengthXY;
	    zPlane3[2]= 0.5*theWorldLengthXY;
	    zPlane3[3]= theWorldLengthXY;
  G4double* rInner3 = new G4double[numZPlanes];
            rInner3[0]=0./2.;
	    rInner3[1]=0.1*theWorldLengthXY/2.;
	    rInner3[2]=0.2*theWorldLengthXY/2.;
	    rInner3[3]=0.4*theWorldLengthXY/2.;
  G4double* rOuter3  = new G4double[numZPlanes];
            rOuter3[0]=0.2*theWorldLengthXY/2.;
	    rOuter3[1]=0.4*theWorldLengthXY/2.;
	    rOuter3[2]=0.6*theWorldLengthXY/2.;
	    rOuter3[3]=0.9*theWorldLengthXY/2.;

  G4Polyhedra* msol = (G4Polyhedra*)theParentSolids[0];
  G4PolyhedraHistorical* origparamMother = msol->GetOriginalParameters();
  G4double rMax = origparamMother->Rmax[0] - origparamMother->Rmin[0];
  msol = (G4Polyhedra*)theParentSolids[1];
  G4double phiMax =  msol->GetEndPhi() - msol->GetStartPhi();
  msol = (G4Polyhedra*)theParentSolids[2];
  origparamMother = msol->GetOriginalParameters();
  G4double zMax = origparamMother->Z_values[origparamMother->Num_z_planes-1] - origparamMother->Z_values[0];

  theWidths.push_back( rMax / theNDiv );
  theWidths.push_back( phiMax / theNDiv );
  theWidths.push_back( zMax / theNDiv );

  theChildSolids.push_back( new G4Polyhedra("child_1", theStartPhi, theDeltaPhi,
                            numSides, numZPlanes, zPlane1, rInner1, rOuter1 ) );
  theChildSolids.push_back( new G4Polyhedra("child_2", theStartPhi, theWidths[0],
                            numSides, numZPlanes, zPlane2, rInner2, rOuter2 ) );
  theChildSolids.push_back( new G4Polyhedra("child_3", theStartPhi, theDeltaPhi,
                            numSides, numZPlanes, zPlane3, rInner3, rOuter3 ) );
}

