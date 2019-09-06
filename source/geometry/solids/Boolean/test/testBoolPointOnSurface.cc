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
// Test for G4BooleanSolid::GetPointOnSurface()
//
// History:
// 2016.05.14 E.Tcherniaev - created
//
//
// --------------------------------------------------------------------

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "ApproxEqual.hh"

#include "G4Box.hh"
#include "G4ScaledSolid.hh"
#include "G4ReflectedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4Timer.hh"
#include "Randomize.hh"
#include "G4Transform3D.hh"


/////////////////////////////////////////////////////////////
//
//  Create boolean solid

G4VSolid* create()
{
  // Create horizontal parallepiped
  G4VSolid* pA    = new G4Box("pA", 6, 2, std::sqrt(3.));
  G4VSolid* pboxA = new G4DisplacedSolid("Box_A", pA, G4TranslateX3D(1));

  // Create vertical parallepiped
  G4VSolid* pB    = new G4Box("pB", 5, 1, 2);
  G4VSolid* pboxB = new G4DisplacedSolid("Box_B", pB, G4RotateY3D(90*deg));

  // Create main box
  G4Box*            tmp0  = new G4Box("tmp0", 4, 3, 8);
  G4ScaledSolid*    tmp1  = new G4ScaledSolid(    "tmp1",  tmp0, G4Scale3D(0.5));
  G4ReflectedSolid* tmp2  = new G4ReflectedSolid( "tmp2",  tmp1, G4ReflectX3D(-1));
  G4VSolid*         pboxC = new G4ScaledSolid(    "Box C", tmp2, G4Scale3D(2));

  // Intersect parallepipeds (incline vertical parallepiped) 
  G4BooleanSolid* pAB  = new G4IntersectionSolid("A x B", pboxA, pboxB, G4RotateY3D(30*deg));

  // Substruct the result of intersection from main box
  G4BooleanSolid* pABC = new G4SubtractionSolid("C - (A x B)", pboxC, pAB, G4TranslateX3D(-1));

  return pABC;
}

/////////////////////////////////////////////////////////////
//
//  Collect statistics

void statistics(const G4ThreeVector& point, G4int action)
{
  const G4double s3  =  std::sqrt(3.);
  const G4double ss3 = 2*s3;

  static G4int i[12];
  G4double x = point.x(), y = point.y(), z = point.z();

  // Inilialization, set all counts to zero
  //
  if (action == 0) { 
    for (G4int k=0; k<12; k++) i[k] = 0; 
  }
 
  // Next point, find face where it is placed and increment corresponding count
  //
  if (action == 1) {
    // Check faces of the main box
         if ((-8 < x && x < 0) && (-3 < y && y < 3) && ApproxEqual(z,-8)) i[1]++;
    else if ((-8 < x && x < 0) && (-3 < y && y < 3) && ApproxEqual(z,+8)) i[2]++;
    else if ((-8 < x && x < 0) && ApproxEqual(y,-3) && (-8 < z && z < 8)) i[3]++;
    else if ((-8 < x && x < 0) && ApproxEqual(y,+3) && (-8 < z && z < 8)) i[4]++;
    else if (ApproxEqual(x,-8) && (-3 < y && y < 3) && (-8 < z && z < 8)) i[5]++;
    else if (ApproxEqual(x, 0) && (-3 < y && y < 3) && (-8 < z && z < 8)) i[6]++;
    // Check faces of the extraction
    else if ((-4/s3-2 < x && x < 0) && ApproxEqual(y,-1) && (-2 < z && z < 2))  i[7]++;
    else if ((-4/s3-2 < x && x < 0) && ApproxEqual(y,+1) && (-2 < z && z < 2))  i[8]++;
    else if ((-4/s3   < x && x < 0) && (-1 < y && y < 1) && ApproxEqual(z, s3)) i[ 9]++;
    else if ((-4/s3-2 < x && x < 0) && (-1 < y && y < 1) && ApproxEqual(z,-s3)) i[10]++;
    else if ((-4/s3-2 < x && x < -4/s3) && (-1 < y && y < 1)  && (-s3 < z && z < s3) &&
	     ApproxEqual((x+4/s3)/(z-s3),std::tan(30*deg))) i[11]++;
    else { // Unrecognized point
      i[0]++;
    }
  } 

  // Print statistics
  //
  if (action > 1) {
    G4double a[12], p[12], q[12];
    a[ 1] = 8*6;          // face_01 (-8<x< 0) (-3<y<3) (    z=-8)  area= 48 (8x6)
    a[ 2] = a[1];         // face_02 (-8<x< 0) (-3<y<3) (    z= 8)  area= 48 (8x6)
    a[ 3] = 8*16;         // face_03 (-8<x< 0) (  y=-3) (-8 <z< 8)  area=128 (8x16)
    a[ 4] = a[3];         // face_04 (-8<x< 0) (  y= 3) (-8 <z< 8)  area=128 (8x16)
    a[ 5] = 6*16;         // face_05 (   x=-8) (-3<y<3) (-8 <z< 8)  area= 96 (6x16)
    a[ 6] = a[5]-ss3*2;   // face_06 (   x= 0) (-3<y<3) (-8 <z< 8)  area= 96 (6x20)-(4xsqrt(3))
    a[ 7] = (4/s3+1)*ss3; // face_07           (  y=-1) (-2 <z< 2)
    a[ 8] = a[7];         // face_08           (  y= 1) (-2 <z< 2)
    a[ 9] = (4/s3  )*2;   // face_10 (-(4:sqrt(3))  <x<0) (-1<y<1) (z= 2)
    a[10] = (4/s3+2)*2;   // face_11 (-(4:sqrt(3)+2)<x<0) (-1<y<1) (z=-2)
    a[11] = 4*2;
    a[ 0] = a[1]+a[2]+a[3]+a[4]+a[5]+a[6]+a[7]+a[8]+a[9]+a[10]+a[11];
    
    for (G4int k=1; k<12; k++) { 
      p[k] = int(a[k]*10000./a[0])  /100.;
      q[k] = int(i[k]* 100./action);
    }

    G4cout << " Face  1: expected  " << p[ 1] << "%, real  " << q[ 1]<< "% (" << i[ 1] << ")\n"
           << " Face  2: expected  " << p[ 2] << "%, real  " << q[ 2]<< "% (" << i[ 2] << ")\n"
           << " Face  3: expected "  << p[ 3] << "%, real "  << q[ 3]<< "% (" << i[ 3] << ")\n"
           << " Face  4: expected "  << p[ 4] << "%, real "  << q[ 4]<< "% (" << i[ 4] << ")\n"
           << " Face  5: expected "  << p[ 5] << "%, real "  << q[ 5]<< "% (" << i[ 5] << ")\n"
           << " Face  6: expected "  << p[ 6] << "%, real "  << q[ 6]<< "% (" << i[ 6] << ")\n"
           << " Face  7: expected  " << p[ 7] << "%, real  " << q[ 7]<< "% (" << i[ 7] << ")\n"
           << " Face  8: expected  " << p[ 8] << "%, real  " << q[ 8]<< "% (" << i[ 8] << ")\n"
           << " Face  9: expected  " << p[ 9] << "%, real  " << q[ 9]<< "% (" << i[ 9] << ")\n"
           << " Face 10: expected  " << p[10] << "%, real  " << q[10]<< "% (" << i[10] << ")\n"
           << " Face 11: expected  " << p[11] << "%, real  " << q[11]<< "% (" << i[11] << ")\n"
           << "________\n"
           << " ????   : expected  0.00%, real  " << i[0]*100/action << "% (" << i[ 0] << ")" << G4endl;
    if (i[0] == 0) {
      G4cout << "\nTest is OK" << G4endl;
    }else{
      G4cout << "\n!!! Attention: there are " << i[0] << " unrecognized points !!!" << G4endl;
    }
  }
}

/////////////////////////////////////////////////////////////
//
//  Main

int main() 
{ 
  G4Timer time;
  G4VSolid* solid = create();

  G4int np = 1000000;
  G4ThreeVector p;
  statistics(p, 0);  

  G4cout << "*\n"
         << "* Test for G4BooleanSolid::GetPointOnSurface()\n"
         << "*\n"
         << "* Generate " << np << " random points on the surface of a boolean solid\n"
         << "* constructed by subtraction of the intersection of two parallepipeds\n"
         << "* from a box. The resulting polyhedron has 11 faces.\n"
         << "*\n"
         << "Distribution of generated points by faces:" << G4endl; 

  time.Start();
  for (G4int i=0; i<np; i++) {
    p = solid->GetPointOnSurface();
    statistics(p, 1);
  }  
  time.Stop();

  // Print statistics
  statistics(p, np);
  G4cout << "Time taken: " << time.GetRealElapsed() << " seconds\n " << G4endl;
  return 0;
}
