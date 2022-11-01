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
// History:
//
// 2016.10.12 E.Tcherniaev - initial version
// 
// --------------------------------------------------------------------
#include <assert.h>

#include "globals.hh"
#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"

#include "G4GeomTools.hh"

using namespace CLHEP;

void CheckTriangelArea();
void CheckQuadArea();
void CheckPolygonArea();
void CheckPointInTriangle();
void CheckPointInPolygon();
void CheckIsConvex();
void CheckTriangulatePolygon();
void CheckRemoveRedundantVertices();
void CheckDiskExtent();
void CheckEllipsePerimeter();
void CheckEllipticConeLateralArea();

void CheckTriangelAreaNormal();
void CheckQuadAreaNormal();
void CheckPolygonAreaNormal();
void CheckClosestPointOnSegment();
void CheckDistancePointSegment();
void CheckClosestPointOnTriangle();
void CheckSphereExtent();

int main()
{

  //===============================================
  //  2D Utilities
  //-----------------------------------------------

  G4cout << "\n====== Test 2D utilities ======\n" << G4endl;

  G4cout << "G4GeomTools::TriangelArea()" << G4endl;
  CheckTriangelArea();

  G4cout << "G4GeomTools::QuadArea()" << G4endl;
  CheckQuadArea();

  G4cout << "G4GeomTools::PolygonArea()" << G4endl;
  CheckPolygonArea();

  G4cout << "G4GeomTools::PointInTriangle()" << G4endl;
  CheckPointInTriangle();

  G4cout << "G4GeomTools::PointInPolygon()" << G4endl;
  CheckPointInPolygon();

  G4cout << "G4GeomTools::IsConvex()" << G4endl;
  CheckIsConvex();

  G4cout << "G4GeomTools::TriangulatePolygon()" << G4endl;
  CheckTriangulatePolygon();

  G4cout << "G4GeomTools::RemoveRedundantVertices()" << G4endl;
  CheckRemoveRedundantVertices();

  G4cout << "G4GeomTools::DiskExtent()" << G4endl;
  CheckDiskExtent();

  G4cout << "G4GeomTools::EllipsePerimeter()" << G4endl;
  CheckEllipsePerimeter();

  G4cout << "G4GeomTools::EllipticConeLateralArea()" << G4endl;
  CheckEllipticConeLateralArea();

  //===============================================
  //  3D Utilities
  //-----------------------------------------------

  G4cout << "\n====== Test 3D utilities ======\n" << G4endl;

  G4cout << "G4GeomTools::TriangelAreaNormal()" << G4endl;
  CheckTriangelAreaNormal();

  G4cout << "G4GeomTools::QuadAreaNormal()" << G4endl;
  CheckQuadAreaNormal();

  G4cout << "G4GeomTools::PolygonAreaNormal()" << G4endl;
  CheckPolygonAreaNormal();

  G4cout << "G4GeomTools::ClosestPointOnSegment()" << G4endl;
  CheckClosestPointOnSegment();

  G4cout << "G4GeomTools::DistancePointSegment()" << G4endl;
  CheckDistancePointSegment();
  
  G4cout << "G4GeomTools::ClosestPointOnTriangle()" << G4endl;
  CheckClosestPointOnTriangle();

  G4cout << "G4GeomTools::SphereExtent()" << G4endl;
  CheckSphereExtent();

  G4cout << "\nTest G4GeomTools is OK!\n" << G4endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////
//
void CheckTriangelArea()
{
  G4TwoVector A(2,1), B(5,1), C(5,5);
  assert(G4GeomTools::TriangleArea(A.x(),A.y(),B.x(),B.y(),C.x(),C.y()) == 6.);
  assert(G4GeomTools::TriangleArea(C, A, B) == 6.);
}

///////////////////////////////////////////////////////////////////////
//
void CheckQuadArea()
{
  G4TwoVector A(2,1), B(5,1), C(5,5), D(4,2);
  assert(G4GeomTools::QuadArea(A, B, C, D) == 3.5);
}

///////////////////////////////////////////////////////////////////////
//
void CheckPolygonArea()
{
  G4TwoVector A(2,1), B(5,1), C(5,5), D(4,2);
  G4TwoVectorList polygon;

  polygon.push_back(D);
  polygon.push_back(C);
  polygon.push_back(B);
  polygon.push_back(A);
  assert(G4GeomTools::PolygonArea(polygon) == -3.5);
}

///////////////////////////////////////////////////////////////////////
//
void CheckPointInTriangle()
{
  G4TwoVector A(2,1), B(5,1), C(5,5), D(4,2);

  assert(G4GeomTools::PointInTriangle(A.x(),A.y(),B.x(),B.y(),C.x(),C.y(), 4.0,2.0));
  assert(G4GeomTools::PointInTriangle(C.x(),C.y(),B.x(),B.y(),A.x(),A.y(), 3.5,3.0));
  assert(G4GeomTools::PointInTriangle(A.x(),A.y(),B.x(),B.y(),C.x(),C.y(), 5.0,1.0));

  assert(!G4GeomTools::PointInTriangle(C.x(),C.y(),B.x(),B.y(),A.x(),A.y(), 5.0,0.5));
  assert(!G4GeomTools::PointInTriangle(A.x(),A.y(),B.x(),B.y(),C.x(),C.y(), 5.5,3.0));
  assert(!G4GeomTools::PointInTriangle(C.x(),C.y(),B.x(),B.y(),A.x(),A.y(), 2.0,2.0));

  assert(G4GeomTools::PointInTriangle(A,B,C,D));
  assert(!G4GeomTools::PointInTriangle(C,B,A,G4TwoVector(4,4)));
}

///////////////////////////////////////////////////////////////////////
//
void CheckPointInPolygon()
{
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-3,-1)); //   
  polygon.push_back(G4TwoVector(-3, 3)); //   |\       / \       /|
  polygon.push_back(G4TwoVector(-2, 1)); // . | \  .  / . \  .  / | .
  polygon.push_back(G4TwoVector(-1, 1)); // . | .\_._/  .  \_._/. | .
  polygon.push_back(G4TwoVector( 0, 3)); //   |                   |
  polygon.push_back(G4TwoVector( 1, 1)); // . |         .         | .
  polygon.push_back(G4TwoVector( 2, 1)); //   |___________________|
  polygon.push_back(G4TwoVector( 3, 3)); //
  polygon.push_back(G4TwoVector( 3,-1)); // .           .           .

  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-4.0, 2.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-1.5, 2.0),polygon));
  assert( G4GeomTools::PointInPolygon(G4TwoVector( 0.0, 2.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 1.5, 2.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-4.0, 2.0),polygon));

  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-4.0, 1.0),polygon));
  assert( G4GeomTools::PointInPolygon(G4TwoVector(-2.5, 1.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-1.5, 1.0),polygon));
  assert( G4GeomTools::PointInPolygon(G4TwoVector( 0.0, 1.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 1.5, 1.0),polygon));
  assert( G4GeomTools::PointInPolygon(G4TwoVector( 2.5, 1.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 4.0, 1.0),polygon));

  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-4.0, 0.0),polygon));
  assert( G4GeomTools::PointInPolygon(G4TwoVector( 0.0, 0.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 4.0, 0.0),polygon));

  assert(!G4GeomTools::PointInPolygon(G4TwoVector(-4.0,-2.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 0.0,-2.0),polygon));
  assert(!G4GeomTools::PointInPolygon(G4TwoVector( 4.0,-2.0),polygon));
}

///////////////////////////////////////////////////////////////////////
//
void CheckIsConvex()
{
  G4TwoVectorList pp;
  pp.push_back(G4TwoVector(0,0)); assert(!G4GeomTools::IsConvex(pp));
  pp.push_back(G4TwoVector(5,5)); assert(!G4GeomTools::IsConvex(pp));
  pp.push_back(G4TwoVector(2,0)); assert( G4GeomTools::IsConvex(pp));
  pp.push_back(G4TwoVector(2,0)); assert(!G4GeomTools::IsConvex(pp));
  pp[3].set(1,0);                 assert(!G4GeomTools::IsConvex(pp));
  pp[3].set(1,-1);                assert( G4GeomTools::IsConvex(pp));
  pp[3].set(5,-5);                assert(!G4GeomTools::IsConvex(pp));
}

///////////////////////////////////////////////////////////////////////
//
void CheckTriangulatePolygon()
{
  G4TwoVectorList a, result;

  a.push_back( G4TwoVector(0,0));
  a.push_back( G4TwoVector(0,5));
  a.push_back( G4TwoVector(1,5));

  a.push_back( G4TwoVector(1,1));
  a.push_back( G4TwoVector(2,1));
  a.push_back( G4TwoVector(2,5));
  a.push_back( G4TwoVector(3,5));

  a.push_back( G4TwoVector(3,1));
  a.push_back( G4TwoVector(4,1));
  a.push_back( G4TwoVector(4,5));
  a.push_back( G4TwoVector(5,5));

  a.push_back( G4TwoVector(5,1));
  a.push_back( G4TwoVector(6,1));
  a.push_back( G4TwoVector(6,5));
  a.push_back( G4TwoVector(7,5));

  a.push_back( G4TwoVector(7,0));
  
  G4GeomTools::TriangulatePolygon(a,result);

  // Triangle 1 => (0,0) (3,1) (4,1)
  assert(result[0] == G4TwoVector(0,0));
  assert(result[1] == G4TwoVector(3,1));
  assert(result[2] == G4TwoVector(4,1));

  // Triangle 2 => (4,1) (5,1) (0,0)
  // Triangle 3 => (0,0) (2,1) (3,1)
  // Triangle 4 => (5,1) (6,1) (0,0)
  // Triangle 5 => (0,0) (1,1) (2,1)
  // Triangle 6 => (6,1) (7,0) (0,0)
  // Triangle 7 => (0,0) (0,5) (1,1)
  // Triangle 8 => (2,1) (2,5) (3,1)
  assert(result[21] == G4TwoVector(2,1));
  assert(result[22] == G4TwoVector(2,5));
  assert(result[23] == G4TwoVector(3,1));

  // Triangle  9 => (4,1) (4,5) (5,1)
  // Triangle 10 => (6,1) (6,5) (7,0)
  // Triangle 11 => (0,5) (1,5) (1,1)
  // Triangle 12 => (2,5) (3,5) (3,1)
  // Triangle 13 => (4,5) (5,5) (5,1)
  // Triangle 14 => (6,5) (7,5) (7,0)
  assert(result[39] == G4TwoVector(6,5));
  assert(result[40] == G4TwoVector(7,5));
  assert(result[41] == G4TwoVector(7,0));
}

///////////////////////////////////////////////////////////////////////
//
void CheckRemoveRedundantVertices()
{
  std::vector<G4int> iout;
  std::vector<G4TwoVector> polygon;

  polygon.resize(0);
  polygon.push_back(G4TwoVector(0,0));
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 0);
  assert(iout.size() == 1);
  assert(iout[0] == 0);

  polygon.resize(0);
  polygon.push_back(G4TwoVector(0,0));
  polygon.push_back(G4TwoVector(1,1));
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 0);
  assert(iout.size() == 2);
  for (G4int i=0; i<2; ++i) assert(iout[i] == i);
  polygon.resize(0);
  polygon.push_back(G4TwoVector(0,0));
  polygon.push_back(G4TwoVector(1,1));
  polygon.push_back(G4TwoVector(0,0));
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 0);
  assert(iout.size() == 3);
  for (G4int i=0; i<2; ++i) assert(iout[i] == i);

  polygon.resize(0);
  polygon.push_back(G4TwoVector(0,-1));          // 3
  polygon.push_back(G4TwoVector(0,0));           //  +---+---+ 2
  polygon.push_back(G4TwoVector(1,0));           //  |   |1
  polygon.push_back(G4TwoVector(-1,0));          //  +---+
  polygon.push_back(G4TwoVector(-1,-1));         // 4     0
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 4);
  assert(iout.size() == 1);
  assert(iout[0] == 2);

  polygon.resize(0);
  polygon.push_back(G4TwoVector(0,-1));          // 4    3
  polygon.push_back(G4TwoVector(0,0));           //  +---++---+ 2
  polygon.push_back(G4TwoVector(1,0));           //  |   |1
  polygon.push_back(G4TwoVector(0,0));           //  +---+ 
  polygon.push_back(G4TwoVector(-1,0));          // 5     0
  polygon.push_back(G4TwoVector(-1,-1));         //
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 4);
  assert(iout.size() == 2);
  assert(iout[0] == 2);
  assert(iout[1] == 3);
  
  G4double tolerance = 1.E-10;
  polygon.resize(0);
  polygon.push_back(G4TwoVector(1,1));
  polygon.push_back(G4TwoVector(1000-tolerance,1000+tolerance));
  polygon.push_back(G4TwoVector(10000,10000));
  G4GeomTools::G4GeomTools::RemoveRedundantVertices(polygon,iout);
  assert(polygon.size() == 3);
  assert(iout.size() == 0);
}

///////////////////////////////////////////////////////////////////////
//
void CheckDiskExtent()
{
  G4TwoVector pmin, pmax;
  G4double rmin = 2, rmax = 3;

  G4GeomTools::DiskExtent(rmin,rmax, 150*deg,370*deg, pmin,pmax);
  assert(pmin == G4TwoVector(-rmax,-rmax));
  assert(pmax == G4TwoVector( rmax, rmax));

  for (G4int istart=-350; istart <= 350; istart += 10)
  {
    for (G4int idelta=10; idelta <= 360; idelta += 10)
    {
      G4double sphi = istart*deg;
      G4double dphi = idelta*deg;
      G4double ephi = sphi + dphi;
      G4double sinStart = std::sin(sphi), cosStart = std::cos(sphi);  
      G4double sinEnd   = std::sin(ephi), cosEnd   = std::cos(ephi);  

      G4GeomTools::DiskExtent(rmin,rmax, sphi,dphi, pmin,pmax);

      G4double xmin = std::min(std::min(std::min(rmin*cosStart,rmin*cosEnd),rmax*cosStart),rmax*cosEnd);
      G4double xmax = std::max(std::max(std::max(rmin*cosStart,rmin*cosEnd),rmax*cosStart),rmax*cosEnd);
      G4double ymin = std::min(std::min(std::min(rmin*sinStart,rmin*sinEnd),rmax*sinStart),rmax*sinEnd);
      G4double ymax = std::max(std::max(std::max(rmin*sinStart,rmin*sinEnd),rmax*sinStart),rmax*sinEnd);

      // check phi = 0
      if ((sphi < 0.0 && 0.0 < ephi) ||
          (sphi  < twopi &&  twopi < ephi))     { xmax =  rmax; }

      // check phi = 90
      if ((sphi < -1.5*pi && -1.5*pi < ephi) ||
          (sphi <  0.5*pi &&  0.5*pi < ephi) ||
          (sphi <  2.5*pi &&  2.5*pi < ephi))   { ymax =  rmax; }

      // check phi = 180
      if ((sphi <    -pi &&    -pi < ephi) ||
          (sphi <     pi &&     pi < ephi) ||
          (sphi < 3.0*pi && 3.0*pi < ephi))     { xmin = -rmax; }

      // check phi = 270
      if ((sphi < -0.5*pi && -0.5*pi < ephi) ||
          (sphi <  1.5*pi &&  1.5*pi < ephi) ||
          (sphi <  3.5*pi &&  3.5*pi < ephi))   { ymin = -rmax; }

      assert(pmin == G4TwoVector(xmin,ymin));
      assert(pmax == G4TwoVector(xmax,ymax));
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
void CheckEllipsePerimeter()
{
  const G4double eps = 1e-15; // relative precision
  assert(std::abs(G4GeomTools::EllipsePerimeter( 2, 2    ) - 4.*pi) < 4*eps);
  assert(std::abs(G4GeomTools::EllipsePerimeter( 2, 2+eps) - 4.*pi) < 4*eps);
  assert(std::abs(G4GeomTools::EllipsePerimeter( 2, 0    ) - 8.   ) < 8*eps);
  assert(std::abs(G4GeomTools::EllipsePerimeter( 2, 4*eps) - 8.   ) < 8*eps);

  const G4int ntab = 12;
  G4double tab[ntab][3] = {
    { 1, 8, 32.7449566001955006 },
    { 2, 8, 34.3136871006273338 },
    { 3, 8, 36.3668627831685498 },
    { 4, 8, 38.7537928821907016 },
    { 5, 8, 41.3862760722288243 },
    { 6, 8, 44.2069843214190072 },
    { 7, 8, 47.1762642420057432 },

    { 2,    3,   15.86543958929059   },
    { 4,    5,   28.361667888974484  },
    { 1,   10,   40.6397418010089595 },
    { 1,  100,  400.109832972265224  },
    { 1, 1000, 4000.01558810468759   },
  };

  const G4int npow = 9;
  G4double a, b, c;
  for (G4int i = 0; i < ntab; ++i) {
    a = tab[i][0];
    b = tab[i][1];
    c = tab[i][2];
    for (G4int n = 0; n < npow; ++n) {
      G4double k = std::pow(10,n); // 10^n
      assert(std::abs(G4GeomTools::EllipsePerimeter(a*k,b*k) - c*k) < eps*c*k);
      assert(std::abs(G4GeomTools::EllipsePerimeter(a/k,b/k) - c/k) < eps*c/k);
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
void CheckEllipticConeLateralArea()
{
  const G4double eps = 1e-15;
  assert(std::abs(G4GeomTools::EllipticConeLateralArea(2,2,3) - pi*2*std::hypot(2,3)) < eps);
  assert(std::abs(G4GeomTools::EllipticConeLateralArea(2,3,0) - pi*2*3              ) < eps);
  assert(std::abs(G4GeomTools::EllipticConeLateralArea(2,3,4) - 36.978409022454123  ) < eps);
}

///////////////////////////////////////////////////////////////////////
//
void CheckTriangelAreaNormal()
{
  G4ThreeVector A(2,1,1), B(5,1,1), C(5,5,1);
  assert(G4GeomTools::TriangleAreaNormal(A,B,C) == G4ThreeVector(0,0,6));
}

///////////////////////////////////////////////////////////////////////
//
void CheckQuadAreaNormal()
{
  G4ThreeVector A(2,1,1), B(5,1,1), C(5,5,1), D(4,2,1);
  assert(G4GeomTools::QuadAreaNormal(A,B,C,D) == G4ThreeVector(0,0,3.5));
}

///////////////////////////////////////////////////////////////////////
//
void CheckPolygonAreaNormal()
{
  G4ThreeVectorList polygon3d;
  G4ThreeVector A(2,1,1), B(5,1,1), C(5,5,1), D(4,2,1);
  polygon3d.push_back(D);
  polygon3d.push_back(C);
  polygon3d.push_back(B);
  polygon3d.push_back(A);
  assert(G4GeomTools::PolygonAreaNormal(polygon3d) == G4ThreeVector(0,0,-3.5));
}

///////////////////////////////////////////////////////////////////////
//
void CheckClosestPointOnSegment()
{
  G4ThreeVector A(3,2,1), B(1,3,4), P; 
  P = G4GeomTools::ClosestPointOnSegment(G4ThreeVector(0,0,-1),A,B);
  assert(P == A);
  P = G4GeomTools::ClosestPointOnSegment(G4ThreeVector(0,1,5),A,B);
  assert(P == B);
  P = G4GeomTools::ClosestPointOnSegment(G4ThreeVector(1,0.5,2.5),A,B);
  assert(P == 0.5*(A+B));
}

///////////////////////////////////////////////////////////////////////
//
void CheckDistancePointSegment()
{
  const G4double eps = 1e-15;
  G4ThreeVector A(0,1,0);
  G4ThreeVector B(1,1,0);
  G4double dist;

  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(0,0,0),A,B);
  assert (std::abs(dist - 1.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(-1,1,0),A,B);
  assert (std::abs(dist - 1.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(1,2,0),A,B);
  assert (std::abs(dist - 1.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(2,1,0),A,B);
  assert (std::abs(dist - 1.) < eps);

  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(0.5,0,0),A,B);
  assert (std::abs(dist - 1.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(0.8,3,0),A,B);
  assert (std::abs(dist - 2.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(0.2,1,1.5),A,B);
  assert (std::abs(dist - 1.5) < eps);

  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(0,1,0),A,B);
  assert (std::abs(dist - 0.) < eps);
  dist  = G4GeomTools::DistancePointSegment(G4ThreeVector(1,1,0),A,B);
  assert (std::abs(dist - 0.) < eps);
}

///////////////////////////////////////////////////////////////////////
//
void CheckClosestPointOnTriangle()
{
  G4ThreeVector A(0,0,0), B(5,0,2), C(-3,4,2), P; 

  // region 0
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-8,-16,20),A,B,C);
  assert(P == G4ThreeVector(0,0,0));
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-8+0.5,-16+1,20+1),A,B,C);
  assert(P == G4ThreeVector(0.5,1,1));
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-8+1,-16+2,20+2),A,B,C);
  assert(P == G4ThreeVector(1,2,2));

  // region 1
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-8+2,-16+4,20+3),A,B,C);
  assert(P == G4ThreeVector(1,2,2));
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(6,1,3),A,B,C);
  assert(P == B);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-3,5,3),A,B,C);
  assert(P == C);

  // region 2
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-4,5,2),A,B,C);
  assert(P == C);

  // region 3
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-4,4,2),A,B,C);
  assert(P == C);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-20,1,-8),A,B,C);
  assert(P == C);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-3.5,0.5,1),A,B,C);
  assert(P == G4ThreeVector(-1.5,2,1));

  // region 4
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-1,-1,0),A,B,C);
  assert(P == A);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(-20,-1,-8),A,B,C);
  assert(P == C);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(9,-12,-8),A,B,C);
  assert(P == B);

  // region 5
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(5,-1,2),A,B,C);
  assert(P == B);
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(2.5,-1,1),A,B,C);
  assert(P == G4ThreeVector(2.5,0,1));

  // region 6
  P = G4GeomTools::ClosestPointOnTriangle(G4ThreeVector(8,-1,2),A,B,C);
  assert(P == B);
}

///////////////////////////////////////////////////////////////////////
//
void CheckSphereExtent()
{
  const G4double eps = 1e-14;
  G4ThreeVector p3min, p3max;
  G4double rmin = 2, rmax = 3;
  for (G4int sthe = 0; sthe <= 170; sthe += 10)
  {
    for (G4int dthe = 10; dthe <= 180-sthe; dthe += 10)
    {
      for (G4int sphi = -350; sphi <= 350; sphi += 10)
      {
        for (G4int dphi = 10; dphi <= 360; dphi += 10)
        {
          G4double xmin = 999, ymin = 999, zmin = 999;
          G4double xmax =-999, ymax =-999, zmax =-999;
          for (G4int i = sthe; i <= sthe+dthe; i += 10)
          {
            for (G4int k = sphi; k <= sphi+dphi; k += 10)
            {
              G4double x = std::sin(i*deg)*std::cos(k*deg);
              G4double y = std::sin(i*deg)*std::sin(k*deg);
              G4double z = std::cos(i*deg);
              if (x*rmin < xmin) xmin = x*rmin;
              if (x*rmin > xmax) xmax = x*rmin;
              if (y*rmin < ymin) ymin = y*rmin;
              if (y*rmin > ymax) ymax = y*rmin;
              if (z*rmin < zmin) zmin = z*rmin;
              if (z*rmin > zmax) zmax = z*rmin;

              if (x*rmax < xmin) xmin = x*rmax;
              if (x*rmax > xmax) xmax = x*rmax;
              if (y*rmax < ymin) ymin = y*rmax;
              if (y*rmax > ymax) ymax = y*rmax;
              if (z*rmax < zmin) zmin = z*rmax;
              if (z*rmax > zmax) zmax = z*rmax;
            }
          }
	  G4GeomTools::SphereExtent(rmin,rmax,
                                    sthe*deg,dthe*deg,
                                    sphi*deg,dphi*deg,
                                    p3min,p3max);
          assert(std::abs(p3min.x()-xmin) < eps);
          assert(std::abs(p3min.y()-ymin) < eps);
          assert(std::abs(p3min.z()-zmin) < eps);
          assert(std::abs(p3max.x()-xmax) < eps);
          assert(std::abs(p3max.y()-ymax) < eps);
          assert(std::abs(p3max.z()-zmax) < eps);
        }
      }
    }
  }
}
