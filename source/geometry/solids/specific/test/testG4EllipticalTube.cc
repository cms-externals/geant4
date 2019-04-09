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
// Unit test file for G4EllipticalTube
//
// Author: Evgueni Tcherniaev (evgueni.tcherniaev@cern.ch)
//

// ensure asserts are compiled in
#undef NDEBUG

#include <assert.h>
#include <cmath>
#include <iomanip>

#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"
#include "globals.hh"

#include "G4EllipticalTube.hh"
#include "G4GeomTools.hh"

#include "G4ScaledSolid.hh"
#include "G4Tubs.hh"
#include "G4Timer.hh"

// Function to print out an error message, can be used instead of assert().
// Example of use:
//   if (dist != 0) ERROR("( 5) DistanceToIn", pnt, dir, dist, 0);
//
void ERROR(const G4String func,
           const G4ThreeVector &p,
           const G4ThreeVector &v,
           double dist,
           double expected)
{
  G4cout << std::setprecision(10) << "ERROR! " << func << "(" << p << "," << v << ") = "
         << dist << " instead of " << expected << G4endl;
}

G4bool ApproxEqual(double a, const double b) {
  return std::abs(a - b) < 1e-10;
}

G4bool ApproxEqual(const G4ThreeVector &a, const G4ThreeVector &b) {
  return (a - b).mag2() < 1e-20;
}

bool TestEllipticalTube() {
  G4double kTolerance =
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // Check surfce area and volume
  //
  std::cout << "\n=== Check Get/Set, StreamInfo, GetSurfaceArea, GetCubicVolume, BoundingLimits" << std::endl;

  // Check Get
  G4EllipticalTube tube("Test_Elliptical_Tube", 1., 2., 3.);
  assert(tube.GetDx() == 1.);
  assert(tube.GetDy() == 2.);
  assert(tube.GetDz() == 3.);

  // Check Set
  double a = 4., b = 3., z = 5.;
  tube.SetDx(a);
  tube.SetDy(b);
  tube.SetDz(z);
  assert(tube.GetDx() == a);
  assert(tube.GetDy() == b);
  assert(tube.GetDz() == z);

  // Print information on Elliptical tube
  tube.StreamInfo(std::cout);

  // Check GetSurfaceArea
  double area = tube.GetSurfaceArea();
  std::cout << "SufaceArea : " << area << std::endl;
  double sbase = CLHEP::pi * a * b;
  double sside = 2. * z * G4GeomTools::EllipsePerimeter(a, b);
  assert(area == 2. * sbase + sside);

  // Check GetSurfaceArea
  double vol = tube.GetCubicVolume();
  std::cout << "CubicVolume : " << vol << std::endl;
  assert(vol == 2. * z * sbase);

  // Check GetBoundingLimits
  G4ThreeVector bmin, bmax;
  tube.BoundingLimits(bmin, bmax);
  std::cout << "BoundingLimits : " << bmin << ", " << bmax << std::endl;
  assert(bmax == G4ThreeVector(a, b, z));
  assert(bmin == -bmax);

  // Check Inside()
  //
  std::cout << "=== Check Inside()" << std::endl;
  int izmax    = 20;                       // number of steps along Z
  int iphimax  = 36;                       // number of steps along Phi
  int iscmax   = 1000;                     // number of scale factors
  double dz    = 2. * z / izmax;           // step along Z
  double dphi  = twopi / iphimax;          // step along Phi
  double dsc   = 1. / iscmax;              // scale factor increment
  double delta = 0.999 * 0.5 * kTolerance; // shift within surface
  double error = kTolerance / 1000.;       // calculation error tolerance

  // Check inside points
  for (int iz = 0; iz <= izmax; ++iz) {
    double curz = iz * dz - z;
    for (int iphi = 0; iphi < iphimax; ++iphi) {
      double phi  = iphi * dphi;
      double curx = a * std::cos(phi);
      double cury = b * std::sin(phi);
      for (double isc = 0; isc < iscmax; ++isc) {
        double scale = isc * dsc;
	assert(tube.Inside(G4ThreeVector(curx, cury, curz) * scale) == kInside);
      }
    }
  }

  // Check points on surface
  for (int iphi = 0; iphi < iphimax; ++iphi) {
    double phi  = iphi * dphi;
    double curx = (a + delta) * std::cos(phi);
    double cury = (b + delta) * std::sin(phi);
    for (double isc = 0; isc <= iscmax; ++isc) {
      double scale = isc * dsc;

      // base at -Z
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, -z)) == kSurface);
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, -z - delta)) == kSurface);
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, -z + delta)) == kSurface);

      // base at +Z
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, z)) == kSurface);
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, z - delta)) == kSurface);
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, z + delta)) == kSurface);
    }
  }

  // lateral surface
  for (int iphi = 0; iphi < iphimax; ++iphi) {
    double phi  = iphi * dphi;
    double curx = a * std::cos(phi);
    double cury = b * std::sin(phi);
    for (int iz = 0; iz <= izmax; ++iz) {
      double curz = iz * dz - z;
      assert(tube.Inside(G4ThreeVector(curx, cury, curz)) == kSurface);
    }
  }
  for (int iphi = 0; iphi < iphimax; ++iphi) {
    double phi  = iphi * dphi;
    double curx = (a + delta) * std::cos(phi);
    double cury = (b + delta) * std::sin(phi);
    for (int iz = 0; iz <= izmax; ++iz) {
      double curz = iz * dz - z;
      assert(tube.Inside(G4ThreeVector(curx, cury, curz)) == kSurface);
    }
  }
  for (int iphi = 0; iphi < iphimax; ++iphi) {
    double phi  = iphi * dphi;
    double curx = (a - delta) * std::cos(phi);
    double cury = (b - delta) * std::sin(phi);
    for (int iz = 0; iz <= izmax; ++iz) {
      double curz = iz * dz - z;
      assert(tube.Inside(G4ThreeVector(curx, cury, curz)) == kSurface);
    }
  }

  // Check outside points
  for (int iphi = 0; iphi < iphimax; ++iphi) {
    double phi  = iphi * dphi;
    double curx = (a + kTolerance) * std::cos(phi);
    double cury = (b + kTolerance) * std::sin(phi);
    for (double isc = 0; isc <= iscmax; ++isc) {
      double scale = isc * dsc;
      // near base at -Z
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, -z - kTolerance)) == kOutside);

      // near base at +Z
      assert(tube.Inside(G4ThreeVector(curx * scale, cury * scale, z + kTolerance)) == kOutside);
    }
  }
  for (int iz = 0; iz <= izmax; ++iz) {
    double curz = iz * dz - z;
    for (int iphi = 0; iphi < iphimax; ++iphi) {
      double phi  = iphi * dphi;
      double curx = a * std::cos(phi);
      double cury = b * std::sin(phi);
      assert(tube.Inside(G4ThreeVector(curx, cury, curz) * 1.001) == kOutside);
    }
  }

  // Check SurfaceNormal()
  //
  std::cout << "=== Check SurfaceNormal()" << std::endl;
  G4ThreeVector normal;

  // points on elliptical section
  double scaleX = tube.GetDx();
  double scaleY = tube.GetDy();
  for (int i = 0; i < 360; ++i) {
    double phi = i * deg;
    double nx = std::cos(phi);
    double ny = std::sin(phi);
    double px = nx * scaleX;
    double py = ny * scaleY;
    normal = tube.SurfaceNormal(G4ThreeVector(px, py, 0));
    assert(ApproxEqual(normal,G4ThreeVector(nx * scaleY, ny * scaleX, 0).unit()));
  }

  // points near axes
  for (int i = 0; i < 21; ++i) {
    double curx = tube.GetDx() + 0.1 * kTolerance * (i - 10);
    normal = tube.SurfaceNormal(G4ThreeVector(curx, 0, 0));
    assert(ApproxEqual(normal, G4ThreeVector(1, 0, 0)));
    normal = tube.SurfaceNormal(G4ThreeVector(-curx, 0, 0));
    assert(ApproxEqual(normal, G4ThreeVector(-1, 0, 0)));

    double cury = tube.GetDy() + 0.1 * kTolerance * (i - 10);
    normal = tube.SurfaceNormal(G4ThreeVector(0, cury, 0));
    assert(ApproxEqual(normal, G4ThreeVector(0, 1, 0)));
    normal = tube.SurfaceNormal(G4ThreeVector(0, -cury, 0));
    assert(ApproxEqual(normal, G4ThreeVector(0, -1, 0)));

    double curz = tube.GetDz() + 0.1 * kTolerance * (i - 10);
    normal = tube.SurfaceNormal(G4ThreeVector(0, 0, curz));
    assert(normal == G4ThreeVector(0, 0, 1));
    normal = tube.SurfaceNormal(G4ThreeVector(0, 0, -curz));
    assert(normal == G4ThreeVector(0, 0, -1));
  }

  // point on edge
  normal = tube.SurfaceNormal(G4ThreeVector(tube.GetDx(), 0, tube.GetDz()));
  assert(normal == G4ThreeVector(1, 0, 1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(0, tube.GetDy(), tube.GetDz()));
  assert(normal == G4ThreeVector(0, 1, 1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(-tube.GetDx(), 0, tube.GetDz()));
  assert(normal == G4ThreeVector(-1, 0, 1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(0, -tube.GetDy(), tube.GetDz()));
  assert(normal == G4ThreeVector(0, -1, 1).unit());

  normal = tube.SurfaceNormal(G4ThreeVector(tube.GetDx(), 0, -tube.GetDz()));
  assert(normal == G4ThreeVector(1, 0, -1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(0, tube.GetDy(), -tube.GetDz()));
  assert(normal == G4ThreeVector(0, 1, -1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(-tube.GetDx(), 0, -tube.GetDz()));
  assert(normal == G4ThreeVector(-1, 0, -1).unit());
  normal = tube.SurfaceNormal(G4ThreeVector(0, -tube.GetDy(), -tube.GetDz()));
  assert(normal == G4ThreeVector(0, -1, -1).unit());

  // special case of point on Z-axis
  assert(tube.SurfaceNormal(G4ThreeVector(0, 0, 0)) == G4ThreeVector(0, 0, 1));
  assert(tube.SurfaceNormal(G4ThreeVector(0, 0, -kTolerance)) == G4ThreeVector(0, 0, -1));

  // Check SafetyToIn()
  //
  std::cout << "=== Check SafetyToIn()" << std::endl;
  double safety;
  assert(a >= b);

  izmax    = 20;                    // number of steps along Z
  iphimax  = 36;                    // number of steps along Phi
  iscmax   = 1000;                  // number of scale factors
  dz    = 2. * z / izmax;           // step along Z
  dphi  = twopi / iphimax;          // step along Phi
  dsc   = 1. / iscmax;              // scale factor increment
  delta = 0.999 * 0.5 * kTolerance; // shift within surface
  error = kTolerance / 1000.;       // calculation error tolerance

  for (int iz = 1; iz < izmax; ++iz) {
    double curz = iz * dz - z;
    for (int iphi = 0; iphi < iphimax; ++iphi) {
      double phi = iphi * dphi;

      // Check outside point
      double curx = 2. * a * std::cos(phi);
      double cury = 2. * b * std::sin(phi);
      safety = tube.DistanceToIn(G4ThreeVector(curx, cury, curz));
      assert(safety > 0);
      assert(safety >= b - error);

      // Check points on surface, within tolerance at both side of the boundary
      curx = (a + delta) * std::cos(phi);
      cury = (b - delta) * std::sin(phi);
      safety = tube.DistanceToIn(G4ThreeVector(curx, cury, curz));
      assert(safety >= 0);
      assert(safety <= delta);

      // Check inside points
      curx = 0.999 * a * std::cos(phi);
      cury = 0.999 * b * std::sin(phi);
      safety = tube.DistanceToIn(G4ThreeVector(curx, cury, curz));
      assert(safety == 0);
    }
  }

  // points around lateral surface
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, 2 * b, 0)) == b);
  assert(tube.DistanceToIn(G4ThreeVector(2 * a, 0, 0)) >= b);

  assert(tube.DistanceToIn(G4ThreeVector(a, 0, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(a - delta, 0, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(a + delta, 0, 0)) <= delta);

  assert(tube.DistanceToIn(G4ThreeVector(0, b, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, b - delta, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, b + delta, 0)) <= delta);

  assert(tube.DistanceToIn(G4ThreeVector(a / 2, b / 2, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(a / 2, 0, 0)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, b / 2, 0)) == 0);

  // points around the bases
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, -z - 1)) == 1);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, z + 2)) == 2);

  assert(tube.DistanceToIn(G4ThreeVector(0, 0, -z)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, -z - delta)) <= delta);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, -z + delta)) == 0);

  assert(tube.DistanceToIn(G4ThreeVector(0, 0, z)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, z - delta)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, z + delta)) <= delta);

  assert(tube.DistanceToIn(G4ThreeVector(0, 0, -z + 1)) == 0);
  assert(tube.DistanceToIn(G4ThreeVector(0, 0, z - 2)) == 0);

  assert(tube.DistanceToIn(G4ThreeVector(-3 * a, 0, z + 10)) >= 10);
  assert(tube.DistanceToIn(G4ThreeVector(0, -3 * b, z + 10)) >= 10);
  assert(tube.DistanceToIn(G4ThreeVector(-3 * a, 0, -z - 1)) >= 2 * b);
  assert(tube.DistanceToIn(G4ThreeVector(0, -3 * b, -z - 2)) >= 2 * b);

  // Check SafetyToOut()
  //
  std::cout << "=== Check SafetyToOut()" << std::endl;

  izmax    = 20;                    // number of steps along Z
  iphimax  = 36;                    // number of steps along Phi
  iscmax   = 1000;                  // number of scale factors
  dz    = 2. * z / izmax;           // step along Z (0.5)
  dphi  = twopi / iphimax;          // step along Phi
  dsc   = 1. / iscmax;              // scale factor increment
  delta = 0.999 * 0.5 * kTolerance; // shift within surface
  error = kTolerance / 1000.;       // calculation error tolerance

  for (int iz = 1; iz < izmax; ++iz) {
    double curz = iz * dz - z;
    for (int iphi = 0; iphi < iphimax; ++iphi) {
      double phi = iphi * dphi;

      // Check inside point
      double curx = 0.9 * a * std::cos(phi);
      double cury = 0.9 * b * std::sin(phi);
      safety = tube.DistanceToOut(G4ThreeVector(curx, cury, curz));
      assert(safety > 0);
      assert(safety >= (1 - 0.9) * b - error);

      // Check points on surface, within tolerance at both side of the boundary
      curx = (a + delta) * std::cos(phi);
      cury = (b - delta) * std::sin(phi);
      safety = tube.DistanceToOut(G4ThreeVector(curx, cury, curz));
      assert(safety >= 0);
      assert(safety <= delta);

      // Check outside points
      curx = 1.001 * a * std::cos(phi);
      cury = 1.001 * b * std::sin(phi);
      safety = tube.DistanceToOut(G4ThreeVector(curx, cury, curz));
      assert(safety == 0);
    }
  }

  // points around lateral surface
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, 0)) == b);
  assert(tube.DistanceToOut(G4ThreeVector(0, b / 2, 0)) == b / 2);
  assert(tube.DistanceToOut(G4ThreeVector(a / 2, 0, 0)) >= b / 2);

  assert(tube.DistanceToOut(G4ThreeVector(a, 0, 0)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(a - delta, 0, 0)) <= delta);
  assert(tube.DistanceToOut(G4ThreeVector(a + delta, 0, 0)) == 0);

  assert(tube.DistanceToOut(G4ThreeVector(0, b, 0)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, b - delta, 0)) <= delta);
  assert(tube.DistanceToOut(G4ThreeVector(0, b + delta, 0)) == 0);

  assert(tube.DistanceToOut(G4ThreeVector(2 * a, 2 * b, 0)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(-a, -b, z / 2)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(2 * a, 0, 0)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, 2 * b, 0)) == 0);

  // points around the bases
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, -z + 1)) == 1);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, z - 2)) == 2);

  assert(tube.DistanceToOut(G4ThreeVector(0, 0, -z)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, -z - delta)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, -z + delta)) <= delta);

  assert(tube.DistanceToOut(G4ThreeVector(0, 0, z)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, z - delta)) <= delta);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, z + delta)) == 0);

  assert(tube.DistanceToOut(G4ThreeVector(0, 0, -z - 1)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, 0, z + 2)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(a, 0, -z - 1)) == 0);
  assert(tube.DistanceToOut(G4ThreeVector(0, b, z + 2)) == 0);

  // Check DistanceToIn()
  //
  std::cout << "=== Check DistanceToIn()" << std::endl;
  tube.SetDx(a = 4.);
  tube.SetDy(b = 3.);
  tube.SetDz(z = 5.);
  double Rsph = std::sqrt(a * a + b * b + z * z) + 0.1; // surrounding sphere
  double dist, del = kTolerance / 3.;
  G4ThreeVector pnt, dir;

  // set coordinates for points in grid
  static const int np = 11;
  double xxx[np] = {-a - 1, -a - del, -a, -a + del, -1, 0, 1, a - del, a, a + del, a + 1};
  double yyy[np] = {-b - 1, -b - del, -b, -b + del, -1, 0, 1, b - del, b, b + del, b + 1};
  double zzz[np] = {-z - 1, -z - del, -z, -z + del, -1, 0, 1, z - del, z, z + del, z + 1};

  // check directions parallel to Z axis
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        pnt.set(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector rho(pnt.x(), pnt.y(), 0);

        dist = tube.DistanceToIn(pnt, dir = G4ThreeVector(0, 0, 1)); // +Z direction
        if (tube.Inside(pnt) == kInside) { // point is Inside
           assert(dist == 0);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          if (tube.Inside(rho) != kInside) {
             assert(dist == kInfinity);
          } else {
            if (pnt.z() > 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == 0);
            }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          if (tube.Inside(rho) != kInside) {
             assert(dist == kInfinity);
          } else {
            if (pnt.z() > 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == (-pnt.z() - z));
            }
          }
        }

        dist = tube.DistanceToIn(pnt, dir = G4ThreeVector(0, 0, -1)); // -Z direction
        if (tube.Inside(pnt) == kInside) { // point is Inside
          assert(dist == 0);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          if (tube.Inside(rho) != kInside) {
             assert(dist == kInfinity);
          } else {
            if (pnt.z() < 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == 0);
            }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          if (tube.Inside(rho) != kInside) {
             assert(dist == kInfinity);
          } else {
            if (pnt.z() < 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == (pnt.z() - z));
            }
          }
        }
      }
    }
  }

  // check directions parallel to X axis
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        pnt.set(xxx[ix], yyy[iy], zzz[iz]);

        G4ThreeVector z00(0, 0, pnt.z());
        G4ThreeVector y00(0, pnt.y(), 0);
        double tmp  = (1 - pnt.y() / b) * (1 + pnt.y() / b);
        double intx = (tmp < 0) ? 0 : std::sqrt(tmp) * a;

        dist = tube.DistanceToIn(pnt, dir = G4ThreeVector(1, 0, 0)); // +X direction
        if (tube.Inside(pnt) == kInside) {  // point is Inside
          assert(dist == 0);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          if (tube.Inside(z00) != kInside) {
            assert(dist == kInfinity);
          } else {
            if (pnt.x() >= 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == 0);
            }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          if (tube.Inside(z00) != kInside) {
	    assert(dist == kInfinity);
          } else {
            if (pnt.x() >= 0 || tube.Inside(y00) != kInside) {
              assert(dist == kInfinity);
            } else {
              assert(ApproxEqual(dist, (-pnt.x() - intx)));
            }
          }
        }

        dist = tube.DistanceToIn(pnt, dir = G4ThreeVector(-1, 0, 0)); // -X direction
        if (tube.Inside(pnt) == kInside) { // point is Inside
          assert(dist == 0);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          if (tube.Inside(z00) != kInside) {
            assert(dist == kInfinity);
          } else {
            if (pnt.x() <= 0) {
              assert(dist == kInfinity);
            } else {
              assert(dist == 0);
            }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          if (tube.Inside(z00) != kInside) {
	    assert(dist == kInfinity);
          } else {
            if (pnt.x() <= 0 || tube.Inside(y00) != kInside) {
              assert(dist == kInfinity);
            } else {
              assert(ApproxEqual(dist, (pnt.x() - intx)));
            }
          }
        }
      }
    }
  }

  // Check inside points ("wrong side")
  for (int ith = 5; ith < 180; ith += 30) {
    for (int iph = 5; iph < 360; iph += 30) {
      double theta = ith * deg;
      double phi   = iph * deg;
      double vx    = std::sin(theta) * std::cos(phi);
      double vy    = std::sin(theta) * std::sin(phi);
      double vz    = std::cos(theta);
      dir.set(vx, vy, vz);
      dist = tube.DistanceToIn(G4ThreeVector(0, 0, 0), dir);
      assert(dist == 0);
    }
  }

  // Check outside points
  for (int ith = 5; ith < 180; ith += 30) {
    for (int iph = 5; iph < 360; iph += 30) {
      double theta = ith * deg;
      double phi   = iph * deg;
      double vx    = std::sin(theta) * std::cos(phi);
      double vy    = std::sin(theta) * std::sin(phi);
      double vz    = std::cos(theta);
      dir.set(vx, vy, vz); // direction
      pnt = Rsph * dir;

      G4ThreeVector axis(1, 0, 0);
      if (std::abs(vy) < std::abs(vx) && std::abs(vy) < std::abs(vz)) axis.set(0, 1, 0);
      if (std::abs(vz) < std::abs(vx) && std::abs(vz) < std::abs(vy)) axis.set(0, 0, 1);
      G4ThreeVector ort = (axis.cross(dir)).unit(); // orthogonal direction

      dist = tube.DistanceToIn(pnt, dir);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, ort);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, -ort);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, -dir);
      assert(dist > 0 && dist < Rsph);
      assert(tube.Inside(pnt - dist * dir) == kSurface);
    }
  }

  // Check "very far" outside points
  double offset = 1.e+6; // 1km
  for (int ith = 0; ith <= 180; ith += 10) {
    for (int iph = 0; iph < 360; iph += 10) {
      double theta = ith * deg;
      double phi   = iph * deg;
      double vx    = std::sin(theta) * std::cos(phi);
      double vy    = std::sin(theta) * std::sin(phi);
      double vz    = std::cos(theta);
      dir.set(vx, vy, vz); // direction
      pnt = (offset + Rsph) * dir;

      G4ThreeVector axis(1, 0, 0);
      if (std::abs(vy) < std::abs(vx) && std::abs(vy) < std::abs(vz)) axis.set(0, 1, 0);
      if (std::abs(vz) < std::abs(vx) && std::abs(vz) < std::abs(vy)) axis.set(0, 0, 1);
      G4ThreeVector ort = (axis.cross(dir)).unit(); // orthogonal direction

      dist = tube.DistanceToIn(pnt, dir);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, ort);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, -ort);
      assert(dist == kInfinity);

      dist = tube.DistanceToIn(pnt, -dir);
      assert(tube.Inside(pnt - dist * dir) == kSurface);
    }
  }

  // Check DistanceToOut()
  //
  std::cout << "=== Check DistanceToOut()" << std::endl;
  G4bool validNorm = false;
  G4ThreeVector norm(0, 0, 0);


  // check directions parallel to Z axis
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        validNorm = false;
        pnt.set(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector z00(0, 0, pnt.z());

        validNorm = false;
        dist = tube.DistanceToOut(pnt, dir = G4ThreeVector(0, 0, 1), true, &validNorm, &norm); // +Z direction
        if (tube.Inside(pnt) == kInside) {  // point is Inside
           assert(dist == (z - pnt.z()));
           assert(validNorm == true);
           assert(norm == dir);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          assert(validNorm == true);
          assert(norm == dir);
          if (tube.Inside(z00) == kInside) {
            assert(dist == (z - pnt.z()));
          } else {
            if (pnt.z() > 0) assert (dist == 0);
            if (pnt.z() < 0) { assert(dist == (z - pnt.z())); }
          }
	}
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          assert(dist == 0);
	}

        validNorm = false;
        dist = tube.DistanceToOut(pnt, dir = G4ThreeVector(0, 0, -1), true, &validNorm, &norm); // -Z direction
        if (tube.Inside(pnt) == kInside) { // point is Inside
          assert(dist == (z + pnt.z()));
          assert(validNorm == true);
          assert(norm == dir);
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          assert(validNorm == true);
          assert(norm == dir);
          if (tube.Inside(z00) == kInside) {
            assert(dist == (z + pnt.z()));
          } else {
            if (pnt.z() < 0) assert (dist == 0);
            if (pnt.z() > 0) { assert(dist == (z + pnt.z())); }
	  }
	}
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          assert(dist == 0);
	}
      }
    }
  }

  // check directions parallel to X axis
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        pnt.set(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector rho(pnt.x(), pnt.y(), 0);

        G4ThreeVector z00(0, 0, pnt.z());
        G4ThreeVector y00(0, pnt.y(), 0);
        double tmp  = (1 - (yyy[iy] / b)) * (1 + (yyy[iy] / b));
        double intx = (tmp < 0) ? 0 : std::sqrt(tmp) * a;

        validNorm = false;
        dist = tube.DistanceToOut(pnt, dir = G4ThreeVector(1, 0, 0), true, &validNorm, &norm); // +X direction
        if (tube.Inside(pnt) == kInside) {  // point is Inside
          assert(ApproxEqual(dist, intx - pnt.x()));
          assert(validNorm == true);
          assert(ApproxEqual(norm, tube.SurfaceNormal(G4ThreeVector(pnt.x() + dist, pnt.y(), 0))));
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          assert(validNorm == true);
          assert(ApproxEqual(norm, tube.SurfaceNormal(G4ThreeVector(pnt.x() + dist, pnt.y(), 0))));
          if (tube.Inside(rho) == kInside) {
	     assert(ApproxEqual(dist, intx - pnt.x()));
          } else {
	    if (pnt.x() > 0) {
              assert(dist == 0);
              assert(ApproxEqual(norm, tube.SurfaceNormal(G4ThreeVector(pnt.x(), pnt.y(), 0))));
            } else {
	      assert(ApproxEqual(dist, intx - pnt.x()));
	    }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          assert(dist == 0);
	}

        validNorm = false;
        dist = tube.DistanceToOut(pnt, dir = G4ThreeVector(-1, 0, 0), true, &validNorm, &norm); // -X direction
        if (tube.Inside(pnt) == kInside) {  // point is Inside
          assert(ApproxEqual(dist, intx + pnt.x()));
          assert(validNorm == true);
          assert(ApproxEqual(norm, tube.SurfaceNormal(G4ThreeVector(pnt.x() - dist, pnt.y(), 0))));
        }
        if (tube.Inside(pnt) == kSurface) { // point is on Surface
          assert(validNorm == true);
          assert(ApproxEqual(norm, tube.SurfaceNormal(G4ThreeVector(pnt.x() - dist, pnt.y(), 0))));
          if (tube.Inside(rho) == kInside) {
	     assert(ApproxEqual(dist, intx + pnt.x()));
          } else {
	    if (pnt.x() < 0) {
              assert(dist == 0);
            } else {
              assert(ApproxEqual(dist, intx + pnt.x()));
	    }
          }
        }
        if (tube.Inside(pnt) == kOutside) { // point is Outside
          assert(dist == 0);
	}
      }
    }
  }

  // Check inside points
  for (int ith = 5; ith < 180; ith += 30) {
    for (int iph = 5; iph < 360; iph += 30) {
      double theta = ith * deg;
      double phi   = iph * deg;
      double vx    = std::sin(theta) * std::cos(phi);
      double vy    = std::sin(theta) * std::sin(phi);
      double vz    = std::cos(theta);
      dir.set(vx, vy, vz);
      dist = tube.DistanceToOut(G4ThreeVector(0, 0, 0), dir);
      assert(dist > 0 && dist < Rsph);
      assert(tube.Inside(dist * dir) == kSurface);
    }
  }

  // Check outside points ("wrong side")
  for (int ith = 5; ith < 180; ith += 30) {
    for (int iph = 5; iph < 360; iph += 30) {
      double theta = ith * deg;
      double phi   = iph * deg;
      double vx    = std::sin(theta) * std::cos(phi);
      double vy    = std::sin(theta) * std::sin(phi);
      double vz    = std::cos(theta);
      dir.set(vx, vy, vz);

      dist = tube.DistanceToOut(Rsph * dir, dir);
      // assert(dist == 0); // convention

      dist = tube.DistanceToOut(Rsph * dir, -dir);
      // assert(dist == 0); // convention
    }
  }

  // Check SamplePointOnSurface()
  //
  std::cout << "=== Check GetPointOnSurface()" << std::endl;

  G4Timer time;
  clock_t t;

  t = clock();
  int nzneg = 0, nzpos = 0, nside = 0, nfactor = 10000,
      ntot = 4 * area * nfactor;
  for (int i = 0; i < ntot; i++) {
    G4ThreeVector rndPoint = tube.GetPointOnSurface();
    assert(tube.Inside(rndPoint) == kSurface);
    if (rndPoint.x() < 0 || rndPoint.y() < 0)
      continue;
    if (rndPoint.z() == -z)
      ++nzneg;
    else if (rndPoint.z() == z)
      ++nzpos;
    else
      ++nside;
  }
  t = clock() - t;

  std::cout << "szneg,sside,szpos = " << sbase << ", \t" << sside << ", \t"
            << sbase << std::endl;
  std::cout << "nzneg,nside,nzpos = " << nzneg << ", \t" << nside << ", \t"
            << nzpos << std::endl;
  G4cout <<"Time for " << ntot << " points : "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;

  assert(std::abs(nzneg - sbase * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nside - sside * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nzpos - sbase * nfactor) < 2. * std::sqrt(ntot));
  return true;
}

int main(int argc, char *argv[]) {
  assert(TestEllipticalTube());
  std::cout << "\n    All tests for G4EllipticalTube passed!\n" << std::endl;

  return 0;
}
