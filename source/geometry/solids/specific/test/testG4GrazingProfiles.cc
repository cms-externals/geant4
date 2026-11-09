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

// Finite cone/cylinder joins, including the shared GenericPolycone kernel.
#undef NDEBUG
#include "G4Polycone.hh"
#include "G4GenericPolycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4PhysicalConstants.hh"
#include <cassert>
#include <cmath>

// Convex end caps must participate in classification just outside their Z.
void CheckGenericNavigation()
{
  for (int sections : {1, 100})
  {
    std::vector<double> r, z;
    for (int i = 0; i <= sections; ++i)
    {
      const double zi = -227. + 454.*i/sections;
      r.push_back(211. + 10.*(1.-zi*zi/(227.*227.)));
      z.push_back(zi);
    }
    r.insert(r.end(), {26., 26.});
    z.insert(z.end(), {227., -227.});
    for (double span : {twopi, 0.5*pi, 1.5*pi})
    {
      G4GenericPolycone tube("generic caps", 0., span, r.size(), r.data(), z.data());
      const double phi = 0.5*span;
      for (double end : {-1., 1.})
      for (double offset : {-4e-9, -.25e-9, 0., .25e-9, 4e-9})
      {
        G4ThreeVector p(120.*std::cos(phi), 120.*std::sin(phi), end*(227.+offset));
        const EInside expected = offset > .5e-9 ? kOutside
                               : offset < -.5e-9 ? kInside : kSurface;
        assert(tube.Inside(p) == expected);
        const double safety = expected == kOutside ? tube.DistanceToIn(p)
                                                   : tube.DistanceToOut(p);
        assert(safety >= 0. && safety <= std::abs(offset)+1e-12);
      }
    }
  }
  // Recorded CMS distant entry: cancellation previously displaced the hit.
  double cr[] = {1197.2, 1275.5, 1275.5, 285., 285., 285.};
  double cz[] = {3066.7, 3145., 3205., 3205., 3145., 3066.7};
  G4GenericPolycone cone("distant generic", 0., twopi, 6, cr, cz);
  G4ThreeVector p(896242520.058671, 101003171.42452301, -901911626.68683863);
  G4ThreeVector v(-.70265881960858034, -.079187014252830187, .70710678118654657);
  assert(std::abs(cone.DistanceToIn(p, v)-1275500000.) < 1e-4);

  // Non-monotonic Z: two radial intervals joined by the bottom of a U.
  double ur[] = {10., 30., 30., 25., 25., 15., 15., 10.};
  double uz[] = {-20., -20., 20., 20., -10., -10., 20., 20.};
  G4GenericPolycone overhang("overhang", 0., twopi, 8, ur, uz);
  assert(overhang.Inside({12., 0., 0.}) == kInside);
  assert(overhang.Inside({27., 0., 0.}) == kInside);
  assert(overhang.Inside({20., 0., 0.}) == kOutside);
  assert(overhang.Inside({20., 0., -15.}) == kInside);
  assert(std::abs(overhang.DistanceToIn({20., 0., 0.}, {1., 0., 0.})-5.) < 1e-9);
  assert(std::abs(overhang.DistanceToIn({1e9, 0., 0.}, {-1., 0., 0.})-(1e9-30.)) < 1e-6);
  assert(overhang.DistanceToIn({1e9, 0., 21.}, {-1., 0., 0.}) == kInfinity);
}

int main()
{
  CheckGenericNavigation();
  double z[] = {-.1, 608.551409512083, 932.300731880565, 1490.1};
  double ri[] = {0, 0, 0, 0}, ro[] = {265, 265, 219.5, 219.5};
  G4Polycone poly("poly", 0, twopi, 4, z, ri, ro);
  double r[] = {0, 265, 265, 219.5, 219.5, 0};
  double gz[] = {z[0], z[0], z[1], z[2], z[3], z[3]};
  G4GenericPolycone generic("generic", 0, twopi, 6, r, gz);
  for (G4VSolid* solid : {static_cast<G4VSolid*>(&poly), static_cast<G4VSolid*>(&generic)})
  for (double offset : {-.25e-9, 0., .25e-9})
  for (double slope : {-1e-12, 0., 1e-12})
  {
    G4ThreeVector p(219.5+offset, 0, 1367.7376507925883);
    G4ThreeVector v(slope, 0, -1);
    auto d = solid->DistanceToIn(p, v);
    assert(d >= 0. && d < 436.);
    assert(solid->Inside(p+d*v) != kOutside);
  }
  // A shallow cone end must not report a hit beyond its end plane.
  double cz[] = {0., 1000.}, cr[] = {100., 99.};
  G4Polycone cap("cap", 0., twopi, 2, cz, ri, cr);
  G4ThreeVector p(99.-.25e-9, 0, 1100), v(0, 0, -1);
  auto d = cap.DistanceToIn(p, v);
  assert(std::abs(d-100.) < 1e-9);
  assert(cap.Inside(p+d*v) != kOutside);
  std::vector<G4TwoVector> vertices = {
    {-1335.33, -2760.58},
    {-1335.33, -1154.7},
    {-1285.36, -1034.06},
    {-1017.44, -1145.04},
    {-990.65, -1080.37},
    {-694.43, -1203.06},
    {-467.05, -1203.06},
    {147.95, -587.13},
    {147.7, -423.3},
    {91.24, 222.13},
    {-574.73, 733.15},
    {-140.51, 1781.5},
    {-158.99, 1789.16},
    {-100.4, 1930.6},
    {-100.4, 2599.44},
    {521.14, 2599.44},
    {521.14, 2506.22},
    {924.48, 2359.42},
    {1024.48, 2359.42},
    {1065.9, 2459.42},
    {1534.62, 2459.42},
    {1534.62, 1659.41},
    {1212.31, 1659.43},
    {846.05, 72.99},
    {464.66, -847.76},
    {464.66, -2760.58},
  };
  G4ExtrudedSolid plate("plate", vertices, 41.3, {0, 0}, 1., {0, 0}, 1.);
  assert(plate.DistanceToOut(
    {-994.70717823467737, -1090.1638677274759, 29.108635903941469},
    {.38271768332605155, .92386534455489777, 0}) == 0.);
}
