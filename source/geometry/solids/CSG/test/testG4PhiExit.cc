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

// Regression for exits from angular tolerance bands and rounded tangencies.
#undef NDEBUG
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4CutTubs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <cassert>
#include <cmath>
#include <random>

void Check(G4VSolid& solid, double start, double span)
{
  for (int end : {0, 1})
  {
    double phi = start + end*span;
    G4ThreeVector radial(std::cos(phi), std::sin(phi), 0);
    G4ThreeVector normal = (end ? 1. : -1.)
      * G4ThreeVector(-std::sin(phi), std::cos(phi), 0);
    for (double offset : {-4e-9, 0., 4e-9})
    {
      auto p = 1200*radial + offset*normal;
      assert(solid.Inside(p) == kSurface);
      for (double slope : {1., 1e-8, 1e-12})
      {
        auto v = (radial + slope*normal).unit();
        bool valid = false;
        G4ThreeVector n;
        double distance = solid.DistanceToOut(p, v, true, &valid, &n);
        double expected = offset < 0 ? -offset / v.dot(normal) : 0.;
        // At the shallowest slope the radial surface is reached first.
        if (expected > 200.) continue;
        assert(std::abs(distance - expected) < 1e-4);
        assert(distance == solid.DistanceToOut(p, v));
        assert(valid == (span <= pi));
        if (valid) assert((n-normal).mag() < 1e-12);
      }
      assert(solid.DistanceToIn(p, (radial-normal).unit()) < 1e-7);
      assert(solid.DistanceToOut(p, -normal) > 1.);
      assert(solid.DistanceToOut(p, radial) > 100.);
    }
  }
}

// Oblique caps clip a long tube: their plane distances are independent of
// the radial/phi exit, including when the ray is parallel to either cap.
void CheckCuts(G4CutTubs& solid, double start, double span)
{
  G4Tubs reference("reference", 100., 1450., 1000., start, span);
  const G4ThreeVector low(.05, .03, -1), high(-.04, .02, 1);
  std::mt19937 engine(81723);
  std::uniform_real_distribution<double> uniform(0., 1.);
  for (int i = 0; i < 4096; ++i)
  {
    double phi = start + span*(.02 + .96*uniform(engine));
    double radius = 300. + 1000.*uniform(engine);
    G4ThreeVector p(radius*std::cos(phi), radius*std::sin(phi),
                    200.*uniform(engine)-100.);
    G4ThreeVector v(2*uniform(engine)-1, 2*uniform(engine)-1,
                    2*uniform(engine)-1);
    if (i%4 == 0) v = low.cross(v);
    if (i%4 == 1) v = high.cross(v);
    v = v.unit();
    assert(solid.Inside(p) == kInside);
    double expected = reference.DistanceToOut(p, v);
    for (const auto& normal : {low, high})
      if (normal.dot(v) > 0.)
        expected = std::min(expected, (250.-normal.dot(p))/normal.dot(v));
    bool valid = false;
    G4ThreeVector n;
    double distance = solid.DistanceToOut(p, v, true, &valid, &n);
    assert(std::abs(distance-expected) < 1e-7);
    assert(distance == solid.DistanceToOut(p, v));
    assert(solid.DistanceToOut(p) <= expected + 1e-9);
    if (valid) assert(std::abs(n.mag()-1.) < 1e-12 && n.dot(v) >= -1e-12);
  }
}

int main()
{
  for (double start : {0., -5*degree, 75*degree, 266*degree})
  for (double span : {20*degree, pi, 270*degree, 359.99*degree})
  {
    G4Tubs tube("tube", 100., 1450., 250., start, span);
    G4Cons cone("cone", 100., 1450., 200., 1450., 250., start, span);
    G4CutTubs flat("flat", 100., 1450., 250., start, span,
                   G4ThreeVector(0, 0, -1), G4ThreeVector(0, 0, 1));
    G4CutTubs oblique("oblique", 100., 1450., 250., start, span,
                      G4ThreeVector(.05, .03, -1).unit(),
                      G4ThreeVector(-.04, .02, 1).unit());
    Check(flat, start, span);
    Check(oblique, start, span);
    CheckCuts(oblique, start, span);
    Check(tube, start, span);
    Check(cone, start, span);
    G4Cons tapered("tapered", 100., 1450., 200., 1500., 250., start, span);
    Check(tapered, start, span);
  }
  G4Tubs distant("distant", 0., 15.000000372529, 1367.2500313152, 0., twopi);
  G4ThreeVector far(-2718368898.9392972, -296582541.36111861, 623.29061689710352);
  G4ThreeVector toward(.99410086734936876, .10845951104090816, 0.);
  assert(std::abs(distant.DistanceToIn(far, toward)-2734500062.6303997) < 1e-4);
  assert(distant.DistanceToIn({1e10, 15.001, 0}, {-1, 0, 0}) == kInfinity);
  G4Cons pad("pad", 1081.35, 1091.95, 1081.35, 1091.95, 142.5,
             -5.536333547024*degree, 11.072667094048*degree);
  assert(std::abs(pad.DistanceToIn(
    {-1088666540.4637313, -84602627.711171284, -90.47063659876585},
    {.99699401855690595, .077478558077393719, 0.})-1091950000.) < 1e-4);
  G4Tubs wedge("wedge", 0., 1135., 315.,
                255.001338739853*degree, 29.9973225202933*degree);
  G4ThreeVector edge(-126.96146625221488, -473.87093113401528, -98.050613247323781);
  G4ThreeVector inward(.25879647577093079, .9659318734406408, 0.);
  assert((edge+wedge.DistanceToOut(edge, inward)*inward).perp() < 1e-7);
  G4Tubs ring("ring", 100., 110., 10., 0., twopi);
  G4CutTubs cutRing("cutRing", 100., 110., 100., 0., twopi,
                    G4ThreeVector(.05, 0, -1).unit(), G4ThreeVector(0, .05, 1).unit());
  for (double offset : {0., .25e-9})
  {
    assert(ring.DistanceToIn({100.+offset, 0, 0}, {0, 1, 0}) == 0.);
    assert(cutRing.DistanceToIn({100.+offset, 0, 0}, {0, 1, 0}) == 0.);
  }
  G4CutTubs full("full", 100., 1450., 250., 0., twopi,
                 G4ThreeVector(.05, .03, -1).unit(),
                 G4ThreeVector(-.04, .02, 1).unit());
  CheckCuts(full, 0., twopi);
  G4Cons cone("tangent", 124.95, 1450., 1020., 1450., 224.75,
              -5*degree, 10*degree);
  G4ThreeVector p(1449.6736655984491, 30.761392529530831, 65.163274140795693);
  G4ThreeVector v(0.021214753468641951, -0.99977494179203386, 0);
  assert(full.Inside(p) == kSurface);
  assert(full.DistanceToOut(p, v) < 1e-7);
  assert(cone.Inside(p) == kSurface);
  bool valid = false;
  G4ThreeVector n;
  assert(cone.DistanceToOut(p, v, true, &valid, &n) == 0.);
  assert(valid && std::abs(n.mag()-1.) < 1e-12);
}
