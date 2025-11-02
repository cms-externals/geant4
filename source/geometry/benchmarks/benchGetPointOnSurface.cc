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
#include <iostream>
#include <iomanip>

// polyhedra
#include "G4Box.hh"              // 1
#include "G4Para.hh"             // 2
#include "G4Trd.hh"              // 3
#include "G4Trap.hh"             // 4
#include "G4Tet.hh"              // 5
#include "G4Polyhedra.hh"        // 6
#include "G4ExtrudedSolid.hh"    // 7
#include "G4TessellatedSolid.hh" // 8
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

// second order surfaces
#include "G4GenericTrap.hh"      // 1
#include "G4Orb.hh"              // 2
#include "G4Sphere.hh"           // 3
#include "G4Ellipsoid.hh"        // 4
#include "G4Tubs.hh"             // 5
#include "G4CutTubs.hh"          // 6
#include "G4TwistedTubs.hh"      // 7
#include "G4EllipticalTube.hh"   // 8
#include "G4Cons.hh"             // 9
#include "G4EllipticalCone.hh"   // 10
#include "G4Polycone.hh"         // 11
#include "G4GenericPolycone.hh"  // 12
#include "G4Hype.hh"             // 13
#include "G4Paraboloid.hh"       // 14

// torus and helicoids
#include "G4Torus.hh"            // 1
#include "G4TwistedBox.hh"       // 2
#include "G4TwistedTrd.hh"       // 3
#include "G4TwistedTrap.hh"      // 4

// booleans
#include "G4UnionSolid.hh"        // 1
#include "G4SubtractionSolid.hh"  // 2
#include "G4IntersectionSolid.hh" // 3
#include "G4MultiUnion.hh"        // 4

#include "G4Version.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Timer.hh"

////////////////////////////////////////////////////////////////////////
//
// Mesure perfomance of GetPointOnSurface() for given solid
//
void Check_GetPointOnSurface(G4String text, G4VSolid* solid)
{
  clock_t time = clock();
  G4ThreeVector p;
  for (G4int i = 0; i < 1000000; ++i) { p = solid->GetPointOnSurface(); }
  if (p.mag2() >= 0) time = clock() - time;

  G4cout << text << " "    << std::setw(7) << (G4double)time/CLOCKS_PER_SEC*1000 << " ms"
         << " \tarea = "   << std::setw(7) << solid->GetSurfaceArea()/100.       << " cm^2"
         << " \tvolume = " << std::setw(7) << solid->GetCubicVolume()/1000.      << " cm^3"
         <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
//
// Main: construct solids and mesure the performance of GetPointOnSurface()
//
int main()
{
  // Construct G4Box
  G4double dx = 0., dy = 0., dz = 0.;
  G4VSolid* box = new G4Box("box", dx=50, dy=50, dz=100);

  // Construct G4Para
  G4double alpha = 0., theta = 0., phi = 0.;
  G4VSolid* para = new G4Para("para", dx=50, dy=50, dz=100, alpha=10*deg, theta=20*deg, phi=30*deg);

  // Construct G4Trd
  G4double dx1 = 0., dx2 = 0., dy1 = 0., dy2 = 0.;
  G4VSolid* trd = new G4Trd("trd", dx1=50, dx2=40, dy1=50, dy2=60, dz=100);

  // Construct G4Trap
  G4double dx3 = 0., dx4 = 0.;
  G4VSolid* trap = new G4Trap("trap", dz=100, theta=20*deg, phi=30*deg,
                              dy1=50, dx1=60, dx2=40, 10*deg,
                              dy2=50, dx3=60, dx4=40, 10*deg);

  // Construct G4Tubs
  G4double rmin = 0., rmax = 0., sphi = 0., dphi = 0.;
  G4VSolid* tubs = new G4Tubs("tubs", rmin=50, rmax=100, dz=100, sphi=30*deg, dphi=120*deg);

  // Construct G4Cons
  G4double rmin1 = 0., rmax1 = 0., rmin2 = 0., rmax2 = 0.;
  G4VSolid* cons = new G4Cons("cons", rmin1=50, rmax1=100, rmin2=60, rmax2=110, dz=100, sphi=30*deg, dphi=120*deg);

  // Construct G4CutTubs
  G4ThreeVector n1 = G4ThreeVector(-0.1,-0.1,-1.).unit(), n2 = G4ThreeVector(0.1, 0.1, 1.).unit();
  G4VSolid* cuttubs = new G4CutTubs("cuttubs", rmin=50, rmax=100, dz=100, sphi=30*deg, dphi=120*deg, n1, n2);

  // Construct G4TwistedTubs
  G4double zbot = 0., ztop = 0., twist = 0.;
  G4VSolid* twistedtubs = new G4TwistedTubs("twistedtubs", twist=40*deg, rmin=50, rmax=120, zbot=-80, ztop=120, dphi=60*deg);

  // Construct G4Sphere
  G4double stheta = 0., dtheta = 0.;
  G4VSolid* sphere = new G4Sphere("sphere", rmin=60, rmax=120, sphi=30*deg, dphi=120*deg, stheta=30*deg, dtheta=120*deg);

  // Construct G4Orb
  G4double radius = 0.;
  G4VSolid* orb = new G4Orb("orb", radius=80);

  // Construct G4Torus
  G4double rtor = 0.;
  G4VSolid* torus = new G4Torus("torus", rmin=20, rmax=40, rtor=80, sphi=30*deg, dphi=180*deg);

  // Construct G4Tet
  G4ThreeVector p0(0., 0., 100.), p1(-80.,-90.,-100.), p2(110.,-30.,-110.), p3(0., 100.,-120.);
  G4VSolid* tet = new G4Tet("tet", p0, p1, p2, p3);

  // Tessellated solid
  std::vector<G4ThreeVector> vt(8);
  vt[0] = G4ThreeVector(-40,-80,-110);
  vt[1] = G4ThreeVector(-40, 80,-110);
  vt[2] = G4ThreeVector( 80, 40,-110);
  vt[3] = G4ThreeVector( 80,-40,-110);
  vt[4] = G4ThreeVector(-80,-40, 110);
  vt[5] = G4ThreeVector(-80, 40, 110);
  vt[6] = G4ThreeVector(-20, 20, 110);
  vt[7] = G4ThreeVector(-20,-20, 110);
  G4TessellatedSolid* tessellatedsolid = new G4TessellatedSolid("tessellatedsolid");
  tessellatedsolid->AddFacet(new G4QuadrangularFacet(vt[0], vt[1], vt[2], vt[3], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4QuadrangularFacet(vt[4], vt[7], vt[6], vt[5], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4QuadrangularFacet(vt[1], vt[0], vt[4], vt[5], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4QuadrangularFacet(vt[3], vt[2], vt[6], vt[7], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4TriangularFacet(vt[2], vt[1], vt[5], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4TriangularFacet(vt[2], vt[5], vt[6], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4TriangularFacet(vt[3], vt[7], vt[4], ABSOLUTE));
  tessellatedsolid->AddFacet(new G4TriangularFacet(vt[3], vt[4], vt[0], ABSOLUTE));
  tessellatedsolid->SetSolidClosed(true);

  // Construct G4ExtrudedSolid (hexagon)
  /*
  G4int nsect = 6;
  G4double rhex = 60.;
  G4double angle = 360*deg/nsect;
  std::vector<G4TwoVector> hexagon(nsect);
  for (G4int i = 0; i < nsect; ++i) { hexagon[i].set(rhex*std::cos(i*angle), rhex*std::sin(i*angle)); }
  G4VSolid* extruded = new G4ExtrudedSolid("extruded", hexagon, dz=100., G4TwoVector(0.,0.), 1., G4TwoVector(0.,0.), 1.);
  */

  // Construct G4ExtrudedSolid (half hexagon)
  G4int nsect = 3;
  std::vector<G4TwoVector> hexagon(2*nsect + 2);
  rmax = 80;
  rmin = 40;
  sphi = 0*deg;
  dphi = 180*deg;
  G4double angle = dphi/nsect;
  for (G4int i = 0; i < nsect + 1; ++i)
  {
    hexagon[i].set(rmax*std::cos(sphi + i*angle), rmax*std::sin(sphi + i*angle));
    hexagon[2*nsect + 1 - i].set(rmin*std::cos(sphi + i*angle), rmin*std::sin(sphi + i*angle));
  }
  G4VSolid* extruded = new G4ExtrudedSolid("extruded", hexagon, dz=100., G4TwoVector(0.,0.), 1., G4TwoVector(0.,0.), 1.);

  // Construct G4Polyhedra (half hexagon)
  G4double nrz = 4;
  G4double rr[] = { rmin, rmax, rmax, rmin }, zz[] = { -dz, -dz, dz, dz };
  G4VSolid* polyhedra = new G4Polyhedra("polyhedra", sphi=0*deg, dphi=180*deg, nsect=3, nrz=4, rr, zz);

  // Construct G4Polycone
  G4int nz = 4;
  G4double zplane[] = { -50,   0,   0, 30 };
  G4double rinner[] = {   0,  32,  32, 30 };
  G4double router[] = { 110, 120,  60, 60 };
  G4VSolid* polycone = new G4Polycone("polycone", sphi=60*deg, dphi=240*deg, nz=4, zplane, rinner, router);

  // Construct G4GenericPolycone
  G4double rrgen[] = { 40, 120, 110, 60, 80, 0 }, zzgen[] = { -50, -50, 0, 0, 30, 30 };
  G4VSolid* genericpolycone = new G4GenericPolycone("genericpolycone", sphi=60*deg, dphi=240*deg, nrz=6, rrgen, zzgen);

  // Construct G4GenericTrap
  std::vector<G4TwoVector> gentrap(8);
  gentrap[0].set(-70,-60);
  gentrap[1].set(-70, 60);
  gentrap[2].set( 70, 60);
  gentrap[3].set( 70,-60);
  gentrap[4].set(-40, 50);
  gentrap[5].set( 40, 50);
  gentrap[6].set( 40,-50);
  gentrap[7].set(-40,-50);
  G4VSolid* generictrap = new G4GenericTrap("generictrap", dz=100, gentrap);

  // Construct G4Ellipsoid
  G4double zmin = 0., zmax = 0.;
  G4VSolid* ellipsoid = new G4Ellipsoid("ellipsoid", dx=80, dy=100, dz=120, zmin=-75, zmax=95);

  // Construct G4EllipticalTube
  G4VSolid* ellipticaltube = new G4EllipticalTube("ellipticaltube", dx=50, dy=70, dz=100);

  // Construct G4EllipticalCone
  G4double xa = 0., ya = 0., zh = 0., zcut = 0.;
  G4VSolid* ellipticalcone = new G4EllipticalCone("ellipticalcone", xa=0.3, ya=0.4, zh=200, zcut=100);

  // Construct G4Hype
  G4double rin = 0., rout = 0., stin = 0., stout = 0.;
  G4VSolid* hype = new G4Hype("hype", rin=20, rout=50, stin=20*deg, stout=30*deg, dz=100);

  // Construct G4Paraboloid
  G4VSolid* paraboloid = new G4Paraboloid("paraboloid", dz=100, rmin=30, rmax=80);

  // Construct G4TwistedBox
  G4VSolid* twistedbox = new G4TwistedBox("twistedbox", twist=40*deg, dx=50, dy=50, dz=100);

  // Construct G4TwistedTrd
  G4VSolid* twistedtrd = new G4TwistedTrd("twistedtrd", dx1=50, dx2=40, dy1=50, dy2=60, dz=100, twist=40*deg);

  // Construct G4TwistedTrap
  G4VSolid* twistedtrap = new G4TwistedTrap("twistedtrap",
                                            twist=40*deg,  dz=100, theta=20*deg, phi=30*deg,
                                            dy1=50, dx1=60, dx2=40,
                                            dy2=50, dx3=60, dx4=40, 10*deg);

  // Construct G4UnionSolid
  G4VSolid* unionsolid = new G4UnionSolid("unionsolid", box, orb);

  // Construct G4SubtractionSolid
  G4VSolid* subtractionsolid = new G4SubtractionSolid("subtractionsolid", box, orb);

  // Construct G4IntersectionSolid
  G4VSolid* intersectionsolid = new G4IntersectionSolid("intersectionsolid", box, orb);

  // Construct G4MultiUnion (six balls)
  /*
  G4VSolid* ball = new G4Orb("ball", radius=50);
  G4MultiUnion* multiunion = new G4MultiUnion("multiunion");
  G4double dist = 1.0*radius;
  G4TranslateX3D xpos( dist);
  G4TranslateX3D xneg(-dist);
  G4TranslateY3D ypos( dist);
  G4TranslateY3D yneg(-dist);
  G4TranslateZ3D zpos( dist);
  G4TranslateZ3D zneg(-dist);
  multiunion->AddNode(*ball, xpos);
  multiunion->AddNode(*ball, xneg);
  multiunion->AddNode(*ball, ypos);
  multiunion->AddNode(*ball, yneg);
  multiunion->AddNode(*ball, zpos);
  multiunion->AddNode(*ball, zneg);
  multiunion->Voxelize();
  */

  // Construct G4MultiUnion (box + orb)
  G4MultiUnion* multiunion = new G4MultiUnion("multiunion");
  G4Transform3D identity;
  multiunion->AddNode(*box, identity);
  multiunion->AddNode(*orb, identity);
  multiunion->Voxelize();

  // Special: construct Disk with 21 holes
  G4VSolid* hole = new G4Tubs("hole", rmin=0, rmax=15, dz=21, sphi=0*deg, dphi=360*deg);
  G4VSolid* disk = new G4Tubs("disk", rmin=0, rmax=100, dz=20, sphi=0*deg, dphi=360*deg);
  for (G4int i = 0; i < 12; ++i)
  {
    disk = new G4SubtractionSolid("disk", disk, hole, G4RotateZ3D(i*30*deg)*G4TranslateX3D(rmax - 20.));
  }
  for (G4int i = 0; i < 8; ++i)
  {
    disk = new G4SubtractionSolid("disk", disk, hole, G4RotateZ3D(i*45*deg)*G4TranslateX3D(rmax - 55.));
  }
  disk = new G4SubtractionSolid("disk", disk, hole, G4Transform3D());

  // Mesure the performance and print out the result
  std::cout << std::endl;
  std::cout << "   " << std::string(G4Version.begin() + 7, G4Version.end() - 1) << G4Date << std::endl;
  std::cout << "   Sampling 1.000.000 random points" << std::endl;
  std::cout << "     === CSG solids ===" << std::endl;
  Check_GetPointOnSurface(" 1. G4Box               :", box);
  Check_GetPointOnSurface(" 2. G4Para              :", para);
  Check_GetPointOnSurface(" 3. G4Trd               :", trd);
  Check_GetPointOnSurface(" 4. G4Trap              :", trap);
  Check_GetPointOnSurface(" 5. G4Tubs              :", tubs);
  Check_GetPointOnSurface(" 6. G4CutTubs           :", cuttubs);
  Check_GetPointOnSurface(" 7. G4Cons              :", cons);
  Check_GetPointOnSurface(" 8. G4Orb               :", orb);
  Check_GetPointOnSurface(" 9. G4Sphere            :", sphere);
  Check_GetPointOnSurface("10. G4Torus             :", torus);
  std::cout << "     === Specific solids ===" << std::endl;
  Check_GetPointOnSurface("11. G4Tet               :", tet);
  Check_GetPointOnSurface("12. G4TessellatedSolid  :", tessellatedsolid);
  Check_GetPointOnSurface("13. G4ExtrudedSolid     :", extruded);
  Check_GetPointOnSurface("14. G4Polyhedra         :", polyhedra);
  Check_GetPointOnSurface("15. G4Polycone          :", polycone);
  Check_GetPointOnSurface("16. G4GenericPolycone   :", genericpolycone);
  Check_GetPointOnSurface("17. G4GenericTrap       :", generictrap);
  Check_GetPointOnSurface("18. G4Ellipsoid         :", ellipsoid);
  Check_GetPointOnSurface("19. G4EllipticalTube    :", ellipticaltube);
  Check_GetPointOnSurface("20. G4EllipticalCone    :", ellipticalcone);
  Check_GetPointOnSurface("21. G4Hype              :", hype);
  Check_GetPointOnSurface("22. G4Paraboloid        :", paraboloid);
  Check_GetPointOnSurface("23. G4TwistedTubs       :", twistedtubs);
  Check_GetPointOnSurface("24. G4TwistedBox        :", twistedbox);
  Check_GetPointOnSurface("25. G4TwistedTrd        :", twistedtrd);
  Check_GetPointOnSurface("26. G4TwistedTrap       :", twistedtrap);
  std::cout << "     === Boolean solids ===" << std::endl;
  Check_GetPointOnSurface("27. G4MuiltiUnion       :", multiunion);
  Check_GetPointOnSurface("28. G4UnionSolid        :", unionsolid);
  Check_GetPointOnSurface("29. G4SubtractionSolid  :", subtractionsolid);
  Check_GetPointOnSurface("30. G4IntersectionSolid :", intersectionsolid);
  std::cout << "     === Composite Boolean solid ===" << std::endl;
  Check_GetPointOnSurface("Disk with 21 holes      :", disk);
  std::cout << std::endl;

  return 0;
}
