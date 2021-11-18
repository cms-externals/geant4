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
// --------------------------------------------------------------------
// by Evgueni Tcherniaev to compare different MC algorithms for
// surface area estimation

#include <assert.h>
#include "Randomize.hh"
#include "G4VoxelLimits.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4ScaledSolid.hh"
#include "G4EllipticalCone.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Timer.hh"

using namespace CLHEP;

G4double error(G4double real, G4double expected)
{
  return (G4int(10000.*(real - expected)/expected)) / 100.;
};

void BenchEstimateArea(const G4VSolid* solid, G4double ell, G4double expected);
G4double EstimateArea_Safety(const G4VSolid* solid, G4int nStat, G4double ell);
G4double EstimateArea_Distance(const G4VSolid* solid, G4int nStat, G4double ell);
void TestScaledTrd(G4double ell);
void TestScaledCone(G4double ell);
void TestUnionSolid(G4double ell);
void TestSubtractionSolid(G4double ell);
void TestIntersectionSolid(G4double ell);

/////////////////////////////////////////////////////////////////////////////////
//
void TestScaledTrd(G4double ell)
{
  G4VSolid* trd1 = new G4Trd("Trd",
                             20*cm, 10*cm, // Dx1, Dx2
                             20*cm, 10*cm, // Dy1, Dy2
                             100*cm);      // Dz
  G4VSolid* trd2 = new G4Trd("Trd",
                             60*cm, 30*cm, // Dx1, Dx2
                             40*cm, 20*cm, // Dy1, Dy2
                             20*cm);       // Dz

  G4double sx = 3, sy = 2, sz = .2;
  G4ScaledSolid* solid = new G4ScaledSolid("Scaled Trd", trd1, G4Scale3D(sx,sy,sz));

  G4double area = trd2->GetSurfaceArea();
  std::cout << "\n***"
            << "\n*** Test G4Scaled Trd - Scale factors (" << sx << ", " << sy << ", " << sz << ")"
            << " *** Surface area : " << area
            << "\n***" <<std::endl;
  std::cout << "G4ScaledSolid::GetSurfaceArea() = " << solid->GetSurfaceArea() << std::endl;

  BenchEstimateArea(solid, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestScaledCone(G4double ell)
{
  G4VSolid* cone = new G4Cons("Cone",
                              0, 5*cm, // Rmax at -Dz
                              0, 1*cm, // Rmax at +Dz
                              5*cm,    // Dz
                              0, 360*deg);
  G4VSolid* econe = new G4EllipticalCone("Elliptical cone",
                                         2,      // xSemiAxis
                                         1,      // ySemiAxis
                                         15*cm,  // zheight
		                         10*cm); // zTopCut

  G4double sx = 10, sy = 5, sz = 2;
  G4ScaledSolid* solid = new G4ScaledSolid("Scaled cone", cone, G4Scale3D(sx,sy,sz));

  G4double area = econe->GetSurfaceArea();
  std::cout << "\n***"
            << "\n*** Test G4Scaled Cone - Scale factors (" << sx << ", " << sy << ", " << sz << ")"
            << " *** Surface area : " << area
            << "\n***" <<std::endl;
  std::cout << "G4ScaledSolid::GetSurfaceArea() = " << solid->GetSurfaceArea() << std::endl;

  BenchEstimateArea(solid, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestUnionSolid(G4double ell)
{
  G4double rbig = 10  *cm, zbig = 10*cm;
  G4double rsml = 2.5 *cm, zsml = 17*cm;

  G4VSolid* BigTube = new G4Tubs("Big Tube",   0, rbig, zbig, 0*deg, 360*deg);
  G4VSolid* SmlTube = new G4Tubs("Small Tube", 0, rsml, zsml, 0*deg, 360*deg);

  G4double sbig  = 2*zbig*twopi*rbig + 2*pi*rbig*rbig;
  assert(sbig == BigTube->GetSurfaceArea());
  G4double shole = 2*(zsml-zbig)*twopi*rsml;

  G4double dd = 6*cm;
  G4VSolid* body1 = new G4UnionSolid("Union1", BigTube, SmlTube, G4Translate3D( 0., 0., 0.));
  G4VSolid* body2 = new G4UnionSolid("Union2", body1,   SmlTube, G4Translate3D( dd, 0., 0.));
  G4VSolid* body3 = new G4UnionSolid("Union3", body2,   SmlTube, G4Translate3D(-dd, 0., 0.));
  G4VSolid* body4 = new G4UnionSolid("Union4", body3,   SmlTube, G4Translate3D( 0., dd, 0.));
  G4VSolid* body5 = new G4UnionSolid("Union5", body4,   SmlTube, G4Translate3D( 0.,-dd, 0.));

  G4double area = sbig + 5*shole;
  std::cout << "\n***"
            << "\n*** Test G4UnionSolid - Disk plus 5 cylinders *** Surface area : " << area
            << "\n***" << std::endl;
  std::cout << "Disk area         : " << sbig << std::endl;
  std::cout << "Increment area    : " << shole << std::endl;
  std::cout << "Real surface area : " << area << std::endl;

  BenchEstimateArea(body5, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestSubtractionSolid(G4double ell)
{
  G4double rbig = 10  *cm, zbig = 10*cm;
  G4double rsml = 2.5 *cm, zsml = 12*cm;

  G4VSolid* BigTube = new G4Tubs("Big Tube",   0, rbig, zbig, 0*deg, 360*deg);
  G4VSolid* SmlTube = new G4Tubs("Small Tube", 0, rsml, zsml, 0*deg, 360*deg);

  G4double sbig  = 2*zbig*twopi*rbig + 2*pi*rbig*rbig;
  assert(sbig == BigTube->GetSurfaceArea());
  G4double shole = 2*zbig*twopi*rsml - 2*pi*rsml*rsml;

  G4double dd = 6*cm;
  G4VSolid* body1 = new G4SubtractionSolid("SubSolid1", BigTube, SmlTube, G4Translate3D( 0., 0., 0.));
  G4VSolid* body2 = new G4SubtractionSolid("SubSolid2", body1,   SmlTube, G4Translate3D( dd, 0., 0.));
  G4VSolid* body3 = new G4SubtractionSolid("SubSolid3", body2,   SmlTube, G4Translate3D(-dd, 0., 0.));
  G4VSolid* body4 = new G4SubtractionSolid("SubSolid4", body3,   SmlTube, G4Translate3D( 0., dd, 0.));
  G4VSolid* body5 = new G4SubtractionSolid("SubSolid5", body4,   SmlTube, G4Translate3D( 0.,-dd, 0.));

  G4double area = sbig + 5*shole;
  std::cout << "\n***"
            << "\n*** Test G4SubtractionSolid - Disk with 5 holes *** Surface area : " << area
            << "\n***" << std::endl;
  std::cout << "Disk area         : " << sbig << std::endl;
  std::cout << "Increment area    : " << shole << std::endl;
  std::cout << "Real surface area : " << area << std::endl;

  BenchEstimateArea(body5, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestSubtractionNull(G4double ell)
{
  G4double dbig = 10*cm;
  G4double dsml =  1*cm;

  G4VSolid* BigBox = new G4Box("Big Box",  dbig, dbig, dbig);
  G4VSolid* SmlBox = new G4Box("Small Box",  dsml, dsml, dsml);

  G4VSolid* solid = new G4SubtractionSolid("Null", SmlBox, BigBox, G4Translate3D( 0., 0., 0.));
  
  G4double area  = 0;
  std::cout << "\n***"
            << "\n*** Test G4SubtractionSolid - Null object *** Surface area: " << area
            << "\n***" << std::endl;
  std::cout << "Big box area      : " << BigBox->GetSurfaceArea() << std::endl;
  std::cout << "Small box area    : " << SmlBox->GetSurfaceArea() << std::endl;

  BenchEstimateArea(solid, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestSubtractionDisk(G4double ell)
{
  static const G4double TopPMTArrayXY[247][2] =
       {{ 726.44,       -47.61 },
        { 714.01,       -142.03 },
        { 689.37,       -234.01 },
        { 652.92,       -321.99 },
        { 605.31,       -404.45 },
        { 547.34,       -480 },
        { 480,  -547.34 },
        { 404.46,       -605.31 },
        { 321.99,       -652.92 },
        { 234.01,       -689.36 },
        { 142.03,       -714.01 },
        { 47.61,        -726.44 },
        { -47.61,       -726.44 },
        { -142.03,      -714.01 },
        { -234.01,      -689.36 },
        { -321.99,      -652.92 },
        { -404.45,      -605.31 },
        { -480, -547.34 },
        { -547.34,      -480 },
        { -605.31,      -404.45 },
        { -652.92,      -321.99 },
        { -689.36,      -234.01 },
        { -714.01,      -142.03 },
        { -726.44,      -47.61 },
        { -726.44,      47.61 },
        { -714.01,      142.03 },
        { -689.36,      234.01 },
        { -652.92,      321.99 },
        { -605.31,      404.46 },
        { -547.34,      480 },
        { -480, 547.34 },
        { -404.45,      605.31 },
        { -321.99,      652.92 },
        { -234.01,      689.37 },
        { -142.03,      714.01 },
        { -47.61,       726.44 },
        { 47.61,        726.44 },
        { 142.03,       714.01 },
        { 234.01,       689.37 },
        { 321.99,       652.92 },
        { 404.46,       605.31 },
        { 480,  547.34 },
        { 547.34,       480 },
        { 605.31,       404.46 },
        { 652.92,       321.99 },
        { 689.37,       234.01 },
        { 714.01,       142.03 },
        { 726.44,       47.61 },
        { 624.61,       -82.23 },
        { 383.52,       -499.81 },
        { 241.09,       -582.04 },
        { -241.09,      -582.04 },
        { -383.52,      -499.81 },
        { -624.61,      -82.23 },
        { -624.61,      82.23 },
        { -383.52,      499.81 },
        { -241.09,      582.04 },
        { 241.09,       582.04 },
        { 383.52,       499.81 },
        { 624.61,       82.23 },
        { 600.52,       -248.74 },
        { 515.68,       -395.69 },
        { 325,  -562.92 },
        { 84.84,        -644.44 },
        { -84.84,       -644.44 },
        { -325, -562.92 },
        { -515.68,      -395.69 },
        { -600.52,      -248.74 },
        { -650, 0 },
        { -600.52,      248.74 },
        { -515.68,      395.7 },
        { -325, 562.92 },
        { -84.84,       644.44 },
        { 84.84,        644.44 },
        { 325,  562.92 },
        { 515.68,       395.7 },
        { 600.52,       248.74 },
        { 650,  0 },
        { 0,    644.32 },
        { 46.5, 563.78 },
        { -46.5,        563.78 },
        { 279,  483.24 },
        { 186,  483.24 },
        { 93,   483.24 },
        { 0,    483.24 },
        { -93,  483.24 },
        { -186, 483.24 },
        { -279, 483.24 },
        { 325.5,        402.7 },
        { 232.5,        402.7 },
        { 139.5,        402.7 },
        { 46.5, 402.7 },
        { -46.5,        402.7 },
        { -139.5,       402.7 },
        { -232.5,       402.7 },
        { -325.5,       402.7 },
        { 558,  322.16 },
        { 465,  322.16 },
        { 372,  322.16 },
        { 279,  322.16 },
        { 186,  322.16 },
        { 93,   322.16 },
        { 0,    322.16 },
        { -93,  322.16 },
        { -186, 322.16 },
        { -279, 322.16 },
        { -372, 322.16 },
        { -465, 322.16 },
        { -558, 322.16 },
        { 511.5,        241.62 },
        { 418.5,        241.62 },
        { 325.5,        241.62 },
        { 232.5,        241.62 },
        { 139.5,        241.62 },
        { 46.5, 241.62 },
        { -46.5,        241.62 },
        { -139.5,       241.62 },
        { -232.5,       241.62 },
        { -325.5,       241.62 },
        { -418.5,       241.62 },
        { -511.5,       241.62 },
        { 465,  161.08 },
        { 372,  161.08 },
        { 279,  161.08 },
        { 186,  161.08 },
        { 93,   161.08 },
        { 0,    161.08 },
        { -93,  161.08 },
        { -186, 161.08 },
        { -279, 161.08 },
        { -372, 161.08 },
        { -465, 161.08 },
        { 511.5,        80.54 },
        { 418.5,        80.54 },
        { 325.5,        80.54 },
        { 232.5,        80.54 },
        { 139.5,        80.54 },
        { 46.5, 80.54 },
        { -46.5,        80.54 },
        { -139.5,       80.54 },
        { -232.5,       80.54 },
        { -325.5,       80.54 },
        { -418.5,       80.54 },
        { -511.5,       80.54 },
        { 558,  0 },
        { 465,  0 },
        { 372,  0 },
        { 279,  0 },
        { 186,  0 },
        { 93,   0 },
        { 0,    0 },
        { -93,  0 },
        { -186, 0 },
        { -279, 0 },
        { -372, 0 },
        { -465, 0 },
        { -558, 0 },
        { 511.5,        -80.54 },
        { 418.5,        -80.54 },
        { 325.5,        -80.54 },
        { 232.5,        -80.54 },
        { 139.5,        -80.54 },
        { 46.5, -80.54 },
        { -46.5,        -80.54 },
        { -139.5,       -80.54 },
        { -232.5,       -80.54 },
        { -325.5,       -80.54 },
        { -418.5,       -80.54 },
        { -511.5,       -80.54 },
        { 465,  -161.08 },
        { 372,  -161.08 },
        { 279,  -161.08 },
        { 186,  -161.08 },
        { 93,   -161.08 },
        { 0,    -161.08 },
        { -93,  -161.08 },
        { -186, -161.08 },
        { -279, -161.08 },
        { -372, -161.08 },
        { -465, -161.08 },
        { 511.5,        -241.62 },
        { 418.5,        -241.62 },
        { 325.5,        -241.62 },
        { 232.5,        -241.62 },
        { 139.5,        -241.62 },
        { 46.5, -241.62 },
        { -46.5,        -241.62 },
        { -139.5,       -241.62 },
        { -232.5,       -241.62 },
        { -325.5,       -241.62 },
        { -418.5,       -241.62 },
        { -511.5,       -241.62 },
        { 558,  -322.16 },
        { 465,  -322.16 },
        { 372,  -322.16 },
        { 279,  -322.16 },
        { 186,  -322.16 },
        { 93,   -322.16 },
        { 0,    -322.16 },
        { -93,  -322.16 },
        { -186, -322.16 },
        { -279, -322.16 },
        { -372, -322.16 },
        { -465, -322.16 },
        { -558, -322.16 },
        { 325.5,        -402.7 },
        { 232.5,        -402.7 },
        { 139.5,        -402.7 },
        { 46.5, -402.7 },
        { -46.5,        -402.7 },
        { -139.5,       -402.7 },
        { -232.5,       -402.7 },
        { -325.5,       -402.7 },
        { 279,  -483.24 },
        { 186,  -483.24 },
        { 93,   -483.24 },
        { 0,    -483.24 },
        { -93,  -483.24 },
        { -186, -483.24 },
        { -279, -483.24 },
        { 46.5, -563.78 },
        { -46.5,        -563.78 },
        { 0,    -644.32 },
        { 549.51,       165.51 },
        { 549.51,       -165.51 },
        { -549.51,      165.51 },
        { -549.51,      -165.51 },
        { 634.02,       169.94 },
        { 634.02,       -169.94 },
        { -634.02,      169.94 },
        { -634.02,      -169.94 },
        { 418.09,       393.13 },
        { 418.09,       -393.13 },
        { -418.09,      393.13 },
        { -418.09,      -393.13 },
        { 464.18,       464.11 },
        { 464.18,       -464.11 },
        { -464.18,      464.11 },
        { -464.18,      -464.11 },
        { 131.42,       558.64 },
        { 131.42,       -558.64 },
        { -131.42,      558.64 },
        { -131.42,      -558.64 },
        { 169.84,       634.04 },
        { 169.84,       -634.04 },
        { -169.84,      634.04 },
        { -169.84,      -634.04 }};

  G4int numPMTs = 247;
  G4double ptfeLinerThickness = 2.*mm;
  G4double outerFlatRadius = 77.9*cm;
  G4double innerFlatRadius = outerFlatRadius - 1.*cm;
  G4double ptfeHoleRadius = 3.2*cm;
  //G4double ptfeGaseousHeight = 10*mm;
  //G4double ptfeGaseousConeHeight = 48*mm;
  G4double TopPMTArrayXYScaling = 1.0;

  // Titanium plates and PTFE liners
  //
  //  Create solid plates for holders and liner, as well as the holes that
  //  get punched into them. Note that the holes are twice the thickness
  //  of the plates and liner so that we guarantee sufficient overlap.

  G4Tubs *ptfeLiner_solid1 = new G4Tubs("ptfePlate_solid1", 0,
                                        innerFlatRadius-2.*mm,
                                        ptfeLinerThickness/2, 0.*deg, 360.*deg );

  G4Tubs *ptfeLiner_solid2 = new G4Tubs( "ptfeLiner_solid2", 0,
                                         ptfeHoleRadius,
                                         ptfeLinerThickness, 0.*deg, 360.*deg );

  // A dummy volume to start the array union solid:
  G4Tubs *dummy_solid = new G4Tubs( "dummy_solid", 0,
                                    2*mm, 0.5*mm, 0.*deg, 360.*deg );


  // Add together cylinders in shape of array to subtract from plate
  G4double disp_x = TopPMTArrayXY[0][0] * TopPMTArrayXYScaling; 
  G4double disp_y = TopPMTArrayXY[0][1] * TopPMTArrayXYScaling; 
  G4UnionSolid *ptfe_array_solid = new G4UnionSolid("ptfe_array_solid",
                                                    dummy_solid, ptfeLiner_solid2,
                                                    0, G4ThreeVector(disp_x, disp_y, 0) );

  for( G4int i=1; i<numPMTs; i++ )
  {
    disp_x = TopPMTArrayXY[i][0] * TopPMTArrayXYScaling; 
    disp_y = TopPMTArrayXY[i][1] * TopPMTArrayXYScaling; 
      
    ptfe_array_solid = new G4UnionSolid( "ptfe_array_solid",
                                         ptfe_array_solid, ptfeLiner_solid2,
                                         0, G4ThreeVector(disp_x, disp_y, 0) );
  }

  // Subtract whole array from PTFE plate
  G4SubtractionSolid*
  ptfeLiner_solid3 = new G4SubtractionSolid( "ptfeLiner_solid3",
                                         ptfeLiner_solid1, ptfe_array_solid,
                                         0, G4ThreeVector(0, 0, 0) );
  /*
  // Subtract together cylinders from plate
  G4double disp_x = TopPMTArrayXY[0][0] * TopPMTArrayXYScaling; 
  G4double disp_y = TopPMTArrayXY[0][1] * TopPMTArrayXYScaling; 
  G4SubtractionSolid *ptfeLiner_solid3 = new G4SubtractionSolid("ptfeLiner_solid3",
                                                    ptfeLiner_solid1, ptfeLiner_solid2,
                                                    0, G4ThreeVector(disp_x, disp_y, 0) );

  for( G4int i=1; i<numPMTs; i++ )
  {
    disp_x = TopPMTArrayXY[i][0] * TopPMTArrayXYScaling; 
    disp_y = TopPMTArrayXY[i][1] * TopPMTArrayXYScaling; 
      
    ptfeLiner_solid3 = new G4SubtractionSolid( "ptfeLiner_solid3",
                                         ptfeLiner_solid3, ptfeLiner_solid2,
                                         0, G4ThreeVector(disp_x, disp_y, 0) );
  }
  */
  G4double area =
    twopi*(innerFlatRadius-2.*mm)*(innerFlatRadius-2.*mm) +
    twopi*(innerFlatRadius-2.*mm)*ptfeLinerThickness -
    247*twopi*ptfeHoleRadius*ptfeHoleRadius + 247*twopi*ptfeHoleRadius*ptfeLinerThickness;

  std::cout << "\n***"
            << "\n*** Test G4SubtractionSolid - Disk with 247 holes *** Surface area: " << area
            << "\n***" << std::endl;

  BenchEstimateArea(ptfeLiner_solid3, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
void TestIntersectionSolid(G4double ell)
{
  G4double dx = 10*cm, dy1 = 5*cm, dy2 = 8*cm, dz = 50*cm;
  G4VSolid* box1 = new G4Box("Box1", dx, dy1, dz);
  G4VSolid* box2 = new G4Box("Box2", dx, dy2, dz);

  G4double ang = 60*deg;
  G4VSolid* solid = new G4IntersectionSolid("IntSolid", box1, box2,G4RotateY3D(ang));

  G4double side = 2*dx/std::sin(ang);
  G4double area = 4*dx*side + 8*side*std::min(dy1,dy2);
  std::cout << "\n***"
            << "\n*** Test G4IntersectionSolid - two intersecting bars *** Surface area : " << area
            << "\n***" << std::endl; 
  BenchEstimateArea(solid, ell, area);
}

/////////////////////////////////////////////////////////////////////////////////
//
int main()
{
  G4double ell = -1;
  TestScaledTrd(ell);
  TestScaledCone(ell);
  TestUnionSolid(ell);
  TestSubtractionSolid(ell);
  TestIntersectionSolid(ell);
  TestSubtractionNull(ell);

  // Tests for a Disk with 247 holes consumes quite a bit of time:
  // to run three algorithms with statistics 10^5 & 10^6 points takes ~20 mins,
  // statistics 10^7 points was never tried
  // 
  //TestSubtractionDisk(ell);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Benchmarking different algorithms for estimation of surface area
//
void BenchEstimateArea(const G4VSolid* solid, G4double ell, G4double area)
{
  clock_t t;
  G4double a1e5 = 0;
  G4double a1e6 = 0;
  G4double a1e7 = 0;

  // Test EstimateArea_Safety()
  //
  std::cout << "\n=== Old algorithm based on Safety (replaced in Geant4 10.05)" << std::endl;
  t = clock();
  a1e5 = EstimateArea_Safety(solid,100000,-1);
  a1e6 = EstimateArea_Safety(solid,1000000,-1);
  a1e7 = EstimateArea_Safety(solid,10000000,-1);
  std::cout << "N = 10^5   Estimated Area = " << a1e5 << " (" << error(a1e5,area) << "%) " << std::endl;
  std::cout << "N = 10^6   Estimated Area = " << a1e6 << " (" << error(a1e6,area) << "%) " << std::endl;
  std::cout << "N = 10^7   Estimated Area = " << a1e7 << " (" << error(a1e7,area) << "%) " << std::endl;
  t = clock() - t;
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;

  // Test EstimateArea_Distance()
  //
  std::cout << "\n=== Algorithm based on Distance along SurfaceNormal(p)" << std::endl;
  t = clock();
  a1e5 = EstimateArea_Distance(solid,100000,-1);
  a1e6 = EstimateArea_Distance(solid,1000000,-1);
  a1e7 = EstimateArea_Distance(solid,10000000,-1);
  std::cout << "N = 10^5   Estimated Area = " << a1e5 << " (" << error(a1e5,area) << "%) " << std::endl;
  std::cout << "N = 10^6   Estimated Area = " << a1e6 << " (" << error(a1e6,area) << "%) " << std::endl;
  std::cout << "N = 10^7   Estimated Area = " << a1e7 << " (" << error(a1e7,area) << "%) " << std::endl;
  t = clock() - t;
  G4cout << "   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;

  // Test G4VSolid::EstimateSurfaceArea()
  //
  std::cout << "\n=== G4VSolid::EstimateSurfaceArea()" << std::endl;
  t = clock();
  a1e5 = solid->EstimateSurfaceArea(100000,-1);
  a1e6 = solid->EstimateSurfaceArea(1000000,-1);
  a1e7 = solid->EstimateSurfaceArea(10000000,-1);
  std::cout << "N = 10^5   Estimated Area = " << a1e5 << " (" << error(a1e5,area) << "%) " << std::endl;
  std::cout << "N = 10^6   Estimated Area = " << a1e6 << " (" << error(a1e6,area) << "%) " << std::endl;
  std::cout << "N = 10^7   Estimated Area = " << a1e7 << " (" << error(a1e7,area) << "%) " << std::endl;
  t = clock() - t;
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Estimate surface area using Safety(). Based on the implementation by
// Mikhail Kosov that was in G4VSolid before Geant4 10.05.
// The algorithm does not work properly when safety is underestimated,
// for example, in case of Scaled solids.
//
G4double EstimateArea_Safety(const G4VSolid* solid,
                             G4int nStat, G4double ell)
{
  G4int inside=0;
  G4double px,py,pz,minX,maxX,minY,maxY,minZ,maxZ,surf;
  G4ThreeVector p;
  EInside in;

  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  solid->CalculateExtent(kXAxis,limit,origin,minX,maxX);
  solid->CalculateExtent(kYAxis,limit,origin,minY,maxY);
  solid->CalculateExtent(kZAxis,limit,origin,minZ,maxZ);

  // limits

  if(nStat < 100) { nStat = 100; }

  G4double dX=maxX-minX;
  G4double dY=maxY-minY;
  G4double dZ=maxZ-minZ;
  if(ell<=0.)          // Automatic definition of skin thickness
  {
    G4double minval=dX;
    if(dY<dX) { minval=dY; }
    if(dZ<minval) { minval=dZ; }
    ell=.01*minval;
  }

  G4double dd=2*ell;
  minX-=ell; minY-=ell; minZ-=ell; dX+=dd; dY+=dd; dZ+=dd;

  for(G4int i = 0; i < nStat; i++ )
  {
    px = minX+dX*G4UniformRand();
    py = minY+dY*G4UniformRand();
    pz = minZ+dZ*G4UniformRand();
    p  = G4ThreeVector(px,py,pz);
    in = solid->Inside(p);
    if(in != kOutside)
    {
      if  (solid->DistanceToOut(p)<ell) { inside++; }
    }
    else if(solid->DistanceToIn(p)<ell) { inside++; }
  }
  // @@ The conformal correction can be upgraded
  surf = dX*dY*dZ*inside/dd/nStat;
  return surf;
}

//////////////////////////////////////////////////////////////////////////
//
// Estimation of the surface area based on calculation of a distance to the
// surface along the normal. It works fine provided that SurfaceNormal(p)
// always retuns a correct normal, i.e. normal to the nearest surface.
// Unfortunately at present time (Geant4 10.05) it is not the case. To
// resolve the issue current G4VSolid::EstimateSurfaceArea() uses more
// complicated algorithm, that can be replaced with the code below after
// SurfaceNormal(p) will be fixed in all G4Solids.
//
G4double EstimateArea_Distance(const G4VSolid* solid,
                               G4int nstat, G4double ell)
{
  G4ThreeVector bmin, bmax;
  solid->BoundingLimits(bmin, bmax);

  G4double dX = bmax.x() - bmin.x();
  G4double dY = bmax.y() - bmin.y();
  G4double dZ = bmax.z() - bmin.z();
  
  // Define statistics and skin thickness
  //
  G4int npoints = (nstat < 1000) ? 1000 : nstat;
  G4double eps = (ell > 0) ? ell : 0.01 * std::min(std::min(dX, dY), dZ);

  G4double minX = bmin.x() - eps;
  G4double minY = bmin.y() - eps;
  G4double minZ = bmin.z() - eps;

  G4double dd = 2. * eps;
  dX += dd;
  dY += dd;
  dZ += dd;

  // Calculate surface area
  //
  G4int icount = 0;
  for(G4int i = 0; i < npoints; ++i)
  {
    G4double px = minX + dX*G4UniformRand();
    G4double py = minY + dY*G4UniformRand();
    G4double pz = minZ + dZ*G4UniformRand();
    G4ThreeVector p  = G4ThreeVector(px, py, pz);
    EInside in = solid->Inside(p);
    G4double dist = 0;
    if (in == kInside)
    {
      if (solid->DistanceToOut(p) >= eps) continue;
      dist = solid->DistanceToOut(p, solid->SurfaceNormal(p));
    }
    else if (in == kOutside)
    {
      if (solid->DistanceToIn(p) >= eps) continue;
      dist = solid->DistanceToIn(p, -(solid->SurfaceNormal(p)));
    }
    if (dist < eps) icount++;
  }
  return dX*dY*dZ*icount/npoints/dd;
}
