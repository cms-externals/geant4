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
 
// testG4Trd
//             Ensure asserts are compiled in

#undef NDEBUG
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

///////////////////////////////////////////////////////////////////
//
// Dave's auxiliary function

const G4String OutputInside(const EInside a)
{
  switch(a) 
  {
    case kInside : return "kInside"  ; 
    case kOutside: return "kOutside" ;
    case kSurface: return "kSurface" ;
  }
  return "????";
}


G4bool testG4Trd()
{
  EInside inside ;
  G4ThreeVector pzero(0,0,0);
  G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40);
  G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
  G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40);

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxy(-1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/std::sqrt(2.0),-1/std::sqrt(2.0),0);
    G4ThreeVector vxmy(1/std::sqrt(2.0),-1/std::sqrt(2.0),0);

    G4double Dist, vol, volCheck;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4Trd trd1("Test Box #1",20,20,30,30,40);
    G4Trd trd2("Test Trd",10,30,20,40,40);
  G4Trd trd3("BABAR Trd",0.14999999999999999,0.14999999999999999, 
                           24.707000000000001, 24.707000000000001, 
	                   22.699999999999999) ;


// Check name
    assert(trd1.GetName()=="Test Box #1");
    assert(trd2.GetName()=="Test Trd");

// check cubic volume

    vol = trd1.GetCubicVolume();
    volCheck = 8*20*30*40;
    assert(ApproxEqual(vol,volCheck));

// Check Inside

    assert(trd1.Inside(pzero)==kInside);
    assert(trd1.Inside(pbigz)==kOutside);
    assert(trd1.Inside(ponxside)==kSurface);
    assert(trd1.Inside(ponyside)==kSurface);
    assert(trd1.Inside(ponzside)==kSurface);

    inside = trd1.Inside(G4ThreeVector(20,30,40)) ;
    //  G4cout << "trd1.Inside((20,30,40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(-20,30,40)) ;
    // G4cout << "trd1.Inside((-20,30,40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(20,-30,40)) ;
    //  G4cout << "trd1.Inside((20,-30,40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(20,30,-40)) ;
    // G4cout << "trd1.Inside((20,30,-40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(20,30,0)) ;
    // G4cout << "trd1.Inside((20,30,0)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(0,30,40)) ;
    // G4cout << "trd1.Inside((0,30,40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(20,0,40)) ;
    // G4cout << "trd1.Inside((20,0,40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);

    inside = trd1.Inside(G4ThreeVector(-20,-30,-40)) ;
    // G4cout << "trd1.Inside((-20,-30,-40)) = " << OutputInside(inside) << G4endl ;
    assert(inside == kSurface);





    assert(trd2.Inside(pzero)==kInside);
    assert(trd2.Inside(pbigz)==kOutside);
    assert(trd2.Inside(ponxside)==kSurface);
    assert(trd2.Inside(ponyside)==kSurface);
    assert(trd2.Inside(ponzside)==kSurface);

// Check Surface Normal
    G4ThreeVector normal;

    normal=trd1.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));
    normal=trd1.SurfaceNormal(ponmxside);
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));
    normal=trd1.SurfaceNormal(ponyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));
    normal=trd1.SurfaceNormal(ponmyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));
    normal=trd1.SurfaceNormal(ponzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=trd1.SurfaceNormal(ponmzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));
    normal=trd1.SurfaceNormal(ponzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=trd1.SurfaceNormal(ponmzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));


    // Normals on Edges

    G4ThreeVector edgeXY(    20.0,  30., 0.0); 
    G4ThreeVector edgemXmY( -20.0, -30., 0.0); 
    G4ThreeVector edgeXmY(   20.0, -30., 0.0); 
    G4ThreeVector edgemXY(  -20.0,  30., 0.0); 
    G4ThreeVector edgeXZ(    20.0, 0.0, 40.0); 
    G4ThreeVector edgemXmZ( -20.0, 0.0, -40.0); 
    G4ThreeVector edgeXmZ(   20.0, 0.0, -40.0); 
    G4ThreeVector edgemXZ(  -20.0, 0.0, 40.0); 
    G4ThreeVector edgeYZ(    0.0,  30.0,  40.0); 
    G4ThreeVector edgemYmZ(  0.0, -30.0, -40.0); 
    G4ThreeVector edgeYmZ(   0.0,  30.0, -40.0); 
    G4ThreeVector edgemYZ(   0.0, -30.0,  40.0); 

    G4double invSqrt2 = 1.0 / std::sqrt( 2.0); 
    G4double invSqrt3 = 1.0 / std::sqrt( 3.0); 

    normal= trd1.SurfaceNormal( edgeXY ); 
    assert(ApproxEqual( normal, G4ThreeVector( invSqrt2, invSqrt2, 0.0) )); 

    // G4cout << " Normal at " << edgeXY << " is " << normal 
    //    << " Expected is " << G4ThreeVector( invSqrt2, invSqrt2, 0.0) << G4endl;     

    normal= trd1.SurfaceNormal( edgemXmY ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt2, -invSqrt2, 0.0) )); 
    normal= trd1.SurfaceNormal( edgeXmY ); 
    assert(ApproxEqual( normal, G4ThreeVector( invSqrt2, -invSqrt2, 0.0) )); 
    normal= trd1.SurfaceNormal( edgemXY ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt2, invSqrt2, 0.0) )); 

    normal= trd1.SurfaceNormal( edgeXZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt2, 0.0, invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgemXmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt2, 0.0, -invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgeXmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt2, 0.0, -invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgemXZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt2, 0.0, invSqrt2) )); 

    normal= trd1.SurfaceNormal( edgeYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( 0.0,  invSqrt2,  invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgemYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( 0.0, -invSqrt2, -invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgeYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( 0.0,  invSqrt2, -invSqrt2) )); 
    normal= trd1.SurfaceNormal( edgemYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( 0.0, -invSqrt2,  invSqrt2) )); 

    // Normals on corners

    G4ThreeVector cornerXYZ(    20.0,  30., 40.0); 
    G4ThreeVector cornermXYZ(  -20.0,  30., 40.0); 
    G4ThreeVector cornerXmYZ(   20.0, -30., 40.0); 
    G4ThreeVector cornermXmYZ( -20.0, -30., 40.0); 
    G4ThreeVector cornerXYmZ(    20.0,  30., -40.0); 
    G4ThreeVector cornermXYmZ(  -20.0,  30., -40.0); 
    G4ThreeVector cornerXmYmZ(   20.0, -30., -40.0); 
    G4ThreeVector cornermXmYmZ( -20.0, -30., -40.0); 
 
    normal= trd1.SurfaceNormal( cornerXYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt3,  invSqrt3, invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornermXYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt3,  invSqrt3, invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornerXmYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt3, -invSqrt3, invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornermXmYZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt3, -invSqrt3, invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornerXYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt3,  invSqrt3, -invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornermXYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt3,  invSqrt3, -invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornerXmYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector(  invSqrt3, -invSqrt3, -invSqrt3) )); 
    normal= trd1.SurfaceNormal( cornermXmYmZ ); 
    assert(ApproxEqual( normal, G4ThreeVector( -invSqrt3, -invSqrt3, -invSqrt3) )); 






    double cosa = 4/std::sqrt(17.), sina = 1/std::sqrt(17.), tanga = 1.0/4.0 ;
    
    normal=trd2.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,G4ThreeVector(cosa,0,-sina)));
    normal=trd2.SurfaceNormal(ponmxside);
    assert(ApproxEqual(normal,G4ThreeVector(-cosa,0,-sina)));
    normal=trd2.SurfaceNormal(ponyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,cosa,-sina)));
    normal=trd2.SurfaceNormal(ponmyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,-cosa,-sina)));
    normal=trd2.SurfaceNormal(ponzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=trd2.SurfaceNormal(ponmzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));
    normal=trd2.SurfaceNormal(ponzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=trd2.SurfaceNormal(ponmzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1))); // (0,cosa,-sina) ?

// DistanceToOut(P)

    Dist=trd1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,20));
    Dist=trd1.DistanceToOut(vx);
    assert(ApproxEqual(Dist,19));
    Dist=trd1.DistanceToOut(vy);
    assert(ApproxEqual(Dist,20));
    Dist=trd1.DistanceToOut(vz);
    assert(ApproxEqual(Dist,20));

    Dist=trd2.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,20*cosa));
    Dist=trd2.DistanceToOut(vx);
    assert(ApproxEqual(Dist,19*cosa));
    Dist=trd2.DistanceToOut(vy);
    assert(ApproxEqual(Dist,20*cosa));
    Dist=trd2.DistanceToOut(vz);
    assert(ApproxEqual(Dist,20*cosa+sina));


// DistanceToOut(P,V)

    Dist=trd1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(norm,vmx)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,vy)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,vmy)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vmz)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,std::sqrt(800.))&&*pgoodNorm);

    Dist=trd1.DistanceToOut(ponxside,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(ponmxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmx)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(ponyside,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vy)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(ponmyside,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmy)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(ponzside,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=trd1.DistanceToOut(ponmzside,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

    Dist=trd2.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(*pNorm,G4ThreeVector(cosa,0,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(norm,G4ThreeVector(-cosa,0,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,G4ThreeVector(0,cosa,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,G4ThreeVector(0,-cosa,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vmz)&&*pgoodNorm);
    Dist=trd2.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,std::sqrt(800.))&&*pgoodNorm);

    Dist=trd2.DistanceToOut(ponxside,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(*pNorm,G4ThreeVector(cosa,0,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(ponmxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,G4ThreeVector(-cosa,0,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(ponyside,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,G4ThreeVector(0,cosa,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(ponmyside,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,G4ThreeVector(0,-cosa,-sina))&&*pgoodNorm);
    Dist=trd2.DistanceToOut(ponzside,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=trd2.DistanceToOut(ponmzside,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmz)&&*pgoodNorm);


//DistanceToIn(P)

    Dist=trd1.DistanceToIn(pbigx);
    assert(ApproxEqual(Dist,80));
    Dist=trd1.DistanceToIn(pbigmx);
    assert(ApproxEqual(Dist,80));
    Dist=trd1.DistanceToIn(pbigy);
    assert(ApproxEqual(Dist,70));
    Dist=trd1.DistanceToIn(pbigmy);
    assert(ApproxEqual(Dist,70));
    Dist=trd1.DistanceToIn(pbigz);
    assert(ApproxEqual(Dist,60));
    Dist=trd1.DistanceToIn(pbigmz);
    assert(ApproxEqual(Dist,60));

    Dist=trd2.DistanceToIn(pbigx);
    assert(ApproxEqual(Dist,80*cosa));
    Dist=trd2.DistanceToIn(pbigmx);
    assert(ApproxEqual(Dist,80*cosa));
    Dist=trd2.DistanceToIn(pbigy);
    assert(ApproxEqual(Dist,70*cosa));
    Dist=trd2.DistanceToIn(pbigmy);
    assert(ApproxEqual(Dist,70*cosa));
    Dist=trd2.DistanceToIn(pbigz);
    assert(ApproxEqual(Dist,60));
    Dist=trd2.DistanceToIn(pbigmz);
    assert(ApproxEqual(Dist,60));


// DistanceToIn(P,V)

    Dist=trd1.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(Dist,80));
    Dist=trd1.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(Dist,80));
    Dist=trd1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,70));
    Dist=trd1.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(Dist,70));
    Dist=trd1.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(Dist,60));
    Dist=trd1.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,60));
    Dist=trd1.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=trd1.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(Dist,kInfinity));

    Dist=trd2.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(Dist,80));
    Dist=trd2.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(Dist,80));
    Dist=trd2.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,70));
    Dist=trd2.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(Dist,70));
    Dist=trd2.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(Dist,60));
    Dist=trd2.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,60));
    Dist=trd2.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=trd2.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(Dist,kInfinity));

    Dist=trd3.DistanceToIn(G4ThreeVector(  0.15000000000000185,
                                         -22.048743592955137,
                                           2.4268539333219472),
                           G4ThreeVector(-0.76165597579890043,
                                          0.64364445891356026,
                                         -0.074515708658524193)) ;

    //    G4cout<<"BABAR trd distance = "<<Dist<<G4endl ;
    assert(ApproxEqual(Dist,0.0));

   // return-value = 2.4415531753644804e-15




// CalculateExtent

    G4VoxelLimits limit;		// Unlimited
    G4AffineTransform origin(pzero);
    G4double min,max;

    assert(trd1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));
    assert(trd1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd2.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));
    assert(trd2.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));


    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);

    assert(trd1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-120)&&ApproxEqual(max,-80));   
    assert(trd1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-140)&&ApproxEqual(max,-80));
    assert(trd1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-160)&&ApproxEqual(max,-80));

    assert(trd2.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-130)&&ApproxEqual(max,-70));   
    assert(trd2.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-150)&&ApproxEqual(max,-70));
    assert(trd2.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-160)&&ApproxEqual(max,-80));


    G4RotationMatrix r90Z;
    r90Z.rotateZ(halfpi);
    G4AffineTransform tRotZ(r90Z,pzero);

    assert(trd1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));
    assert(trd1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));
    assert(trd2.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd2.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));


// Check that clipped away

    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-50);

    assert(!trd1.CalculateExtent(kXAxis,xClip,origin,min,max));

    assert(!trd2.CalculateExtent(kXAxis,xClip,origin,min,max));


// Assert clipped to volume

    G4VoxelLimits allClip;
    allClip.AddLimit(kXAxis,-5,+5);
    allClip.AddLimit(kYAxis,-5,+5);
    allClip.AddLimit(kZAxis,-5,+5);
    
    G4RotationMatrix genRot;
    genRot.rotateX(pi/6);
    genRot.rotateY(pi/6);
    genRot.rotateZ(pi/6);
    
    G4AffineTransform tGen(genRot,vx);

    assert(trd1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(trd2.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd2.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd2.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));


    G4VoxelLimits buggyClip2;
    buggyClip2.AddLimit(kXAxis,5,15);

    assert(trd1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd2.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));
    assert(trd2.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));


    buggyClip2.AddLimit(kYAxis,5,15);

    assert(trd1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd2.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(trd2.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));


    G4VoxelLimits buggyClip1;
    buggyClip1.AddLimit(kXAxis,-5,+5);

    assert(trd1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd2.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));
    assert(trd2.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));


    buggyClip1.AddLimit(kYAxis,-5,+5);

    assert(trd1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(trd2.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd2.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(trd2.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    
    G4VoxelLimits newvoxlim;
    newvoxlim.AddLimit(kZAxis,-5,+5);
    
    assert(trd1.CalculateExtent(kXAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));
    assert(trd1.CalculateExtent(kYAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(trd1.CalculateExtent(kZAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(trd2.CalculateExtent(kXAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-(20+5*tanga))&&ApproxEqual(max,20+5*tanga));
    assert(trd2.CalculateExtent(kYAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-(30+5*tanga))&&ApproxEqual(max,30+5*tanga));
    assert(trd2.CalculateExtent(kZAxis,newvoxlim,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    G4VoxelLimits nonsymvox;
    nonsymvox.AddLimit(kZAxis,5,15);
    
    assert(trd2.CalculateExtent(kXAxis,nonsymvox,origin,min,max));
    assert(ApproxEqual(min,-(20+15*tanga))&&ApproxEqual(max,20+15*tanga));
    assert(trd2.CalculateExtent(kYAxis,nonsymvox,origin,min,max));
    assert(ApproxEqual(min,-(30+15*tanga))&&ApproxEqual(max,30+15*tanga));
    assert(trd2.CalculateExtent(kZAxis,nonsymvox,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Trd());
    return 0;
}





