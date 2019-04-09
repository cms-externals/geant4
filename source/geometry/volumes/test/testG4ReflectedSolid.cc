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
// Test for G4ReflectedSolid class 
//
// 23.07.01 V.Grichine
// 10.05.02 V.Grichine, bug fixed in CalculateExtent: box like algorithm proposed by
//                      J. Apostolakis 
// 




#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4Normal3D.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4Transform3D.hh"
#include "G4VoxelLimits.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4ReflectedSolid.hh"
#include "G4DisplacedSolid.hh"


const G4String OutputInside(const EInside a)
{
	switch(a) 
        {
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}


int main()
{
    G4ThreeVector pzero(0,0,0), p;
    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40),
                   ponb2x(10,0,0),ponb2y(0,10,0),ponb2z(0,0,10),
                   ponb2mx(-10,0,0),ponb2my(0,-10,0),ponb2mz(0,0,-10);
    G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
    G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40),
                  ponb2zy(0,5,10),ponb2mzy(0,5,-10) ;

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxy(-1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/std::sqrt(2.0),-1/std::sqrt(2.0),0);
    G4ThreeVector vxmy(1/std::sqrt(2.0),-1/std::sqrt(2.0),0);
    G4ThreeVector vxmz(1/std::sqrt(2.0),0,-1/std::sqrt(2.0));

    G4double dist;
    G4int i;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4RotationMatrix identity, xRot ;
    G4Transform3D moveX50(identity,G4ThreeVector(50.0, 0.0, 0.0));    
    G4Transform3D moveY50(identity,G4ThreeVector(0.0, 50.0, 0.0)); 
   
// NOTE: xRot = rotation such that x axis->y axis & y axis->-x axis

    xRot.rotateZ(-pi*0.5) ;

    G4ReflectX3D reflectX ;

    G4Point3D newPointN = reflectX*G4Normal3D(1,1,1) ;
    G4cout<<"x component of normal newPoint = "<<newPointN.x()<<G4endl ;

    G4Point3D newPoint = reflectX*G4Point3D(1,1,1) ;
    G4cout<<"x component of point newPoint = "<<newPoint.x()<<G4endl ;

    //  G4Transform3D transfomX = G4ReflectX3D() ;
    //  G4Transform3D transform(xRot,pzero) ;

    G4Box b1("Test Box #1",20,30,40);
    G4Box b2("Test Box #2",10,10,10);
    G4Box b3("Test Box #3",10,50,10);
    G4Box b4("Test Box #4",50,50,50);

    G4DisplacedSolid db4("db4",&b4,moveX50);
    G4ReflectedSolid rdb4("rdb4",&db4,reflectX);

    G4DisplacedSolid dyb4("db4",&b4,moveY50);
    G4ReflectedSolid rdyb4("rdb4",&dyb4,reflectX);



    G4Tubs t1("Solid Tube #1",0,50,50,0,360);
    G4Tubs t2("Hole Tube #2",45,50,50,0,360);
    G4Tubs t3("Solid cutted Tube #3",0,50,50,0,pi/2.0);

    G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*pi),
	   c2("Full Cone",0,50,0,100,50,0,2*pi) ;
    G4Sphere s1("s1",0.0,50.0,0.0,0.5*pi,0.0,pi);    



    G4IntersectionSolid b1Ib2("b1Intersectionb2",&b1,&b2),
                        t1Ib2("t1Intersectionb2",&t1,&b2),
                        c2Ib2("c2Intersectionb2",&c2,&b2) ;

    // passRotT3 should be in 2 octant, while actiRotT3 in 4 one

    G4ReflectedSolid passRotT3("passRotT3",&t3,reflectX) ;
    G4ReflectedSolid rs1("rs1",&s1,reflectX) ;

    //    G4ReflectedSolid actiRotT3("actiRotT3",&t3,transform) ;
    //    G4ReflectedSolid actiRotB1("actiRotB3",&b1,transform) ;

    G4ThreeVector pRmaxPlus(50,1,0) ;
   
    dist = passRotT3.DistanceToIn(pRmaxPlus,vmx) ;
    assert(ApproxEqual(dist,50));
    //  cout<<"passRotT3.DistanceToIn(pRmaxPlus,vmx) = "<<dist<<G4endl ;

    //    dist = actiRotT3.DistanceToIn(pRmaxPlus,vmx) ;
    //   cout<<"actiRotT3.DistanceToIn(pRmaxPlus,vmx) = "<<dist<<G4endl ;

// Check Inside

    assert(passRotT3.Inside(pzero)==kSurface);
    assert(passRotT3.Inside(pbigz)==kOutside);
    assert(passRotT3.Inside(ponb2x)==kOutside);
    assert(passRotT3.Inside(ponb2y)==kSurface);
    assert(passRotT3.Inside(ponb2z)==kSurface);

    //    assert(t3.Inside(pzero)==kSurface);
    //    assert(actiRotT3.Inside(pzero)==kSurface);

    //  assert(actiRotT3.Inside(pbigz)==kOutside);
    //  assert(actiRotT3.Inside(ponb2x)==kSurface);
    //   assert(actiRotT3.Inside(ponb2y)==kOutside);

    //    assert(actiRotT3.Inside(ponb2z)==kSurface);

    EInside in = rdb4.Inside(G4ThreeVector(-50.0,0.0,0.0));
    G4cout<<"rdb4.Inside(G4ThreeVector(-50.0,0.0,0.0) = "<<OutputInside(in)<<G4endl;

    
// Check Surface Normal

    G4ThreeVector normal;

 //   normal=actiRotT3.SurfaceNormal(G4ThreeVector(20,0,0));
 //   assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=passRotT3.SurfaceNormal(G4ThreeVector(-15,0,0));
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));

    normal=passRotT3.SurfaceNormal(G4ThreeVector(0,10,0));
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    // (0,-1,0) is transformed to (0,1,0) ??
    //    normal=actiRotT3.SurfaceNormal(ponb2my);
    //    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

  //  normal=actiRotT3.SurfaceNormal(G4ThreeVector(1,-1,50));
  //  assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

 //   normal=actiRotT3.SurfaceNormal(G4ThreeVector(1,-1,-50));
 //   assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    normal=passRotT3.SurfaceNormal(G4ThreeVector(-15,1,50));
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

    normal=passRotT3.SurfaceNormal(G4ThreeVector(-15,1,-50));
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));


// DistanceToIn(P,V)

    dist=passRotT3.DistanceToIn(G4ThreeVector(100,10,0),vmx);
    assert(ApproxEqual(dist,100));

    //  dist=actiRotT3.DistanceToIn(G4ThreeVector(-100,-10,0),vx);
    //  assert(ApproxEqual(dist,100));

    // dist=actiRotT3.DistanceToIn(G4ThreeVector(10,100,0),vmy);
    //  assert(ApproxEqual(dist,100));

    dist=passRotT3.DistanceToIn(G4ThreeVector(-20,-100,0),vy);
    assert(ApproxEqual(dist,100));

    dist=passRotT3.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(dist,50));

    dist=passRotT3.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(dist,50));

    dist=passRotT3.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    //    dist=actiRotT3.DistanceToIn(pbigmx,vxy);
    //    assert(ApproxEqual(dist,kInfinity));


//DistanceToIn(P)

    //   dist=actiRotT3.DistanceToIn(G4ThreeVector(10,1,0));
    //  assert(ApproxEqual(dist,1));

    //  dist=actiRotT3.DistanceToIn(ponb2x);
    //  assert(ApproxEqual(dist,0));

    dist=passRotT3.DistanceToIn(G4ThreeVector(10,0,0));
    assert(ApproxEqual(dist,10));

    dist=passRotT3.DistanceToIn(G4ThreeVector(0,10,0));
    assert(ApproxEqual(dist,0));

// DistanceToOut(P,V)

    //    dist=actiRotT3.DistanceToOut(G4ThreeVector(1,-1,0),vy,
    //                         calcNorm,pgoodNorm,pNorm);
    //  assert(ApproxEqual(dist,1)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);

    //  dist=actiRotT3.DistanceToOut(G4ThreeVector(1,-1,0),vxmy,
    //                             calcNorm,pgoodNorm,pNorm);
    // assert(ApproxEqual(dist,50-std::sqrt(2.0))&&ApproxEqual(norm,vxmy)&&*pgoodNorm);

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,10,0),vx,
                                 calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vx)&&*pgoodNorm);

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,10,0),vmy,
                                 calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,10,0),vz,
                                 calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,10,0),vmz,
                                 calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

//  cout<<"b1.DistanceToOut(ponxside,vy) = "<<dist<<G4endl;
//  assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);

// DistanceToOut(P)

    // dist=actiRotT3.DistanceToOut(pzero);
    //  assert(ApproxEqual(dist,0));

    //  dist=actiRotT3.DistanceToOut(vx);
    //  assert(ApproxEqual(dist,0));

    //  dist=actiRotT3.DistanceToOut(G4ThreeVector(1,-1,0));
    //  assert(ApproxEqual(dist,1));

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,20,45));
    assert(ApproxEqual(dist,5));

    dist=passRotT3.DistanceToOut(G4ThreeVector(-2,20,0));
    assert(ApproxEqual(dist,2));

    dist=passRotT3.DistanceToOut(G4ThreeVector(-20,2,0));
    assert(ApproxEqual(dist,2));

  // Point on surface
    G4cout<<G4endl;
    G4cout<<"Point on surface of 10x10x10 box shifted -10 along x-axis:"<<G4endl<<G4endl;
    for(i=0;i<10;i++)
    {
      p = rdb4.GetPointOnSurface();
      G4cout<<p.x()<<"\t"<<p.y()<<"\t"<<p.z()<<G4endl;
    }
    G4cout<<G4endl;





// CalculateExtent

    G4VoxelLimits limit ;		// Unlimited
    G4RotationMatrix noRot ;
    G4AffineTransform origin ;
    G4double min,max;
    G4bool calcExt;
    G4cout<<G4endl;
    G4cout<<"passRotT3.CalculateExtent: X, Y, Z:"<<G4endl<<G4endl;

    calcExt=passRotT3.CalculateExtent(kXAxis,limit,origin,min,max);
    G4cout<<"minX = "<<min<<"\t"<<"maxX = "<<max<<G4endl;
    //  assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(passRotT3.CalculateExtent(kYAxis,limit,origin,min,max));
    G4cout<<"minY = "<<min<<"\t"<<"maxY = "<<max<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(passRotT3.CalculateExtent(kZAxis,limit,origin,min,max));
    G4cout<<"minZ = "<<min<<"\t"<<"maxZ = "<<max<<G4endl<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));


    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);

    assert(passRotT3.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    G4cout<<"passRotT3.CE(kXAxis,limit,tPosOnly,min = "
          <<min<<"; max = "<<max<<G4endl ;
    //  assert(ApproxEqual(min,-100)&&ApproxEqual(max,-50));

    assert(passRotT3.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    G4cout<<"passRotT3.CE(kYAxis,limit,tPosOnly,min = "
          <<min<<"; max = "<<max<<G4endl ;
    //assert(ApproxEqual(min,-160)&&ApproxEqual(max,-60));

    assert(passRotT3.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    G4cout<<"passRotT3.CE(kZAxis,limit,tPosOnly,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl ;
    //  assert(ApproxEqual(min,-170)&&ApproxEqual(max,-70));


    G4RotationMatrix r90Z;
    r90Z.rotateZ(pi/2);
    G4AffineTransform tRotZ(r90Z,pzero);

    assert(passRotT3.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    G4cout<<"passRotT3.CE(kXAxis,limit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl ;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(passRotT3.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    G4cout<<"passRotT3.CE(kYAxis,limit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl ;
    //    assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(passRotT3.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    G4cout<<"passRotT3.CE(kZAxis,limit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl ;
    //  assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    G4cout<<"rdb4.CalculateExtent: X, Y, Z:"<<G4endl<<G4endl;

    assert(rdb4.CalculateExtent(kXAxis,limit,origin,min,max));
    G4cout<<"minX = "<<min<<"\t"<<"maxX = "<<max<<G4endl;
    //  assert(ApproxEqual(min,-50)&&ApproxEqual(max,0));

    assert(rdb4.CalculateExtent(kYAxis,limit,origin,min,max));
    G4cout<<"minY = "<<min<<"\t"<<"maxY = "<<max<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(rdb4.CalculateExtent(kZAxis,limit,origin,min,max));
    G4cout<<"minZ = "<<min<<"\t"<<"maxZ = "<<max<<G4endl<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));









    /* *********************************************************



    G4cout<<"rs1.CalculateExtent: X, Y, Z:"<<G4endl<<G4endl;

    assert(rs1.CalculateExtent(kXAxis,limit,origin,min,max));
    G4cout<<"minX = "<<min<<"\t"<<"maxX = "<<max<<G4endl;
    //  assert(ApproxEqual(min,-50)&&ApproxEqual(max,0));

    assert(rs1.CalculateExtent(kYAxis,limit,origin,min,max));
    G4cout<<"minY = "<<min<<"\t"<<"maxY = "<<max<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

    assert(rs1.CalculateExtent(kZAxis,limit,origin,min,max));
    G4cout<<"minZ = "<<min<<"\t"<<"maxZ = "<<max<<G4endl<<G4endl;
    // assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));







// Check that clipped away

    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-50);
    assert(!passRotT3.CalculateExtent(kXAxis,xClip,origin,min,max));

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

    assert(passRotT3.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(passRotT3.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(passRotT3.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    G4VoxelLimits buggyClip2;
    buggyClip2.AddLimit(kXAxis,5,15);

    assert(passRotT3.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(passRotT3.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(passRotT3.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    buggyClip2.AddLimit(kYAxis,5,15);

    assert(passRotT3.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(passRotT3.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(passRotT3.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4VoxelLimits buggyClip1;
    buggyClip1.AddLimit(kXAxis,-5,+5);

    assert(passRotT3.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(passRotT3.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(passRotT3.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    buggyClip1.AddLimit(kYAxis,-5,+5);

    assert(passRotT3.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(passRotT3.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(passRotT3.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    ****************************************************** */


  return 0 ;
}














