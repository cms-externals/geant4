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

// testG4Hype
//
//  Test file for class G4Hype [NOT thorough]
//
//             Ensure asserts are compiled in

#undef NDEBUG
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Hype.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool testG4Hype()
{
    G4ThreeVector pzero(0,0,0);

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector ponxside(50,0,0);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxy(-1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/std::sqrt(2.0),-1/std::sqrt(2.0),0);
    G4ThreeVector vxmy(1/std::sqrt(2.0),-1/std::sqrt(2.0),0);

    G4double Dist;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    double iR1=0;
    double oR1=50;
    double iR2=45;
    double oR2=50;
    double noStereo=0;
    double yesStereo=0.3;
    double len=50;

    G4Hype t1("Solid Hype #1",iR1,oR1,noStereo,yesStereo,len);    
    G4Hype t2("Hole  Hype #2",iR2,oR2,yesStereo,yesStereo,len);

    G4ThreeVector Spoint ;
    G4double dist ;
    for ( int i = 0 ; i < 10 ; i++ ) {
      //  G4cout << "Event " << i << G4endl << G4endl ;
      Spoint = t1.GetPointOnSurface() ;
      dist = t1.DistanceToIn(Spoint,-Spoint/Spoint.mag()) ;
      G4cout << "Spoint " << Spoint << " " <<  dist << G4endl ;
    }

    t1.GetPointOnSurface() ;
    t2.GetPointOnSurface() ; 


// Check name
    assert(t1.GetName()=="Solid Hype #1");

// Check Inside
    assert(t1.Inside(pzero)==kInside);
    assert(t1.Inside(pbigx)==kOutside);

// Check Surface Normal
    G4ThreeVector normal;

    normal=t1.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,vx));

// DistanceToOut(P)
    Dist=t1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,50));

// DistanceToOut(P,V)
    Dist=t1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vx)&&(!*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmx)&&(!*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vy)&&(!*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmy)&&(!*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&(*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&(*pgoodNorm));
    Dist=t1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&(!*pgoodNorm));

    Dist=t2.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=t2.DistanceToOut(pzero,vxy) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(ponxside,vmx,calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=t2.DistanceToOut(ponxside,vmx) = "<<Dist<<G4endl;
    Dist=t2.DistanceToOut(ponxside,vmxmy,calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=t2.DistanceToOut(ponxside,vmxmy) = "<<Dist<<G4endl;
    Dist=t2.DistanceToOut(ponxside,vz,calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=t2.DistanceToOut(ponxside,vz) = "<<Dist<<G4endl;


//DistanceToIn(P)
    Dist=t1.DistanceToIn(pbigx);
    assert(Dist <= 50*(1.0+DBL_EPSILON));
    Dist=t1.DistanceToIn(pbigmx);
    assert(Dist <= 50*(1.0+DBL_EPSILON));
    Dist=t1.DistanceToIn(pbigy);
    assert(Dist <= 50*(1.0+DBL_EPSILON));
    Dist=t1.DistanceToIn(pbigmy);
    assert(Dist <= 50*(1.0+DBL_EPSILON));
    Dist=t1.DistanceToIn(pbigz);
    assert(Dist <= 50*(1.0+DBL_EPSILON));
    Dist=t1.DistanceToIn(pbigmz);
    assert(Dist <= 50*(1.0+DBL_EPSILON));

// DistanceToIn(P,V)
    Dist=t1.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));

// CalculateExtent
    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    assert(t1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    
    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);
    assert(t1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(min<=-150&&max>=-50);
    assert(t1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(min<=-160&&max>=-60);
    assert(t1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);

    G4RotationMatrix r90Z;
    r90Z.rotateZ(halfpi);
    G4AffineTransform tRotZ(r90Z,pzero);
    assert(t1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);

// Check that clipped away
    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-60);
    assert(!t1.CalculateExtent(kXAxis,xClip,origin,min,max));

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
    assert(t1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(t1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(t1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);


// Test z clipping ok
    for (G4double zTest=-100;zTest<100;zTest+=9)
	{
	    G4VoxelLimits zTestClip;
	    zTestClip.AddLimit(kZAxis,-kInfinity,zTest);
	    if (zTest<-50)
		{
		    assert(!t1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		}
	    else
		{
		    assert(t1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    G4double testMin=-50;
		    G4double testMax=(zTest<50) ? zTest : 50;
		    assert (ApproxEqual(min,testMin)
			    &&ApproxEqual(max,testMax));
		}
	}

// Test y clipping ok

    // Hype end Outer Radius -> clipping scale
    double eOR=std::sqrt(std::tan(yesStereo)*std::tan(yesStereo)*len*len+oR1*oR1);

    for (G4double xTest=-100;xTest<100;xTest+=9)
	{
	    G4VoxelLimits xTestClip;
	    xTestClip.AddLimit(kXAxis,-kInfinity,xTest);
	    if (xTest<-eOR)
		{
		    assert(!t1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
		}
	    else
		{
		   assert(t1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
// Calc max y coordinate
		   G4double testMax=(xTest<0) ? std::sqrt(eOR*eOR-xTest*xTest) : eOR;
		   assert ((min < -testMax) && (max > testMax));
		}
	}

// Test x clipping ok
    for (G4double yTest=-100;yTest<100;yTest+=9)
	{
	    G4VoxelLimits yTestClip;
	    yTestClip.AddLimit(kYAxis,-kInfinity,yTest);
	    if (yTest<-eOR)
		{
		    assert(!t1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
		}
	    else
		{
		   assert(t1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
// Calc max y coordinate
		   G4double testMax=(yTest<0) ? std::sqrt(eOR*eOR-yTest*yTest) : eOR;
		   assert ((min < -testMax) && (max > testMax));
		}
	}


    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Hype());
    return 0;
}

