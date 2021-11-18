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

// testG4Torus
//
//  Test file for class G4Torus [NOT thorough]
//
//   History
// 30.10.96     V.Grichine First version for first G4Torus implementation

#include "G4ios.hh"
#undef NDEBUG
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Torus.hh"
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
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}


G4bool testG4Torus()
{
   G4int i ;
   G4double Rtor = 100 ;
   G4double Rmax = Rtor*0.9 ;
   G4double Rmin = Rtor*0.1 ;

//G4double z = atof ( argv[2] );
   G4double x; 
   G4double z;
 
   G4double Dist, dist, vol, volCheck;
   EInside side;
   G4ThreeVector *pNorm,norm;
   G4bool *pgoodNorm,goodNorm,calcNorm=true;

   G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

   pNorm=&norm;
   pgoodNorm=&goodNorm;
   
   G4ThreeVector pzero(0,0,0);

   G4ThreeVector pbigx(240,0,0),pbigy(0,240,0),pbigz(0,0,240);
   G4ThreeVector pbigmx(-240,0,0),pbigmy(0,-240,0),pbigmz(0,0,-240);

   G4ThreeVector ponrmax(190,0,0);
   G4ThreeVector ponrmin(0,110,0);
   G4ThreeVector ponrtor(0,100,0);
   G4ThreeVector ponphi1(100/std::sqrt(2.),100/std::sqrt(2.),0) ;
   G4ThreeVector ponphi2(-100/std::sqrt(2.),100/std::sqrt(2.),0) ;
   G4ThreeVector ponphi12(120/std::sqrt(2.),120/std::sqrt(2.),0) ;
   G4ThreeVector ponphi22(-120/std::sqrt(2.),120/std::sqrt(2.),0) ;
   G4ThreeVector ponphi23(-120/std::sqrt(2.)+0.5,120/std::sqrt(2.),0) ;
    

   G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
   G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
   G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);
   G4ThreeVector vmxy(-1/std::sqrt(2.0),1/std::sqrt(2.0),0);
   G4ThreeVector vmxmy(-1/std::sqrt(2.0),-1/std::sqrt(2.0),0);
   G4ThreeVector vxmy(1/std::sqrt(2.0),-1/std::sqrt(2.0),0);

   G4ThreeVector pstart((Rtor+Rmax)/std::sqrt(2.0),(Rtor+Rmax)/std::sqrt(2.0),0) ;
   G4ThreeVector vdirect(1/std::sqrt(2.0),-1/std::sqrt(2.0),0) ;
   
   G4ThreeVector pother(110,0,0);
   vdirect = vdirect.unit() ;
   G4ThreeVector p1;
   G4ThreeVector v1(1,0,0);
   v1 = v1.unit();


// Check torus roots
   
   G4Torus t1("Solid Torus #1",0,Rmax,Rtor,0,twopi);
   G4Torus t2("Hole cutted Torus #2",Rmin,Rmax,Rtor,pi/4,halfpi);
   G4Torus tn2("tn2",Rmin,Rmax,Rtor,halfpi,halfpi);
   G4Torus tn3("tn3",Rmin,Rmax,Rtor,halfpi,3*halfpi);



   G4Torus t3("Hole cutted Torus #3",4*Rmin,Rmax,Rtor,halfpi-pi/24,pi/12);
   G4Torus t4("Solid Torus #4",0,Rtor - 2.e3*kCarTolerance,Rtor,0,twopi);
   G4Torus t5("Solid cutted Torus #5",0,Rtor - 2.e3*kCarTolerance,Rtor,pi/4,halfpi);
   
   G4Torus * aTub = new G4Torus("Ring1", 0*cm, 10*cm, 
                                         1*m, 0*deg, 360*deg ); 
   G4Torus t6("t6",100*mm, 150*mm, 200*mm,0*degree,
                                       59.99999999999999*degree);

    G4Torus* clad =
      new G4Torus("clad",0.,1.*cm,10.*cm,0.*deg,180.*deg);    // external

    G4Torus* core =
      new G4Torus("core",0.,0.5*cm,10.*cm,0.*deg,180.*deg); // internal


   G4cout.precision(20);

   G4ThreeVector p1t6( 60.73813233071262,  
                      -27.28494547459707,
                       37.47827539879173);

   G4ThreeVector vt6( 0.3059312222729116,
		      0.8329513862588347,
		     -0.461083588265824);

   G4ThreeVector p2t6(70.75950555416668,
		      -3.552713678800501e-15,
                      22.37458414788935       );


  //   num = t1.TorusRoots(Ri,pstart,vdirect) ;
  //   num = t1.TorusRoots(Ri,pother,vx) ;
  //    G4Torus t2("Hole Torus #2",45,50,50,0,360);
  // Check name
  // assert(t1.GetName()=="Solid Torus #1");

  // check cubic volume

   vol  = t1.GetCubicVolume();
   volCheck = twopi*pi*Rtor*(Rmax*Rmax);
   assert(ApproxEqual(vol,volCheck));
    
   vol  = t2.GetCubicVolume();
   volCheck = halfpi*pi*Rtor*(Rmax*Rmax-Rmin*Rmin);
   assert(ApproxEqual(vol,volCheck));

    // Check Inside
    assert(t1.Inside(pzero)==kOutside);
    assert(t1.Inside(pbigx)==kOutside);
    assert(t1.Inside(ponrmax)==kSurface);
    assert(t2.Inside(ponrmin)==kSurface);
    assert(t2.Inside(pbigx)==kOutside);
    assert(t2.Inside(pbigy)==kOutside);
    
    assert(t2.Inside(ponphi1)==kOutside);
    assert(t2.Inside(ponphi2)==kOutside);
    assert(t2.Inside(ponphi12)==kSurface);
    assert(t2.Inside(ponphi22)==kSurface);

    side=t6.Inside(p1t6);    
    // G4cout<<"t6.Inside(p1t6) = "<<OutputInside(side)<<G4endl;
    side=t6.Inside(p2t6);    
    // G4cout<<"t6.Inside(p2t6) = "<<OutputInside(side)<<G4endl;

// Check Surface Normal
    G4ThreeVector normal;
    G4double p2=1./std::sqrt(2.); // ,p3=1./std::sqrt(3.);


    normal=t1.SurfaceNormal(ponrmax);
    assert(ApproxEqual(normal,vx));
    normal=t1.SurfaceNormal(G4ThreeVector(0.,190.,0.));
    assert(ApproxEqual(normal,vy));
    normal=tn2.SurfaceNormal(G4ThreeVector(0.,Rtor+Rmax,0.));
    assert(ApproxEqual(normal,G4ThreeVector(p2,p2,0.)));
    normal=tn2.SurfaceNormal(G4ThreeVector(0.,Rtor+Rmin,0.));
    assert(ApproxEqual(normal,G4ThreeVector(p2,-p2,0.)));
    normal=tn2.SurfaceNormal(G4ThreeVector(0.,Rtor-Rmin,0.));
    assert(ApproxEqual(normal,G4ThreeVector(p2,p2,0.)));
    normal=tn2.SurfaceNormal(G4ThreeVector(0.,Rtor-Rmax,0.));
    assert(ApproxEqual(normal,G4ThreeVector(p2,-p2,0.)));
    normal=tn3.SurfaceNormal(G4ThreeVector(Rtor,0.,Rmax));
    assert(ApproxEqual(normal,G4ThreeVector(0.,p2,p2)));
    normal=tn3.SurfaceNormal(G4ThreeVector(0.,Rtor,Rmax));
    assert(ApproxEqual(normal,G4ThreeVector(p2,0.,p2)));
    
    normal=t2.SurfaceNormal(ponrmin);
    assert(ApproxEqual(normal,vmy));
    normal=t2.SurfaceNormal(ponphi1);
    assert(ApproxEqual(normal,vxmy));
    normal=t2.SurfaceNormal(ponphi2);
    assert(ApproxEqual(normal,vmxmy));

// DistanceToOut(P)
    Dist=t1.DistanceToOut(ponrmin);
    assert(ApproxEqual(Dist,80));
    Dist=t1.DistanceToOut(ponrmax);
    assert(ApproxEqual(Dist,0));
    /* // later: why it was introduced, while they are outside (see above)
    Dist=t2.DistanceToOut(ponphi1);
    assert(ApproxEqual(Dist,0));
    Dist=t2.DistanceToOut(ponphi2);
    assert(ApproxEqual(Dist,0));
    */
// DistanceToOut(P,V)
    Dist=t1.DistanceToOut(ponrmax,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vx));
    Dist=t1.DistanceToOut(ponphi1,vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"t1.DistanceToOut(ponphi1,vz,...) = "<<Dist<<G4endl ;
    assert(ApproxEqual(Dist,90)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vz));
    
    Dist=t1.DistanceToOut(ponrmin,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,80)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vy));
    Dist=t1.DistanceToOut(ponrmin,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,100) && !*pgoodNorm);
//    Dist=t1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
//    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&*pgoodNorm);
//    Dist=t1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
//    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
//    Dist=t1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
//    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);


    Dist=t2.DistanceToOut(ponphi12,vxmy,calcNorm,pgoodNorm,pNorm);
//    G4cout<<"Dist=t2.DistanceToOut(ponphi12,vxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vxmy));
    Dist=t2.DistanceToOut(ponphi22,vmxmy,calcNorm,pgoodNorm,pNorm);
//    G4cout<<"Dist=t2.DistanceToOut(ponphi22,vmxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vmxmy));
    Dist=t2.DistanceToOut(ponphi23,vmxmy,calcNorm,pgoodNorm,pNorm);
//    G4cout<<"Dist=t2.DistanceToOut(ponphi23,vmxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.353553)&&*pgoodNorm&&ApproxEqual(pNorm->unit(),vmxmy));
  

// Check for Distance to Out ( start from an internal point )

  for ( i=0; i<12; i++ ) 
  {
    x = -1050;
    z = G4double(i)/10;
    p1 = G4ThreeVector(x,0,z);
    //    G4cout << p1 << " - " << v1 << G4endl;

    Dist = aTub->DistanceToIn (p1,v1); 
    //    G4cout << "Distance to in dir: " << Dist ;

    Dist = aTub->DistanceToOut (p1,v1) ;
    //    G4cout << "   Distance to out dir: " << Dist << G4endl ;

    //    G4cout << "Distance to in : " << aTub->DistanceToIn (p1);
    //    G4cout << "   Distance to out : " << aTub->DistanceToOut (p1)
    //           << G4endl;
    //    G4cout << "   Inside : " << aTub->Inside (p1);
    //    G4cout << G4endl;
  }



//DistanceToIn(P)
    Dist=t1.DistanceToIn(pbigx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigz);
//    G4cout<<"Dist=t1.DistanceToIn(pbigz) = "<<Dist<<G4endl;
//    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmz);
//    G4cout<<"Dist=t1.DistanceToIn(pbigmz) = "<<Dist<<G4endl;
//    assert(ApproxEqual(Dist,50));

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
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t1.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t1.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t1.DistanceToIn(ponrmax,vx);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t1.DistanceToIn(ponrmax,vmx);
    assert(ApproxEqual(Dist,0));

    G4ThreeVector vnew(1,0,0) ;
    vnew.rotateZ(pi/4-5*1e-9) ;    // old test: check pzero with vxy
    Dist=t2.DistanceToIn(pzero,vnew);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t2.DistanceToIn(pzero,vy);
    assert(ApproxEqual(Dist,10));
    Dist=t2.DistanceToIn(ponphi12,vy);
    assert(ApproxEqual(Dist,0));
    Dist=t2.DistanceToIn(ponphi12,vmy);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t2.DistanceToIn(ponphi1,vy);
//    G4cout<<"Dist=t2.DistanceToIn(ponphi1,vy) = "<<Dist<<G4endl;  // about 13
    Dist=t2.DistanceToIn(ponrmin,vy);
    assert(ApproxEqual(Dist,0));
    Dist=t2.DistanceToIn(ponrmin,vmy);
    assert(ApproxEqual(Dist,20));

    Dist=t3.DistanceToIn(ponrtor,vy);
    assert(ApproxEqual(Dist,40));
    Dist=t3.DistanceToIn(ponrtor,vmy);
    assert(ApproxEqual(Dist,40));
    Dist=t3.DistanceToIn(ponrtor,vz);
    assert(ApproxEqual(Dist,40));
    Dist=t3.DistanceToIn(ponrtor,vmz);
    assert(ApproxEqual(Dist,40));
    Dist=t3.DistanceToIn(ponrtor,vx);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=t3.DistanceToIn(ponrtor,vmx);
    assert(ApproxEqual(Dist,kInfinity));

    Dist=t6.DistanceToIn(p1t6,vt6);
    // G4cout<<"t6.DistanceToIn(p1t6,vt6) = "<<Dist<<G4endl;

    // Bug 810

    G4ThreeVector pTmp(0.,0.,0.);

    dist = clad->DistanceToIn(pTmp,vy);   
    pTmp += dist*vy;
    G4cout<<"pTmpX = "<<pTmp.x()<<";  pTmpY = "<<pTmp.y()<<";  pTmpZ = "<<pTmp.z()<<G4endl;
    side=core->Inside(pTmp);    
    G4cout<<"core->Inside(pTmp) = "<<OutputInside(side)<<G4endl;
    side=clad->Inside(pTmp);    
    G4cout<<"clad->Inside(pTmp) = "<<OutputInside(side)<<G4endl;

    dist = core->DistanceToIn(pTmp,vy);   
    pTmp += dist*vy;
    G4cout<<"pTmpX = "<<pTmp.x()<<";  pTmpY = "<<pTmp.y()<<";  pTmpZ = "<<pTmp.z()<<G4endl;
    side=core->Inside(pTmp);    
    G4cout<<"core->Inside(pTmp) = "<<OutputInside(side)<<G4endl;
    side=clad->Inside(pTmp);    
    G4cout<<"clad->Inside(pTmp) = "<<OutputInside(side)<<G4endl;

    dist = core->DistanceToOut(pTmp,vy,calcNorm,pgoodNorm,pNorm);   
    pTmp += dist*vy;
    G4cout<<"pTmpX = "<<pTmp.x()<<";  pTmpY = "<<pTmp.y()<<";  pTmpZ = "<<pTmp.z()<<G4endl;
    side=core->Inside(pTmp);    
    G4cout<<"core->Inside(pTmp) = "<<OutputInside(side)<<G4endl;
    side=clad->Inside(pTmp);    
    G4cout<<"clad->Inside(pTmp) = "<<OutputInside(side)<<G4endl;

    dist = clad->DistanceToOut(pTmp,vy,calcNorm,pgoodNorm,pNorm);   
    pTmp += dist*vy;
    G4cout<<"pTmpX = "<<pTmp.x()<<";  pTmpY = "<<pTmp.y()<<";  pTmpZ = "<<pTmp.z()<<G4endl;
    side=core->Inside(pTmp);    
    G4cout<<"core->Inside(pTmp) = "<<OutputInside(side)<<G4endl;
    side=clad->Inside(pTmp);    
    G4cout<<"clad->Inside(pTmp) = "<<OutputInside(side)<<G4endl;

// Check for Distance to In ( start from an external point )

   for ( i=0; i<12; i++ ) 
   {
     x = -1200;
     z = G4double(i)/10;
     p1 = G4ThreeVector(x,0,z);
     //     G4cout << p1 << " - " << v1 << G4endl;

     Dist = aTub->DistanceToIn (p1,v1) ;
     //     G4cout << "Distance to in dir: " << Dist ;

     Dist = aTub->DistanceToOut (p1,v1) ;
     //     G4cout << "   Distance to out dir: " << Dist << G4endl ;

     // G4cout << "Distance to in : " << aTub->DistanceToIn (p1);
     //  G4cout << "   Distance to out : " << aTub->DistanceToOut (p1)
     //    << G4endl;
     //     G4cout << "   Inside : " << aTub->Inside (p1);
     //     G4cout << G4endl;
   }
    
    
// CalculateExtent

    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    assert(t1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-190&&max>=190);
    assert(t1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-190&&max>=190);
    assert(t1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-90&&max>=90);
    
    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);
    assert(t1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(min<=-290&&max>=-90);
    assert(t1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(min<=-300&&max>=-100);
    assert(t1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(min<=-210&&max>=-30);

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
    xClip.AddLimit(kXAxis,-300,-200);
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
    assert(t4.CalculateExtent(kXAxis,allClip,tGen,min,max));
//    G4cout<<"min = "<<min<<"    max = "<<max<<G4endl;
    assert(min<=-5&&max>=5);
    assert(t4.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(t4.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);

    assert(t5.CalculateExtent(kXAxis,allClip,tGen,min,max));
    //  G4cout<<"min = "<<min<<"    max = "<<max<<G4endl;
    assert(t5.CalculateExtent(kYAxis,allClip,tGen,min,max));
    //  G4cout<<"min = "<<min<<"    max = "<<max<<G4endl;
    assert(t5.CalculateExtent(kZAxis,allClip,tGen,min,max));
    //  G4cout<<"min = "<<min<<"    max = "<<max<<G4endl;

    t1.CalculateExtent(kZAxis,allClip,tGen,min,max);
    //  G4cout<<"min = "<<min<<"    max = "<<max<<G4endl;

// Test z clipping ok
    for (G4double zTest=-200;zTest<200;zTest+=9)
	{
	    G4VoxelLimits zTestClip;
	    zTestClip.AddLimit(kZAxis,-kInfinity,zTest);
	    if (zTest<-100)
		{
		    assert(!t4.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		}
	    else
		{
		    assert(t4.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    G4double testMin=-100+2e3*kCarTolerance;
		    G4double testMax=(zTest<100-2e3*kCarTolerance) ? zTest : 100-2e3*kCarTolerance;
		    assert (ApproxEqual(min,testMin)
			    &&ApproxEqual(max,testMax));
		}
	}


    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Torus());
    return 0;
}

