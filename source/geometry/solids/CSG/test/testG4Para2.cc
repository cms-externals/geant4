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
// Test of G4Para
// Includes all/most of the tests Done for a box

#undef NDEBUG
#include <assert.h>
#include <cmath>
#include "G4ios.hh"

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Para.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"



//#include "G4ios.hh"
//#include "globals.hh"
//#include "G4Para.hh"

#define	DELTA 0.0001

// Returns false if actual is within wanted+/- DELTA
//         true if error
G4bool OutRange(G4double actual,G4double wanted)
{
    G4bool rng=false;
    if (actual<wanted-DELTA||actual>wanted+DELTA) rng=true;
    return rng;
}
G4bool OutRange(G4ThreeVector actual,G4ThreeVector wanted)
{
    G4bool rng=false;
    if (OutRange(actual.x(),wanted.x())
	||OutRange(actual.y(),wanted.y())
	||OutRange(actual.z(),wanted.z())  ) rng=true;
    return rng;
}

int main(void)
{
    G4double Dist;
    G4ThreeVector pzero(0,0,0),px(10,0,0),py(0,10,0),pz(0,0,10);
    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);
    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40);
    G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
    G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40);
    G4RotationMatrix runit;
    G4RotationMatrix r90X,r90Y,r90Z,r45X,r30Y;
    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1,1,0);
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    r90X.rotateX(halfpi);
    r90Y.rotateY(halfpi);
    r90Z.rotateZ(halfpi);
    r45X.rotateX(pi/4);
    r30Y.rotateY(pi/6);

    vxy=vxy.unit();

    G4Para p1("Box",20,30,40,0,0,0),
	p2("2",50,50,50,pi/6,0,0),
	p3("3",50,50,50,0,pi/6,0),
	p4("4",50,50,50,0,0,pi/6),
	p5("5",50,50,50,0,pi/6,pi/6),	
	p6("6",50,50,50,pi/6,pi/6,pi/6);	

    G4cout << "Name:"<< p1.GetName()
	 << " ID=" <<G4endl;

    G4cout << "Checking G4Para::Get Alpha/Theta/Phi ...\n";
    if (p1.GetAlpha()!=0.0){
       const G4String Err1("Error1: p1 - alpha incorrect");
       G4cout << Err1 << G4endl; // "Error 1: p1 - alpha" << G4endl;
       G4cerr << Err1 << G4endl; // "Error 1: p1 - alpha" << G4endl;
    }
    if ( ! ApproxEqual( p2.GetAlpha() , pi/6.0 ) ){
       const G4String Err2("Error 2: p2 - Get alpha incorrect");
       G4cout << Err2 << G4endl;
       G4cerr << Err2 << G4endl;
    }
    if ( ! ApproxEqual( p3.GetTheta() , pi/6.0 ) ){
       const G4String Err3("Error 3: p3 - Get theta incorrect.");
       G4cout << Err3 << G4endl;
       G4cerr << Err3 << G4endl;
    }
    // p4.SetThetaAndPhi(0,pi/6);
    if ( false ) { //  ( ! ApproxEqual( p4.GetPhi() , pi/6.0) ){  // Test fails!
                   //  Reason: theta=0 phi=finite is same as theta=0 phi=0
       const G4String Err4("Error 4: p4 - Get Phi incorrect. ");
       G4cout << Err4 << " value = " << p4.GetPhi() << " vs expected = " << pi/6.0 << G4endl;
       G4cerr << Err4 << G4endl;
    }
    if ( ! ApproxEqual( p5.GetTheta() , pi/6.0 ) ){
       const G4String Err5a("Error 5/a: p5 - Get theta incorrect.");
       G4cout << Err5a << G4endl;
       G4cerr << Err5a << G4endl;
    }
    if ( ! ApproxEqual( p5.GetPhi() , pi/6.0) ){
       const G4String Err5b("Error 5/b: p5 - Get Phi incorrect. ");
       G4cout << Err5b << " value = " << p5.GetPhi() << " vs expected = " << pi/6.0 << G4endl;
       G4cerr << Err5b << G4endl;
    }
    if ( ! ApproxEqual( p6.GetAlpha() , pi/6.0 ) ){
       const G4String Err6a("Error6/a: p6 - alpha incorrect");
       G4cout << Err6a << G4endl;
       G4cerr << Err6a << G4endl;
    }
    if ( ! ApproxEqual( p6.GetTheta() , pi/6.0 ) ){
       const G4String Err6b("Error 6/b: p6 - Get theta incorrect.");
       G4cout << Err6b << G4endl;
       G4cerr << Err6b << G4endl;
    }
    if ( ! ApproxEqual( p6.GetPhi() , pi/6.0) ){
       const G4String Err6c("Error 6/c: p6 - Get Phi incorrect. ");
       G4cout << Err6c << " value = " << p6.GetPhi() << " vs expected = " << pi/6.0 << G4endl;
       G4cerr << Err6c << G4endl;
    }

    G4cout << "Checking G4Para::Inside...\n";
    if (p1.Inside(pzero)!=kInside){
	G4cout << "Error A" << G4endl;
	G4cerr << "Error A" << G4endl;
    }
    if (p1.Inside(pbigz)!=kOutside) {
	G4cout << "Error B" << G4endl;
	G4cerr << "Error B" << G4endl;
    }
    if (p1.Inside(ponxside)!=kSurface){
	G4cout << "Error C" << G4endl;
	G4cerr << "Error C" << G4endl;
    }
    if (p1.Inside(ponyside)!=kSurface){
	G4cout << "Error D" << G4endl;
	G4cerr << "Error D" << G4endl;
    }
    if (p1.Inside(ponzside)!=kSurface){
	G4cout << "Error E" << G4endl;
	G4cerr << "Error E" << G4endl;
    }


    G4cout << "Checking G4Para::SurfaceNormal...\n";
    norm=p1.SurfaceNormal(ponxside);
    if (OutRange(norm,G4ThreeVector(1,0,0)))
	G4cout << "Error A " << norm << G4endl;
    norm=p1.SurfaceNormal(ponmxside);
    if (OutRange(norm,G4ThreeVector(-1,0,0)))
	G4cout << "Error B " << norm << G4endl;
    norm=p1.SurfaceNormal(ponyside);
    if (OutRange(norm,G4ThreeVector(0,1,0)))
	G4cout << "Error C " << norm << G4endl;
    norm=p1.SurfaceNormal(ponmyside);
    if (OutRange(norm,G4ThreeVector(0,-1,0)))
	G4cout << "Error D " << norm << G4endl;
    norm=p1.SurfaceNormal(ponzside);
    if (OutRange(norm,G4ThreeVector(0,0,1)))
	G4cout << "Error E " << norm << G4endl;
    norm=p1.SurfaceNormal(ponmzside);
    if (OutRange(norm,G4ThreeVector(0,0,-1)))
	G4cout << "Error F " << norm << G4endl;
    norm=p1.SurfaceNormal(ponzsidey);
    if (OutRange(norm,G4ThreeVector(0,0,1)))
	G4cout << "Error G " << norm << G4endl;
    norm=p1.SurfaceNormal(ponmzsidey);
    if (OutRange(norm,G4ThreeVector(0,0,-1)))
	G4cout << "Error H " << norm << G4endl;


    G4cout << "Checking G4Para::DistanceToOut(P)...\n";
    Dist=p1.DistanceToOut(pzero);
    if (OutRange(Dist,20))
	G4cout << "Error A1 " << Dist << G4endl;
    Dist=p2.DistanceToOut(pzero);
    if (OutRange(Dist,50*std::cos(pi/6)))
	G4cout << "Error A2 " << Dist << G4endl;
     Dist=p3.DistanceToOut(pzero);
    if (OutRange(Dist,50*std::cos(pi/6)))
	G4cout << "Error A3 " << Dist << G4endl;
    Dist=p5.DistanceToOut(pzero);
    if (OutRange(Dist,2*50/std::sqrt(5.)))
	G4cout << "Error A4 " << Dist << G4endl;

     Dist=p1.DistanceToOut(px);
    if (OutRange(Dist,10))
	G4cout << "Error B " << Dist << G4endl;
    Dist=p1.DistanceToOut(py);
    if (OutRange(Dist,20))
	G4cout << "Error C " << Dist << G4endl;
    Dist=p1.DistanceToOut(pz);
    if (OutRange(Dist,20))
	G4cout << "Error D " << Dist << G4endl;




    G4cout << "Checking G4Para::DistanceToOut(P,V)...\n";

    Dist=p1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,20)||OutRange(*pNorm,vx)||!*pgoodNorm)
	G4cout << "Error A " << Dist << G4endl;
    Dist=p1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,20)||OutRange(norm,vmx)||!*pgoodNorm)
 	G4cout << "Error B " << Dist << G4endl;
    Dist=p1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,30)||OutRange(norm,vy)||!*pgoodNorm)
 	G4cout << "Error C " << Dist << G4endl;
    Dist=p1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,30)||OutRange(norm,vmy)||!*pgoodNorm)
 	G4cout << "Error D " << Dist << G4endl;
     Dist=p1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,40)||OutRange(norm,vz)||!*pgoodNorm)
 	G4cout << "Error E " << Dist << G4endl;
     Dist=p1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,40)||OutRange(norm,vmz)||!*pgoodNorm)
 	G4cout << "Error F " << Dist << G4endl;
    Dist=p1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,std::sqrt(800.))||!*pgoodNorm)
 	G4cout << "Error F " << Dist << G4endl;

    Dist=p1.DistanceToOut(ponxside,vx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(*pNorm,vx)||!*pgoodNorm)
	G4cout << "Error A2 " << Dist << G4endl;
    Dist=p1.DistanceToOut(ponmxside,vmx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmx)||!*pgoodNorm)
 	G4cout << "Error B2 " << Dist << G4endl;
    Dist=p1.DistanceToOut(ponyside,vy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vy)||!*pgoodNorm)
 	G4cout << "Error C2 " << Dist << G4endl;
    Dist=p1.DistanceToOut(ponmyside,vmy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmy)||!*pgoodNorm)
 	G4cout << "Error D2 " << Dist << G4endl;
     Dist=p1.DistanceToOut(ponzside,vz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vz)||!*pgoodNorm)
 	G4cout << "Error E2 " << Dist << G4endl;
    Dist=p1.DistanceToOut(ponmzside,vmz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmz)||!*pgoodNorm)
 	G4cout << "Error F2 " << Dist << G4endl;
  
    G4cout << "Checking G4Para::DistanceToIn(P)...\n";
    Dist=p1.DistanceToIn(pbigx);
    if (OutRange(Dist,80))
 	G4cout << "Error A " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmx);
    if (OutRange(Dist,80))
 	G4cout << "Error B " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigy);
    if (OutRange(Dist,70))
 	G4cout << "Error C " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmy);
    if (OutRange(Dist,70))
 	G4cout << "Error D " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigz);
    if (OutRange(Dist,60))
 	G4cout << "Error E " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmz);
    if (OutRange(Dist,60))
 	G4cout << "Error F " << Dist << G4endl;

    Dist=p3.DistanceToIn(pbigx);
    if (OutRange(Dist,50*std::cos(pi/6)))
	G4cout << "Error G1 " << Dist <<G4endl;
    Dist=p3.DistanceToIn(pbigy);
    if (OutRange(Dist,50))
	G4cout << "Error G2 " << Dist <<G4endl;

    G4cout << "Checking G4Para::DistanceToIn(P,V)...\n";

    Dist=p1.DistanceToIn(pbigx,vmx);
    if (OutRange(Dist,80))
 	G4cout << "Error A " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmx,vx);
    if (OutRange(Dist,80))
 	G4cout << "Error B " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigy,vmy);
    if (OutRange(Dist,70))
 	G4cout << "Error C " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmy,vy);
    if (OutRange(Dist,70))
 	G4cout << "Error D " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigz,vmz);
    if (OutRange(Dist,60))
 	G4cout << "Error E " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmz,vz);
    if (OutRange(Dist,60))
 	G4cout << "Error F " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigx,vxy);
    if (OutRange(Dist,kInfinity))
 	G4cout << "Error G " << Dist << G4endl;
    Dist=p1.DistanceToIn(pbigmx,vxy);
    if (OutRange(Dist,kInfinity))
 	G4cout << "Error H " << Dist << G4endl;

    return 0;    
}




