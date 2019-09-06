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

// testG4Tet
//
//  Test file for class G4Tet [NOT thorough]
//
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
#include "G4Tet.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool testG4Tet()
{
    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pnt1(10.,0.,0.),pnt2(5.0,10.,0), pnt3(5.,5.,10.);

    G4ThreeVector norm;

    G4bool  goodTet;
    G4Tet   t1( "Solid Tet #1", pzero, pnt1, pnt2, pnt3, &goodTet); 

// Check name
    assert(t1.GetName()=="Solid Tet #1");

    G4ThreeVector pntA( 1.0 , 1.0 , 1.0 ); 
    G4ThreeVector pntB( 1.5 , 0.5 , 1.0 );  
    // G4ThreeVector pntBbis= 0.1 * (pnt1-pzero) + 0.1 * (pnt3-pzero); 
    G4ThreeVector pntBr023= (1.0/3.0) * (pzero + pnt2 + pnt3); 
    G4ThreeVector pntC( 0.0,  5.0 , 1.5 );  

// Check Inside
    assert(t1.Inside(pntA)==kInside);
    assert(t1.Inside(pntB)==kSurface);
    assert(t1.Inside(pntBr023)==kSurface);
    assert(t1.Inside(pntC)==kOutside);

// Check Surface Normal
    G4ThreeVector normal;
    G4ThreeVector pntOnBotSurf012( 5.0, 5.0, 0.0); 
    G4ThreeVector vmz(0,0,-1.0); 

    normal=t1.SurfaceNormal(pntOnBotSurf012);
    assert(ApproxEqual(normal,vmz));


    return true;   
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Tet());
    return 0;
}

