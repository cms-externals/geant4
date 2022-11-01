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

// testG4TwistedBox
//
//  Test file for class G4TwistedBox
//
//             Ensure asserts are compiled in

#undef NDEBUG
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4TwistedBox.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"


G4bool testG4TwistedBox()
{

    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pout ( 200*cm,200*cm,200*cm) ;
    G4ThreeVector dir = pzero-pout ;
    
    dir *= 1./dir.mag();

    G4TwistedBox t1("Solid Twisted Box #1",
			  20*deg ,    // twisted angle
			  10*cm,       // half x-length
			  15*cm,       // half y-length
			  20*cm) ;      // half z-length

    G4double dist = t1.DistanceToIn(pout,dir) ;
    G4cout << "distance = " << dist << G4endl ;

// Check Inside

    G4cout << "Test Inside(p) " << G4endl << G4endl ;

    G4cout << "pzero : " << t1.Inside(pzero) << G4endl ;
    //assert(t1.Inside(pzero)==kInside) ;

    // testing the volume

    G4double volume = t1.GetCubicVolume() ;
    G4cout << G4endl << G4endl ;
    G4cout << "Twisted Box volume = " << volume / cm / cm /cm << " cm^3" << G4endl ;

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedBox());
    return 0;
}


