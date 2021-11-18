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
// class ExDivTesterTrd Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterTrd.hh"
#include "G4Trd.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterTrd::
ExDivTesterTrd( PVType& pvtype, PlaceType& postype,
                std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, postype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kXAxis );
  theAxis.push_back( kYAxis );
  theAxis.push_back( kZAxis );
}
 
//--------------------------------------------------------------------------
void ExDivTesterTrd::GenerateScanPoints()
{
  GenerateScanPointsAsBox();
}

//--------------------------------------------------------------------------
void ExDivTesterTrd::BuildParentSolids()
{
  theParentSolids.push_back( new G4Trd("parent_1", theWorldLengthXY, theWorldLengthXY,
                             theWorldLengthXY*0.5, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Trd("parent_2", theWorldLengthXY*0.5, theWorldLengthXY,
                             theWorldLengthXY, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Trd("parent_3", theWorldLengthXY*0.5, theWorldLengthXY,
                             theWorldLengthXY*0.5, theWorldLengthXY, theWorldLengthXY) );
}

//--------------------------------------------------------------------------
void ExDivTesterTrd::BuildChildrenSolids()
{
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );

  theChildSolids.push_back( new G4Trd("child_1", theWorldLengthXY*0.2,
                            theWorldLengthXY*0.2, theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY) );
  theChildSolids.push_back( new G4Trd("child_2", theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.2,
                            theWorldLengthXY*0.2, theWorldLengthXY) );
  theChildSolids.push_back( new G4Trd("child_3", theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.2) );
}

