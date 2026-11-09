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
// class ExDivTesterPara Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterPara.hh"

#include "G4PVPlacement.hh"
#include "G4Para.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <fstream>

//--------------------------------------------------------------------------
ExDivTesterPara::ExDivTesterPara(PVType& pvtype, PlaceType& postype,
                                 std::vector<G4String>& extraPars)
  : ExVDivTester(pvtype, postype, extraPars)
{
  //----- Get the axis of division
  theAxis.push_back(kXAxis);
  theAxis.push_back(kYAxis);
  theAxis.push_back(kZAxis);
}

//--------------------------------------------------------------------------
void ExDivTesterPara::GenerateScanPoints()
{
  GenerateScanPointsAsBox();
}

//--------------------------------------------------------------------------
void ExDivTesterPara::BuildParentSolids()
{
  // Cartesian replicas use axis-aligned planes to locate their copies and
  // therefore cannot represent a skewed para.  Divisions use a dedicated
  // parameterisation and retain the general skewed geometry.
  G4double alpha = (thePVType == pvReplica) ? 0. : 30. * deg;
  G4double theta = (thePVType == pvReplica) ? 0. : 45. * deg;
  G4double phi = (thePVType == pvReplica) ? 0. : 60. * deg;
  theParentSolids.push_back(new G4Para("parent_1", theWorldLengthXY * 0.4, theWorldLengthXY * 0.5,
                                       theWorldLengthXY * 0.6, alpha, theta, phi));
  theParentSolids.push_back(new G4Para("parent_2", theWorldLengthXY * 0.4, theWorldLengthXY * 0.5,
                                       theWorldLengthXY * 0.6, alpha, theta, phi));
  theParentSolids.push_back(new G4Para("parent_3", theWorldLengthXY * 0.4, theWorldLengthXY * 0.5,
                                       theWorldLengthXY * 0.6, alpha, theta, phi));
}

//--------------------------------------------------------------------------
void ExDivTesterPara::BuildChildrenSolids()
{
  G4double alpha = (thePVType == pvReplica) ? 0. : 30. * deg;
  G4double theta = (thePVType == pvReplica) ? 0. : 45. * deg;
  G4double phi = (thePVType == pvReplica) ? 0. : 60. * deg;
  theChildSolids.push_back(new G4Para("child_1", theWorldLengthXY * 0.4 / theNDiv,
                                      theWorldLengthXY * 0.5, theWorldLengthXY * 0.6, alpha, theta,
                                      phi));
  theChildSolids.push_back(new G4Para("child_2", theWorldLengthXY * 0.4,
                                      theWorldLengthXY * 0.5 / theNDiv, theWorldLengthXY * 0.6,
                                      alpha, theta, phi));
  theChildSolids.push_back(new G4Para("child_3", theWorldLengthXY * 0.4, theWorldLengthXY * 0.5,
                                      theWorldLengthXY * 0.6 / theNDiv, alpha, theta, phi));
  theWidths.push_back(2 * ((G4Para*)theParentSolids[0])->GetXHalfLength() / theNDiv);
  theWidths.push_back(2 * ((G4Para*)theParentSolids[1])->GetYHalfLength() / theNDiv);
  theWidths.push_back(2 * ((G4Para*)theParentSolids[2])->GetZHalfLength() / theNDiv);
}
