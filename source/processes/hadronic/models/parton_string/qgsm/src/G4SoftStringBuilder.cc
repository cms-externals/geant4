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
#include "G4SoftStringBuilder.hh"
//#define debug_G4SoftStringBuilder

//***************************************************************************************************

G4SoftStringBuilder::G4SoftStringBuilder()
{
}

G4SoftStringBuilder::~G4SoftStringBuilder()
{
}

//***************************************************************************************************

G4ExcitedString* G4SoftStringBuilder::BuildString(G4PartonPair * aPair)       
{
  #ifdef debug_G4SoftStringBuilder
  G4cout << "BuildSoft " << aPair->GetParton1()->GetPDGcode() << " "
			 << aPair->GetParton1()->Get4Momentum() << " "
			 << aPair->GetParton1()->GetX() << " " << G4endl;
  G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " " 
			<< aPair->GetParton2()->Get4Momentum() << " "
			<< aPair->GetParton2()->GetX() << " " << G4endl;
  #endif

  return new G4ExcitedString(aPair->GetParton1(), aPair->GetParton2(), aPair->GetDirection());
}

//***********************************************************************************************

