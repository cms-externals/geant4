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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4DecayStrongResonances
//
// Modified:
// 02.11.2010 V.Ivanchenko moved constructor and destructor to source
// 07.27.2011 M.Kelsey -- Use new decay utility to process input list

#include "G4DecayStrongResonances.hh"

#include "G4DecayKineticTracks.hh"

G4ReactionProductVector* G4DecayStrongResonances::Propagate(G4KineticTrackVector* theSecondaries,
                                                            G4V3DNucleus*)
{
  G4DecayKineticTracks decay(theSecondaries);  // Changes input list in situ

  // translate to ReactionProducts
  G4ReactionProductVector* theResult = new G4ReactionProductVector;
  if (nullptr == theSecondaries || theSecondaries->empty())
  {
    return theResult;
  }

  theResult->reserve(theSecondaries->size());
  for (auto& aSecondary : *theSecondaries)
  {
    if (nullptr == aSecondary) continue;  // Skip null pointers
    auto it = new G4ReactionProduct();

    it->SetDefinition(aSecondary->GetDefinition());
    it->SetMass(aSecondary->GetDefinition()->GetPDGMass());
    it->SetTotalEnergy(aSecondary->Get4Momentum().t());
    it->SetMomentum(aSecondary->Get4Momentum().vect());
    it->SetCreatorModelID(aSecondary->GetCreatorModelID());
    it->SetParentResonanceDef(aSecondary->GetParentResonanceDef());
    it->SetParentResonanceID(aSecondary->GetParentResonanceID());
    delete aSecondary;
    theResult->push_back(it);
  }
  delete theSecondaries;

  return theResult;
}
