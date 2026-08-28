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

#include "Tst69RunAction.hh"

#include "G4AblaInterface.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLXXInterface.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4Run.hh"
#include "Randomize.hh"

#include "Tst69INCLXXTallyAnalysis.hh"

#include <vector>

Tst69RunAction::Tst69RunAction(const char* const physList)
{
  theTally = new Tst69INCLXXTallyAnalysis(physList);
}

Tst69RunAction::~Tst69RunAction()
{
  delete theTally;
}

void Tst69RunAction::BeginOfRunAction(const G4Run*)
{
  if (std::getenv("TEST69_USE_ABLA"))
  {
    std::vector<G4HadronicInteraction*> interactions =
      G4HadronicInteractionRegistry::Instance()->FindAllModels(
        G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName());
    for (std::vector<G4HadronicInteraction*>::const_iterator iInter = interactions.begin(),
                                                             e = interactions.end();
         iInter != e; ++iInter)
    {
      G4INCLXXInterface* theINCLInterface = static_cast<G4INCLXXInterface*>(*iInter);
      if (theINCLInterface)
      {
        G4HadronicInteraction* interaction =
          G4HadronicInteractionRegistry::Instance()->FindModel("ABLA");
        G4AblaInterface* theAblaInterface = static_cast<G4AblaInterface*>(interaction);
        if (!theAblaInterface) theAblaInterface = new G4AblaInterface;
        G4cout << "Coupling INCLXX to ABLA" << G4endl;
        theINCLInterface->SetDeExcitation(theAblaInterface);
      }
    }
  }

  // set the INCL++ tally object if necessary
  G4INCLXXInterfaceStore* theStore = G4INCLXXInterfaceStore::GetInstance();
  G4cout << "Activating tally for INCLXX" << G4endl;
  theStore->SetTally(theTally);
  theTally->Open();
}

void Tst69RunAction::EndOfRunAction(const G4Run*)
{
  theTally->Close();
}
