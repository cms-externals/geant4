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
#include "TimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ITTrackHolder.hh"
#include "G4MoleculeTable.hh"
#include "G4MoleculeCounter.hh"
#include "G4Scheduler.hh"
#include "G4MolecularConfiguration.hh"

TimeStepAction::TimeStepAction() :
    G4UserTimeStepAction()
{
  /**
   * Give to G4ITStepManager the user defined time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   * Those time steps are used for the chemistry of G4DNA
   */

  AddTimeStep(1 * picosecond, 0.1 * picosecond);
  AddTimeStep(10 * picosecond, 1 * picosecond);
  AddTimeStep(100 * picosecond, 3 * picosecond);
  AddTimeStep(1000 * picosecond, 10 * picosecond);
  AddTimeStep(10000 * picosecond, 100 * picosecond);

}

TimeStepAction::~TimeStepAction()
{
}

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
    G4UserTimeStepAction(other)
{
}

TimeStepAction& TimeStepAction::operator=(const TimeStepAction& rhs)
{
  if(this == &rhs) return *this;
  return *this;
}

void TimeStepAction::UserPostTimeStepAction()
{
  G4ConfigurationIterator speciesIt = G4MoleculeTable::Instance()
      ->GetConfigurationIterator();

  speciesIt.reset();
  double time = G4Scheduler::Instance()->GetGlobalTime();

  while(speciesIt())
  {
    const G4MolecularConfiguration* speciesType = speciesIt.value();

    auto weak_counter = G4MoleculeCounterManager::Instance()->GetMoleculeCounter(0);
    auto counter = std::dynamic_pointer_cast<G4MoleculeCounter>(weak_counter.lock());
    if (counter == nullptr) {
      G4Exception("ScoreSpecies::EndOfEvent", "BAD_REFERENCE", FatalException,
                  "The molecule counter could not be received!");
    }
    if(!counter->IsReactantIgnored(speciesType->GetDefinition()))
    {

    G4TrackList* trackList =
        G4ITTrackHolder::Instance()->GetMainList(speciesType->GetMoleculeID());

    auto index = counter->BuildSimpleIndex(speciesType);

      auto MoleculeIndex = dynamic_cast<G4MoleculeCounterIndex*>(index.get());
      if (!MoleculeIndex)
      {
        G4Exception("TimeStepAction::UserPostTimeStepAction",
                    "INVALID_CAST",
                    JustWarning,
                    "Failed to cast index to G4MoleculeCounterIndex.");
        return;
      }

      int number = counter->GetNbMoleculesAtTime(*MoleculeIndex,
                                                                    time);

    int trackListSize = 0;

    if(trackList)
    {
      trackListSize = (int) trackList->size();
    }
    if(trackListSize != number)
    {
      G4cout << "Time asked = " << G4BestUnit(time, "Time") << G4endl;
      counter->Dump();
      G4ExceptionDescription errMsg;
      errMsg << "There seems to be a problem related to G4MoleculeCounter" << G4endl
             << " Species = " <<  speciesType->GetName()
             << " n tracks = " << trackListSize
             << " numberFromCounter = " << number
             << G4endl;
      G4Exception("TimeStepAction::UserPostTimeStepAction",
                  "CHECK_COUNTER_CONSISTENCY",
                  FatalException,
                  errMsg);
    }
    }
  };
}

void TimeStepAction::EndProcessing()
{
  G4cout << " --> TimeStepAction::EndProcessing" << G4endl;
}
