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
// Author: Christian Velten (2026)

#ifndef G4VPRIMITIVEMOLECULESCORER_HH
#define G4VPRIMITIVEMOLECULESCORER_HH

#include "G4EventManager.hh"
#include "G4MoleculeCounterManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "G4VMoleculeCounterInternalBase.hh"
#include "G4VPrimitiveScorer.hh"

class G4HCofThisEvent;

//------------------------------------------------------------------------------

// T    : G4VMoleculeCounter | G4VMoleculeReactionCounter
// U    : G4VMoleculeCounter::G4VMoleculeCounterIndex |
//        G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex
// V... : List of 'simple' types like G4String, G4int

template<class T, typename U, typename... V>
class G4VPrimitiveMoleculeScorer : public G4VPrimitiveScorer, public G4UImessenger
{
    static_assert(std::is_base_of<G4VMoleculeCounterInternalBase, T>::value,
                  "T must be derived from G4VMoleculeCounter or G4VMoleculeReactionCounter! "
                  "No forward declaration is allowed.");
    static_assert(
      std::is_base_of<G4VMoleculeCounterInternalBase::G4VMoleculeCounterIndexInterface, U>::value,
      "U must be derived from G4VMoleculeCounter::G4VMoleculeCounterIndex or "
      "from G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex! "
      "No forward declaration is allowed.");

  public:
    G4VPrimitiveMoleculeScorer(G4String scorerName, G4String counterName, G4int depth = 0);
    ~G4VPrimitiveMoleculeScorer() override = default;

    void AbsorbResultsFromWorkerScorer(G4VPrimitiveMoleculeScorer*);

    virtual std::tuple<V...> ConvertCounterIndexToKey(const U&) const = 0;
    virtual void WriteWithAnalysisManager() const = 0;

  private:
    G4String fCounterName;

  protected:
    G4int fNbOfScoredEvents{0};
    std::set<G4double> fTimesToRecord{};
    std::map<std::tuple<V...>, InnerCounterMapType> fCountPerIndexPerTime;

  public:
    G4String GetCounterName() const;
    const std::map<std::tuple<V...>, InnerCounterMapType>& GetCountMap() const;

    void AddTimeToRecord(G4double);
    void ClearTimesToRecord();
    G4int GetNumberOfRecordedEvents() const;

    // G4VPrimitiveScorer
  public:
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    void EndOfEvent(G4HCofThisEvent*) override;
    void clear() override;

    // G4UImessenger
  public:
    void SetNewValue(G4UIcommand*, G4String) override;

  protected:
    const char* ScorerCommand(G4String) const;
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fCmdAddTimeToRecord;
    std::unique_ptr<G4UIcmdWithAnInteger> fCmdSetNbOfTimeBins;
};

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
G4VPrimitiveMoleculeScorer<T, U, V...>::G4VPrimitiveMoleculeScorer(G4String scorerName,
                                                                     G4String counterName,
                                                                     G4int depth)
  : G4VPrimitiveScorer(scorerName, depth), fCounterName(counterName)
{
  fCmdAddTimeToRecord =
    std::make_unique<G4UIcmdWithADoubleAndUnit>(ScorerCommand("addTimeToRecord"), this);
  fCmdSetNbOfTimeBins =
    std::make_unique<G4UIcmdWithAnInteger>(ScorerCommand("setNbOfTimeBins"), this);
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
void G4VPrimitiveMoleculeScorer<T, U, V...>::AbsorbResultsFromWorkerScorer(
  G4VPrimitiveMoleculeScorer* worker)
{
  // NB: If called from G4Run::Merge, this ends up being called on the worker threads
  //     with `this == master` and `worker == worker`.
  if (worker == nullptr)
  {
    G4Exception("G4VPrimitiveMoleculeScorer<T, U, V...>::AbsorbResultsFromWorkerScorer",
                "BAD_REFERENCE", FatalErrorInArgument, "The method was called with a nullptr!");
  }
  if (worker == this)
  {
    G4Exception("G4VPrimitiveMoleculeScorer<T, U, V...>::AbsorbResultsFromWorkerScorer",
                "WRONG_IMPL", FatalException, "WORKER == this. Cannot absorb from self!");
  }
  if (G4MoleculeCounterManager::Instance()->GetAccumulateCounterIntoMaster())
  {
    G4ExceptionDescription description;
    description
      << "The `" << GetName() << "` scorer is asked to merge results from workers, while "
      << "the G4MoleculeCounterManager is also set to accumulate counter results into master! "
      << "This creates an inconsistency and risks counting molecules more than once!";
    G4Exception("G4VPrimitiveMoleculeScorer<T, U, V...>::AbsorbResultsFromWorkerScorer", "MOLMAN",
                JustWarning, description);
  }

  fNbOfScoredEvents += worker->fNbOfScoredEvents;

  for (auto const& worker_it : worker->GetCountMap())
  {
    auto [master_it, indexIsNew] =
      fCountPerIndexPerTime.emplace(worker_it.first, worker_it.second);

    if (!indexIsNew)
    {  // worker_it.first already exists, emplacedPair.first == it
      for (auto const& inner_worker_it : worker_it.second)
      {
        auto inner_master_it = master_it->second.find(inner_worker_it.first);
        inner_master_it->second += inner_worker_it.second;
      }
    }
  }

  worker->clear();
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
G4bool G4VPrimitiveMoleculeScorer<T, U, V...>::ProcessHits(G4Step*, G4TouchableHistory*)
{
  return true;
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
void G4VPrimitiveMoleculeScorer<T, U, V...>::EndOfEvent(G4HCofThisEvent*)
{
  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) return;

  T* counter = nullptr;

  auto moleculeCounters = G4MoleculeCounterManager::Instance()->GetMoleculeCounters(fCounterName);
  if (moleculeCounters.rbegin() != moleculeCounters.rend())
  {
    counter = const_cast<T*>(dynamic_cast<const T*>(*moleculeCounters.begin()));
    if (counter == nullptr)
      G4Exception("G4VPrimitiveMoleculeScorer<T>::EndOfEvent", "SCOREMOLCOUNT", FatalException,
                  "Molecule counter has wrong type!");
  }
  if (counter == nullptr)  // fCounterName is not a molecule counter
  {
    auto reactionCounters = G4MoleculeCounterManager::Instance()->GetMoleculeReactionCounters(fCounterName);
    if (reactionCounters.rbegin() != reactionCounters.rend()) {
      counter = const_cast<T*>(dynamic_cast<const T*>(*reactionCounters.begin()));
      if (counter == nullptr)
        G4Exception("G4VPrimitiveMoleculeScorer<T>::EndOfEvent", "SCOREMOLCOUNT", FatalException,
                    "Molecule reaction counter has wrong type!");
    }
  }
  if (counter == nullptr)  // fCounterName is also not a reaction counter
  {
    G4Exception("G4VPrimitiveMoleculeScorer<T>::EndOfEvent", "SCOREMOLCOUNT", FatalException,
                "No molecule counter with given name found!");
  }

  for (auto const& [idx, _] : counter->GetCounterMap())
  {
    auto key = ConvertCounterIndexToKey(idx);
    for (auto const& time : fTimesToRecord)
    {
      G4int count = counter->GetCountAtIndexAndTime(&idx, time);
      if (count > -1) fCountPerIndexPerTime[key][time] += count;
    }
  }

  ++fNbOfScoredEvents;
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
void G4VPrimitiveMoleculeScorer<T, U, V...>::clear()
{
  fCountPerIndexPerTime.clear();
  fNbOfScoredEvents = 0;

  G4VPrimitiveScorer::clear();
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
inline G4String G4VPrimitiveMoleculeScorer<T, U, V...>::GetCounterName() const
{
  return fCounterName;
}

template<typename T, typename U, typename... V>
inline const std::map<std::tuple<V...>, InnerCounterMapType>&
G4VPrimitiveMoleculeScorer<T, U, V...>::GetCountMap() const
{
  return fCountPerIndexPerTime;
}

template<typename T, typename U, typename... V>
inline void G4VPrimitiveMoleculeScorer<T, U, V...>::AddTimeToRecord(G4double time)
{
  fTimesToRecord.insert(time);
}

template<typename T, typename U, typename... V>
inline void G4VPrimitiveMoleculeScorer<T, U, V...>::ClearTimesToRecord()
{
  fTimesToRecord.clear();
}

template<typename T, typename U, typename... V>
inline G4int G4VPrimitiveMoleculeScorer<T, U, V...>::GetNumberOfRecordedEvents() const
{
  return fNbOfScoredEvents;
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename... V>
void G4VPrimitiveMoleculeScorer<T, U, V...>::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fCmdAddTimeToRecord.get())
  {
    G4double cmdTime = fCmdAddTimeToRecord->GetNewDoubleValue(newValue);
    AddTimeToRecord(cmdTime);
  }
  if (command == fCmdSetNbOfTimeBins.get())
  {
    ClearTimesToRecord();
    G4int cmdBins = fCmdSetNbOfTimeBins->GetNewIntValue(newValue);
    G4double timeMin = 1 * ps;
    G4double timeMax = G4Scheduler::Instance()->GetEndTime();
    G4double timeLogMin = std::log10(timeMin);
    G4double timeLogMax = std::log10(timeMax);
    for (G4int i = 0; i < cmdBins; i++)
    {
      AddTimeToRecord(std::pow(10, timeLogMin + i * (timeLogMax - timeLogMin) / (cmdBins - 1)));
    }
  }
}

template<typename T, typename U, typename... V>
const char* G4VPrimitiveMoleculeScorer<T, U, V...>::ScorerCommand(G4String command) const
{
  command = "/scorer/" + GetName() + "/" + command;
  return command.c_str();
}

//------------------------------------------------------------------------------

#endif
