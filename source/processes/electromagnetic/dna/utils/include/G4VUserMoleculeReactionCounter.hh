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
// Author: Christian Velten (2025)

#ifndef G4VUSERMOLECULEREACTIONCOUNTER_HH
#define G4VUSERMOLECULEREACTIONCOUNTER_HH 1

#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounterTemplates.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4VMoleculeReactionCounter.hh"

//------------------------------------------------------------------------------

template<class TIndex>
class G4VUserMoleculeReactionCounter : public G4VMoleculeReactionCounter
{
    static_assert(std::is_base_of<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex, TIndex>::value,
                  "TIndex must be derived from G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex! "
                  "No forward declaration is allowed.");

  public:
    G4VUserMoleculeReactionCounter();
    G4VUserMoleculeReactionCounter(G4String, MoleculeReactionCounterType = MoleculeReactionCounterType::Basic);
    virtual ~G4VUserMoleculeReactionCounter() = default;

  public:
    void Initialize() final;
    void InitializeUser() override = 0;
    void ResetCounter() override;
    void Dump() const override;
    void DumpCounterMapIndices() const override;

	void AbsorbCounter(std::weak_ptr<G4VMoleculeCounterInternalBase>) override;
    //no needed
    //std::shared_ptr<G4VMoleculeReactionCounterIndex> BuildIndex(const G4Track*, const G4Track*, const G4DNAMolecularReactionData*) const override = 0;
    std::shared_ptr<G4VMoleculeReactionCounterIndex> BuildSimpleIndex(const G4DNAMolecularReactionData*) const override = 0;

    void RecordReaction(std::shared_ptr<G4VMoleculeReactionCounterIndex>, G4double, G4int = 1) override;

    std::set<const G4DNAMolecularReactionData*> GetRecordedReactions() const override;
    std::set<G4double> GetRecordedTimes() const override;

  protected:
    static std::shared_ptr<TIndex> ConvertBaseIndexPointer(std::weak_ptr<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex>);

    std::map<TIndex, InnerCounterMapType> fCounterMap{};

  public:
    const std::map<TIndex, InnerCounterMapType>& GetCounterMap() const { return fCounterMap; }
    std::vector<TIndex> GetMapIndices() const;

    virtual G4int GetNbReactionsAtTime(const TIndex&, G4double);

    //-SEARCH-----------------------------------------------------------------------
  protected:
    struct Search
    {
        Search() : fLowerBoundSet(false) {}
        typename std::map<TIndex, InnerCounterMapType>::const_iterator fLastIndexSearched;
        InnerCounterMapType::const_iterator fLowerBoundTime;
        G4bool fLowerBoundSet;
    };
    std::unique_ptr<Search> fpLastSearch{};
    G4bool SearchIndexUpdated(const TIndex&);
    G4int SearchUpperBoundTime(G4double, G4bool);
};

//------------------------------------------------------------------------------

// #include "G4VUserMoleculeReactionCounter.icc"

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeReactionCounter<T>::G4VUserMoleculeReactionCounter() : G4VMoleculeReactionCounter()
{}

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeReactionCounter<T>::G4VUserMoleculeReactionCounter(G4String name, MoleculeReactionCounterType type)
  : G4VMoleculeReactionCounter(name, type)
{}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::Initialize()
{
  InitializeUser();
  fIsInitialized = true;
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4int G4VUserMoleculeReactionCounter<TIndex>::GetNbReactionsAtTime(const TIndex& index, G4double time)
{
  G4bool sameIndex = !SearchIndexUpdated(index);
  return SearchUpperBoundTime(time, sameIndex);
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::shared_ptr<TIndex> G4VUserMoleculeReactionCounter<TIndex>::ConvertBaseIndexPointer(
  std::weak_ptr<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex> pIndex)
{
  if (pIndex.expired()) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type " << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << "!\n"
           << "Because the weak_ptr has expired since it has been passed, indicating something "
              "very wrong!"
           << G4endl;
    G4Exception(
      G4String("G4VUserMoleculeReactionCounter<" + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::ConvertBaseIndexPointer"),
      "BAD_REFERENCE", FatalException, errMsg);
  }

  std::shared_ptr<TIndex> mapIndex = std::dynamic_pointer_cast<TIndex>(pIndex.lock());
  if (mapIndex == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type " << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << "!" << G4endl;
    G4Exception(
      G4String("G4VUserMoleculeReactionCounter<" + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::ConvertBaseIndexPointer"),
      "BAD_REFERENCE", FatalException, errMsg);
  }
  return mapIndex;
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::RecordReaction(
  std::shared_ptr<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex> pIndex,
  G4double time, G4int number)
{
  std::shared_ptr<TIndex> mapIndex = ConvertBaseIndexPointer(pIndex);

  if (IsTimeAboveUpperBound(time)) {
    if (fVerbose > 3) {
      G4cout << "G4VUserMoleculeReactionCounter<"
             << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
             << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
             << " at time : " << G4BestUnit(time, "Time") << G4endl;
      G4cout << ":: [IsTimeAboveUpperBound] Skipping since IsTimeAboveUpperBound == true" << G4endl;
    }
    return;
  }
  else if (IsTimeBelowLowerBound(time)) {
    if (fVerbose > 3) {
      G4cout << "G4VUserMoleculeReactionCounter<"
             << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
             << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
             << " at time : " << G4BestUnit(time, "Time") << G4endl;
      G4cout << ":: [IsTimeBelowLowerBound] Skipping since IsTimeBelowLowerBound == true" << G4endl;
    }
    return;
  }

  if (fVerbose > 2) {
    G4cout << "G4VUserMoleculeReactionCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">(" << GetName()
           << ")::RecordReaction : " << mapIndex->GetReactionData()->GetReactionID()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  auto [it, indexIsNew] = fCounterMap.emplace(*mapIndex, InnerCounterMapType{fTimeComparer});

  if (indexIsNew || it->second.empty()) {
    it->second.emplace(time, number);
    // it->second[time] = number;
  }
  else {
	  if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent())
         // can only do consistency check if the counters are cleared before each event
    {
      auto end = it->second.rbegin();

      if ((end->first <= time
              || fabs(end->first - time) <= fTimeComparer.GetPrecisionAtTime(time)))
      // Case 1 = new time comes after last recorded data
      // Case 2 = new time is about the same as the last recorded one
      {
        // it->second[time] = end->second + number;
        auto [it_time, _] = it->second.emplace(time, end->second);
        it_time->second += number;
      }
      else {
        G4ExceptionDescription errMsg;
        errMsg << "Time of reaction " << mapIndex->GetReactionData()->GetReactionID() << " is "
               << G4BestUnit(time, "Time") << "while the global time is "
               << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
               << "(last counter time: " << G4BestUnit(end->first, "Time") << ")" << G4endl;
        G4Exception(G4String("G4VUserMoleculeReactionCounter<"
                             + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                             + ">::RecordReaction"),
                    "TIME_DONT_MATCH", FatalException, errMsg);
      }
    }
    else {
      // since counters are not cleared after chemical run (i.e., after event)
      // there will already be numbers in the map, so...
      // (1) find the closest time
      // (2) emplace entry using closest value as init + number
      // (3) add number to all "future" entries as well
      auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(it->second, time);
      auto [it_new, _] = it->second.emplace(time, it_closest->second);
      do {
        it_new->second += number;
      } while (++it_new != it->second.end());
    }
  }
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::vector<TIndex> G4VUserMoleculeReactionCounter<TIndex>::GetMapIndices() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4VUserMoleculeReactionCounter::GetMapIndices" << G4endl;
  }
  return G4::MoleculeCounter::GetMapIndices(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
std::set<const G4DNAMolecularReactionData*> G4VUserMoleculeReactionCounter<T>::GetRecordedReactions() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4VUserMoleculeReactionCounter::GetRecordedReactions" << G4endl;
  }
  std::set<const G4DNAMolecularReactionData*> output{};
  for (const auto& it : fCounterMap) {
    output.insert(it.first.GetReactionData());
  }
  return output;
}

//------------------------------------------------------------------------------

template<typename T>
std::set<G4double> G4VUserMoleculeReactionCounter<T>::GetRecordedTimes() const
{
  return G4::MoleculeCounter::GetRecordedTimes<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeReactionCounter<T>::Dump() const
{
  DumpCounterMapIndices();
  G4::MoleculeCounter::DumpCounterMapContents<T>(fCounterMap);
}

template<typename T>
void G4VUserMoleculeReactionCounter<T>::DumpCounterMapIndices() const
{
  G4::MoleculeCounter::DumpCounterMapIndices<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeReactionCounter<T>::ResetCounter()
{
  if (fVerbose > 1) {
    G4cout << "G4VUserMoleculeReactionCounter<" << G4::MoleculeCounter::GetTemplateTypeName<T>() << ">(" << GetName()
           << ")::ResetCounter" << G4endl;
  }
  fCounterMap.clear();
  fpLastSearch.reset();
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4bool G4VUserMoleculeReactionCounter<TIndex>::SearchIndexUpdated(const TIndex& index)
{
  if (fpLastSearch == nullptr) {
    fpLastSearch = std::make_unique<Search>();
  }
  else {
    if (fpLastSearch->fLowerBoundSet && !(fpLastSearch->fLastIndexSearched->first < index)
        && !(index < fpLastSearch->fLastIndexSearched->first))
    {
      return true;
    }
  }

  auto mol_it = fCounterMap.find(index);
  fpLastSearch->fLastIndexSearched = mol_it;

  if (mol_it != fCounterMap.end()) {
    fpLastSearch->fLowerBoundTime = fpLastSearch->fLastIndexSearched->second.end();
    fpLastSearch->fLowerBoundSet = true;
  }
  else {
    fpLastSearch->fLowerBoundSet = false;
  }

  return false;
}

//------------------------------------------------------------------------------

template<typename T>
G4int G4VUserMoleculeReactionCounter<T>::SearchUpperBoundTime(G4double time, G4bool sameIndex)
{
  auto mol_it = fpLastSearch->fLastIndexSearched;
  if (mol_it == fCounterMap.end()) {
    return 0;
  }

  InnerCounterMapType const& timeMap = mol_it->second;
  if (timeMap.empty()) {
    return 0;
  }

  if (sameIndex) {
    if (fpLastSearch->fLowerBoundSet && fpLastSearch->fLowerBoundTime != timeMap.end()) {
      if (fpLastSearch->fLowerBoundTime->first < time) {
        auto upperToLast = fpLastSearch->fLowerBoundTime;
        upperToLast++;

        if (upperToLast == timeMap.end()) {
          return fpLastSearch->fLowerBoundTime->second;
        }

        if (upperToLast->first > time) {
          return fpLastSearch->fLowerBoundTime->second;
        }
      }
    }
  }

  auto up_time_it = timeMap.upper_bound(time);

  if (up_time_it == timeMap.end()) {
    auto last_time = timeMap.rbegin();
    return last_time->second;
  }
  if (up_time_it == timeMap.begin()) {
    return 0;
  }

  up_time_it--;

  fpLastSearch->fLowerBoundTime = up_time_it;
  fpLastSearch->fLowerBoundSet = true;

  return fpLastSearch->fLowerBoundTime->second;
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeReactionCounter<TIndex>::AbsorbCounter(
  std::weak_ptr<G4VMoleculeCounterInternalBase> wpCounter)
{
  if (wpCounter.expired()) {
	G4ExceptionDescription errMsg;
	errMsg << "Could not cast the pointer to type G4VUserMoleculeReactionCounter<"
		   << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
		   << "Because the weak_ptr has expired since it has been passed, indicating something "
			  "very wrong!"
		   << G4endl;
	G4Exception(G4String("G4VUserMoleculeReactionCounter<"
						 + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
				"BAD_REFERENCE", FatalException, errMsg);
  }

  std::shared_ptr<G4VUserMoleculeReactionCounter<TIndex>> pCounter =
	std::dynamic_pointer_cast<G4VUserMoleculeReactionCounter<TIndex>>(wpCounter.lock());

  if (pCounter == nullptr) {
	G4ExceptionDescription errMsg;
	errMsg << "Could not cast the pointer to type G4VUserMoleculeReactionCounter<"
		   << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
		   << "Because the objects aren't of the same type!" << G4endl;
	G4Exception(G4String("G4VUserMoleculeReactionCounter<"
						 + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
				"BAD_REFERENCE", FatalException, errMsg);
  }

  if (pCounter->GetType() != GetType()) {
	G4ExceptionDescription errMsg;
	errMsg << "You are trying to absorb a counter with different type!" << G4endl;
	G4Exception(G4String("G4VUserMoleculeReactionCounter<"
						 + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
				"TYPE_DIFF", JustWarning, errMsg);
  }

  for (auto const& worker_it : pCounter->GetCounterMap()) {
	auto [master_it, indexIsNew] =
	  fCounterMap.emplace(worker_it.first, InnerCounterMapType{fTimeComparer});

	G4int currentNumber = 0, previousNumber = 0;
	for (auto const& [time, number] : worker_it.second) {
	  currentNumber = number - previousNumber;
	  previousNumber = number;

	  if (master_it->second.empty()) {
		master_it->second.emplace(time, currentNumber);
	  }
	  else {  // at least one element exists, so we can try to find the closest key
		auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(master_it->second, time);
		auto [it, _] = master_it->second.emplace(time, it_closest->second);
		do {
		  it->second += currentNumber;
		} while (++it != master_it->second.end());
	  }
	}
  }
}

//------------------------------------------------------------------------------

#endif
