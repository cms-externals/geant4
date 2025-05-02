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

#ifndef G4VUSERMOLECULECOUNTER_HH
#define G4VUSERMOLECULECOUNTER_HH 1

#include "G4MoleculeCounterManager.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeCounterTemplates.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "G4VMoleculeCounter.hh"

//------------------------------------------------------------------------------

template<class TIndex>
class G4VUserMoleculeCounter : public G4VMoleculeCounter
{
    static_assert(std::is_base_of<G4VMoleculeCounter::G4VMoleculeCounterIndex, TIndex>::value,
                  "TIndex must be derived from G4VMoleculeCounter::G4VMoleculeCounterIndex! "
                  "No forward declaration is allowed.");

  public:
    G4VUserMoleculeCounter();
    G4VUserMoleculeCounter(G4String, MoleculeCounterType = MoleculeCounterType::Other);
    ~G4VUserMoleculeCounter() override = default;

  public:
    void Initialize() final;
    void InitializeUser() override = 0;
    void ResetCounter() override;
    void Dump() const override;
    void DumpCounterMapIndices() const override;

    void AbsorbCounter(std::weak_ptr<G4VMoleculeCounterInternalBase>) override;

    std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
    BuildIndex(const G4Track*) const override = 0;
    std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
    BuildIndex(const G4Track*, const G4StepPoint*) const override = 0;
    std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
    BuildSimpleIndex(const G4MolecularConfiguration*) const override = 0;

    void AddMolecule(std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>, G4double,
                     G4int = 1) override;
    void RemoveMolecule(std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>, G4double,
                        G4int = 1) override;

    std::set<const G4MolecularConfiguration*> GetRecordedMolecules() const override;
    std::set<G4double> GetRecordedTimes() const override;

    void SchedulerFinalizedTracking() override;

  protected:
    static std::shared_ptr<TIndex>
      ConvertBaseIndexPointer(std::weak_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>);

    std::map<TIndex, InnerCounterMapType> fCounterMap{};
    std::map<TIndex, G4int> fShadowCounterMap{};

  public:
    const std::map<TIndex, InnerCounterMapType>& GetCounterMap() const { return fCounterMap; }
    std::vector<TIndex> GetMapIndices() const;

    virtual G4int GetNbMoleculesAtTime(const TIndex&, G4double);

    virtual G4bool IsReactantIgnored(const G4MoleculeDefinition* molecule);
    virtual G4bool IsReactantIgnored(const G4MolecularConfiguration* reactant);
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

// #include "G4VUserMoleculeCounter.icc"

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeCounter<T>::G4VUserMoleculeCounter() : G4VMoleculeCounter()
{}

//------------------------------------------------------------------------------

template<typename T>
G4VUserMoleculeCounter<T>::G4VUserMoleculeCounter(G4String name, MoleculeCounterType type)
  : G4VMoleculeCounter(name, type)
{}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeCounter<TIndex>::Initialize()
{
  InitializeUser();
  fIsInitialized = true;
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4int G4VUserMoleculeCounter<TIndex>::GetNbMoleculesAtTime(const TIndex& index, G4double time)
{
  G4bool sameIndex = !SearchIndexUpdated(index);
  return SearchUpperBoundTime(time, sameIndex);
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::shared_ptr<TIndex> G4VUserMoleculeCounter<TIndex>::ConvertBaseIndexPointer(
  std::weak_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex> pIndex)
{
  if (pIndex.expired()) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type "
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << "!\n"
           << "Because the weak_ptr has expired since it has been passed, indicating something "
              "very wrong!"
           << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                         + ">::ConvertBaseIndexPointer"),
                "BAD_REFERENCE", FatalException, errMsg);
  }

  std::shared_ptr<TIndex> mapIndex = std::dynamic_pointer_cast<TIndex>(pIndex.lock());
  if (mapIndex == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type "
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << "!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                         + ">::ConvertBaseIndexPointer"),
                "BAD_REFERENCE", FatalException, errMsg);
  }
  return mapIndex;
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeCounter<TIndex>::AddMolecule(
  std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex> pIndex, G4double time, G4int number)
{
  std::shared_ptr<TIndex> mapIndex = ConvertBaseIndexPointer(pIndex);

  if (G4::MoleculeCounter::Contains(fIgnoredMolecules, mapIndex->GetMolecule()->GetDefinition())
      || G4::MoleculeCounter::Contains(fIgnoredReactants, mapIndex->GetMolecule()))
  {
    return;
  }

  if (fCheckTimeIsConsistentWithScheduler && G4Scheduler::Instance()->IsRunning()
      && fabs(time - G4Scheduler::Instance()->GetGlobalTime())
           > G4Scheduler::Instance()->GetTimeTolerance())
  {
    G4ExceptionDescription errMsg;
    errMsg << "Time of species " << mapIndex->GetMolecule()->GetName() << " is "
           << G4BestUnit(time, "Time") << "while the global time is "
           << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time") << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AddMolecule"),
                "TIME_DONT_MATCH", FatalException, errMsg);
  }

  if (IsTimeAboveUpperBound(time)) {
    if (fVerbose > 3) {
      G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
             << ">(" << GetName() << ")::AddMolecule : " << mapIndex->GetMolecule()->GetName()
             << " at time : " << G4BestUnit(time, "Time") << G4endl;
      G4cout << ":: [IsTimeAboveUpperBound] Skipping since IsTimeAboveUpperBound == true" << G4endl;
    }
    return;
  }
  else if (IsTimeBelowLowerBound(time)) {
    // put into shadow counter
    auto [it, indexIsNew] = fShadowCounterMap.emplace(*mapIndex, number);
    if (!indexIsNew) it->second += number;
    if (fVerbose > 3) {
      G4cout << ":: [IsTimeBelowLowerBound] Adding " << mapIndex->GetInfo()
             << " shadow count: " << it->second - number << " + " << number << G4endl;
    }
    return;
  }

  if (fVerbose > 1) {
    G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
           << ">(" << GetName() << ")::AddMolecule : " << mapIndex->GetMolecule()->GetName()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  // within time bounds && not-ignored molecule
  // -> continue

  auto it_shadow = fShadowCounterMap.find(*mapIndex);
  auto [it, indexIsNew] = fCounterMap.emplace(*mapIndex, InnerCounterMapType{fTimeComparer});

  if (it_shadow != fShadowCounterMap.end()) {
    // entry found in shadow counter
    if (indexIsNew) {
      // mapIndex is new, initialize with shadow counter
      it->second[fActiveLowerBound] = it_shadow->second;
    }
    else {
      // mapIndex existed, we need to add the shadow count to it
      // this happens if we are in a subsequent event and have just crossed over the lower activity
      // bound the counter has then already entries from the previous event
      InnerCounterMapType::iterator it_time;
      G4bool timeIsNew;
      std::tie(it_time, timeIsNew) = it->second.emplace(fActiveLowerBound, 0);
      do {
        it_time->second += it_shadow->second;
      } while (++it_time != it->second.end());
    }
    // either way, remove the shadow count
    fShadowCounterMap.erase(it_shadow);
  }
  //  else if (indexIsNew)
  //  {
  //	  it->second.emplace(time, 0);
  //  }

  //  if (indexIsNew || it->second.empty()) {
  //    // Use shadow count, if it exists, then remove shadow count
  //    auto it_shadow = fShadowCounterMap.find(*mapIndex);
  //    if (it_shadow != fShadowCounterMap.end()) {
  //      it->second[time] = it_shadow->second + number;
  //      if (fVerbose > 3) {
  //        G4cout << ":: Initializing " << mapIndex->GetInfo()
  //               << " from shadow count n = " << it_shadow->second << " + " << number
  //               << " at time: " << G4BestUnit(time, "Time") << G4endl;
  //      }
  //      fShadowCounterMap.erase(it_shadow);
  //    }
  //    else {
  //      // it->second[time] = number;
  //      auto [it_time, timeIsNew] = it->second.emplace(time, number);
  // #ifdef G4VERBOSE
  //      if (fVerbose > 1 && !timeIsNew) {
  //        G4cout << "G4VUserMoleculeCounter<" <<
  //        G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
  //               << ">(" << GetName() << ")::AddMolecule : Initialized new map but time is not
  //               new?!"
  //               << G4endl;
  //      }
  // #endif
  //    }
  //  }

  // map for index existed (and was not empty)
  if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent())
  // can only do consistency check if the counters are cleared for each event (= chem run)
  {
    auto end = it->second.rbegin();
    auto init_n = end == it->second.rend() ? 0 : end->second;

    auto [it_time, timeIsNew] = it->second.emplace(time, init_n);
    it_time->second += number;

    if (fCheckRecordedTimesAreConsistent
        && !(end->first <= time
             || fabs(end->first - time) <= fTimeComparer.GetPrecisionAtTime(time)))
    // Case 1 = new time comes after last recorded data
    // Case 2 = new time is about the same as the last recorded one
    {
      G4ExceptionDescription errMsg;
      errMsg << "Time of species " << mapIndex->GetMolecule()->GetName() << " is "
             << G4BestUnit(time, "Time") << "while the global time is "
             << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
             << "(last counter time: " << G4BestUnit(end->first, "Time") << ")" << G4endl;
      G4Exception(G4String("G4VUserMoleculeCounter<"
                           + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AddMolecule"),
                  "TIME_DONT_MATCH", FatalException, errMsg);
    }
  }
  else  // counters are not (automatically) reset by manager
  {
    // since counters are not cleared after chemical run (i.e., after event)
    // there will already be numbers in the map, so...
    // (1) find the closest time
    // (2) emplace entry using closest value as init + number
    // (3) add number to all "future" entries as well
    if (it->second.empty()) {
      it->second.emplace(time, number);
    }
    else {  // at least one element exists, so we can try to find the closest key
      auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(it->second, time);
      auto [it_time, _] = it->second.emplace(time, it_closest->second);
      do {
        it_time->second += number;
      } while (++it_time != it->second.end());
    }
  }
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeCounter<TIndex>::RemoveMolecule(
  std::shared_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex> pIndex, G4double time, G4int number)
{
  std::shared_ptr<TIndex> mapIndex = ConvertBaseIndexPointer(pIndex);

  if (G4::MoleculeCounter::Contains(fIgnoredMolecules, mapIndex->GetMolecule()->GetDefinition())
      || G4::MoleculeCounter::Contains(fIgnoredReactants, mapIndex->GetMolecule()))
  {
    return;
  }

  if (fCheckTimeIsConsistentWithScheduler && G4Scheduler::Instance()->IsRunning()
      && fabs(time - G4Scheduler::Instance()->GetGlobalTime())
           > G4Scheduler::Instance()->GetTimeTolerance())
  {
    G4ExceptionDescription errMsg;
    errMsg << "Time of species " << mapIndex->GetMolecule()->GetName() << " is "
           << G4BestUnit(time, "Time") << "while the global time is "
           << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time") << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                         + ">::RemoveMolecule"),
                "TIME_DONT_MATCH", FatalException, errMsg);
  }

  if (IsTimeBelowLowerBound(time)) {
    auto it = fShadowCounterMap.find(*mapIndex);
    if (it == fShadowCounterMap.end()) {
      G4ExceptionDescription errMsg;
      errMsg << "There was no " << mapIndex->GetMolecule()->GetName()
             << " recorded at the time or even before the time asked" << G4endl;
      G4Exception(G4String("G4VUserMoleculeCounter<"
                           + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                           + ">::RemoveMolecule"),
                  "", FatalErrorInArgument, errMsg);
    }
    if (fVerbose > 3) {
      G4cout << ":: [IsTimeBelowLowerBound] Removing " << mapIndex->GetInfo()
             << " shadow count: " << it->second << " - " << number << G4endl;
    }
    it->second -= number;
    return;
  }
  else if (IsTimeAboveUpperBound(time)) {
    // if the "active" counter was not filled, add the shadow counter to it at the lower bound
    // only needed for remove since remove will always be called at the end when the molecule is
    // destroyed
    auto [it, indexIsNew] = fCounterMap.emplace(*mapIndex, InnerCounterMapType{fTimeComparer});
    if (indexIsNew || it->second.empty()) {
      auto it_shadow = fShadowCounterMap.find(*mapIndex);
      if (it_shadow != fShadowCounterMap.end()) {
        it->second[fActiveLowerBound] = it_shadow->second;
        if (fVerbose > 3) {
          G4cout << ":: [IsTimeAboveUpperBound] Set " << mapIndex->GetInfo()
                 << " Map[ActiveLowerBound] with shadow count:" << it_shadow->second << G4endl;
        }
        fShadowCounterMap.erase(it_shadow);
      }
      else if (fVerbose > 3) {
        G4cout << ":: [IsTimeAboveUpperBound] Not updating with shadow count since"
                  " no shadow count was found!"
               << G4endl;
      }
    }
    else if (fVerbose > 3) {
      G4cout << ":: [IsTimeAboveUpperBound] Not updating with shadow count"
                " since it already exists "
             << G4endl;
    }
    return;
  }

  if (fVerbose > 2) {
    G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
           << ">(" << GetName() << ")::RemoveMolecule : " << mapIndex->GetMolecule()->GetName()
           << " at time : " << G4BestUnit(time, "Time") << G4endl;
  }

  // within time bounds && not-ignored molecule
  // -> continue

  auto it_shadow = fShadowCounterMap.find(*mapIndex);
  auto it = fCounterMap.find(*mapIndex);

  if (it_shadow != fShadowCounterMap.end()) {
    // entry found in shadow counter
    if (it == fCounterMap.end()) {
      // no mapIndex found, initialize with shadow counter
      G4bool indexIsNew = false;
      std::tie(it, indexIsNew) = fCounterMap.emplace(*mapIndex, InnerCounterMapType{fTimeComparer});
      if (!indexIsNew) {
        G4ExceptionDescription errMsg;
        errMsg << "We tried to emplace the index after it was found to not exist, but now it "
                  "says it existed!?"
               << G4endl;
        G4Exception(G4String("G4VUserMoleculeCounter<"
                             + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                             + ">::RemoveMolecule"),
                    "NONSENSICAL", FatalErrorInArgument, errMsg);
      }
      it->second[fActiveLowerBound] = it_shadow->second;
    }
    else {  // it != fCounterMap.end()
      // mapIndex exists, we need to add the shadow count to it
      // this happens if we are in a subsequent event and have just crossed over the lower activity
      // bound the counter has then already entries from the previous event
      auto [it_time, _] = it->second.emplace(fActiveLowerBound, 0);
      do {
        it_time->second += it_shadow->second;
      } while (++it_time != it->second.end());
    }

    // either way, remove the shadow count
    fShadowCounterMap.erase(it_shadow);

    //	  if (it_shadow != fShadowCounterMap.end()) {
    //      G4bool indexIsNew = false;
    //      // auto [it_, indexIsNew] = fCounterMap.emplace(*mapIndex,
    //      // InnerCounterMapType{fTimeComparer});
    //      std::tie(it, indexIsNew) = fCounterMap.emplace(*mapIndex,
    //      InnerCounterMapType{fTimeComparer}); if (!indexIsNew) {
    //        G4ExceptionDescription errMsg;
    //        errMsg << "We tried to emplace the index after it was found to not exist, but now it "
    //                  "says it existed!?"
    //               << G4endl;
    //        G4Exception(G4String("G4VUserMoleculeCounter<"
    //                             + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
    //                             + ">::RemoveMolecule"),
    //                    "NONSENSICAL", FatalErrorInArgument, errMsg);
    //      }
    //      it->second[fActiveLowerBound] = it_shadow->second;
    //      // it = it_;
    //      if (fVerbose > 3) {
    //        G4cout << ":: Initialize " << mapIndex->GetInfo()
    //               << " Map[ActiveLowerBound] with shadow count:" << it_shadow->second << G4endl;
    //      }
    //      fShadowCounterMap.erase(it_shadow);
    //    }
  }

  InnerCounterMapType& nbMolPerTime = it->second;
  InnerCounterMapType::iterator it_time;
  G4bool isNewTime = false;
  G4double oldTime = 0;

  if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent())
  {
    auto end = nbMolPerTime.rbegin();  // get last entry
    oldTime = end->first;

    // CHECK: no molecules have been recorded for this index
    if (end == nbMolPerTime.rend()) {
      if (fVerbose > 2) {
        mapIndex->GetMolecule()->PrintState();
        Dump();
      }
      G4ExceptionDescription errMsg;
      errMsg << "There was no " << mapIndex->GetMolecule()->GetName()
             << " recorded at the time or even before the time asked" << G4endl;
      G4Exception("G4VUserMoleculeCounter::RemoveMolecule", "", FatalErrorInArgument, errMsg);
    }
    // CHECK: current time is less (by more than counter precision) than the most recently
    // recorded index
    if (fCheckRecordedTimesAreConsistent
        && time - end->first < -fTimeComparer.GetPrecisionAtTime(time))
    {
      if (fVerbose > 2) {
        mapIndex->GetMolecule()->PrintState();
        Dump();
      }
      G4ExceptionDescription errMsg;
      errMsg << "Is time going back?? " << mapIndex->GetMolecule()->GetName()
             << " is being removed at time " << G4BestUnit(time, "Time")
             << "while last recorded time was " << G4BestUnit(end->first, "Time") << ".";
      G4Exception("G4VUserMoleculeCounter::RemoveMolecule", "RETURN_TO_THE_FUTUR",
                  FatalErrorInArgument, errMsg);
    }

    it_time = std::prev(end.base());

    std::tie(it_time, isNewTime) = nbMolPerTime.emplace(time, it_time->second);
    // auto oldNumber = it_time->second;
    it_time->second -= number;
    // if (time > 0.001)
    //   G4cout << "(t=" << time << ")=" << number << " | old_it(t=" << it_time->first
    //          << ") = " << oldNumber << " | timeIsNew=" << isNewTime
    //          << " | new_it(t=" << it_time->first << ") = " << it_time->second << G4endl;
  }
  else {
    // since counters are not cleared after chemical run (i.e., after event)
    // there will already be numbers in the map, so...
    // (1) find the closest time
    // (2) emplace entry using closest value as init - number
    // (3) remove number from all "future" entries as well
    auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(nbMolPerTime, time);
    std::tie(it_time, isNewTime) = nbMolPerTime.emplace(time, it_closest->second);
    auto _it = it_time;
    do {
      _it->second -= number;
    } while (++_it != it->second.end());
  }

  // Check that count at new time is >= 0
  // This currently throws tons of errors for non-basic counters.
  //    auto it_time = nbMolPerTime.find(time);
  if (it_time == nbMolPerTime.end() || it_time->second < 0) {
    if (fVerbose > 2) Dump();
    G4ExceptionDescription errMsg;
    errMsg << "After removal of " << number << " species of " << mapIndex->GetMolecule()->GetName()
           << " the final number at time " << G4BestUnit(time, "Time")
           << " is less than zero and so not valid."
           << "\nIndex was :" << mapIndex->GetInfo() << "\nGlobal time is "
           << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
           << "\nPrevious selected time is " << G4BestUnit(oldTime, "Time") << G4endl;
    if (fNegativeCountsAreFatal) {
      G4Exception(G4String("G4VUserMoleculeCounter<"
                           + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                           + ">::RemoveMolecule"),
                  "N_INF_0", FatalException, errMsg);
    }
    else if (fVerbose > 0) {
      G4Exception(G4String("G4VUserMoleculeCounter<"
                           + G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
                           + ">::RemoveMolecule"),
                  "N_INF_0", JustWarning, errMsg);
    }
  }
}

//------------------------------------------------------------------------------

template<typename TIndex>
void G4VUserMoleculeCounter<TIndex>::SchedulerFinalizedTracking()
{
  // Add record to fCounterMap for each fShadowCounterMap index unless they exist already
  for (auto& it_shadow : fShadowCounterMap) {
    auto [it, indexIsNew] =
      fCounterMap.emplace(it_shadow.first, InnerCounterMapType{fTimeComparer});
    if (indexIsNew || it->second.empty()) {
      it->second[fActiveLowerBound] = it_shadow.second;
      if (fVerbose > 3) {
        G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
               << ">(" << GetName() << ")::SchedulerEndedTracking : " << "setting map index '"
               << it_shadow.first.GetInfo() << "' from shadow counter to n = " << it_shadow.second
               << G4endl;
      }
    }
    else if (fVerbose > 2) {
      G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<TIndex>()
             << ">(" << GetName() << ")::SchedulerEndedTracking : "
             << "encountered dangling shadow counter iterator for index '"
             << it_shadow.first.GetInfo() << "'" << G4endl;
    }
  }
  fShadowCounterMap.clear();
}

//------------------------------------------------------------------------------

template<typename TIndex>
std::vector<TIndex> G4VUserMoleculeCounter<TIndex>::GetMapIndices() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4VUserMoleculeCounter::GetMapIndices" << G4endl;
  }
  return G4::MoleculeCounter::GetMapIndices(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
std::set<const G4MolecularConfiguration*> G4VUserMoleculeCounter<T>::GetRecordedMolecules() const
{
  if (fVerbose > 2) {
    G4cout << "Entering in G4MoleculeCounter::RecordMolecules" << G4endl;
  }
  std::set<const G4MolecularConfiguration*> output{};
  for (const auto& it : fCounterMap) {
    output.insert(it.first.GetMolecule());
  }
  return output;
}

//------------------------------------------------------------------------------

template<typename T>
std::set<G4double> G4VUserMoleculeCounter<T>::GetRecordedTimes() const
{
  return G4::MoleculeCounter::GetRecordedTimes<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeCounter<T>::Dump() const
{
  DumpCounterMapIndices();
  G4::MoleculeCounter::DumpCounterMapContents<T>(fCounterMap);
}

template<typename T>
void G4VUserMoleculeCounter<T>::DumpCounterMapIndices() const
{
  G4::MoleculeCounter::DumpCounterMapIndices<T>(fCounterMap);
}

//------------------------------------------------------------------------------

template<typename T>
void G4VUserMoleculeCounter<T>::ResetCounter()
{
  if (fVerbose > 1) {
    G4cout << "G4VUserMoleculeCounter<" << G4::MoleculeCounter::GetTemplateTypeName<T>() << ">("
           << GetName() << ")::ResetCounter" << G4endl;
  }
  fCounterMap.clear();
  fpLastSearch.reset(nullptr);
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4bool G4VUserMoleculeCounter<TIndex>::SearchIndexUpdated(const TIndex& index)
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
G4int G4VUserMoleculeCounter<T>::SearchUpperBoundTime(G4double time, G4bool sameIndex)
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
void G4VUserMoleculeCounter<TIndex>::AbsorbCounter(
  std::weak_ptr<G4VMoleculeCounterInternalBase> wpCounter)
{
  if (wpCounter.expired()) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type G4VUserMoleculeCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
           << "Because the weak_ptr has expired since it has been passed, indicating something "
              "very wrong!"
           << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
                "BAD_REFERENCE", FatalException, errMsg);
  }

  std::shared_ptr<G4VUserMoleculeCounter<TIndex>> pCounter =
    std::dynamic_pointer_cast<G4VUserMoleculeCounter<TIndex>>(wpCounter.lock());

  if (pCounter == nullptr) {
    G4ExceptionDescription errMsg;
    errMsg << "Could not cast the pointer to type G4VUserMoleculeCounter<"
           << G4::MoleculeCounter::GetTemplateTypeName<TIndex>() << ">!\n"
           << "Because the objects aren't of the same type!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
                         + G4::MoleculeCounter::GetTemplateTypeName<TIndex>() + ">::AbsorbCounter"),
                "BAD_REFERENCE", FatalException, errMsg);
  }

  if (pCounter->GetType() != GetType()) {
    G4ExceptionDescription errMsg;
    errMsg << "You are trying to absorb a counter with different type!" << G4endl;
    G4Exception(G4String("G4VUserMoleculeCounter<"
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
template<typename TIndex>
G4bool G4VUserMoleculeCounter<TIndex>::IsReactantIgnored(const G4MoleculeDefinition* molecule){
  return G4::MoleculeCounter::Contains(fIgnoredMolecules, molecule);
}

//------------------------------------------------------------------------------

template<typename TIndex>
G4bool G4VUserMoleculeCounter<TIndex>::IsReactantIgnored(const G4MolecularConfiguration* reactant){
    return G4::MoleculeCounter::Contains(fIgnoredReactants, reactant);
}

#endif
