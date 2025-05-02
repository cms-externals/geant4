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

#ifndef G4MOLECULECOUNTERMANAGER_HH
#define G4MOLECULECOUNTERMANAGER_HH 1

#include "G4Exception.hh"
#include "G4Threading.hh"
#include "G4Types.hh"
#include "G4VMoleculeCounter.hh"
#include "G4VMoleculeReactionCounter.hh"
#include "G4ios.hh"
#include <functional>
#include <map>
#include <vector>

class G4Event;
class G4Run;
class G4Step;
class G4StepPoint;
class G4Track;
class G4MoleculeCounterManagerMessenger;

class G4MoleculeCounterManager final
{
  private:
    struct Private
    {
        explicit Private() = default;
    };

  public:
    G4MoleculeCounterManager(Private);
    virtual ~G4MoleculeCounterManager();

  private:
    G4bool fInstancesRegistered{false};
    void RegisterInstance();

  public:
    static std::shared_ptr<G4MoleculeCounterManager> Instance();

    void Initialize();

    //
    // Management

    G4int RegisterCounter(std::unique_ptr<G4VMoleculeCounter>);
    G4int RegisterCounter(std::unique_ptr<G4VMoleculeReactionCounter>);
    void DeregisterAllCounters();

  private:
    template<typename T>
    G4int RegisterCounter(std::map<G4int, std::shared_ptr<T>>&, std::unique_ptr<T>,
                          std::function<G4int()>);
    void InitializeMaster();
    void InitializeWorker();

  public:
    // methods to receive forwarding from G4DNAChemistryManager
    void BeginOfEventAction(const G4Event*);
    void BeginOfRunAction(const G4Run*);
    void EndOfEventAction(const G4Event*);
    void EndOfRunAction(const G4Run*);
    // Dumping counters
    void DumpMasterCounters() const;
    void DumpWorkerCounters() const;
	// Accumulation of counters into master
	void AbsorbWorkerManagerCounters(const G4MoleculeCounterManager* = nullptr);

    //
    // Calls to Molecule Counters
  public:
    //[[deprecated("This should only be used for IRT and may be replaced as well.")]]
    void AddMoleculeWithoutTrack(const G4MolecularConfiguration*, G4double, G4int = 1);
    //[[deprecated("This should only be used for IRT and may be replaced as well.")]]
    void RemoveMoleculeWithoutTrack(const G4MolecularConfiguration*, G4double, G4int = 1);

    void AddMolecule(const G4Track*, G4double, G4int = 1);
    void RemoveMolecule(const G4Track*, G4double, G4int = 1);

    void AddMolecule(const G4Track*, const G4StepPoint*, G4double, G4int = 1);
    void RemoveMolecule(const G4Track*, const G4StepPoint*, G4double, G4int = 1);

    void ActivateCounterAtTimes(G4int, G4double, G4double, G4bool = true, G4bool = true);

    //
    // Calls to Molecule Reaction Counters

    void RecordReaction(const G4DNAMolecularReactionData*, G4double, G4int = 1);
//    void RecordReaction(const G4Track*, const G4Track*, const G4DNAMolecularReactionData*, G4double,
//                        G4int = 1);

    void ActivateReactionCounterAtTimes(G4int, G4double, G4double, G4bool = true, G4bool = true);

    //
    // Other Calls to both types or only one (e.g., BroadcastIgnoreMolecule)

    void ResetCounters();

    void NotifyOfStep(const G4Step*);
    void NotifyOfFinalize();

    // These will be broadcast to all counters registered at time of call
    void BroadcastIgnoreMolecule(const G4MoleculeDefinition*);
    void BroadcastIgnoreReactant(const G4MolecularConfiguration*);
    void BroadcastRegisterAllMoleculesAndReactants();

  private:
    static std::shared_ptr<G4MoleculeCounterManager> fpMasterInstance;
    static std::vector<std::weak_ptr<G4MoleculeCounterManager>> fWorkerInstances;
    G4ThreadLocalStatic std::shared_ptr<G4MoleculeCounterManager> fpInstance;

    std::unique_ptr<G4MoleculeCounterManagerMessenger> fpMessenger;

    G4int fVerbosity;
    G4bool fIsInitialized;
    G4bool fIsActive;
    static std::atomic<G4bool> fResetCountersBeforeEvent;
    static std::atomic<G4bool> fResetCountersBeforeRun;
    static std::atomic<G4bool> fResetMasterCounterWithWorkers;
    static std::atomic<G4bool> fAccumulateCounterIntoMaster;

    std::map<G4int, std::shared_ptr<G4VMoleculeCounter>> fCounters{};
    std::map<G4int, std::shared_ptr<G4VMoleculeReactionCounter>> fReactionCounters{};

  public:
    G4bool GetIsActive() const { return fIsActive; }
    void SetIsActive(G4bool flag) { fIsActive = flag; }

    G4int GetVerbosity() const { return fVerbosity; }
    void SetVerbosity(G4int v) { fVerbosity = v; }

    G4bool GetResetCountersBeforeEvent() const;
    void SetResetCountersBeforeEvent(G4bool = true);

    G4bool GetResetCountersBeforeRun() const;
    void SetResetCountersBeforeRun(G4bool = true);

    G4bool GetResetMasterCounterWithWorkers() const;
    void SetResetMasterCounterWithWorkers(G4bool = true);

    G4bool GetAccumulateCounterIntoMaster() const;
    void SetAccumulateCounterIntoMaster(G4bool = true);

    std::vector<std::weak_ptr<G4VMoleculeCounter>> GetMoleculeCounters();
    std::vector<std::weak_ptr<G4VMoleculeCounter>> GetMoleculeCounters(G4String);
    std::weak_ptr<G4VMoleculeCounter> GetMoleculeCounter(G4int);

    std::vector<std::weak_ptr<G4VMoleculeReactionCounter>> GetMoleculeReactionCounters();
    std::vector<std::weak_ptr<G4VMoleculeReactionCounter>> GetMoleculeReactionCounters(G4String);
    std::weak_ptr<G4VMoleculeReactionCounter> GetMoleculeReactionCounter(G4int);

  private:
    G4ThreadLocalStatic std::atomic<G4bool> fBeginOfEventTriggered;
    static std::atomic<G4bool> fBeginOfRunTriggered;
};

template<typename T>
G4int G4MoleculeCounterManager::RegisterCounter(std::map<G4int, std::shared_ptr<T>>& map,
                                                std::unique_ptr<T> counter,
                                                std::function<G4int()> idProvider)
{
  // this template allows for more than one counter type to be registered in the same way without
  // duplicating code; see public methods for G4VMoleculeCounter and G4VMoleculeReactionCounter

  if (fVerbosity > 0) {
    G4cout << "G4MoleculeCounterManager::RegisterCounter ("
           << (G4Threading::IsMasterThread() ? "master" : "worker") << ")" << G4endl;
  }

  if (counter->GetManagedId() >= 0) {
    G4ExceptionDescription description;
    description << "Trying to add a counter whose id was already altered to be non-negative!\n";
    description << "  Id: " << counter->GetManagedId() << "\n";
    description << "Name: " << counter->GetName();
    G4Exception("G4MoleculeCounterManager::RegisterCounter", "MOLMAN002", FatalErrorInArgument,
                description);
  }

  std::shared_ptr<T> sp = std::move(counter);
  // Set managed Id
  sp->SetManagedId(idProvider());

  map.emplace(sp->GetManagedId(), sp);

  return sp->GetManagedId();
}

#endif  // G4MOLECULECOUNTERMANAGER_HH
