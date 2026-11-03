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
// G4WorkerThread implementation
//
// Authors: X.Dong, A.Dotti, 2013
// --------------------------------------------------------------------

#include "G4WorkerThread.hh"

#include "G4GeometryWorkspace.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4MTRunManager.hh"
#include "G4ParticlesWorkspace.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicsListWorkspace.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SolidsWorkspace.hh"
#include "G4WorkerRunManager.hh"

// --------------------------------------------------------------------
void G4WorkerThread::SetThreadId(G4int tid)
{
  threadId = tid;
}

// --------------------------------------------------------------------
G4int G4WorkerThread::GetThreadId() const
{
  return threadId;
}

// --------------------------------------------------------------------
void G4WorkerThread::SetNumberThreads(G4int nw)
{
  numThreads = nw;
}

// --------------------------------------------------------------------
G4int G4WorkerThread::GetNumberThreads() const
{
  return numThreads;
}

// --------------------------------------------------------------------
void G4WorkerThread::BuildGeometryAndPhysicsVector()
{
  // Initialise all split classes
  // with copy of data from master thread

  G4GeometryWorkspace::GetPool()->CreateAndUseWorkspace();
  G4SolidsWorkspace::GetPool()->CreateAndUseWorkspace();
  G4ParticlesWorkspace::GetPool()->CreateAndUseWorkspace();
  G4PhysicsListWorkspace::GetPool()->CreateAndUseWorkspace();
}

// --------------------------------------------------------------------
void G4WorkerThread::DestroyGeometryAndPhysicsVector()
{
  // Clear all split classes

  G4GeometryWorkspace::GetPool()->CleanUpAndDestroyAllWorkspaces();
  G4SolidsWorkspace::GetPool()->CleanUpAndDestroyAllWorkspaces();
  G4ParticlesWorkspace::GetPool()->CleanUpAndDestroyAllWorkspaces();
  G4PhysicsListWorkspace::GetPool()->CleanUpAndDestroyAllWorkspaces();
}

// --------------------------------------------------------------------
void G4WorkerThread::UpdateGeometryAndPhysicsVectorFromMaster()
{
  // =================================================
  // Step-0: keep sensitive detector and field manager
  // =================================================
  // First remember SD and Field associated with worker
  // in order to re-use it
  // (note that all the stuff after this will reset SD and Field)
  using LV2SDFM = std::map<G4LogicalVolume*, std::pair<G4VSensitiveDetector*, G4FieldManager*>>;
  LV2SDFM lvmap;

  using R2FSM = std::map<G4Region*, std::pair<G4FastSimulationManager*, G4UserSteppingAction*>>;
  R2FSM rgnmap;

  G4LogicalVolumeStore* mLogVolStore = G4LogicalVolumeStore::GetInstance();
  for (auto lv : *mLogVolStore)
  {
    // Volumes added between runs do not yet have worker-specific SD/FM
    // associations. Before preserving existing associations, check that
    // their instance IDs fit within this worker's current split-data array.
    //
    // The store already contains the master's current geometry, whereas
    // this worker still has the split-data array from its previous update.
    // Adding or rebuilding volumes creates new instance IDs; deleting old
    // volumes does not reclaim their IDs. When the master array grows,
    // existing worker arrays retain their previous capacity until Step-2.
    // GetSensitiveDetector() and GetFieldManager() index the worker array
    // directly, so even a valid master SD/FM pointer does not make these
    // getters safe for a new volume whose ID is beyond that capacity.
    //
    // Skip such volumes here: there are no worker associations to preserve
    // for them yet. Step-2 creates the new workspace, and the subsequent
    // worker geometry initialisation calls ConstructSDandField(). We must
    // save existing associations before resetting the workspace, otherwise
    // the worker-specific SD and field manager pointers would be lost.
    //
    if (!lv->HasSubInstance())
    {
      continue;
    }

    // New volumes within the old array's spare capacity have null worker
    // SD/FM pointers: G4GeomSplitter value-initialises these entries before
    // copying them to workers, so no associations will be saved for them.
    //
    // Preserve existing worker associations only when the corresponding
    // master association is non-null. This retains the assumption that
    // ConstructSDandField() sets up SD/FM associations on the master as
    // well as the workers; worker-only associations are not preserved.
    // These master-pointer checks govern preservation, whereas the
    // capacity check above ensures the worker getters can be called safely.
    //
    G4VSensitiveDetector* sd = nullptr;
    G4FieldManager* fmgr = nullptr;
    if (lv->GetMasterSensitiveDetector() != nullptr)
    {
      sd = lv->GetSensitiveDetector();
    }
    if (lv->GetMasterFieldManager() != nullptr)
    {
      fmgr = lv->GetFieldManager();
    }
    if (sd != nullptr || fmgr != nullptr)
    {
      lvmap[lv] = std::make_pair(sd, fmgr);
    }
  }
  G4RegionStore* mRegStore = G4RegionStore::GetInstance();
  for (auto reg : *mRegStore)
  {
    // Regions added since the previous run have the same allocation issue.
    if (!reg->HasSubInstance())
    {
      continue;
    }
    G4FastSimulationManager* fsm = reg->GetFastSimulationManager();
    G4UserSteppingAction* usa = reg->GetRegionalSteppingAction();
    if (fsm != nullptr || usa != nullptr)
    {
      rgnmap[reg] = std::make_pair(fsm, usa);
    }
  }

  //===========================
  // Step-1: Clean the workspace
  //===========================
  G4GeometryWorkspace* geomWorkspace = G4GeometryWorkspace::GetPool()->GetWorkspace();
  geomWorkspace->DestroyWorkspace();
  G4SolidsWorkspace* solidWorkspace = G4SolidsWorkspace::GetPool()->GetWorkspace();
  solidWorkspace->DestroyWorkspace();

  //===========================
  // Step-2: Re-create and initialize workspace
  //===========================
  geomWorkspace->InitialiseWorkspace();
  solidWorkspace->InitialiseWorkspace();

  //===================================================
  // Step-4: Restore sensitive detector and field manager
  //===================================================
  for (const auto& it : lvmap)
  {
    G4LogicalVolume* lv = it.first;
    G4VSensitiveDetector* sd = (it.second).first;
    G4FieldManager* fmgr = (it.second).second;
    if (fmgr != nullptr)  // What should be the second parameter?
    {  // We use always false for MT mode
      lv->SetFieldManager(fmgr, false);
    }
    if (sd != nullptr)
    {
      lv->SetSensitiveDetector(sd);
    }
  }
  for (const auto& it3 : rgnmap)
  {
    G4Region* reg = it3.first;
    G4FastSimulationManager* fsm = (it3.second).first;
    if (fsm != nullptr) reg->SetFastSimulationManager(fsm);
    G4UserSteppingAction* usa = (it3.second).second;
    if (usa != nullptr) reg->SetRegionalSteppingAction(usa);
  }
}

// --------------------------------------------------------------------
void G4WorkerThread::SetPinAffinity(G4int affinity) const
{
  if (affinity == 0) return;

#if !defined(WIN32)
  G4cout << "AFFINITY SET" << G4endl;
  // Assign this thread to cpus in a round robin way
  G4int offset = affinity;
  G4int cpuindex = 0;
  if (std::abs(offset) > G4Threading::G4GetNumberOfCores())
  {
    G4Exception("G4WorkerThread::SetPinAffinity()", "Run0100", JustWarning,
                "Cannot set thread affinity, affinity parameter larger than "
                "number of cores");
    return;
  }
  if (offset > 0)  // Start assigning affinity to given CPU
  {
    --offset;
    cpuindex = (GetThreadId() + offset) % G4Threading::G4GetNumberOfCores();
    // Round robin
  }
  else  // Exclude the given CPU
  {
    offset *= -1;
    --offset;
    G4int myidx = GetThreadId() % (G4Threading::G4GetNumberOfCores() - 1);
    cpuindex = myidx + static_cast<G4int>(myidx >= offset);
  }
  G4cout << "Setting affinity to:" << cpuindex << G4endl;

#  if defined(G4MULTITHREADED)
  // Avoid compilation warning in C90 standard w/o MT
  G4NativeThread t = pthread_self();
#  else
  G4NativeThread t;
#  endif
  G4bool success = G4Threading::G4SetPinAffinity(cpuindex, t);
  if (!success)
  {
    G4Exception("G4MTRunManagerKernel::StarThread()", "Run0101", JustWarning,
                "Cannot set thread affinity.");
  }
#endif
}
