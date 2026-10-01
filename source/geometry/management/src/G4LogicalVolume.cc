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
// G4LogicalVolume implementation
//
// 11.07.95 P.Kent: Initial version
// 04.08.97 P.M.DeFreitas: Added methods for parameterised simulation
// 09.11.98 M.Verderi, J.Apostolakis: Added BiasWeight member and accessors
// 12.02.99 S.Giani: Added user defined optimisation quality
// 18.04.01 G.Cosmo: Migrated to STL vector
// 17.05.02 G.Cosmo: Added IsToOptimise() method and related flag
// 24.09.02 G.Cosmo: Added flags and accessors for region cuts handling
// 12.11.04 G.Cosmo: Added GetMass() method for computing mass of the tree
// 01.03.05 G.Santin: Added flag for optional propagation of GetMass()
// 15.01.13 G.Cosmo, A.Dotti: Modified for thread-safety for MT
// --------------------------------------------------------------------

#include "G4LogicalVolume.hh"

#include "G4GeomTaskHelper.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4VPVParameterisation.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <algorithm>
#include <set>
#include <vector>

// This new field helps to use the class G4LVManager
//
G4LVManager G4LogicalVolume::subInstanceManager;

// Helps to use the class G4LVMutableManager, for logical volumes marked
// mutable (used inside a parameterised/division-type placement).
//
G4LVMutableManager G4LogicalVolume::subInstanceManagerMutable;

// Thread-local data fields are accessed via the (out-of-line) helper GetLVData()

namespace
{
  thread_local G4bool tlsInMassTask = false;

  class ScopedMassTaskFlag
  {
    public:

      ScopedMassTaskFlag() : fPrevious(tlsInMassTask) { tlsInMassTask = true; }
      ~ScopedMassTaskFlag() { tlsInMassTask = fPrevious; }

      ScopedMassTaskFlag(const ScopedMassTaskFlag&) = delete;
      ScopedMassTaskFlag& operator=(const ScopedMassTaskFlag&) = delete;

    private:

      G4bool fPrevious = false;
  };

  template <typename F>
  G4bool ParallelFor(std::size_t count, F&& task)
  {
#ifdef G4MULTITHREADED
    if (count < 2 || tlsInMassTask
     || G4GeometryManager::GetInstance()->IsSerialGetMassForced())
    {
      return false;
    }

    return G4GeomTaskHelper::ParallelForWithWorkerContext(count, 1,
      [&](std::size_t begin, std::size_t end, std::size_t)
      {
        for (std::size_t i = begin; i < end; ++i)
        {
          task(i);
        }
      },
      [](std::size_t, const auto& run)
      {
        ScopedMassTaskFlag massTaskFlag;
        G4GeomTaskHelper::ScopedGeometryWorkspaces workspaces;
        run();
      });
#else
    (void)count;
    (void)task;
    return false;
#endif
  }

#ifdef G4MULTITHREADED
  void WarmCubicVolumeCache(G4LogicalVolume* logical,
                            G4bool propagate,
                            std::set<G4VSolid*>& warmedSolids)
  {
    if (logical == nullptr)
    {
      return;
    }

    auto* solid = logical->GetSolid();
    if (solid != nullptr && warmedSolids.insert(solid).second)
    {
      solid->GetCubicVolume();
    }

    auto noDaughters = logical->GetNoDaughters();
    for (std::size_t idaughter = 0; idaughter < noDaughters; ++idaughter)
    {
      auto* physDaughter = logical->GetDaughter((G4int)idaughter);
      auto* logDaughter = physDaughter->GetLogicalVolume();
      auto* physParam = physDaughter->GetParameterisation();
      if (physParam != nullptr)
      {
        for (auto icopy = 0; icopy < physDaughter->GetMultiplicity(); ++icopy)
        {
          auto* daughterSolid = physParam->ComputeSolid(icopy, physDaughter);
          if (daughterSolid != nullptr)
          {
            daughterSolid->ComputeDimensions(physParam, icopy, physDaughter);
            daughterSolid->GetCubicVolume();
          }
        }
      }
      if (propagate)
      {
        WarmCubicVolumeCache(logDaughter, true, warmedSolids);
      }
      else if (physParam == nullptr)
      {
        auto* daughterSolid = logDaughter->GetSolid();
        if (daughterSolid != nullptr
         && warmedSolids.insert(daughterSolid).second)
        {
          daughterSolid->GetCubicVolume();
        }
      }
    }
  }
#endif

}  // namespace

// ********************************************************************
// Constructor - sets member data and adds to logical Store,
//               voxel pointer for optimisation set to 0 by default.
//               Initialises daughter vector to 0 length.
// ********************************************************************
//
G4LogicalVolume::G4LogicalVolume(G4VSolid* pSolid, G4Material* pMaterial, const G4String& name,
                                 G4FieldManager* pFieldMgr, G4VSensitiveDetector* pSDetector,
                                 G4UserLimits* pULimits, G4bool optimise)
  : fDaughters(0, (G4VPhysicalVolume*)nullptr), fDaughtersVolumeType(kNormal), fOptimise(optimise)
{
  // Initialize 'Shadow'/master pointers - for use in copying to workers
  //
  fSolid = pSolid;
  fSensitiveDetector = pSDetector;
  fFieldManager = pFieldMgr;

  instanceID = subInstanceManager.CreateSubInstance();
  AssignFieldManager(pFieldMgr);

  ResetMass();
  fCutsCouple = nullptr;

  SetSolid(pSolid);
  SetMaterial(pMaterial);
  SetName(name);
  SetSensitiveDetector(pSDetector);
  SetUserLimits(pULimits);

  //
  // Add to store
  //
  G4LogicalVolumeStore::Register(this);
}

// ********************************************************************
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
// ********************************************************************
//
G4LogicalVolume::G4LogicalVolume(__void__&) : fDaughters(0, (G4VPhysicalVolume*)nullptr), fName("")
{
  instanceID = subInstanceManager.CreateSubInstance();

  SetSensitiveDetector(nullptr);
  SetFieldManager(nullptr, false);

  ResetMass();
  fCutsCouple = nullptr;

  // Add to store
  //
  G4LogicalVolumeStore::Register(this);
}

// ********************************************************************
// Destructor - Removes itself from solid Store
// NOTE: Not virtual
// ********************************************************************
//
G4LogicalVolume::~G4LogicalVolume()
{
  if (!fLock && fRootRegion)  // De-register root region first if not locked
  {  // and flagged as root logical-volume
    fRegion->RemoveRootLogicalVolume(this, true);
  }
  G4LogicalVolumeStore::DeRegister(this);
}

// ********************************************************************
// SetName - Set volume name and notify store of the change
// ********************************************************************
//
void G4LogicalVolume::SetName(const G4String& pName)
{
  fName = pName;
  G4LogicalVolumeStore::GetInstance()->SetMapValid(false);
}

// ********************************************************************
// InitialiseWorker
//
// This method is similar to the constructor. It is used by each worker
// thread to achieve the same effect as that of the master thread exept
// to register the new created instance. This method is invoked explicitly.
// It does not create a new G4LogicalVolume instance. It only assign the value
// for the fields encapsulated by the class G4LVData.
// ********************************************************************
//
void G4LogicalVolume::InitialiseWorker(G4LogicalVolume* /*pMasterObject*/, G4VSolid* pSolid,
                                       G4VSensitiveDetector* pSDetector)
{
  subInstanceManager.SlaveCopySubInstanceArray();

  if (!fImmutable)
  {
    // This thread's G4LVMutableData slot was already populated (as a raw
    // copy of the master's values) by subInstanceManagerMutable's own
    // SlaveCopySubInstanceArray(), called from G4GeometryWorkspace.
    // Install pSolid directly (bypassing SetSolid()'s public guard, which
    // forbids worker threads from touching a *shared*, non-mutable
    // volume's solid) so that G4GeometryWorkspace::CloneReplicaSolid()'s
    // genuinely thread-private clone takes effect for this thread.
    SolidSlot() = pSolid;
  }
  // else: solid is shared (same object, same address, for every thread) -
  // nothing to do here.
  SetSensitiveDetector(pSDetector);  //  How this object is available now ?
  AssignFieldManager(fFieldManager);
  // Should be set - but a per-thread copy is not available yet
  // Must not call SetFieldManager(), which propagates FieldMgr

#ifdef CLONE_FIELD_MGR
  // Create a field FieldManager by cloning
  //
  G4FieldManager workerFldMgr = fFieldManager->GetWorkerClone(G4bool * created);
  if (created || (GetFieldManager() != workerFldMgr))
  {
    SetFieldManager(fFieldManager, false);  // which propagates FieldMgr
  }
  else
  {
    // Field manager existed and is equal to current one
    //
    AssignFieldManager(workerFldMgr);
  }
#endif
}

// ********************************************************************
// Clean
// ********************************************************************
//
void G4LogicalVolume::Clean()
{
  subInstanceManager.FreeSlave();
  subInstanceManagerMutable.FreeSlave();
}

// ********************************************************************
// TerminateWorker
//
// This method is similar to the destructor. It is used by each worker
// thread to achieve the partial effect as that of the master thread.
// For G4LogicalVolume instances, nothing more to do here.
// ********************************************************************
//
void G4LogicalVolume::TerminateWorker(G4LogicalVolume* /*pMasterObject*/) {}

// ********************************************************************
// GetSubInstanceManager
// ********************************************************************
//
const G4LVManager& G4LogicalVolume::GetSubInstanceManager()
{
  return subInstanceManager;
}

// ********************************************************************
// GetMutableSubInstanceManager
// ********************************************************************
//
const G4LVMutableManager& G4LogicalVolume::GetMutableSubInstanceManager()
{
  return subInstanceManagerMutable;
}

/**
 * ********************************************************************
 * MarkMaterialMutable
 *
 * Marks this logical volume as needing a per-thread copy of its solid,
 * material, mass and cuts-couple, instead of sharing them across all
 * threads.
 *
 * Notes:
 * - Must be called for a logical volume which is used in any type of
 *   physical that can change any of those attributes!
 * - Called automatically when the volume is used as the daughter
 *   of a G4PVParameterised, G4PVDivision or G4ReplicatedSlice.
 * Idempotent.
 * ********************************************************************
 */
void G4LogicalVolume::MarkMaterialMutable()
{
  if (!fImmutable) { return; }  // Already mutable.

#ifdef G4GEOM_DEBUG
  G4cout <<"G4LogicalVolume::MarkMaterialMutable() called for lv="
      << this->GetName() << "  and its attributes are "
      << " VolumeType " << this->fDaughtersVolumeType
      << " all-of-it= " << this << G4endl;
#endif

  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4ExceptionDescription ed;
    ed << "Attempt to call MarkMaterialMutable() for volume '" << fName
       << "' after the geometry has been closed.\n"
       << "This is only supported during initial detector construction.";
    G4Exception("G4LogicalVolume::MarkMaterialMutable()", "GeomMgt0002",
                FatalException, ed);
  }

  fInstanceIDifMutable = subInstanceManagerMutable.CreateSubInstance();
  fImmutable = false;

  // The values set so far live only in the shared fields; copy them as
  // the new per-thread slot's initial value.
  G4LVMutableData& slot = subInstanceManagerMutable.offset[fInstanceIDifMutable];
  slot.fSolid = fSolid;
  slot.fMaterial = fMaterial;
  slot.fCutsCouple = fCutsCouple;
  slot.fMass = fMass.load(std::memory_order_relaxed);
}

// ********************************************************************
// GetFieldManager
// ********************************************************************
//
G4FieldManager* G4LogicalVolume::GetFieldManager() const
{
  return GetLVData().fFieldManager;
}

// ********************************************************************
// AssignFieldManager
// ********************************************************************
//
void G4LogicalVolume::AssignFieldManager(G4FieldManager* fldMgr)
{
  GetLVData().fFieldManager = fldMgr;
  if (G4Threading::IsMasterThread())
  {
    fFieldManager = fldMgr;
  }
}

// ********************************************************************
// IsExtended
// ********************************************************************
//
G4bool G4LogicalVolume::IsExtended() const
{
  return false;
}

// ********************************************************************
// SetFieldManager
// ********************************************************************
//
void G4LogicalVolume::SetFieldManager(G4FieldManager* pNewFieldMgr, G4bool forceAllDaughters)
{
  AssignFieldManager(pNewFieldMgr);

  auto NoDaughters = GetNoDaughters();
  while ((NoDaughters--) > 0)
  {
    G4LogicalVolume* DaughterLogVol;
    DaughterLogVol = GetDaughter(NoDaughters)->GetLogicalVolume();
    if (forceAllDaughters || (DaughterLogVol->GetFieldManager() == nullptr))
    {
      DaughterLogVol->SetFieldManager(pNewFieldMgr, forceAllDaughters);
    }
  }
}

// ********************************************************************
// AddDaughter
// ********************************************************************
//
void G4LogicalVolume::AddDaughter(G4VPhysicalVolume* pNewDaughter)
{
  EVolume daughterType = pNewDaughter->VolumeType();

  // The type of the navigation needed is determined by the first daughter
  //
  if (fDaughters.empty())
  {
    fDaughtersVolumeType = daughterType;
  }
  else
  {
    // Check consistency of detector description

    // 1. A replica or parameterised volume can have only one daughter
    //
    if (fDaughters[0]->IsReplicated())
    {
      std::ostringstream message;
      message << "ERROR - Attempt to place a volume in a mother volume" << G4endl
              << "        already containing a replicated volume." << G4endl
              << "        A volume can either contain several placements" << G4endl
              << "        or a unique replica or parameterised volume !" << G4endl
              << "           Mother logical volume: " << GetName() << G4endl
              << "           Placing volume: " << pNewDaughter->GetName() << G4endl;
      G4Exception("G4LogicalVolume::AddDaughter()", "GeomMgt0002", FatalException, message,
                  "Replica or parameterised volume must be the only daughter!");
    }
    else
    {
      // 2. Ensure that Placement and External physical volumes do not mix
      //
      if (daughterType != fDaughtersVolumeType)
      {
        std::ostringstream message;
        message << "ERROR - Attempt to place a volume in a mother volume" << G4endl
                << "        already containing a different type of volume." << G4endl
                << "        A volume can either contain" << G4endl
                << "        - one or more placements, OR" << G4endl
                << "        - one or more 'external' type physical volumes." << G4endl
                << "          Mother logical volume: " << GetName() << G4endl
                << "          Volume being placed: " << pNewDaughter->GetName() << G4endl;
        G4Exception("G4LogicalVolume::AddDaughter()", "GeomMgt0002", FatalException, message,
                    "Cannot mix placements and external physical volumes !");
      }
    }
  }

  // Invalidate previous calculation of mass - if any - for all threads
  //
  ResetMass();
  fDaughters.push_back(pNewDaughter);

  G4LogicalVolume* pDaughterLogical = pNewDaughter->GetLogicalVolume();

  // Propagate the Field Manager, if the daughter has no field Manager
  //
  G4FieldManager* pDaughterFieldManager = pDaughterLogical->GetFieldManager();

  // Avoid propagating the fieldManager pointer if null
  // and daughter's one is null as well...
  //
  if ((GetLVData().fFieldManager != nullptr) && (pDaughterFieldManager == nullptr))
  {
    pDaughterLogical->SetFieldManager(GetLVData().fFieldManager, false);
  }
  if (fRegion != nullptr)
  {
    PropagateRegion();
    fRegion->RegionModified(true);
  }
}

// ********************************************************************
// RemoveDaughter
// ********************************************************************
//
void G4LogicalVolume::RemoveDaughter(const G4VPhysicalVolume* p)
{
  for (auto i = fDaughters.cbegin(); i != fDaughters.cend(); ++i)
  {
    if (**i == *p)
    {
      fDaughters.erase(i);
      break;
    }
  }
  if (fRegion != nullptr)
  {
    fRegion->RegionModified(true);
  }
  ResetMass();
}

// ********************************************************************
// ClearDaughters
// ********************************************************************
//
void G4LogicalVolume::ClearDaughters()
{
  fDaughters.erase(fDaughters.cbegin(), fDaughters.cend());
  if (fRegion != nullptr)
  {
    fRegion->RegionModified(true);
  }
  ResetMass();
}

// ********************************************************************
// SetSensitiveDetector
// ********************************************************************
//
void G4LogicalVolume::SetSensitiveDetector(G4VSensitiveDetector* pSDetector)
{
  GetLVData().fSensitiveDetector = pSDetector;
  if (G4Threading::IsMasterThread())  { fSensitiveDetector = pSDetector; }
}

// ********************************************************************
// CheckGeometryConfiguration
// ********************************************************************
//
void G4LogicalVolume::CheckGeometryConfiguration(const char* methodName,
                                                 const G4String& volName,
                                                 G4bool isMutable)
{
  if (!isMutable && !G4Threading::IsMasterThread())
  {
    G4ExceptionDescription ed;
    ed << "Attempt to call G4LogicalVolume::" << methodName
       << "() for volume '" << volName
       << "' from a worker thread.\n"
       << "This logical volume's solid/material/cuts-couple are shared\n"
       << "by all threads (it is not used inside a parameterised or\n"
       << "division-type placement) and may only be modified from the\n"
       << "master thread.";
    G4Exception((G4String("G4LogicalVolume::") + methodName + "()").c_str(),
                "GeomMgt0002", FatalException, ed);
    }
#ifdef G4GEOM_DEBUG
  else if (!isMutable && G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
     G4ExceptionDescription ed;
     ed << "Call to G4LogicalVolume::" << methodName
        << "() for volume '" << volName
        << "' while the geometry is closed.\n"
        << "This logical volume's solid/material/cuts-couple are shared\n"
        << "by all threads.\n"
        << "Ensure that no worker threads are active when making this modification.";
     G4Exception((G4String("G4LogicalVolume::") + methodName + "()").c_str(),
                 "GeomMgt0002", JustWarning, ed);
  }
#endif
}


// The following methods are out-of-line to ensure correctness in DLLs / Windows
//

// ********************************************************************
// GetLVData: 
// ********************************************************************
//
G4LVData& G4LogicalVolume::GetLVData() const
{
  return subInstanceManager.offset[instanceID];
}


// ********************************************************************
// GetSensitiveDetector
// ********************************************************************
//
G4VSensitiveDetector* G4LogicalVolume::GetSensitiveDetector() const
{
  return GetLVData().fSensitiveDetector;
}

// ********************************************************************
// SolidSlot / MaterialSlot / CutsCoupleSlot
// ********************************************************************
//
G4VSolid*& G4LogicalVolume::SolidSlot() const
{
  if (fImmutable)
  {
    return const_cast<G4LogicalVolume*>(this)->fSolid;
  }
  return subInstanceManagerMutable.offset[fInstanceIDifMutable].fSolid;
}

G4Material*& G4LogicalVolume::MaterialSlot() const
{
  if (fImmutable)
  {
    return const_cast<G4LogicalVolume*>(this)->fMaterial;
  }
  return subInstanceManagerMutable.offset[fInstanceIDifMutable].fMaterial;
}

G4MaterialCutsCouple*& G4LogicalVolume::CutsCoupleSlot() const
{
  if (fImmutable)
  {
    return const_cast<G4LogicalVolume*>(this)->fCutsCouple;
  }
  return subInstanceManagerMutable.offset[fInstanceIDifMutable].fCutsCouple;
}

// ********************************************************************
// GetSolidMutable / GetMaterialMutable / GetCutsCoupleMutable /
// GetCachedMassMutable / SetCachedMassMutable
//
// Out-of-line bodies for the mutable branch of GetSolid()/GetMaterial()/
// GetMaterialCutsCouple()/GetCachedMass()/SetCachedMass() (inline in
// G4LogicalVolume.icc for their shared fast path only). See the note
// above SolidSlot().
// ********************************************************************
//
G4VSolid* G4LogicalVolume::GetSolidMutable() const
{
  return SolidSlot();
}

G4Material* G4LogicalVolume::GetMaterialMutable() const
{
  return MaterialSlot();
}

G4MaterialCutsCouple* G4LogicalVolume::GetCutsCoupleMutable() const
{
  return CutsCoupleSlot();
}

G4double G4LogicalVolume::GetCachedMassMutable() const
{
  return subInstanceManagerMutable.offset[fInstanceIDifMutable].fMass;
}

void G4LogicalVolume::SetCachedMassMutable(G4double val)
{
  subInstanceManagerMutable.offset[fInstanceIDifMutable].fMass = val;
}

// ********************************************************************
// SetSolid
// ********************************************************************
//
void G4LogicalVolume::SetSolid(G4VSolid* pSolid)
{
  CheckGeometryConfiguration("SetSolid", fName, IsMaterialMutable());
  SolidSlot() = pSolid;
  // Keep the master-visible shadow (fSolid) in sync for a mutable volume,
  // exactly as SetSensitiveDetector()/AssignFieldManager() already do for
  // their own master-shadow fields: G4GeometryWorkspace::
  // InitialisePhysicalVolumes() seeds every newly-created worker's
  // per-thread slot from GetMasterSolid() (== fSolid), so a master-thread
  // change must land here or new workers see a stale solid.
  if (!fImmutable && G4Threading::IsMasterThread()) { fSolid = pSolid; }
  ResetMass();
}

// ********************************************************************
// SetMaterial
// ********************************************************************
//
void G4LogicalVolume::SetMaterial(G4Material* pMaterial)
{
  CheckGeometryConfiguration("SetMaterial", fName, IsMaterialMutable());
  MaterialSlot() = pMaterial;
  if (!fImmutable && G4Threading::IsMasterThread()) { fMaterial = pMaterial; }
  ResetMass();
}

// ********************************************************************
// UpdateMaterial
// ********************************************************************
//
void G4LogicalVolume::UpdateMaterial(G4Material* pMaterial)
{
  CheckGeometryConfiguration("UpdateMaterial", fName, IsMaterialMutable());
  MaterialSlot() = pMaterial;
  if (fRegion != nullptr) { CutsCoupleSlot() = fRegion->FindCouple(pMaterial); }
  if (!fImmutable && G4Threading::IsMasterThread())
  {
    fMaterial = pMaterial;
    if (fRegion != nullptr) { fCutsCouple = CutsCoupleSlot(); }
  }
  ResetMass();
}

// ********************************************************************
// SetMaterialCutsCouple
//
// NOTE: legitimately re-assigned on the master thread at the start of
// every run, after the geometry is closed (G4RunManagerKernel::
// RunInitialization -> UpdateRegion -> G4ProductionCutsTable::
// ScanAndSetCouple). For a shared (non-mutable) volume, workers never call
// it and must not: the guard is master-thread-only, not closed-geometry.
// ********************************************************************
//
void G4LogicalVolume::SetMaterialCutsCouple(G4MaterialCutsCouple* cuts)
{
  if (fImmutable && !G4Threading::IsMasterThread())
  {
    G4ExceptionDescription ed;
    ed << "Attempt to call SetMaterialCutsCouple() for volume '" << fName
       << "' from a worker thread.\n"
       << "This logical volume's material-cuts-couple is shared by all\n"
       << "threads (it is not used inside a parameterised or division-type\n"
       << "placement) and may only be set from the master thread.";
    G4Exception("G4LogicalVolume::SetMaterialCutsCouple()", "GeomMgt0002",
                FatalException, ed);
  }
  CutsCoupleSlot() = cuts;
  if (!fImmutable && G4Threading::IsMasterThread()) { fCutsCouple = cuts; }
}

// ********************************************************************
// IsAncestor
//
// Finds out if the current logical volume is an ancestor of a given
// physical volume
// ********************************************************************
//
G4bool G4LogicalVolume::IsAncestor(const G4VPhysicalVolume* aVolume) const
{
  G4bool isDaughter = IsDaughter(aVolume);
  if (!isDaughter)
  {
    for (const auto& daughter : fDaughters)
    {
      isDaughter = daughter->GetLogicalVolume()->IsAncestor(aVolume);
      if (isDaughter)
      {
        break;
      }
    }
  }
  return isDaughter;
}

// ********************************************************************
// TotalVolumeEntities
//
// Returns the total number of physical volumes (replicated or placed)
// in the tree represented by the current logical volume.
// ********************************************************************
//
G4int G4LogicalVolume::TotalVolumeEntities() const
{
  G4int vols = 1;
  for (auto physDaughter : fDaughters)
  {
    vols +=
      physDaughter->GetMultiplicity() * physDaughter->GetLogicalVolume()->TotalVolumeEntities();
  }
  return vols;
}

// ********************************************************************
// ComputeDaughterMassContribution
//
// Computes daughter volume mass from its logical volume.
// ********************************************************************
//
G4double G4LogicalVolume::
ComputeDaughterMassContribution(G4VPhysicalVolume* physDaughter,
                                G4double motherDensity, G4bool propagate)
{
  G4LogicalVolume* logDaughter = physDaughter->GetLogicalVolume();
  G4double contribution = 0.0;

  for (auto i = 0; i < physDaughter->GetMultiplicity(); ++i)
  {
    G4VSolid* daughterSolid = nullptr;
    G4Material* daughterMaterial = nullptr;
    G4VPVParameterisation* physParam = physDaughter->GetParameterisation();
    if (physParam != nullptr)
    {
      daughterSolid = physParam->ComputeSolid(i, physDaughter);
      daughterSolid->ComputeDimensions(physParam, i, physDaughter);
      daughterMaterial = physParam->ComputeMaterial(i, physDaughter);
    }
    else
    {
      daughterSolid = logDaughter->GetSolid();
      daughterMaterial = logDaughter->GetMaterial();
    }

    contribution -= daughterSolid->GetCubicVolume() * motherDensity;
    if (propagate)
    {
      contribution += logDaughter->ComputeMassNoCache(true, daughterMaterial);
    }
  }
  return contribution;
}

// ********************************************************************
// ComputeMassNoCache
//
// Computes mass without reading or publishing the cached value. Used for
// recursive daughter contributions, including parameterised material cases.
// ********************************************************************
//
G4double G4LogicalVolume::ComputeMassNoCache(G4bool propagate, G4Material* parMaterial)
{
  G4Material* logMaterial = parMaterial != nullptr ? parMaterial : GetMaterial();
  if (logMaterial == nullptr)
  {
    std::ostringstream message;
    message << "No material associated to the logical volume: " << GetName() << " !"
            << G4endl << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002", FatalException, message);
    return 0.0;
  }
  if (GetSolid() == nullptr)
  {
    std::ostringstream message;
    message << "No solid is associated to the logical volume: " << GetName() << " !"
            << G4endl << "Sorry, cannot compute the mass ...";
    G4Exception("G4LogicalVolume::GetMass()", "GeomMgt0002", FatalException, message);
    return 0.0;
  }

  G4double globalDensity = logMaterial->GetDensity();
  G4double massSum = GetSolid()->GetCubicVolume() * globalDensity;
  auto noDaughters = GetNoDaughters();
  if (noDaughters == 0)
  {
    return massSum;
  }

#ifdef G4MULTITHREADED
  if (!G4GeometryManager::GetInstance()->IsSerialGetMassForced()
   && !tlsInMassTask && noDaughters > 1)
  {
    std::set<G4VSolid*> warmedSolids;
    WarmCubicVolumeCache(this, propagate, warmedSolids);
  }
#endif

  std::vector<G4double> contributions(noDaughters, 0.0);
  auto compute = [&](std::size_t i)
  {
    contributions[i] = ComputeDaughterMassContribution(GetDaughter(i),
                                                       globalDensity, propagate);
  };

  if (!ParallelFor(noDaughters, compute))
  {
    for (std::size_t i = 0; i < noDaughters; ++i)
    {
      compute(i);
    }
  }

  for (const auto contribution : contributions)
  {
    massSum += contribution;
  }
  return massSum;
}

// ********************************************************************
// GetMass
//
// Returns the mass of the logical volume tree computed from the
// estimated geometrical volume of each solid and material associated
// to the logical volume and its daughters.
// NOTE: the computation may require considerable amount of time,
//       depending from the complexity of the geometry tree.
//       The returned value is cached and can be used for successive
//       calls (default), unless recomputation is forced by providing
//       'true' for the boolean argument in input. Computation should
//       be forced if the geometry setup has changed after the previous
//       call. By setting the 'propagate' boolean flag to 'false' the
//       method returns the mass of the present logical volume only
//       (subtracted for the volume occupied by the daughter volumes).
//       The extra argument 'parMaterial' is internally used to
//       consider cases of geometrical parameterisations by material.
// ********************************************************************
//
G4double G4LogicalVolume::GetMass(G4bool forced, G4bool propagate, G4Material* parMaterial)
{
  // Return the cached non-zero value, if not forced
  //
  const G4double cachedMass = GetCachedMass();
  if ((cachedMass != 0.0) && (!forced))
  {
    return cachedMass;
  }

  G4double massSum = ComputeMassNoCache(propagate, parMaterial);
  SetCachedMass(massSum);
  return massSum;
}

// ********************************************************************
// Change the daughters volume type -- checking proposed values
//
//  Undertakes primitive checking, to ensure that only 'legal' changes
//  are made:
//    - any type to 'external'  ( user responsibility )
//    - the type proposed is checked against the deduced type
//       (for potential switch back to 'internal' type.)
//  Returns success (true) or failure (false)
//
G4bool G4LogicalVolume::ChangeDaughtersType(EVolume aType)
{
  G4bool works = false;
  if (aType == kExternal)
  {
    // It is the responsibility of External Navigator to handle types selected
    //
    fDaughtersVolumeType = aType;
    works = true;
  }
  else
  {
    EVolume expectedVType = DeduceDaughtersType();
    works = (expectedVType == aType);
    if (works)
    {
      fDaughtersVolumeType = aType;
    }
  }
  return works;
}

// ********************************************************************
// SetVisAttributes - copy version
// ********************************************************************
//
void G4LogicalVolume::SetVisAttributes(const G4VisAttributes& VA)
{
  if (G4Threading::IsWorkerThread())
  {
    return;
  }
  fVisAttributes = std::make_shared<const G4VisAttributes>(VA);
}

// ********************************************************************
// SetVisAttributes
// ********************************************************************
//
void G4LogicalVolume::SetVisAttributes(const G4VisAttributes* pVA)
{
  if (G4Threading::IsWorkerThread())
  {
    return;
  }
  fVisAttributes = std::shared_ptr<const G4VisAttributes>(pVA, [](const G4VisAttributes*) {});
}
