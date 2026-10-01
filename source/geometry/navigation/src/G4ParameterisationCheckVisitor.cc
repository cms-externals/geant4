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
// G4ParameterisationCheckVisitor class Implementation
//
// Author: John Apostolakis (CERN), August 2026
// --------------------------------------------------------------------

#include "G4ParameterisationCheckVisitor.hh"

#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4PVParameterised.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistoryHandle.hh"
#include "G4TouchableWalker.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VVolumeMaterialScanner.hh"
#include "G4ios.hh"

#include <algorithm>
#include <ostream>
#include <utility>

// ---------------------------------------------------------------------------

G4ParameterisationCheckVisitor::G4ParameterisationCheckVisitor(G4bool verbose) : fVerbose(verbose) {}

// ---------------------------------------------------------------------------
// CheckGeometry
// ---------------------------------------------------------------------------
//
G4int G4ParameterisationCheckVisitor::CheckGeometry(G4VPhysicalVolume* startVol, G4bool repair,
                                                     G4bool verbose)
{
  if (!G4Threading::IsMasterThread())
  {
    return 0;
  }

  const G4bool geometryOpen = !G4GeometryManager::GetInstance()->IsGeometryClosed();
  if (!geometryOpen && repair)
  {
    G4Exception("G4ParameterisationCheckVisitor::CheckGeometry", "GeomNav1001", JustWarning,
               "Geometry is already closed: MarkMaterialMutable() would be illegal here. "
               "Running in diagnose-only mode (repair request ignored).");
    repair = false;
  }

  // --- Snapshot: physical-volume state for every parameterised PV --------
  // The walk changes each parameterised PV's logical volume/copy number/
  // translation/rotation/name via the parameterisation's own hooks; this
  // must be undone afterwards, exactly as G4ParameterisationChecker does.
  struct PVSnapshot
  {
    G4LogicalVolume* fLV = nullptr;
    G4int fCopyNo = 0;
    G4ThreeVector fTranslation;
    G4RotationMatrix* fRotation = nullptr;
    G4String fName;
  };
  std::map<G4VPhysicalVolume*, PVSnapshot> savedPV;
  std::map<G4VPhysicalVolume*, G4LogicalVolume*> declaredLV;

  for (auto* pv : *G4PhysicalVolumeStore::GetInstance())
  {
    if (pv == nullptr || !pv->IsParameterised())
    {
      continue;
    }
    PVSnapshot snap;
    snap.fLV = pv->GetLogicalVolume();
    snap.fCopyNo = pv->GetCopyNo();
    snap.fTranslation = pv->GetTranslation();
    snap.fRotation = pv->GetRotation();
    snap.fName = pv->GetName();
    savedPV[pv] = snap;
    declaredLV[pv] = snap.fLV;
  }

  // --- Snapshot: every logical volume's solid/material -------------------
  // Required because the walk can silently overwrite the shared solid/
  // material of a *non-mutable* logical volume with no diagnostic when run
  // on the master thread (see G4LogicalVolume::CheckGeometryConfiguration()),
  // and because a nested ComputeMaterial() can rebind a physical volume to a
  // logical volume that was never known to be parameterised beforehand.
  struct LVSnapshot
  {
    G4VSolid* fSolid = nullptr;
    G4Material* fMaterial = nullptr;
  };
  std::map<G4LogicalVolume*, LVSnapshot> savedLV;
  for (auto* lv : *G4LogicalVolumeStore::GetInstance())
  {
    if (lv == nullptr)
    {
      continue;
    }
    LVSnapshot snap;
    snap.fSolid = lv->GetSolid();
    snap.fMaterial = lv->GetMaterial();
    savedLV[lv] = snap;
  }

  // --- Walk ----------------------------------------------------------------
  G4ParameterisationCheckVisitor visitor(verbose);
  visitor.SetDeclaredLogicalVolumes(declaredLV);
  G4TouchableWalker walker(visitor);
  walker.Walk(startVol);

  // Finalise per-volume evidence that depends on the whole-walk accumulators.
  for (auto& entry : visitor.fReports)
  {
    G4VPhysicalVolume* pv = entry.first;
    VolumeReport& rep = entry.second;

    auto pcIt = visitor.fParentContexts.find(pv);
    if (pcIt != visitor.fParentContexts.end())
    {
      rep.fDistinctParentContexts = pcIt->second.size();
    }

    auto mbrIt = visitor.fMaterialByReplica.find(pv);
    if (mbrIt != visitor.fMaterialByReplica.end())
    {
      for (const auto& perReplica : mbrIt->second)
      {
        if (perReplica.second.size() > 1)
        {
          rep.fMaterialVariesWithParent = true;
          break;
        }
      }
    }
  }

  // --- Restore ---------------------------------------------------------
  // Logical-volume state first, then physical-volume state - the reverse of
  // the order the walk applied them in.
  for (auto& entry : savedLV)
  {
    G4LogicalVolume* lv = entry.first;
    const LVSnapshot& snap = entry.second;
    if (lv->GetSolid() != snap.fSolid)
    {
      lv->SetSolid(snap.fSolid);
    }
    if (lv->GetMaterial() != snap.fMaterial)
    {
      lv->UpdateMaterial(snap.fMaterial);
    }
  }
  for (auto& entry : savedPV)
  {
    G4VPhysicalVolume* pv = entry.first;
    const PVSnapshot& snap = entry.second;
    pv->SetLogicalVolume(snap.fLV);
    pv->SetCopyNo(snap.fCopyNo);
    pv->SetTranslation(snap.fTranslation);
    pv->SetRotation(snap.fRotation);
    pv->SetName(snap.fName);
  }

  // --- Repair --------------------------------------------------------------
  G4int nRepaired = 0;
  if (repair && geometryOpen)
  {
    for (auto& entry : visitor.fReports)
    {
      G4VPhysicalVolume* pv = entry.first;
      const VolumeReport& rep = entry.second;
      for (auto* lv : rep.fReboundLVs)
      {
        if (!lv->IsMaterialMutable())
        {
          lv->MarkMaterialMutable();
          ++nRepaired;
          if (verbose)
          {
            G4cout << "G4ParameterisationCheckVisitor: marked logical volume '" << lv->GetName()
                   << "' as mutable - reassigned by parameterised physical volume '"
                   << pv->GetName() << "'." << G4endl;
          }
        }
      }
    }
  }

  // --- Report ----------------------------------------------------------
  if (verbose)
  {
    visitor.PrintReport(G4cout);
  }

  const G4int nProblems = visitor.CountProblems();

  if (!geometryOpen && nProblems > 0)
  {
    G4ExceptionDescription ed;
    ed << nProblems
       << " problem(s) found among nested/parameterised volumes, and the "
          "geometry was already closed so they could not be repaired.\n"
          "Call G4ParameterisationCheckVisitor::CheckGeometry() at the end of "
          "DetectorConstruction::Construct(), while the geometry is still "
          "open, so offending logical volumes can be marked mutable "
          "automatically; or call MarkMaterialMutable() on them there "
          "directly.\n";
    visitor.PrintReport(ed);
    G4Exception("G4ParameterisationCheckVisitor::CheckGeometry", "GeomNav1002", FatalException, ed);
  }

  G4cout << "G4ParameterisationCheckVisitor::CheckGeometry: " << nProblems
         << " problem(s) found, " << nRepaired << " logical volume(s) marked mutable."
         << (visitor.fGlobalTruncated ? " Coverage TRUNCATED - this is a sample, not a proof."
                                      : "")
         << G4endl;

  return nProblems;
}

// ---------------------------------------------------------------------------

void G4ParameterisationCheckVisitor::SetDeclaredLogicalVolumes(
  std::map<G4VPhysicalVolume*, G4LogicalVolume*> declared)
{
  fDeclaredLV = std::move(declared);
}

// ---------------------------------------------------------------------------

G4ParameterisationCheckVisitor::VolumeReport&
G4ParameterisationCheckVisitor::GetOrCreateReport(G4VPhysicalVolume* pv)
{
  auto it = fReports.find(pv);
  if (it == fReports.end())
  {
    VolumeReport rep;
    rep.fPV = pv;
    it = fReports.emplace(pv, rep).first;
  }
  return it->second;
}

// ---------------------------------------------------------------------------

void G4ParameterisationCheckVisitor::VisitParameterised(const G4TouchableHistoryHandle& touchable,
                                                        G4PVParameterised* pv, G4int replicaNo,
                                                        G4int depth)
{
  ++fTotalVisits;
  if (fMaxTotalVisits >= 0 && fTotalVisits > fMaxTotalVisits)
  {
    fGlobalTruncated = true;
    return;
  }

  VolumeReport& rep = GetOrCreateReport(pv);
  rep.fInstancesDeclared = pv->GetMultiplicity();
  if (rep.fInstancesVisited == 0)
  {
    rep.fMinDepth = depth;
    rep.fMaxDepth = depth;
  }
  else
  {
    rep.fMinDepth = std::min(rep.fMinDepth, depth);
    rep.fMaxDepth = std::max(rep.fMaxDepth, depth);
  }
  ++rep.fInstancesVisited;

  G4VPVParameterisation* param = pv->GetParameterisation();
  rep.fIsNested = (param != nullptr) && param->IsNested();

  G4LogicalVolume* boundLV = pv->GetLogicalVolume();
  rep.fBoundLVs.insert(boundLV);

  auto declIt = fDeclaredLV.find(pv);
  const G4bool rebound = (declIt != fDeclaredLV.end()) && (boundLV != declIt->second);
  if (rebound)
  {
    rep.fReboundLVs.insert(boundLV);
    if (fVerbose)
    {
      G4cout << "G4ParameterisationCheckVisitor: '" << pv->GetName() << "' replica " << replicaNo
             << " at depth " << depth << " rebinds to logical volume '"
             << (boundLV != nullptr ? boundLV->GetName() : G4String("(null)")) << "'." << G4endl;
    }
  }

  G4Material* mat = (boundLV != nullptr) ? boundLV->GetMaterial() : nullptr;
  if (mat == nullptr)
  {
    // If this instance rebound the logical volume, the walker's SetSolid()/
    // UpdateMaterial() calls (issued before ComputeMaterial() ran) applied
    // to the *previous* logical volume, so reading the new one here is
    // indeterminate - not evidence of a genuinely null material.
    if (!rebound)
    {
      ++rep.fNullMaterialCount;
    }
  }
  else
  {
    rep.fMaterialsSeen.insert(mat);
  }

  const G4int histDepth = touchable->GetHistoryDepth();

  if (rep.fIsNested)
  {
    if (histDepth < 1)
    {
      ++rep.fShallowNestedCount;
    }
    else
    {
      G4VPhysicalVolume* parentPV = touchable->GetVolume(1);
      const G4int parentReplicaNo = touchable->GetReplicaNumber(1);
      fParentContexts[pv].insert(std::make_pair(parentPV, parentReplicaNo));
      if (mat != nullptr)
      {
        fMaterialByReplica[pv][replicaNo].insert(mat);
      }
    }

    // Declared-material cross-check: G4Region::ScanVolumeTree() trusts the
    // scanner's GetMaterial(0..N-1) list exclusively (it never calls
    // ComputeMaterial()) to build production-cut couples, so a material
    // returned here that is absent from that list means silently wrong cuts.
    if (mat != nullptr && param != nullptr)
    {
      G4VVolumeMaterialScanner* scanner = param->GetMaterialScanner();
      if (scanner != nullptr)
      {
        auto declMatIt = fDeclaredMaterials.find(pv);
        if (declMatIt == fDeclaredMaterials.end())
        {
          std::set<G4Material*> declared;
          for (G4int i = 0; i < scanner->GetNumberOfMaterials(); ++i)
          {
            declared.insert(scanner->GetMaterial(i));
          }
          declMatIt = fDeclaredMaterials.emplace(pv, std::move(declared)).first;
        }
        const std::set<G4Material*>& declared = declMatIt->second;
        G4bool found = (declared.find(mat) != declared.end());
        if (!found && mat->GetBaseMaterial() != nullptr)
        {
          found =
            (declared.find(const_cast<G4Material*>(mat->GetBaseMaterial())) != declared.end());
        }
        if (!found)
        {
          rep.fUndeclaredMaterials.insert(mat);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------

void G4ParameterisationCheckVisitor::VisitRegularStructure(const G4TouchableHistoryHandle& touchable,
                                                            G4PVParameterised* pv, G4int nInstances,
                                                            G4int depth)
{
  VolumeReport& rep = GetOrCreateReport(pv);
  rep.fIsRegularStructure = true;
  rep.fInstancesDeclared = nInstances;
  rep.fMinDepth = depth;
  rep.fMaxDepth = depth;

  G4VPVParameterisation* param = pv->GetParameterisation();
  rep.fIsNested = (param != nullptr) && param->IsNested();

  if (fVerbose)
  {
    G4cout << "G4ParameterisationCheckVisitor: '" << pv->GetName() << "' is a regular structure ("
           << nInstances << " instance(s), history depth=" << touchable->GetHistoryDepth()
           << ") - not individually exercised." << G4endl;
  }
}

// ---------------------------------------------------------------------------

G4bool G4ParameterisationCheckVisitor::DescendInto(const G4TouchableHistoryHandle& touchable,
                                                    G4LogicalVolume* lv, G4int /*depth*/)
{
  if (fMaxTotalVisits >= 0 && fTotalVisits >= fMaxTotalVisits)
  {
    fGlobalTruncated = true;
    return false;
  }

  if (fPruneUnparameterisedSubtrees && !SubtreeHasParameterised(lv))
  {
    return false;
  }

  G4VPhysicalVolume* curPV = touchable->GetVolume(0);
  if (curPV == nullptr)
  {
    return true;
  }

  // Sample only a few parent instances of a replicated/parameterised mother:
  // enough to give the C4 nesting-exercised check distinct parent contexts,
  // without paying for every combination down a deep hierarchy.
  if (fMaxParentSamples >= 0 && curPV->GetMultiplicity() > 1)
  {
    const G4int copyNo = touchable->GetReplicaNumber(0);
    const G4int mult = curPV->GetMultiplicity();
    const G4bool withinSampleWindow = (copyNo < fMaxParentSamples) || (copyNo == mult - 1);
    if (!withinSampleWindow)
    {
      fGlobalTruncated = true;
      return false;
    }
  }

  // Once a parameterised PV has been checked enough times, stop recursing
  // into its daughters. ComputeMaterial() etc. still ran for every declared
  // copy in VisitParameterised(); only the deeper recursion is capped.
  if (fMaxInstancesPerVolume >= 0 && curPV->IsParameterised())
  {
    auto it = fReports.find(curPV);
    if (it != fReports.end() && it->second.fInstancesVisited > fMaxInstancesPerVolume)
    {
      it->second.fTruncated = true;
      fGlobalTruncated = true;
      return false;
    }
  }

  return true;
}

// ---------------------------------------------------------------------------

G4bool G4ParameterisationCheckVisitor::SubtreeHasParameterised(G4LogicalVolume* lv)
{
  auto it = fSubtreeHasParamCache.find(lv);
  if (it != fSubtreeHasParamCache.end())
  {
    return it->second;
  }

  fSubtreeHasParamCache[lv] = false;  // guard against (disallowed) cycles

  G4bool found = false;
  const std::size_t nDaughters = lv->GetNoDaughters();
  for (std::size_t i = 0; i < nDaughters && !found; ++i)
  {
    G4VPhysicalVolume* d = lv->GetDaughter(i);
    if (d->IsParameterised())
    {
      found = true;
      break;
    }
    G4LogicalVolume* dLV = d->GetLogicalVolume();
    if (dLV != nullptr && SubtreeHasParameterised(dLV))
    {
      found = true;
    }
  }
  fSubtreeHasParamCache[lv] = found;
  return found;
}

// ---------------------------------------------------------------------------

G4bool G4ParameterisationCheckVisitor::HasProblem(const VolumeReport& rep) const
{
  return !rep.fReboundLVs.empty() || (rep.fNullMaterialCount > 0) ||
         !rep.fUndeclaredMaterials.empty() || (rep.fShallowNestedCount > 0) ||
         (rep.fIsNested && !rep.fIsRegularStructure && (rep.fInstancesVisited > 0) &&
          !rep.fMaterialVariesWithParent);
}

// ---------------------------------------------------------------------------

G4int G4ParameterisationCheckVisitor::CountProblems() const
{
  G4int n = 0;
  for (const auto& entry : fReports)
  {
    if (HasProblem(entry.second))
    {
      ++n;
    }
  }
  return n;
}

// ---------------------------------------------------------------------------

void G4ParameterisationCheckVisitor::PrintReport(std::ostream& os) const
{
  for (const auto& entry : fReports)
  {
    G4VPhysicalVolume* pv = entry.first;
    const VolumeReport& rep = entry.second;
    const G4bool problem = HasProblem(rep);
    if (!rep.fIsNested && !problem)
    {
      continue;  // uninteresting: a plain parameterisation with nothing to report
    }

    os << "G4ParameterisationCheckVisitor: volume '" << pv->GetName() << "'"
       << (rep.fIsNested ? " [nested]" : "")
       << (rep.fIsRegularStructure ? " [regular-structure]" : "")
       << (problem ? " ** PROBLEM **" : "") << "\n";

    if (rep.fIsRegularStructure)
    {
      os << "    regular structure, " << rep.fInstancesDeclared
         << " instance(s) not individually exercised; the flat "
            "ComputeMaterial(copyNo) form is used by G4RegularNavigation and "
            "is covered by G4ParameterisationChecker.\n";
      continue;
    }

    os << "    instances visited=" << rep.fInstancesVisited << "/" << rep.fInstancesDeclared
       << "  depth=[" << rep.fMinDepth << "," << rep.fMaxDepth << "]"
       << "  distinctLVsBound=" << rep.fBoundLVs.size()
       << "  distinctMaterials=" << rep.fMaterialsSeen.size();
    if (rep.fIsNested)
    {
      os << "  distinctParentContexts=" << rep.fDistinctParentContexts
         << "  materialVariesWithParent=" << (rep.fMaterialVariesWithParent ? "yes" : "no");
    }
    if (rep.fTruncated)
    {
      os << "  TRUNCATED";
    }
    os << "\n";

    if (!rep.fReboundLVs.empty())
    {
      os << "    rebinds to " << rep.fReboundLVs.size() << " other logical volume(s):";
      for (auto* lv : rep.fReboundLVs)
      {
        os << " '" << lv->GetName() << "'(mutable=" << (lv->IsMaterialMutable() ? "yes" : "no")
           << ")";
      }
      os << "\n";
    }
    if (rep.fNullMaterialCount > 0)
    {
      os << "    " << rep.fNullMaterialCount << " instance(s) with a null material.\n";
    }
    if (!rep.fUndeclaredMaterials.empty())
    {
      os << "    " << rep.fUndeclaredMaterials.size()
         << " material(s) returned by ComputeMaterial() are absent from the "
            "parameterisation's declared GetMaterial() list:";
      for (auto* mat : rep.fUndeclaredMaterials)
      {
        os << " '" << mat->GetName() << "'";
      }
      os << "\n";
    }
    if (rep.fShallowNestedCount > 0)
    {
      os << "    " << rep.fShallowNestedCount
         << " nested instance(s) visited at history depth 0 (no parent replica to key on).\n";
    }
  }
}

namespace
{
  // Adapts CheckGeometry() to the raw function-pointer signature that
  // G4GeometryManager::RegisterNestedParameterisationChecker() expects.
  G4int NestedParameterisationCheckerEntryPoint(G4VPhysicalVolume* startVol,
                                                G4bool repair, G4bool verbose)
  {
    return G4ParameterisationCheckVisitor::CheckGeometry(startVol, repair, verbose);
  }

  // Registers this module's stronger, nested-aware checker with
  // G4GeometryManager (in the management module, which must not depend on
  // navigation) so that UseNestedParameterisationCheckAtClose(true) can find
  // it, without management ever including a navigation header.
  const G4bool registeredNestedParameterisationChecker = []
  {
    G4GeometryManager::RegisterNestedParameterisationChecker(
      &NestedParameterisationCheckerEntryPoint);
    return true;
  }();
}  // namespace
