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
// G4GeometryMemoryAdvisor class implementation
//
// Author: John Apostolakis (CERN), July 2026
// --------------------------------------------------------------------

#include "G4GeometryMemoryAdvisor.hh"

#include "G4LogicalVolume.hh"
#include "G4LVData.hh"
#include "G4LVMutableData.hh"
#include "G4PVData.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VolumeTreeWalker.hh"
#include "G4Threading.hh"

#include <algorithm>
#include <iomanip>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <typeindex>
#include <vector>

// ---------------------------------------------------------------------------
// Two sets of size figures are used, and both are printed side by side so
// the report is trustworthy on whatever platform it is compiled/run on:
//
//   - "reference" values: measured once via a standalone sizeof() probe
//     against this source tree, Apple/LLVM libc++, arm64, 2026-07-20.
//     Platform/standard-library dependent (the classes contain a
//     std::string, a std::vector and/or a std::shared_ptr, whose sizes
//     differ between e.g. libc++ and libstdc++ and between 32/64-bit ABIs).
//     Kept only as a fixed point of prior experience to compare against.
//   - "this platform" values: plain sizeof() of the real classes, evaluated
//     by the compiler building this file - always exact for whatever
//     platform/standard-library/compiler is in use (portable standard
//     C++17; no allocator- or OS-specific introspection).
//
// All memory-savings figures below are computed from the "this platform"
// values; the reference values are shown only as a labelled comparison
// column. sizeof() necessarily excludes any heap memory a solid or other
// object owns internally (e.g. G4Polycone z-planes, tessellated facets,
// boolean constituents).
// ---------------------------------------------------------------------------

namespace
{
  // -- Reference (Apple/LLVM libc++, arm64, 2026-07-20) ----------------------
  constexpr std::size_t kRefLVSize = 176;            // sizeof(G4LogicalVolume)
  constexpr std::size_t kRefLVDataSize = 16;         // sizeof(G4LVData)
  constexpr std::size_t kRefLVMutableDataSize = 32;  // sizeof(G4LVMutableData)

  // Reference point only: what a single monolithic per-thread G4LVData
  // (SD+field+solid+material+mass+cuts, as in the pre-2026-07-20 default
  // build, before the runtime mutable/shared split existed) used to cost.
  // The legacy classes no longer exist in this source tree, so these two
  // figures cannot be re-measured "this platform" directly; see
  // kThisLegacyLVDataSize below for a same-platform-derived estimate.
  constexpr std::size_t kRefLegacyLVSize = 152;
  constexpr std::size_t kRefLegacyLVDataSize = 48;

  // -- This platform (sizeof(), evaluated by the compiler building this) ----
  constexpr std::size_t kThisLVSize = sizeof(G4LogicalVolume);
  constexpr std::size_t kThisLVDataSize = sizeof(G4LVData);
  constexpr std::size_t kThisLVMutableDataSize = sizeof(G4LVMutableData);
  constexpr std::size_t kThisPVDataSize = sizeof(G4PVData);

  // -- Legacy (pre-2026-07-20) monolithic layout, reconstructed for
  // measurement --------------------------------------------------------------
  //
  // The legacy classes no longer exist in this source tree (superseded by
  // the G4LVData/G4LVMutableData split), so their "this platform" size
  // cannot be obtained from the real classes. The two structs below are
  // private, reporting-only mimics of the pre-split layout - same member
  // types/order as the fields that used to live in the always-per-thread
  // monolithic G4LVData, or directly on G4LogicalVolume before the
  // mutable/shared split existed - using void* in place of any forward-
  // declared pointer member (identical size on every supported platform).
  // Never instantiated as functional objects; used only for sizeof().
  struct G4TempOldLVData
  {
    void* fSensitiveDetector = nullptr;
    void* fFieldManager = nullptr;
    void* fSolid = nullptr;
    void* fMaterial = nullptr;
    void* fCutsCouple = nullptr;
    G4double fMass = 0.0;
  };

  struct G4TempOldLogicalVolume
  {
    std::vector<void*> fDaughters;
    G4String fName;
    void* fUserLimits = nullptr;
    void* fVoxel = nullptr;
    G4double fSmartless = 2.0;
    void* fRegion = nullptr;
    G4double fBiasWeight = 1.0;
    std::shared_ptr<const void> fVisAttributes;
    G4int instanceID = 0;
    EVolume fDaughtersVolumeType = kNormal;
    G4bool fOptimise = true;
    G4bool fRootRegion = false;
    G4bool fLock = false;
    G4bool fIsUnchanged = true;
  };

  // The legacy monolithic G4LVData folded exactly the fields that are now
  // split between G4LVData and G4LVMutableData into one per-thread struct;
  // sizeof(G4TempOldLVData) measures this directly, "this platform".
  constexpr std::size_t kThisLegacyLVDataSize = sizeof(G4TempOldLVData);

  // sizeof(G4TempOldLogicalVolume) measures, "this platform", what
  // G4LogicalVolume itself used to cost once the six fields above are
  // excluded (they lived only in the per-thread monolithic G4LVData, never
  // as direct members of the legacy G4LogicalVolume).
  constexpr std::size_t kThisLegacyLVSize = sizeof(G4TempOldLogicalVolume);

  // Inserts ',' thousands separators into a non-negative integer, e.g.
  // 12000000 -> "12,000,000". Implemented by hand (rather than via
  // std::locale/std::numpunct, whose thousands-separator behaviour is
  // locale- and platform-dependent) to keep the output identical on every
  // supported platform.
  std::string FormatWithCommas(std::size_t value)
  {
    std::string digits = std::to_string(value);
    std::string result;
    int sinceSeparator = 0;
    for (auto it = digits.rbegin(); it != digits.rend(); ++it)
    {
      if (sinceSeparator > 0 && sinceSeparator % 3 == 0) result.push_back(',');
      result.push_back(*it);
      ++sinceSeparator;
    }
    std::reverse(result.begin(), result.end());
    return result;
  }

  constexpr G4int kCompositionLabelField = 54;
  constexpr G4int kCompositionValueField = 14;

  // Prints one "label: value" line of the geometry-composition summary,
  // with the ':' aligned across all lines and the value right-aligned
  // with thousands separators.
  void PrintCompositionLine(std::ostream& os, const std::string& label, std::size_t value,
                            G4int indent = 2)
  {
    std::string full(static_cast<std::size_t>(indent), ' ');
    full += label;
    os << std::left << std::setw(kCompositionLabelField) << full << ": " << std::right
       << std::setw(kCompositionValueField) << FormatWithCommas(value) << "\n";
  }

  void PrintLine(std::ostream& os, const char* label, G4double bytes)
  {
    os << "    " << label << ": " << (bytes / 1.0e6) << " MB\n";
  }

  void PrintSavingsLine(std::ostream& os, const char* label, const char* baselineLabel,
                        G4double bytes, G4double baselineBytes)
  {
    G4double savings = baselineBytes - bytes;
    G4double pct = (baselineBytes > 0.) ? (100. * savings / baselineBytes) : 0.;
    os << "    " << label << ": " << (baselineBytes / 1.0e6) << " MB total, " << (savings / 1.0e6)
       << " MB saved (" << pct << "% of " << baselineLabel << ")\n";
  }

  constexpr G4int kSizeTableLabelField = 24;
  constexpr G4int kSizeTableValueField = 19;

  // Prints the header of the reference-vs-this-platform structure-size table.
  void PrintSizeTableHeader(std::ostream& os)
  {
    os << "  " << std::left << std::setw(kSizeTableLabelField) << "Structure" << std::right
       << std::setw(kSizeTableValueField) << "reference(B)" << std::setw(kSizeTableValueField)
       << "this-platform(B)" << std::setw(kSizeTableValueField) << "Delta(B)" << "\n";
  }

  // Prints one row of the reference-vs-this-platform structure-size table.
  // 'hasRef'=false prints '-' in the reference column and 'n/a' for Delta
  // (structures, such as G4PVData, that have no reference-platform figure).
  void PrintSizeTableRow(std::ostream& os, const char* label, long long refBytes, G4bool hasRef,
                        long long thisBytes)
  {
    os << "  " << std::left << std::setw(kSizeTableLabelField) << label << std::right
       << std::setw(kSizeTableValueField) << (hasRef ? std::to_string(refBytes) : std::string("-"))
       << std::setw(kSizeTableValueField) << thisBytes;
    if (hasRef)
    {
      long long delta = thisBytes - refBytes;
      std::string deltaStr = ((delta >= 0) ? "+" : "") + std::to_string(delta);
      if (delta > 8 || delta < -8) deltaStr += " **";
      os << std::setw(kSizeTableValueField) << deltaStr;
    }
    else
    {
      os << std::setw(kSizeTableValueField) << "n/a";
    }
    os << "\n";
  }

  // -- Physical-volume-type -> {label, sizeof()} registry, for
  // SpaceOfPhysicalVolumes() -------------------------------------------------
  //
  // Keyed by std::type_index(typeid(*pv)), so it correctly distinguishes
  // G4PVParameterised from other kParameterised-reporting types. A concrete
  // physical-volume type not in this map is still counted (see
  // SpaceOfPhysicalVolumes()), just labelled "other" with
  // sizeof(G4VPhysicalVolume) as a lower bound.
  //
  // NOTE: G4PVDivision and G4ReplicatedSlice were deliberately omitted
  // (2026-07-28) to avoid a G4geomdivision module dependency; they report
  // as "other" until this is revisited - see G4navigation History.
  struct PVTypeInfo
  {
    const char* label;
    std::size_t bytes;
  };

  const std::map<std::type_index, PVTypeInfo>& PVTypeRegistry()
  {
    static const std::map<std::type_index, PVTypeInfo> registry = {
      {std::type_index(typeid(G4PVPlacement)), {"G4PVPlacement", sizeof(G4PVPlacement)}},
      {std::type_index(typeid(G4PVReplica)), {"G4PVReplica", sizeof(G4PVReplica)}},
      {std::type_index(typeid(G4PVParameterised)), {"G4PVParameterised", sizeof(G4PVParameterised)}},
    };
    return registry;
  }

  // Returns {label, bytes, isKnown} for the concrete (most-derived) type of
  // '*pv'; unregistered ("other"/user) types report sizeof(G4VPhysicalVolume)
  // as a lower bound (excludes derived-class members) with isKnown=false.
  struct PVSizeResult
  {
    G4String label;
    std::size_t bytes;
    G4bool known;
  };

  PVSizeResult PVSizeOf(const G4VPhysicalVolume& pv)
  {
    const auto& registry = PVTypeRegistry();
    auto it = registry.find(std::type_index(typeid(pv)));
    if (it != registry.end()) return {it->second.label, it->second.bytes, true};
    return {"other", sizeof(G4VPhysicalVolume), false};
  }
}  // namespace

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::Analyse(G4VPhysicalVolume* world)
{
  fLogicalVolumesSeen.clear();
  fMutableLogicalVolumesSeen.clear();
  fPlacementCount = 0;
  fReplicaInstanceCount = 0;
  fParamInstanceCount = 0;
  fUnknownCount = 0;

  if (world == nullptr) return;

  G4VolumeTreeWalker walker(*this);
  walker.Walk(world);
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::VisitPlacement(G4PVPlacement* pv, G4int /*depth*/)
{
  if (pv == nullptr) return;
  ++fPlacementCount;
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::VisitReplica(G4PVReplica* pv, G4int /*depth*/)
{
  if (pv == nullptr) return;
  G4int multiplicity = pv->GetMultiplicity();
  fReplicaInstanceCount += (multiplicity > 0) ? static_cast<std::size_t>(multiplicity) : 1;
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::VisitParameterised(G4PVReplica* pv, G4int /*depth*/)
{
  if (pv == nullptr) return;
  // 'pv' is actually a G4PVParameterised, G4PVDivision or G4ReplicatedSlice
  // (VolumeType() == kParameterised covers all three); G4PVReplica is their
  // real common base, so only members declared there (or in
  // G4VPhysicalVolume) may be called here - see G4VVolumeVisitor.hh.
  G4int multiplicity = pv->GetMultiplicity();
  fParamInstanceCount += (multiplicity > 0) ? static_cast<std::size_t>(multiplicity) : 1;
  fMutableLogicalVolumesSeen.insert(pv->GetLogicalVolume());
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::VisitUnknown(G4VPhysicalVolume* pv, G4int /*depth*/)
{
  if (pv == nullptr) return;
  ++fUnknownCount;
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::EnterLogicalVolume(G4LogicalVolume* lv, G4int /*depth*/)
{
  if (lv == nullptr) return;
  fLogicalVolumesSeen.insert(lv);
}

// ---------------------------------------------------------------------------

std::size_t G4GeometryMemoryAdvisor::GetTotalPhysicalVolumes() const
{
  return fPlacementCount + fReplicaInstanceCount + fParamInstanceCount + fUnknownCount;
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::PrintReport(std::ostream& os, G4int nThreads) const
{
  if (nThreads <= 0)
  {
    nThreads = G4Threading::G4GetNumberOfCores();
    if (nThreads <= 0) nThreads = 1;
  }

  std::size_t nLV = fLogicalVolumesSeen.size();
  std::size_t nMutableLV = fMutableLogicalVolumesSeen.size();
  std::size_t nSharedLV = nLV - nMutableLV;

  os << "\n========== G4GeometryMemoryAdvisor Report ==========\n";
  os << "Geometry composition:\n";
  PrintCompositionLine(os, "Distinct logical volumes", nLV);
  PrintCompositionLine(os, "- mutable (per-thread solid/material/mass/cuts)", nMutableLV, 4);
  PrintCompositionLine(os, "- shared  (across all threads)", nSharedLV, 4);
  PrintCompositionLine(os, "G4PVPlacement instances", fPlacementCount);
  PrintCompositionLine(os, "G4PVReplica instances (all copies)", fReplicaInstanceCount);
  PrintCompositionLine(os, "Parameterised/Division instances", fParamInstanceCount);
  if (fUnknownCount > 0)
  {
    PrintCompositionLine(os, "Other (external) volume instances", fUnknownCount);
  }

  ReportLogicalVolumeFootprint(os, nThreads);

  os << "=============================================================\n";
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::ReportLogicalVolumeFootprint(std::ostream& os, G4int nThreads) const
{
  std::size_t nLV = fLogicalVolumesSeen.size();
  std::size_t nMutableLV = fMutableLogicalVolumesSeen.size();

  os << "\nStructure sizes - reference (arm64/libc++, 2026-07-20) vs. this platform:\n";
  PrintSizeTableHeader(os);
  PrintSizeTableRow(os, "G4LogicalVolume", static_cast<G4long>(kRefLVSize), true,
                    static_cast<G4long>(kThisLVSize));
  PrintSizeTableRow(os, "G4LVData", static_cast<G4long>(kRefLVDataSize), true,
                    static_cast<G4long>(kThisLVDataSize));
  PrintSizeTableRow(os, "G4LVMutableData", static_cast<G4long>(kRefLVMutableDataSize), true,
                    static_cast<G4long>(kThisLVMutableDataSize));
  PrintSizeTableRow(os, "G4PVData", 0, false, static_cast<G4long>(kThisPVDataSize));
  PrintSizeTableRow(os, "Monolithic G4LVData", static_cast<G4long>(kRefLegacyLVDataSize), true,
                    static_cast<G4long>(kThisLegacyLVDataSize));
  PrintSizeTableRow(os, "Legacy G4LogicalVolume", static_cast<G4long>(kRefLegacyLVSize), true,
                    static_cast<G4long>(kThisLegacyLVSize));

  // One-time G4LogicalVolume master storage: allocated ONCE per logical
  // volume (holding the shared solid/material/mass/cuts-couple pointers
  // for non-mutable LVs), NOT duplicated per thread, so it must NOT be
  // added to a per-thread quantity below.
  G4double lvMasterBytes = static_cast<G4double>(nLV) * static_cast<G4double>(kThisLVSize);
  os << "\nOne-time G4LogicalVolume master storage (NOT per-thread, shared\n"
        "across all threads - excluded from the per-thread comparison below,\n"
        "which covers only what is genuinely duplicated per worker thread):\n";
  PrintLine(os, "Master storage                   ", lvMasterBytes);

  // Per-thread REPLICATED-data footprint: the fair, same-platform (all
  // sizes measured in this build) comparison of what genuinely gets
  // duplicated per worker thread. Legacy: every logical volume got one
  // full per-thread copy (solid+material+mass+cuts-couple, folded with
  // SD+field into one struct) regardless of mutability - there was no
  // shared/mutable distinction before 2026-07-20. Current: G4LVData is
  // always per-thread; G4LVMutableData is per-thread ONLY for the
  // nMutableLV logical volumes actually marked mutable.
  G4double currentPerThreadBytes = static_cast<G4double>(nThreads) *
    (static_cast<G4double>(nLV) * static_cast<G4double>(kThisLVDataSize)
     + static_cast<G4double>(nMutableLV) * static_cast<G4double>(kThisLVMutableDataSize));
  G4double legacyPerThreadBytes = static_cast<G4double>(nThreads) *
    static_cast<G4double>(nLV) * static_cast<G4double>(kThisLegacyLVDataSize);

  os << "\nLogical-volume-side per-thread REPLICATED-data footprint for "
     << nThreads << " worker thread(s):\n";
  PrintLine(os, "Current                          ", currentPerThreadBytes);
  PrintSavingsLine(os, "vs. legacy (pre-2026-07-20, always-per-thread) baseline", "legacy baseline",
                    currentPerThreadBytes, legacyPerThreadBytes);

  // Per-thread split-base overhead: the raw cost of the per-thread slot
  // arrays themselves (G4GeomSplitter<T>::offset), independent of the
  // "legacy vs current" comparison above. Each G4GeomSplitter<T> grows its
  // thread-local offset array in +512-slot chunks (see G4GeomSplitter::
  // Reallocate/PushBack), so actual allocated space can exceed nLV/nPV
  // slots by up to 512 entries per manager; that slack plus the
  // thread-local base pointer itself are not counted below.
  G4double splitBaseBytes = static_cast<G4double>(nThreads) *
    (static_cast<G4double>(nLV) * static_cast<G4double>(kThisLVDataSize)
     + static_cast<G4double>(nMutableLV) * static_cast<G4double>(kThisLVMutableDataSize)
     + static_cast<G4double>(fPlacementCount + fReplicaInstanceCount + fParamInstanceCount)
       * static_cast<G4double>(kThisPVDataSize));
  os << "\nPer-thread split-base overhead (G4LVData + G4LVMutableData + G4PVData slot arrays)"
        " for " << nThreads << " worker thread(s):\n";
  PrintLine(os, "Slot arrays (excl. +512-slot slack per manager, this-platform sizes)", splitBaseBytes);
}

// ---------------------------------------------------------------------------

void G4GeometryMemoryAdvisor::SpaceOfPhysicalVolumes(std::ostream& os, G4int nThreads) const
{
  if (nThreads <= 0)
  {
    nThreads = G4Threading::G4GetNumberOfCores();
    if (nThreads <= 0) nThreads = 1;
  }

  os << "\n========== G4GeometryMemoryAdvisor: SpaceOfPhysicalVolumes ==========\n";
  os << "Inventory of physical-volume instances in the global G4PhysicalVolumeStore, by type:\n";

  struct Row
  {
    G4String label;
    std::size_t count = 0;
    std::size_t bytesPerInstance = 0;
    G4bool known = true;
  };

  // Group by the concrete (most-derived) type via typeid(*pv), NOT by
  // VolumeType(), which reports kParameterised for G4PVParameterised,
  // G4PVDivision and G4ReplicatedSlice alike - typeid distinguishes them.
  std::map<G4String, Row> byLabel;
  for (const auto* pv : *G4PhysicalVolumeStore::GetInstance())
  {
    if (pv == nullptr) continue;
    PVSizeResult info = PVSizeOf(*pv);
    Row& row = byLabel[info.label];
    row.label = info.label;
    ++row.count;
    row.bytesPerInstance = info.bytes;
    row.known = info.known;
  }

  std::vector<Row> rows;
  rows.reserve(byLabel.size());
  std::size_t totalCount = 0;
  G4double totalMasterBytes = 0.;
  G4double totalPerThreadBytes = 0.;
  G4bool anyUnknown = false;

  for (const auto& [label, row] : byLabel)
  {
    rows.push_back(row);
    totalCount += row.count;
    totalMasterBytes += static_cast<G4double>(row.count) * static_cast<G4double>(row.bytesPerInstance);
    totalPerThreadBytes += static_cast<G4double>(row.count) * static_cast<G4double>(kThisPVDataSize)
                          * static_cast<G4double>(nThreads);
    if (!row.known) anyUnknown = true;
  }

  std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
    return (a.count * a.bytesPerInstance) > (b.count * b.bytesPerInstance);
  });

  constexpr G4int kLabelField = 22;
  constexpr G4int kValueField = 19;

  os << "  " << std::left << std::setw(kLabelField) << "Type" << std::right
     << std::setw(kValueField) << "Count" << std::setw(kValueField) << "Bytes/instance"
     << std::setw(kValueField) << "Master (MB)" << std::setw(kValueField) << "Extra/inst.(B)"
     << std::setw(kValueField) << "Extra total (MB)" << "\n";
  for (const auto& row : rows)
  {
    G4double masterTotal = static_cast<G4double>(row.count) * static_cast<G4double>(row.bytesPerInstance);
    G4double perThreadTotal = static_cast<G4double>(row.count) * static_cast<G4double>(kThisPVDataSize)
                            * static_cast<G4double>(nThreads);
    G4String bytesStr = std::to_string(row.bytesPerInstance) + (row.known ? "" : "*");
    os << "  " << std::left << std::setw(kLabelField) << row.label << std::right
       << std::setw(kValueField) << FormatWithCommas(row.count) << std::setw(kValueField) << bytesStr
       << std::setw(kValueField) << (masterTotal / 1.0e6) << std::setw(kValueField) << kThisPVDataSize
       << std::setw(kValueField) << (perThreadTotal / 1.0e6) << "\n";
  }
  os << "  --------------------------------------------------------------------------------------\n";
  os << "  TOTAL: " << FormatWithCommas(totalCount) << " physical volumes, " << (totalMasterBytes / 1.0e6)
     << " MB (master), " << (totalPerThreadBytes / 1.0e6) << " MB (per-thread G4PVData, x"
     << nThreads << " threads)\n";
  if (anyUnknown)
  {
    os << "  (* type unknown to the registry, labelled \"other\" - sizeof(G4VPhysicalVolume) used as"
          " a lower bound; add the type to PVTypeRegistry() in G4GeometryMemoryAdvisor.cc for an"
          " exact figure.)\n";
  }
  os << "  Note: every physical volume carries one G4PVData slot per thread\n"
        "  (rotation + translation), but only parameterised/replica/division\n"
        "  placements actually vary it between threads - see G4PVData.hh.\n";

  os << "=============================================================\n";
}
