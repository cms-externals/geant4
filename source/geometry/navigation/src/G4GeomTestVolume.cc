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
// Class G4GeomTestVolume implementation
//
// Author: Gabriele Cosmo (CERN), 22 August 2013
// --------------------------------------------------------------------

#include "G4GeomTestVolume.hh"

#include "G4GeomTaskHelper.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiUnion.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

#include <cstdint>
#include <functional>
#include <queue>
#include <set>
#include <vector>

namespace
{
  struct OverlapCheckEntry
  {
    G4VPhysicalVolume* volume = nullptr;
    G4PVPlacement* placement = nullptr;
    G4bool overlap = false;
    G4bool checked = false;
    G4bool completed = true;
    std::vector<G4String> warnings;
  };

  G4bool NeedsSurfaceWarmup(const G4VSolid* solid)
  {
    const auto& type = solid->GetEntityType();
    return type == "G4BooleanSolid" || type == "G4UnionSolid" || type == "G4SubtractionSolid"
           || type == "G4IntersectionSolid" || type == "G4GenericPolycone" || type == "G4Polycone"
           || type == "G4Polyhedra";
  }

  void WarmSurfaceSampling(G4VSolid* solid, std::set<G4VSolid*>& warmedSolids)
  {
    if (solid->GetEntityType() == "G4MultiUnion")
    {
      auto* multiUnion = static_cast<G4MultiUnion*>(solid);
      for (G4int i = 0; i < multiUnion->GetNumberOfSolids(); ++i)
      {
        WarmSurfaceSampling(multiUnion->GetSolid(i), warmedSolids);
      }
    }
    if (NeedsSurfaceWarmup(solid) && warmedSolids.insert(solid).second)
    {
      G4GeomTaskHelper::SeedQuickRand(G4GeomTaskHelper::NextSeed());
      solid->GetPointOnSurface();
    }
  }
}  // namespace

//
// Constructor
//
G4GeomTestVolume::G4GeomTestVolume(G4VPhysicalVolume* theTarget, G4double theTolerance,
                                   G4int numberOfPoints, G4bool theVerbosity)
  : target(theTarget), tolerance(theTolerance), resolution(numberOfPoints), verbosity(theVerbosity)
{
  ;
}

//
// Destructor
//
G4GeomTestVolume::~G4GeomTestVolume()
{
  ;
}

//
// Get error tolerance
//
G4double G4GeomTestVolume::GetTolerance() const
{
  return tolerance;
}

//
// Set error tolerance
//
void G4GeomTestVolume::SetTolerance(G4double tol)
{
  tolerance = tol;
}

//
// Get number of points to check (resolution)
//
G4int G4GeomTestVolume::GetResolution() const
{
  return resolution;
}

//
// Set number of points to check (resolution)
//
void G4GeomTestVolume::SetResolution(G4int np)
{
  resolution = np;
}

//
// Get verbosity
//
G4bool G4GeomTestVolume::GetVerbosity() const
{
  return verbosity;
}

//
// Set verbosity
//
void G4GeomTestVolume::SetVerbosity(G4bool verb)
{
  verbosity = verb;
}

//
// Get errors reporting threshold
//
G4int G4GeomTestVolume::GetErrorsThreshold() const
{
  return maxErr;
}

//
// Set maximum number of errors to report
//
void G4GeomTestVolume::SetErrorsThreshold(G4int max)
{
  maxErr = max;
}

//
// Test overlap in tree
//
void G4GeomTestVolume::TestOverlapInTree() const
{
  std::queue<G4VPhysicalVolume*> volumes;
  std::set<G4LogicalVolume*> checked;
  std::vector<G4VPhysicalVolume*> checks;

  volumes.push(target);
  while (!volumes.empty())
  {
    G4VPhysicalVolume* current = volumes.front();
    volumes.pop();

    // check overlaps for daughters
    G4LogicalVolume* logical = current->GetLogicalVolume();
    std::size_t ndaughters = logical->GetNoDaughters();
    for (std::size_t i = 0; i < ndaughters; ++i)
    {
      G4VPhysicalVolume* daughter = logical->GetDaughter(i);
      checks.push_back(daughter);
    }

    // append the queue of volumes
    G4LogicalVolume* previousLogical = nullptr;
    for (std::size_t i = 0; i < ndaughters; ++i)
    {
      G4VPhysicalVolume* daughter = logical->GetDaughter(i);
      G4LogicalVolume* daughterLogical = daughter->GetLogicalVolume();
      if (daughterLogical->GetNoDaughters() == 0)
      {
        continue;
      }
      G4bool found = (daughterLogical == previousLogical);
      if (!found)
      {
        found = (checked.find(daughterLogical) != checked.cend());
      }
      if (!found)
      {
        checked.emplace(daughterLogical);
        previousLogical = daughterLogical;
        volumes.push(daughter);
      }
      else
      {
        if (verbosity)
        {
          G4cout << "Checking overlaps in tree of volume " << daughter->GetName() << " ("
                 << daughterLogical->GetSolid()->GetEntityType() << ")"
                 << " is omitted, to avoid duplication" << G4endl;
        }
      }
    }
  }
  CheckVolumes(checks);
}

//
// TestRecursiveOverlap
//
void G4GeomTestVolume::TestRecursiveOverlap(G4int slevel, G4int depth)
{
  std::vector<G4VPhysicalVolume*> checks;
  std::function<void(G4VPhysicalVolume*, G4int, G4int)> visit;
  visit = [&](G4VPhysicalVolume* current, G4int currentLevel, G4int currentDepth) {
    if (currentDepth == 0)
    {
      return;
    }
    if (currentDepth != -1)
    {
      --currentDepth;
    }
    if (currentLevel != 0)
    {
      --currentLevel;
    }

    if (currentLevel == 0)
    {
      checks.push_back(current);
    }

    const G4LogicalVolume* logical = current->GetLogicalVolume();
    auto nDaughter = (G4int)logical->GetNoDaughters();
    for (auto iDaughter = 0; iDaughter < nDaughter; ++iDaughter)
    {
      visit(logical->GetDaughter(iDaughter), currentLevel, currentDepth);
    }
  };
  visit(target, slevel, depth);

  CheckVolumes(checks);
}

//
// CheckVolumes
//
void G4GeomTestVolume::CheckVolumes(const std::vector<G4VPhysicalVolume*>& volumes) const
{
  if (volumes.empty())
  {
    return;
  }

  std::vector<OverlapCheckEntry> entries(volumes.size());
  for (std::size_t i = 0; i < volumes.size(); ++i)
  {
    entries[i].volume = volumes[i];
    entries[i].placement = dynamic_cast<G4PVPlacement*>(volumes[i]);
  }

  const auto forceSerial = G4GeometryManager::GetInstance()->IsSerialOverlapsCheckForced();
  G4GeomTaskHelper::PoolContext taskContext(forceSerial ? 0 : volumes.size(), true);
  G4bool useParallel = !forceSerial && taskContext.IsValid();
  if (useParallel)
  {
    std::set<G4VSolid*> warmedSolids;

    for (const auto& entry : entries)
    {
      if (entry.placement == nullptr)
      {
        continue;
      }

      G4LogicalVolume* mother = entry.placement->GetMotherLogical();
      if (mother == nullptr)
      {
        continue;
      }

      WarmSurfaceSampling(entry.placement->GetLogicalVolume()->GetSolid(), warmedSolids);
      WarmSurfaceSampling(mother->GetSolid(), warmedSolids);
      for (std::size_t i = 0; i < mother->GetNoDaughters(); ++i)
      {
        WarmSurfaceSampling(mother->GetDaughter((G4int)i)->GetLogicalVolume()->GetSolid(),
                            warmedSolids);
      }
    }
  }

  std::uint64_t baseSeed = useParallel ? G4GeomTaskHelper::NextSeed() : 0;

  auto compute = [&](std::size_t i)
  {
    if (entries[i].placement == nullptr)
    {
      return;
    }
    entries[i].overlap = entries[i].placement->ComputeOverlaps(
      resolution, tolerance, maxErr,
      useParallel ? G4GeomTaskHelper::MixSeed(baseSeed + i) : 0,
      entries[i].checked, entries[i].completed, entries[i].warnings);
  };

  if (useParallel)
  {
    G4GeomTaskHelper::ParallelForWithWorkerContext(taskContext, entries.size(), 1,
      [&](std::size_t begin, std::size_t end, std::size_t)
      {
        for (std::size_t i = begin; i < end; ++i)
        {
          compute(i);
        }
      },
      [](std::size_t, const auto& run)
      {
        G4GeomTaskHelper::ScopedGeometryWorkspaces workspaces;
        run();
      });
  }
  else
  {
    for (std::size_t i = 0; i < entries.size(); ++i)
    {
      compute(i);
    }
  }

  for (const auto& entry : entries)
  {
    if (entry.placement != nullptr)
    {
      entry.placement->ReportOverlaps(entry.overlap, entry.checked,
                                      entry.completed, entry.warnings, verbosity);
    }
    else
    {
      entry.volume->CheckOverlaps(resolution, tolerance, verbosity, maxErr);
    }
  }
}
