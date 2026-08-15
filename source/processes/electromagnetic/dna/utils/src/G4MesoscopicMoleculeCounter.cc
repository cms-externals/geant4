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

#include "G4MesoscopicMoleculeCounter.hh"

#include "G4DNAMesh.hh"
#include "G4Molecule.hh"
#include "G4MoleculeCounterTemplates.hh"
#include "G4UnitsTable.hh"

//------------------------------------------------------------------------------

G4String G4MesoscopicMoleculeCounterIndex::GetInfo() const
{
  auto molecule_name = Molecule != nullptr ? Molecule->GetName() : G4String();
  return "Molecule: " + molecule_name
         + " Resolution (nm): " + std::to_string(std::round(resolution / CLHEP::nm))
         + " VoxelLinear: " + std::to_string(linear);
}

//------------------------------------------------------------------------------

G4MesoscopicMoleculeCounter::G4MesoscopicMoleculeCounter()
  : G4VUserMoleculeCounter("", MoleculeCounterType::Mesoscopic)
{}

G4MesoscopicMoleculeCounter::G4MesoscopicMoleculeCounter(G4String name)
  : G4VUserMoleculeCounter(name, MoleculeCounterType::Mesoscopic)
{}

//------------------------------------------------------------------------------

void G4MesoscopicMoleculeCounter::InitializeUser() {}

//------------------------------------------------------------------------------

void G4MesoscopicMoleculeCounter::SetMeshSnapshot(const G4DNAMesh* const mesh, G4double time)
{
  if (mesh == nullptr)
  {
    G4ExceptionDescription errMsg;
    errMsg << "Mesh pointer is nullptr!" << G4endl;
    G4Exception("G4MesoscopicMoleculeCounter::SetMeshSnapshot", "BAD_REFERENCE", FatalException,
                errMsg);
    return;
  }

  if (!IsActiveAtGlobalTime(time)) return;

  auto const& domain = mesh->GetBoundingBox();
  const G4float res = mesh->GetResolution();
  const auto nx = static_cast<G4int>(std::floor((domain.Getxhi() - domain.Getxlo()) / res));
  const auto ny = static_cast<G4int>(std::floor((domain.Getyhi() - domain.Getylo()) / res));

  if (fVerbose > 1)
  {
    G4cout << "G4MesoscopicMoleculeCounter(" << GetName()
           << ")::SetMeshSnapshot at time: " << G4BestUnit(time, "Time")
           << " (Resolution: " << G4BestUnit(res, "Length") << ")" << G4endl;
  }

  for (auto it_mesh = mesh->const_begin(); it_mesh != mesh->const_end(); ++it_mesh)
  {
    auto const& voxelIdx = std::get<0>(*it_mesh);
    auto const& data = std::get<2>(*it_mesh);

    if (data.empty()) continue;

    const G4int linear = voxelIdx.x + nx * voxelIdx.y + nx * ny * voxelIdx.z;

    for (auto const& [mol, count] : data)
    {
      if (IsReactantIgnored(mol) || count == 0) continue;

      auto idx = fRecordMeshData ? G4MesoscopicMoleculeCounterIndex(mol, res, linear)
                                 : G4MesoscopicMoleculeCounterIndex(mol);
      auto [it, indexIsNew] = fCounterMap.emplace(idx, InnerCounterMapType{fTimeComparer});

      // G4DNAMesh stores counts as size_t, for now we cast this to G4int and throw on <0
      auto count_ = static_cast<G4int>(count);
      if (count_ < 0)
      {
        G4ExceptionDescription desc;
        desc << "The count `" << count
             << "` is too large to fit into the G4int used by the G4MesoscopicMoleculeCounter!"
             << G4endl;
        G4Exception("G4MesoscopicMoleculeCounter::SetMeshSnapshot", "OVERFLOW",
                    FatalErrorInArgument, desc);
      }

      // NO CONSISTENCY CHECKS [TBD]
      // there will already be numbers in the map, so...
      // (1) find the closest time
      // (2) emplace entry using closest value as init + number
      // (3) add number to all "future" entries as well
      if (it->second.empty())
      {
        it->second.emplace(time, count_);
      }
      else
      {  // at least one element exists, so we can try to find the closest key
        auto it_closest = G4::MoleculeCounter::FindClosestEntryForKey(it->second, time);
        auto [it_time, _] = it->second.emplace(time, it_closest->second);
        do
        {
          it_time->second += count_;
        } while (++it_time != it->second.end());
      }
    }
  }
}

//------------------------------------------------------------------------------

G4int G4MesoscopicMoleculeCounter::GetNbMoleculesAtTime(const G4MesoscopicMoleculeCounterIndex& idx,
                                                        G4double time) const
{
  if (GetVerbose() > 1)
  {
    auto it = fCounterMap.find(idx);
    if (it == fCounterMap.cend())
    {
      G4Exception("G4MesoscopicMoleculeCounter::GetNbMoleculesAtTime", "MOLCOUNTER", JustWarning,
                  "The index provided was not found in the molecule counter map! "
                  "This will return a zero count.\n\n"
                  "Note that extracting molecule count values from the mesosocpic counters "
                  "is non-trivial -- double check everything you receive as these interfaces "
                  "may change in a future version!");
      return 0;
    }

    std::vector<G4double> recordedTimes = G4::MoleculeCounter::GetMapIndices(it->second);
    auto [tmin_it, tmax_it] = std::minmax_element(recordedTimes.begin(), recordedTimes.end());
    if (time < *tmin_it || time > *tmax_it)
    {
      G4ExceptionDescription desc;
      desc << "The search time " << G4BestUnit(time, "Time") << " is before or after the boundary "
           << "times for this index! The return value is determined by:\n"
           << " * t < " << G4BestUnit(*tmin_it, "Time") << "->0 \n"
           << " * t > " << G4BestUnit(*tmax_it, "Time") << "-> count@t_max\n\n"
           << "Note that extracting molecule count values from the mesosocpic counters "
              "is non-trivial -- double check everything you receive as these interfaces "
              "may change in a future version!"
           << G4endl;
      G4Exception("G4MesoscopicMoleculeCounter::GetNbMoleculesAtTime", "MOLCOUNTER", JustWarning,
                  desc);
    }
  }

  return G4VUserMoleculeCounter<G4MesoscopicMoleculeCounterIndex>::GetNbMoleculesAtTime(idx, time);
}

G4int G4MesoscopicMoleculeCounter::GetNbMoleculesAtTime(Search& search,
                                                        const G4MesoscopicMoleculeCounterIndex& idx,
                                                        G4double time) const
{
  return G4VUserMoleculeCounter::GetNbMoleculesAtTime(search, idx, time);
}

//------------------------------------------------------------------------------
