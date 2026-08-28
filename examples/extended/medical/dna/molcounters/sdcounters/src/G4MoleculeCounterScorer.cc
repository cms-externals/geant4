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
/// \file ScoreBasicMoleculeCounts.cc
/// \brief Implementation of the ScoreBasicMoleculeCounts class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025~26)
//

#include "G4MoleculeCounterScorer.hh"

#include "G4AnalysisManager.hh"

G4MoleculeCounterScorer::G4MoleculeCounterScorer(G4String name, G4String counterName, G4int depth)
  : G4VPrimitiveMoleculeScorer<G4MoleculeCounter, G4MoleculeCounterIndex, G4String>(
      name, counterName, depth)
{}

std::tuple<G4String>
G4MoleculeCounterScorer::ConvertCounterIndexToKey(const G4MoleculeCounterIndex& idx) const
{
  return std::make_tuple(idx.GetMolecule()->GetName());
}

void G4MoleculeCounterScorer::WriteWithAnalysisManager() const
{
  if (G4Threading::IsWorkerThread()) return;

  auto analysisManager = G4AnalysisManager::Instance();

  int fNtupleID = analysisManager->CreateNtuple(GetName() + "_" + GetCounterName(),
                                                GetName() + "_" + GetCounterName());

  analysisManager->CreateNtupleDColumn("Time__ps_");
  analysisManager->CreateNtupleDColumn("Molecule_Count");
  analysisManager->CreateNtupleSColumn("Molecule_Name");
  analysisManager->FinishNtuple();

  for (auto const& [key, time_map] : fCountPerIndexPerTime)
  {
    auto moleculeName = std::get<0>(key);
    G4cout << "(" << moleculeName << "):" << G4endl;

    for (auto const& [time, count] : time_map)
    {
      G4cout << "  t=" << G4BestUnit(time, "Time") << ", n=" << count << G4endl;

      analysisManager->FillNtupleDColumn(fNtupleID, 0, time / ps);
      analysisManager->FillNtupleDColumn(fNtupleID, 1, count);
      analysisManager->FillNtupleSColumn(fNtupleID, 2, moleculeName);
      analysisManager->AddNtupleRow(fNtupleID);
    }
  }
}
