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
//---------------------------------------------------------------------------//
// Micro benchmarks for G4MuonDecayChannel::DecayIt.
//---------------------------------------------------------------------------//

#include "G4MuonDecayChannel.hh"

#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4DecayProducts.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4Positron.hh"

#include "CLHEP/Random/Random.h"
#include <benchmark/benchmark.h>

#include <algorithm>

namespace
{
/** Ensure all particles needed by the decay channel are in the particle table.*/
void EnsureParticleDefinitions()
{
  G4MuonMinus::Definition();
  G4MuonPlus::Definition();
  G4Electron::Definition();
  G4Positron::Definition();
  G4NeutrinoE::Definition();
  G4AntiNeutrinoE::Definition();
  G4NeutrinoMu::Definition();
  G4AntiNeutrinoMu::Definition();
}

static void BM_MuonDecayIt_AtRest(benchmark::State& state)
{
  const auto caseId = state.range(0);
  const char* parentName = nullptr;
  long seed = 0;

  switch (caseId)
  {
    case 0:
      parentName = "mu-";
      seed = 1234567;
      break;
    case 1:
      parentName = "mu+";
      seed = 7654321;
      break;
    default:
      state.SkipWithError("Unknown benchmark case id");
      return;
  }

  EnsureParticleDefinitions();
  CLHEP::HepRandom::setTheSeed(seed);
  state.SetLabel(parentName);

  G4MuonDecayChannel channel(parentName, 1.0);

  for (auto _ : state)
  {
    G4DecayProducts* products = channel.DecayIt(0.0);
    benchmark::DoNotOptimize(products);
    delete products;
  }

  state.SetItemsProcessed(state.iterations());
}

/** Parametrized benchmark over mu-/mu+ */
BENCHMARK(BM_MuonDecayIt_AtRest)
  ->Arg(0)
  ->Arg(1)
  ->Unit(benchmark::kNanosecond)
  ->Repetitions(20)
  ->ReportAggregatesOnly(true)
  ->ComputeStatistics("min",
                      [](const std::vector<double>& values) {
                        return *std::min_element(values.begin(), values.end());
                      })
  ->ComputeStatistics("max", [](const std::vector<double>& values) {
    return *std::max_element(values.begin(), values.end());
  });

}  // namespace

BENCHMARK_MAIN();
