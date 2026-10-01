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
// * regarding  this  software system or assume any liability, express *
// * or implied, for use of this software.  Please see the license in  *
// * the file LICENSE and URL above for the full disclaimer and the    *
// * limitation of liability.                                         *
// *                                                                  *
// * This  code  implementation is the result of the scientific and    *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Benchmark G4Tubs/G4Cons navigation queries for all topology tags.

#include "G4Cons.hh"
#include "G4CutTubs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kDataSize = 8192;
volatile G4double gSink = 0.;

struct Topology
{
  const char* name;
  G4double dphi;
  G4bool hollow;
};

struct Data
{
  std::vector<G4ThreeVector> inside, outside, surface, directionIn, directionOut;
};

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> u(0., 1.);
  const auto z = 2. * u(engine) - 1.;
  const auto phi = twopi * u(engine);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

Data MakeData(G4bool isCons, const Topology& topology)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> u(0., 1.);
  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  const G4double rmin1 = topology.hollow ? 4. * cm : 0.;
  const G4double rmin2 = topology.hollow ? (isCons ? 7. * cm : 4. * cm) : 0.;
  const G4double rmax1 = 12. * cm;
  const G4double rmax2 = isCons ? 20. * cm : 12. * cm;
  const G4double dz = 15. * cm;
  const G4double sphi = -35. * deg;

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto z = (-.7 + 1.4 * u(engine)) * dz;
    const auto t = (z / dz + 1.) * .5;
    const auto rmin = rmin1 + t * (rmin2 - rmin1);
    const auto rmax = rmax1 + t * (rmax2 - rmax1);
    const auto phi = sphi + (.08 + .84 * u(engine)) * topology.dphi;
    const auto radius = rmin + (.15 + .7 * u(engine)) * (rmax - rmin);
    const G4ThreeVector pin(radius * std::cos(phi), radius * std::sin(phi), z);
    const G4ThreeVector psurf(rmax * std::cos(phi), rmax * std::sin(phi), z);
    const G4ThreeVector pout(2.5 * rmax * std::cos(phi), 2.5 * rmax * std::sin(phi), z);
    data.inside.push_back(pin);
    data.surface.push_back(psurf);
    data.outside.push_back(pout);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

template<class Function>
G4double Time(const std::uint64_t calls, Function function)
{
  G4double sum = 0.;
  const auto begin = std::chrono::steady_clock::now();
  for (std::uint64_t i = 0; i < calls; ++i) sum += function(i % kDataSize);
  const auto end = std::chrono::steady_clock::now();
  gSink += sum;
  return std::chrono::duration<G4double, std::milli>(end - begin).count();
}

void PrintTime(const char* label, G4double ms)
{
  std::cout << "  " << std::left << std::setw(34) << label << std::right << std::setw(10) << ms
            << " ms\n";
}

void Run(const char* shape, const Topology& topology, std::uint64_t calls)
{
  const G4bool isCons = std::string(shape) == "G4Cons";
  const G4bool isCutTubs = std::string(shape) == "G4CutTubs";
  const auto rmin1 = topology.hollow ? 4. * cm : 0.;
  const auto rmin2 = topology.hollow ? (isCons ? 7. * cm : 4. * cm) : 0.;
  std::unique_ptr<G4VSolid> solid;
  if (isCons)
    solid = std::make_unique<G4Cons>(topology.name, rmin1, 12. * cm, rmin2, 20. * cm,
                                     15. * cm, -35. * deg, topology.dphi);
  else if (isCutTubs)
    solid = std::make_unique<G4CutTubs>(topology.name, rmin1, 12. * cm, 15. * cm,
                                        -35. * deg, topology.dphi,
                                        G4ThreeVector(.1, .05, -1.).unit(),
                                        G4ThreeVector(-.05, .1, 1.).unit());
  else
    solid = std::make_unique<G4Tubs>(topology.name, rmin1, 12. * cm, 15. * cm,
                                     -35. * deg, topology.dphi);

  std::cout << '\n' << shape << "  topology=" << topology.name << '\n';

  const auto data = MakeData(isCons, topology);
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    if (solid->Inside(data.inside[i]) != kInside || solid->Inside(data.outside[i]) != kOutside
        || solid->Inside(data.surface[i]) != kSurface)
    {
      std::cerr << "invalid generated point for " << shape << '/' << topology.name << '\n';
      std::exit(3);
    }
  }
  const auto emit = [&](const char* query, auto function) {
    PrintTime(query, Time(calls, function));
  };

  // A short untimed pass faults code/data pages in before measurement.
  for (std::size_t i = 0; i < kDataSize; ++i)
    gSink += static_cast<G4double>(solid->Inside(data.inside[i]));
  emit("Inside(p)", [&](auto i) { return G4double(solid->Inside(data.inside[i])); });
  emit("SurfaceNormal(p)", [&](auto i) { return solid->SurfaceNormal(data.surface[i]).x(); });
  emit("DistanceToIn(p)", [&](auto i) { return solid->DistanceToIn(data.outside[i]); });
  emit("DistanceToIn(p,v)", [&](auto i) {
    return solid->DistanceToIn(data.outside[i], data.directionIn[i]);
  });
  emit("DistanceToOut(p)", [&](auto i) { return solid->DistanceToOut(data.inside[i]); });
  emit("DistanceToOut(p,v)", [&](auto i) {
    return solid->DistanceToOut(data.inside[i], data.directionOut[i]);
  });
  emit("DistanceToOut(p,v,n)", [&](auto i) {
    G4bool valid = false;
    G4ThreeVector normal;
    return solid->DistanceToOut(data.inside[i], data.directionOut[i], true, &valid, &normal)
           + normal.x() + valid;
  });
}
}  // namespace

int main(int argc, char** argv)
{
  const std::uint64_t calls = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 1000000;
  if (calls == 0)
  {
    std::cerr << "calls must be positive\n";
    return 2;
  }
  const Topology topologies[] = {{"solid-full", twopi, false},
                                 {"solid-less-pi", 90. * deg, false},
                                 {"solid-pi", pi, false},
                                 {"solid-greater-pi", 270. * deg, false},
                                 {"hollow-full", twopi, true},
                                 {"hollow-less-pi", 90. * deg, true},
                                 {"hollow-pi", pi, true},
                                 {"hollow-greater-pi", 270. * deg, true}};
  std::cout << "*********************************************************************\n";
  std::cout << "* Benchmark for G4Tubs, G4CutTubs and G4Cons                        *\n";
  std::cout << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *\n";
  std::cout << "* Configurations: full/partial phi, zero/nonzero inner radius        *\n";
  std::cout << "*********************************************************************\n";
  std::cout << "\n       Number of calls per method : " << calls << '\n';
  std::cout << "       Number of sampled points   : " << kDataSize << '\n';
  for (const auto& topology : topologies)
  {
    Run("G4Tubs", topology, calls);
    Run("G4Cons", topology, calls);
    Run("G4CutTubs", topology, calls);
  }
  return gSink == -1. ? 1 : 0;
}
