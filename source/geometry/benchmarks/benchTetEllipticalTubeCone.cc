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
// * technical work of the Geant4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Benchmark G4Tet/G4EllipticalTube/G4EllipticalCone navigation queries.
// Analytic, backend-independent samples; balanced surface classes, not area weighted.

#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"
#include "G4Tet.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "globals.hh"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace
{
constexpr std::size_t kDataSize = 8192;
volatile G4double gSink = 0.;

struct Data
{
    std::vector<G4ThreeVector> inside, outside, surface, directionIn, directionOut;
};

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  const auto z = 2. * uniform(engine) - 1.;
  const auto phi = twopi * uniform(engine);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

// Elliptical caps and lateral wall receive equal sample counts. At an apex,
// only the lower cap and lateral wall have a well-defined surface normal.
Data MakeEllipseData(G4double a, G4double b, G4double dz, G4double height = 0.)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  Data data;
  const auto point = [&](G4double z, G4double phi, G4double radius) {
    const auto scale = height == 0. ? 1. : height - z;
    return G4ThreeVector(a * scale * radius * std::cos(phi),
                         b * scale * radius * std::sin(phi), z);
  };
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto z = (.3 + 1.4 * uniform(engine) - 1.) * dz;
    const auto phi = twopi * uniform(engine);
    const auto pin = point(z, phi, .2 + .6 * uniform(engine));
    auto pout = point(z, phi, 1.4);
    if (i % 3 != 0) pout.setZ(i % 3 == 1 ? -1.6 * dz : 1.6 * dz);
    const auto face = i % (height == dz ? 2 : 3);
    const auto end = face == 1 ? -dz : dz;
    const auto surface = face == 0 ? point(z, phi, 1.)
                                  : point(end, phi, .2 + .6 * uniform(engine));
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(surface);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

using Vertices = std::array<G4ThreeVector, 4>;
Data MakeTetData(const Vertices& vertices)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(.2, 1.);
  const auto centre = .25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3]);
  Data data;
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    // Positive barycentric coordinates keep samples away from edges/vertices.
    const auto face = i % 4;
    G4ThreeVector surface;
    G4double total = 0.;
    for (std::size_t j = 0; j < 4; ++j)
    {
      if (j == face) continue;
      const auto weight = uniform(engine);
      surface += weight * vertices[j];
      total += weight;
    }
    surface /= total;
    const auto pin = centre + .55 * (surface - centre);
    const auto pout = centre + 1.5 * (surface - centre);
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(surface);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

void Validate(const G4VSolid& solid, const Data& data, const std::string& shape,
              const std::string& topology)
{
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto inside = solid.Inside(data.inside[i]);
    const auto outside = solid.Inside(data.outside[i]);
    const auto surface = solid.Inside(data.surface[i]);
    const auto checkQueries = (i % 127) == 0;
    const auto normal = checkQueries ? solid.SurfaceNormal(data.surface[i]) : G4ThreeVector();
    const auto safetyIn = checkQueries ? solid.DistanceToIn(data.outside[i]) : 0.;
    const auto distanceIn =
      checkQueries ? solid.DistanceToIn(data.outside[i], data.directionIn[i]) : 0.;
    const auto safetyOut = checkQueries ? solid.DistanceToOut(data.inside[i]) : 0.;
    const auto distanceOut =
      checkQueries ? solid.DistanceToOut(data.inside[i], data.directionOut[i]) : 0.;
    G4bool validOutNormal = false;
    G4ThreeVector outNormal;
    const auto distanceOutWithNormal = checkQueries
                                         ? solid.DistanceToOut(data.inside[i], data.directionOut[i],
                                                               true, &validOutNormal, &outNormal)
                                         : 0.;
    // Exact distances must reach their respective boundaries.
    const auto entry = data.outside[i] + distanceIn * data.directionIn[i];
    const auto exit = data.inside[i] + distanceOut * data.directionOut[i];
    // All three solid families are convex; exit normals must be valid.
    if (inside != kInside || outside != kOutside || surface != kSurface
        || (checkQueries
            && (!validOutNormal || solid.Inside(entry) != kSurface || solid.Inside(exit) != kSurface
                || !std::isfinite(normal.mag2()) || normal.mag2() < .9 || normal.mag2() > 1.1
                || !std::isfinite(safetyIn) || !std::isfinite(distanceIn)
                || !std::isfinite(safetyOut) || !std::isfinite(distanceOut)
                || !std::isfinite(distanceOutWithNormal) || safetyIn < 0. || safetyIn >= kInfinity
                || distanceIn < 0. || distanceIn >= kInfinity || safetyOut < 0.
                || safetyOut >= kInfinity || distanceOut < 0. || distanceOut >= kInfinity
                || std::abs(distanceOutWithNormal - distanceOut) > 1.e-9 * mm
                || !std::isfinite(outNormal.mag2()) || outNormal.mag2() < .9
                || outNormal.mag2() > 1.1 || safetyIn > distanceIn + 1.e-9 * mm
                || safetyOut > distanceOut + 1.e-9 * mm)))
    {
      std::cerr << "invalid generated data at index " << i << " for " << shape << '/' << topology
                << ": Inside=" << inside << ", Outside=" << outside << ", Surface=" << surface
                << ", Normal2=" << normal.mag2() << ", SafetyToIn=" << safetyIn
                << ", DistanceToIn=" << distanceIn << ", SafetyToOut=" << safetyOut
                << ", DistanceToOut=" << distanceOut
                << ", DistanceToOutWithNormal=" << distanceOutWithNormal
                << ", OutNormal2=" << outNormal.mag2() << '\n';
      std::exit(3);
    }
  }
}

template<class Function>
G4double Time(const std::uint64_t calls, Function function)
{
  G4double sum = 0.;
  const auto begin = std::chrono::steady_clock::now();
  for (std::uint64_t i = 0; i < calls; ++i)
    sum += function(i % kDataSize);
  const auto end = std::chrono::steady_clock::now();
  gSink = gSink + sum;
  return std::chrono::duration<G4double, std::milli>(end - begin).count();
}

void PrintTime(const char* label, G4double ms)
{
  G4cout << "  " << std::left << std::setw(34) << label << std::right << std::setw(10) << ms
         << " ms" << G4endl;
}

void RunQueries(const char* implementation, const char* shape, const std::string& topology,
                const G4VSolid& solid, const Data& data, std::uint64_t calls)
{
  G4cout << "\n"
         << shape << "  implementation=" << implementation << "  topology=" << topology << G4endl;

  Validate(solid, data, shape, topology);
  const auto emit = [&](const char* query, auto function) {
    for (std::size_t i = 0; i < kDataSize; ++i)
      gSink = gSink + function(i);
    PrintTime(query, Time(calls, function));
  };

  // A short untimed pass faults code/data pages in before measurement.
  for (std::size_t i = 0; i < kDataSize; ++i)
    gSink = gSink + static_cast<G4double>(solid.Inside(data.inside[i]));
  emit("Inside(p)", [&](auto i) { return G4double(solid.Inside(data.inside[i])); });
  emit("SurfaceNormal(p)", [&](auto i) { return solid.SurfaceNormal(data.surface[i]).x(); });
  emit("DistanceToIn(p)", [&](auto i) { return solid.DistanceToIn(data.outside[i]); });
  emit("DistanceToIn(p,v)",
       [&](auto i) { return solid.DistanceToIn(data.outside[i], data.directionIn[i]); });
  emit("DistanceToOut(p)", [&](auto i) { return solid.DistanceToOut(data.inside[i]); });
  emit("DistanceToOut(p,v)",
       [&](auto i) { return solid.DistanceToOut(data.inside[i], data.directionOut[i]); });
  emit("DistanceToOut(p,v,n)", [&](auto i) {
    G4bool valid = false;
    G4ThreeVector normal;
    return solid.DistanceToOut(data.inside[i], data.directionOut[i], true, &valid, &normal)
           + normal.x() + valid;
  });
}

void RunTet(const std::string& name, const Vertices& vertices, std::uint64_t calls)
{
  const G4Tet solid(name, vertices[0], vertices[1], vertices[2], vertices[3]);
#if defined(G4GEOM_USE_UTET)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Tet", name, solid, MakeTetData(vertices), calls);
}

void RunTube(const char* name, G4double a, G4double b, std::uint64_t calls)
{
  const G4double dz = 60. * mm;
  const G4EllipticalTube solid(name, a, b, dz);
#if defined(G4GEOM_USE_UELLIPTICALTUBE)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4EllipticalTube", name, solid,
             MakeEllipseData(a, b, dz), calls);
}

void RunCone(const std::string& name, G4double a, G4double b, G4double cut,
             std::uint64_t calls)
{
  const G4double height = 100. * mm;
  const G4EllipticalCone solid(name, a, b, height, cut);
#if defined(G4GEOM_USE_UELLIPTICALCONE)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4EllipticalCone", name, solid,
             MakeEllipseData(a, b, std::min(cut, height), height), calls);
}
}  // namespace

int main(int argc, char** argv)
{
  std::uint64_t calls = 1000000;
  if (argc > 2)
  {
    std::cerr << "usage: " << argv[0] << " [positive-number-of-calls]\n";
    return 2;
  }
  if (argc == 2)
  {
    const std::string_view argument(argv[1]);
    const auto result = std::from_chars(argument.begin(), argument.end(), calls);
    if (result.ec != std::errc{} || result.ptr != argument.end() || calls == 0)
    {
      std::cerr << "calls must be a positive integer\n";
      return 2;
    }
  }

  G4cout << "*********************************************************************\n"
         << "* Benchmark for G4Tet, G4EllipticalTube and G4EllipticalCone          *\n"
         << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *\n"
         << "*********************************************************************\n"
         << "\n       Number of calls per method : " << calls << G4endl
         << "       Number of sampled points   : " << kDataSize << G4endl;
  G4cout << std::fixed << std::setprecision(3);

  // Tetrahedra have one nondegenerate boundary topology. Exercise regular,
  // right, general skew/translated and slender forms, in both orientations.
  // Degenerate (zero-volume) tetrahedra are not valid navigation solids.
  const std::array<Vertices, 4> tetrahedra = {{
    {{{50., 50., 50.}, {50., -50., -50.}, {-50., 50., -50.}, {-50., -50., 50.}}},
    {{{0., 0., 0.}, {100., 0., 0.}, {0., 80., 0.}, {0., 0., 120.}}},
    {{{20., -30., 10.}, {130., -10., 25.}, {45., 90., -15.}, {-10., 15., 140.}}},
    {{{0., 0., 0.}, {100., 0., 0.}, {20., 80., 0.}, {40., 30., 2.}}}
  }};
  const char* names[] = {"regular", "right", "general-translated", "slender"};
  for (std::size_t i = 0; i < tetrahedra.size(); ++i)
  {
    auto vertices = tetrahedra[i];
    for (auto& vertex : vertices) vertex *= mm;
    RunTet(std::string(names[i]) + "/original-order", vertices, calls);
    std::swap(vertices[1], vertices[2]);
    RunTet(std::string(names[i]) + "/reversed-order", vertices, calls);
  }

  struct Axes { const char* name; G4double a, b; };
  const Axes axes[] = {{"circular", .6, .6}, {"x-major", .8, .4}, {"y-major", .4, .8}};
  for (const auto& axis : axes)
  {
    RunTube(axis.name, axis.a * 100. * mm, axis.b * 100. * mm, calls);
    RunCone(std::string(axis.name) + "/truncated", axis.a, axis.b, 60. * mm, calls);
    RunCone(std::string(axis.name) + "/apex", axis.a, axis.b, 100. * mm, calls);
    // Cuts above zMax clamp to zMax; cover this supported constructor case.
    RunCone(std::string(axis.name) + "/clamped-apex", axis.a, axis.b, 140. * mm, calls);
  }
  return gSink == -1. ? 1 : 0;
}
