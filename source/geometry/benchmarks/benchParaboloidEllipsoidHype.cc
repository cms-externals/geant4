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
// Benchmark G4Paraboloid/G4Ellipsoid/G4Hype navigation queries.
// Analytic, backend-independent samples; balanced surface classes, not area weighted.

#include "G4Ellipsoid.hh"
#include "G4Hype.hh"
#include "G4Paraboloid.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "globals.hh"

#include <algorithm>
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

// Cross sections have elliptical outer and optional circular inner boundaries.
struct Section
{
    G4double x, y, inner;
};
G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  const auto z = 2. * uniform(engine) - 1.;
  const auto phi = twopi * uniform(engine);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

template<class CrossSection>
Data MakeData(G4double zmin, G4double zmax, G4bool hollow, CrossSection section)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  // 0: outer wall, 1: inner wall, 2/3: lower/upper caps.
  std::vector<int> surfaces{0};
  if (hollow) surfaces.push_back(1);
  for (int cap = 2; cap <= 3; ++cap)
  {
    const auto s = section(cap == 2 ? zmin : zmax);
    if (std::min(s.x, s.y) > s.inner + 1.e-9 * mm) surfaces.push_back(cap);
  }
  Data data;
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto z = zmin + (.15 + .7 * uniform(engine)) * (zmax - zmin);
    const auto phi = twopi * uniform(engine);
    const auto c = std::cos(phi), s = std::sin(phi);
    const auto cross = section(z);
    const auto radial = [&](const Section& a, G4double fraction, G4double height) {
      return G4ThreeVector((a.inner + fraction * (a.x - a.inner)) * c,
                           (a.inner + fraction * (a.y - a.inner)) * s, height);
    };
    const auto pin = radial(cross, .2 + .6 * uniform(engine), z);
    G4ThreeVector pout;
    // Outside samples cover outer wall, hole (where present), and both ends.
    const auto region = i % (hollow ? 4 : 3);
    if (region == 0)
      pout = radial(cross, 1.4, z);
    else if (hollow && region == 3)
      pout = radial(cross, 0., z) * .5;
    else
      pout = G4ThreeVector(pin.x(), pin.y(),
                           region == 1 ? zmin - .3 * (zmax - zmin) : zmax + .3 * (zmax - zmin));
    // Preserve z when sampling the conical/curved hole.
    if (hollow && region == 3) pout.setZ(z);
    const auto face = surfaces[i % surfaces.size()];
    G4ThreeVector psurface;
    if (face < 2)
      psurface = radial(cross, face == 0 ? 1. : 0., z);
    else
    {
      const auto end = face == 2 ? zmin : zmax;
      psurface = radial(section(end), .2 + .6 * uniform(engine), end);
    }
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(psurface);
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
    // Concave surfaces may legitimately return validOutNormal=false.
    (void)validOutNormal;
    if (inside != kInside || outside != kOutside || surface != kSurface
        || (checkQueries
            && (solid.Inside(entry) != kSurface || solid.Inside(exit) != kSurface
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

void RunParaboloid(const char* name, G4double r1, std::uint64_t calls)
{
  const G4double dz = 60. * mm, r2 = 70. * mm;
  const G4Paraboloid solid(name, dz, r1, r2);
#if defined(G4GEOM_USE_UPARABOLOID)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  const auto data = MakeData(-dz, dz, false, [&](G4double z) {
    const auto r = std::sqrt(r1 * r1 + (r2 * r2 - r1 * r1) * (z + dz) / (2. * dz));
    return Section{r, r, 0.};
  });
  RunQueries(implementation, "G4Paraboloid", name, solid, data, calls);
}

void RunEllipsoid(const char* name, G4double a, G4double b, G4double low, G4double high,
                  std::uint64_t calls)
{
  const G4double c = 100. * mm;
  const G4Ellipsoid solid(name, a, b, c, low * c, high * c);
#if defined(G4GEOM_USE_UELLIPSOID)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  const auto data = MakeData(low * c, high * c, false, [&](G4double z) {
    const auto scale = std::sqrt(std::max(0., (1. - z / c) * (1. + z / c)));
    return Section{a * scale, b * scale, 0.};
  });
  RunQueries(implementation, "G4Ellipsoid", name, solid, data, calls);
}

void RunHype(const char* name, G4double ri, G4double ti, G4double to, std::uint64_t calls)
{
  const G4double dz = 60. * mm, ro = 70. * mm;
  const G4Hype solid(name, ri, ro, std::atan(ti), std::atan(to), dz);
#if defined(G4GEOM_USE_UHYPE)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  const auto data = MakeData(-dz, dz, ri > 0. || ti != 0., [&](G4double z) {
    const auto outer = std::sqrt(ro * ro + z * z * to * to);
    return Section{outer, outer, std::sqrt(ri * ri + z * z * ti * ti)};
  });
  RunQueries(implementation, "G4Hype", name, solid, data, calls);
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
         << "* Benchmark for G4Paraboloid, G4Ellipsoid and G4Hype                  *\n"
         << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *\n"
         << "*********************************************************************\n"
         << "\n       Number of calls per method : " << calls << G4endl
         << "       Number of sampled points   : " << kDataSize << G4endl;
  G4cout << std::fixed << std::setprecision(3);

  // Native G4Paraboloid requires 0 <= r1 < r2: no reversed/cylindrical case.
  RunParaboloid("apex", 0., calls);
  RunParaboloid("truncated", 30. * mm, calls);

  // All cap-presence and cut-location cases, including one-sided segments.
  struct Cuts
  {
      const char* name;
      G4double low, high;
  };
  const Cuts cuts[] = {{"uncut", -1., 1.},
                       {"lower-cut", -.5, 1.},
                       {"upper-cut", -1., .6},
                       {"two-cuts-straddling", -.5, .6},
                       {"lower-half", -1., 0.},
                       {"upper-half", 0., 1.},
                       {"lower-cap", -1., -.4},
                       {"upper-cap", .4, 1.},
                       {"two-cuts-from-equator", 0., .6},
                       {"two-cuts-to-equator", -.6, 0.},
                       {"two-cuts-negative", -.8, -.3},
                       {"two-cuts-positive", .3, .8}};
  // Cross cut topologies with triaxial, oblate, prolate and spherical cases.
  struct Axes
  {
      const char* name;
      G4double a, b;
  };
  const Axes axes[] = {{"triaxial", 50. * mm, 70. * mm},
                       {"oblate", 130. * mm, 130. * mm},
                       {"prolate", 60. * mm, 60. * mm},
                       {"sphere", 100. * mm, 100. * mm}};
  for (const auto& axis : axes)
    for (const auto& cut : cuts)
    {
      const auto name = std::string(axis.name) + "/" + cut.name;
      RunEllipsoid(name.c_str(), axis.a, axis.b, cut.low, cut.high, calls);
    }

  // Outer wall: cylindrical/hyperbolic. Inner: absent/cylinder/hyperboloid/cone.
  RunHype("solid-cylinder", 0., 0., 0., calls);
  RunHype("solid-hyperboloid", 0., 0., .5, calls);
  RunHype("cylindrical-shell", 20. * mm, 0., 0., calls);
  RunHype("cylindrical-hole", 20. * mm, 0., .5, calls);
  RunHype("hyperbolic-hole-cylinder", 20. * mm, .3, 0., calls);
  RunHype("equal-stereo-shell", 20. * mm, .5, .5, calls);
  RunHype("outer-stereo-larger", 20. * mm, .2, .5, calls);
  RunHype("inner-stereo-larger", 20. * mm, .5, .2, calls);
  RunHype("conical-hole-cylinder", 0., .3, 0., calls);
  RunHype("conical-hole-hyperboloid", 0., .3, .5, calls);
  // End annuli collapse to rings; avoid singular rims in the normal samples.
  RunHype("closed-ends-hyperbolic-hole", 20. * mm, 1.125, .125, calls);
  RunHype("closed-ends-hyperbolic-hole-cylinder", 20. * mm, std::sqrt(1.25), 0., calls);
  RunHype("closed-ends-conical-hole-hyperboloid", 0., 37. / 24., 25. / 24., calls);
  RunHype("closed-ends-conical-hole", 0., 7. / 6., 0., calls);
  return gSink == -1. ? 1 : 0;
}
