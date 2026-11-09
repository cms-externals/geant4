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
// Benchmark closed tessellations and extruded polygons; see benchTessellatedExtruded.md.
#include "G4ExtrudedSolid.hh"
#include "G4PhysicalConstants.hh"
#include "G4QuadrangularFacet.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4VSolid.hh"

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
constexpr std::size_t kDataSize = 1024;
volatile G4double gSink = 0.;
std::uint64_t gRepeats = 3;
std::string gFilter;

using Face = std::vector<G4ThreeVector>;
using Mesh = std::vector<Face>;
using Section = G4ExtrudedSolid::ZSection;
inline std::vector<G4TwoVector> Polygon(int n, bool concave)
{
  std::vector<G4TwoVector> p;
  for (int i = 0; i < n; ++i)
  {
    const auto a = twopi * i / n;
    const auto r = (concave && i % 2 ? 55. : 100.) * mm;
    p.emplace_back(r * std::cos(a), r * std::sin(a));
  }
  return p;  // Counterclockwise; the extrusion constructor accepts either winding.
}

inline Mesh Loft(const std::vector<G4TwoVector>& p, const std::vector<Section>& s)
{
  Mesh mesh;
  const auto vertex = [&](std::size_t i, std::size_t j) {
    const auto q = s[j].fScale * p[i] + s[j].fOffset;
    return G4ThreeVector(q.x(), q.y(), s[j].fZ);
  };
  const auto centre = [&](std::size_t j) {
    return G4ThreeVector(s[j].fOffset.x(), s[j].fOffset.y(), s[j].fZ);
  };
  // All generated contours are star-shaped about their origin.
  for (std::size_t i = 0; i < p.size(); ++i)
  {
    const auto next = (i + 1) % p.size();
    mesh.push_back({centre(0), vertex(next, 0), vertex(i, 0)});
    mesh.push_back({centre(s.size() - 1), vertex(i, s.size() - 1), vertex(next, s.size() - 1)});
    for (std::size_t j = 0; j + 1 < s.size(); ++j)
      mesh.push_back({vertex(i, j), vertex(next, j), vertex(next, j + 1), vertex(i, j + 1)});
  }
  return mesh;
}

inline Mesh Tessellation(const std::string& topology, int n)
{
  const std::vector<Section> sections = {{-100. * mm, {}, 1.}, {100. * mm, {}, 1.}};
  if (topology == "handle")
  {
    Mesh mesh;
    constexpr int sides = 8;
    const auto point = [&](int i, int j) {
      const auto u = twopi * (i % n) / n, v = twopi * (j % sides) / sides;
      return G4ThreeVector((100. + 30. * std::cos(v)) * std::cos(u) * mm,
                           (100. + 30. * std::cos(v)) * std::sin(u) * mm, 30. * std::sin(v) * mm);
    };
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < sides; ++j)
        mesh.push_back({point(i, j), point(i + 1, j), point(i + 1, j + 1), point(i, j + 1)});
    return mesh;
  }
  auto mesh = Loft(Polygon(n, topology == "concave"), sections);
  if (topology == "cavity" || topology == "disconnected")
  {
    auto extra = mesh;
    for (auto& face : extra)
    {
      for (auto& p : face)
        p = topology == "cavity" ? .4 * p : p + G4ThreeVector(300. * mm, 0., 0.);
      if (topology == "cavity") std::reverse(face.begin(), face.end());
    }
    mesh.insert(mesh.end(), extra.begin(), extra.end());
  }
  return mesh;
}

inline Mesh FacetLayout(const Mesh& input, int layout)
{
  if (layout == 2)
  {
    Mesh mixed;
    for (std::size_t i = 0; i < input.size(); ++i)
    {
      const auto& f = input[i];
      if (f.size() == 4 && i % 2 == 0)
      {
        mixed.push_back({f[0], f[1], f[2]});
        mixed.push_back({f[0], f[2], f[3]});
      }
      else
        mixed.push_back(f);
    }
    return mixed;
  }
  Mesh triangles;
  for (const auto& f : input)
    for (std::size_t j = 1; j + 1 < f.size(); ++j)
      triangles.push_back({f[0], f[j], f[j + 1]});
  if (layout == 0) return triangles;
  Mesh quads;
  // Conforming subdivision: each original triangle becomes three convex quads.
  for (const auto& t : triangles)
  {
    const auto c = (t[0] + t[1] + t[2]) / 3.;
    for (int j = 0; j < 3; ++j)
      quads.push_back({t[j], .5 * (t[j] + t[(j + 1) % 3]), c, .5 * (t[j] + t[(j + 2) % 3])});
  }
  return quads;
}

struct Data
{
    std::vector<G4ThreeVector> inside, outside, surface, directionIn, directionOut, directionNormal;
};

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  const auto z = 2. * uniform(engine) - 1.;
  const auto phi = twopi * uniform(engine);
  const auto r = std::sqrt(1. - z * z);
  return {r * std::cos(phi), r * std::sin(phi), z};
}

Data MakeData(const Mesh& mesh)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  Data data;
  for (auto* points : {&data.inside, &data.outside, &data.surface, &data.directionIn,
                       &data.directionOut, &data.directionNormal})
    points->reserve(kDataSize);
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto& f = mesh[engine() % mesh.size()];
    const auto j = 1 + (i / mesh.size()) % (f.size() - 2);
    const auto a = .2 + .2 * uniform(engine), b = .2 + .2 * uniform(engine);
    const auto p = (1. - a - b) * f[0] + a * f[j] + b * f[j + 1];
    const auto normal = (f[j] - f[0]).cross(f[j + 1] - f[0]).unit();
    // Local probes stay clear of edges even in the highest-resolution mesh.
    const auto altitude =
      (f[j] - f[0]).cross(f[j + 1] - f[0]).mag()
      / std::max({(f[j] - f[0]).mag(), (f[j + 1] - f[0]).mag(), (f[j + 1] - f[j]).mag()});
    const auto depth = .0001 * altitude * (.5 + uniform(engine));
    const auto pin = p - depth * normal, pout = p + depth * normal;
    data.surface.push_back(p);
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.directionNormal.push_back(-normal);
    data.directionIn.push_back((-normal + .15 * UnitVector(engine)).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

// Independent, unaccelerated ray/triangle oracle, used only outside timing.
G4double ReferenceDistance(const Mesh& mesh, const G4ThreeVector& p, const G4ThreeVector& v)
{
  G4double best = kInfinity;
  for (const auto& f : mesh)
    for (std::size_t j = 1; j + 1 < f.size(); ++j)
    {
      const auto e1 = f[j] - f[0], e2 = f[j + 1] - f[0];
      const auto h = v.cross(e2);
      const auto det = e1.dot(h);
      if (std::abs(det) < 1.e-12 * mm2) continue;
      const auto s = p - f[0];
      const auto u = s.dot(h) / det;
      const auto q = s.cross(e1);
      const auto w = v.dot(q) / det;
      const auto t = e2.dot(q) / det;
      if (u >= -1.e-10 && w >= -1.e-10 && u + w <= 1. + 1.e-10 && t > 1.e-9 * mm)
        best = std::min(best, t);
    }
  return best;
}

void Validate(const G4VSolid& solid, const Data& data, const Mesh& mesh, const std::string& shape,
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
    // A concave surface may not provide a supporting exit normal. Obtain the
    // local surface normal separately when the exit output is not usable.
    if (checkQueries && !validOutNormal)
      outNormal =
        solid.SurfaceNormal(data.inside[i] + distanceOutWithNormal * data.directionOut[i]);
    if (inside != kInside || outside != kOutside || surface != kSurface
        || (checkQueries
            && (!std::isfinite(normal.mag2()) || normal.mag2() < .9 || normal.mag2() > 1.1
                || !std::isfinite(safetyIn + safetyOut + distanceIn + distanceOut) || safetyIn < 0.
                || safetyIn >= kInfinity || distanceIn < 0. || distanceIn >= kInfinity
                || safetyOut < 0. || safetyOut >= kInfinity || distanceOut < 0.
                || distanceOut >= kInfinity
                || std::abs(distanceOutWithNormal - distanceOut) > 1.e-9 * mm
                || outNormal.mag2() < .9 || outNormal.mag2() > 1.1
                || std::abs(distanceIn
                            - ReferenceDistance(mesh, data.outside[i], data.directionIn[i]))
                     > 1.e-7 * mm
                || std::abs(distanceOut
                            - ReferenceDistance(mesh, data.inside[i], data.directionOut[i]))
                     > 1.e-7 * mm
                || solid.Inside(data.outside[i] + distanceIn * data.directionIn[i]) != kSurface
                || solid.Inside(data.inside[i] + distanceOut * data.directionOut[i]) != kSurface
                || safetyIn > distanceIn + 1.e-9 * mm || safetyOut > distanceOut + 1.e-9 * mm)))
    {
      std::cerr << "invalid generated data at index " << i << " for " << shape << '/' << topology
                << ": point=" << data.inside[i] << ", Inside=" << inside << ", Outside=" << outside
                << ", Surface=" << surface << ", Normal2=" << normal.mag2()
                << ", SafetyToIn=" << safetyIn << ", DistanceToIn=" << distanceIn
                << ", SafetyToOut=" << safetyOut << ", DistanceToOut=" << distanceOut
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
                const G4VSolid& solid, const Mesh& mesh, std::uint64_t calls)
{
  G4cout << "\n"
         << shape << "  implementation=" << implementation << "  topology=" << topology << G4endl;

  const auto data = MakeData(mesh);
  Validate(solid, data, mesh, shape, topology);
  std::size_t normalIncidenceFailures = 0;
  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto expected = (data.outside[i] - data.surface[i]).mag();
    const auto actual = solid.DistanceToIn(data.outside[i], data.directionNormal[i]);
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1.e-7 * mm)
    {
      if (normalIncidenceFailures == 0)
        G4cout << std::setprecision(17)
               << "  First normal-incidence failure: point=" << data.outside[i]
               << " direction=" << data.directionNormal[i] << " expected_mm=" << expected / mm
               << " actual_mm=" << actual / mm << std::setprecision(6) << G4endl;
      ++normalIncidenceFailures;
    }
  }
  G4cout << "  Normal-incidence entry failures: " << normalIncidenceFailures << " / " << kDataSize
         << G4endl;
  std::vector<G4double> times;
  times.reserve(gRepeats);
  const auto emit = [&](const char* query, auto function) {
    Time(std::min<std::uint64_t>(256, calls), function);
    times.clear();
    for (std::uint64_t repeat = 0; repeat < gRepeats; ++repeat)
      times.push_back(Time(calls, function));
    std::sort(times.begin(), times.end());
    const auto middle = times.size() / 2;
    const auto median = times.size() % 2 ? times[middle] : .5 * (times[middle - 1] + times[middle]);
    PrintTime(query, median);
    G4cout << "    min=" << times.front() << " ms  max=" << times.back() << " ms" << G4endl;
  };

#if defined(G4BENCH_ACCURATE_SAFETY)
  // Native-only build mode for comparing equal-quality surface safeties.
  const auto& faceted = dynamic_cast<const G4TessellatedSolid&>(solid);
  emit("SafetyFromOutside(accurate)",
       [&](auto i) { return faceted.SafetyFromOutside(data.outside[i], true); });
  emit("SafetyFromInside(accurate)",
       [&](auto i) { return faceted.SafetyFromInside(data.inside[i], true); });
  return;
#endif

  emit("Inside(p)", [&](auto i) { return G4double(solid.Inside(data.inside[i])); });
  emit("Inside(outside)", [&](auto i) { return G4double(solid.Inside(data.outside[i])); });
  emit("Inside(surface)", [&](auto i) { return G4double(solid.Inside(data.surface[i])); });
  emit("SurfaceNormal(p)", [&](auto i) { return solid.SurfaceNormal(data.surface[i]).x(); });
  emit("DistanceToIn(p)", [&](auto i) { return solid.DistanceToIn(data.outside[i]); });
  emit("DistanceToIn(p,v)",
       [&](auto i) { return solid.DistanceToIn(data.outside[i], data.directionIn[i]); });
  if (normalIncidenceFailures == 0)
    emit("DistanceToIn(normal-incidence)",
         [&](auto i) { return solid.DistanceToIn(data.outside[i], data.directionNormal[i]); });
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

void RunTessellated(const std::string& topology, int n, int layout, std::uint64_t calls)
{
  const char* layouts[] = {"triangles", "quadrangles", "mixed"};
  const auto name = topology + "-" + layouts[layout];
  if (!gFilter.empty() && ("G4TessellatedSolid/" + name).find(gFilter) == std::string::npos) return;
  const auto mesh = FacetLayout(Tessellation(topology, n), layout);
  G4TessellatedSolid solid(name);
  for (const auto& f : mesh)
  {
    G4VFacet* facet =
      f.size() == 3
        ? static_cast<G4VFacet*>(new G4TriangularFacet(f[0], f[1], f[2], ABSOLUTE))
        : static_cast<G4VFacet*>(new G4QuadrangularFacet(f[0], f[1], f[2], f[3], ABSOLUTE));
    if (!solid.AddFacet(facet))
    {
      delete facet;
      std::cerr << "AddFacet failed\n";
      std::exit(3);
    }
#if defined(G4GEOM_USE_UTESSELLATEDSOLID)
    // The adapter copies vertices; unlike the native solid it does not own this input facet.
    delete facet;
#endif
  }
  solid.SetSolidClosed(true);
#if defined(G4GEOM_USE_UTESSELLATEDSOLID)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4TessellatedSolid",
             name + "  resolution=" + std::to_string(n) + "  facets=" + std::to_string(mesh.size()),
             solid, mesh, calls);
}

void RunExtruded(bool concave, int n, int kind, int nz, std::uint64_t calls)
{
  const char* kinds[] = {"right", "tapered", "offset", "scaled-offset", "multisection"};
  const auto name = std::string(concave ? "concave-" : "convex-") + kinds[kind];
  if (!gFilter.empty() && ("G4ExtrudedSolid/" + name).find(gFilter) == std::string::npos) return;
  auto polygon = Polygon(n, concave);
  std::vector<Section> sections;
  for (int j = 0; j < nz; ++j)
  {
    const auto t = G4double(j) / (nz - 1);
    G4TwoVector offset;
    G4double scale = 1.;
    if (kind == 1 || kind == 3) scale = 1. - .4 * t;
    if (kind == 2 || kind == 3) offset = {30. * t * mm, -20. * t * mm};
    if (kind == 4)
    {
      scale = 1. + .2 * std::sin(3. * pi * t);
      offset = {20. * std::sin(2. * pi * t) * mm, 15. * std::sin(pi * t) * mm};
    }
    sections.emplace_back((200. * t - 100.) * mm, offset, scale);
  }
  G4ExtrudedSolid solid(name, polygon, sections);
#if defined(G4GEOM_USE_UEXTRUDEDSOLID)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  // Geant4's defining facets: quadrangular sides, triangular caps except n=4.
  const auto facets = n * (nz - 1) + (n == 4 ? 2 : 2 * (n - 2));
  RunQueries(implementation, "G4ExtrudedSolid",
             name + "  vertices=" + std::to_string(n) + "  sections=" + std::to_string(nz)
               + "  facets=" + std::to_string(facets),
             solid, Loft(polygon, sections), calls);
}
}  // namespace

int main(int argc, char** argv)
{
  // Sixteen complete cycles through the samples give a practical default
  // for all 116 cases. Longer precision runs can set the call count explicitly.
  std::uint64_t calls = 16 * kDataSize;
  if (argc > 4)
  {
    std::cerr << "usage: " << argv[0] << " [positive-calls [positive-repeats [case-substring]]]\n";
    return 2;
  }
  for (int i = 1; i < std::min(argc, 3); ++i)
  {
    const std::string_view arg(argv[i]);
    auto& value = i == 1 ? calls : gRepeats;
    const auto result = std::from_chars(arg.begin(), arg.end(), value);
    if (result.ec != std::errc{} || result.ptr != arg.end() || value == 0)
    {
      std::cerr << "calls and repeats must be positive integers\n";
      return 2;
    }
  }
  if (argc == 4) gFilter = argv[3];
  G4cout << "*********************************************************************\n"
         << "* Benchmark for G4TessellatedSolid and G4ExtrudedSolid               *\n"
         << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut    *\n"
         << "*********************************************************************\n"
         << "\n       Number of calls per method : " << calls
         << "\n       Number of sampled points   : " << kDataSize
         << "\n       Repetitions (median time)  : " << gRepeats << G4endl;
  G4cout << std::fixed << std::setprecision(6);
  for (const auto* topology : {"convex", "concave", "cavity", "handle", "disconnected"})
    for (int layout = 0; layout < 3; ++layout)
      for (int n : {8, 32, 128, 512})
        RunTessellated(topology, n, layout, calls);
  for (bool concave : {false, true})
  {
    for (int kind = 0; kind < 5; ++kind)
      for (int n : {3, 4, 8, 32, 128, 512})
      {
        if (concave && n < 8) continue;
        RunExtruded(concave, n, kind, kind == 4 ? 4 : 2, calls);
      }
    // Independent Z-section scaling with fixed polygon complexity.
    for (int nz : {8, 32, 128})
      RunExtruded(concave, 32, 4, nz, calls);
  }
  return gSink == -1. ? 1 : 0;
}
