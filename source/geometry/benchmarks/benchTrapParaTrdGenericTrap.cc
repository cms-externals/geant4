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
// Benchmark G4Trap/G4Para/G4Trd/G4GenericTrap navigation queries.

#include "G4GenericTrap.hh"
#include "G4Para.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4TwoVector.hh"
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
#include <utility>
#include <vector>

namespace
{
constexpr std::size_t kDataSize = 8192;
volatile G4double gSink = 0.;

using Face = std::array<G4TwoVector, 4>;

struct Loft
{
    G4double dz;
    Face bottom;
    Face top;
};

struct Data
{
    std::vector<G4ThreeVector> inside, outside, surface, directionIn, directionOut;
};

struct TrapTopology
{
    const char* name;
    G4double dz, theta, phi;
    G4double dy1, dx1, dx2, alpha1;
    G4double dy2, dx3, dx4, alpha2;
};

struct ParaTopology
{
    const char* name;
    G4double dx, dy, dz, alpha, theta, phi;
};

struct TrdTopology
{
    const char* name;
    G4double dx1, dx2, dy1, dy2, dz;
};

struct GenericTrapTopology
{
    const char* name;
    std::vector<G4TwoVector> vertices;
};

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  const auto z = 2. * uniform(engine) - 1.;
  const auto phi = twopi * uniform(engine);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

G4TwoVector Interpolate(const G4TwoVector& a, const G4TwoVector& b, G4double fraction)
{
  return (1. - fraction) * a + fraction * b;
}

G4ThreeVector PointAt(const Loft& loft, std::size_t edge, G4double alongEdge, G4double alongZ)
{
  const auto bottom = Interpolate(loft.bottom[edge], loft.bottom[(edge + 1) % 4], alongEdge);
  const auto top = Interpolate(loft.top[edge], loft.top[(edge + 1) % 4], alongEdge);
  const auto point = Interpolate(bottom, top, alongZ);
  return {point.x(), point.y(), (2. * alongZ - 1.) * loft.dz};
}

G4ThreeVector CentreAt(const Loft& loft, G4double alongZ)
{
  G4TwoVector centre;
  for (std::size_t i = 0; i < 4; ++i)
    centre += Interpolate(loft.bottom[i], loft.top[i], alongZ);
  centre *= .25;
  return {centre.x(), centre.y(), (2. * alongZ - 1.) * loft.dz};
}

G4double FaceArea(const Face& face)
{
  G4double twiceArea = 0.;
  for (std::size_t i = 0; i < 4; ++i)
    twiceArea += face[i].x() * face[(i + 1) % 4].y() - face[i].y() * face[(i + 1) % 4].x();
  return .5 * std::abs(twiceArea);
}

Data MakeData(const Loft& loft)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  std::vector<G4int> surfaces;
  std::vector<G4double> endFaces;
  std::vector<std::size_t> edges;
  if (FaceArea(loft.bottom) > 1.e-12 * mm2)
  {
    surfaces.push_back(-1);
    endFaces.push_back(0.);
  }
  if (FaceArea(loft.top) > 1.e-12 * mm2)
  {
    surfaces.push_back(-2);
    endFaces.push_back(1.);
  }
  for (std::size_t edge = 0; edge < 4; ++edge)
  {
    const auto next = (edge + 1) % 4;
    const auto bottomLength = (loft.bottom[next] - loft.bottom[edge]).mag2();
    const auto topLength = (loft.top[next] - loft.top[edge]).mag2();
    if (std::max(bottomLength, topLength) > 1.e-12 * mm2)
    {
      surfaces.push_back(static_cast<G4int>(edge));
      edges.push_back(edge);
    }
  }

  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto alongZ = .12 + .76 * uniform(engine);
    const auto edge = edges[(i / 2) % edges.size()];
    const auto alongEdge = .12 + .76 * uniform(engine);
    const auto centre = CentreAt(loft, alongZ);
    const auto boundary = PointAt(loft, edge, alongEdge, alongZ);
    G4ThreeVector pin;
    G4ThreeVector pout;
    if ((i % 2) != 0)
    {
      G4ThreeVector entry;
      G4ThreeVector innerCentre;
      if (endFaces.empty())
      {
        entry = boundary;
        const auto innerZ = alongZ < .5 ? alongZ + .2 : alongZ - .2;
        innerCentre = CentreAt(loft, innerZ);
      }
      else
      {
        const auto end = endFaces[(i / 2) % endFaces.size()];
        const auto endCentre = CentreAt(loft, end);
        const auto endBoundary = PointAt(loft, edge, alongEdge, end);
        entry = endCentre + .55 * (endBoundary - endCentre);
        innerCentre = CentreAt(loft, .5);
      }
      pin = entry + .35 * (innerCentre - entry);
      pout = entry + .75 * (entry - pin);
    }
    else
    {
      const auto radialFraction = .15 + .55 * uniform(engine);
      pin = centre + radialFraction * (boundary - centre);
      pout = centre + 2.25 * (boundary - centre);
    }

    const auto surface = surfaces[i % surfaces.size()];
    G4ThreeVector psurface;
    if (surface < 0)
    {
      const auto zFraction = surface == -1 ? 0. : 1.;
      const auto faceCentre = CentreAt(loft, zFraction);
      const auto faceBoundary = PointAt(loft, edge, alongEdge, zFraction);
      psurface = faceCentre + .65 * (faceBoundary - faceCentre);
    }
    else
    {
      psurface = PointAt(loft, static_cast<std::size_t>(surface), alongEdge, alongZ);
    }

    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(psurface);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

Loft MakeTrapLoft(const TrapTopology& topology)
{
  const auto tanTheta = std::tan(topology.theta);
  const auto thetaX = tanTheta * std::cos(topology.phi);
  const auto thetaY = tanTheta * std::sin(topology.phi);
  const auto alpha1 = std::tan(topology.alpha1);
  const auto alpha2 = std::tan(topology.alpha2);
  const auto z = topology.dz;

  Loft loft;
  loft.dz = z;
  loft.bottom = {
    {{-z * thetaX - topology.dy1 * alpha1 - topology.dx1, -z * thetaY - topology.dy1},
     {-z * thetaX - topology.dy1 * alpha1 + topology.dx1, -z * thetaY - topology.dy1},
     {-z * thetaX + topology.dy1 * alpha1 + topology.dx2, -z * thetaY + topology.dy1},
     {-z * thetaX + topology.dy1 * alpha1 - topology.dx2, -z * thetaY + topology.dy1}}};
  loft.top = {{{z * thetaX - topology.dy2 * alpha2 - topology.dx3, z * thetaY - topology.dy2},
               {z * thetaX - topology.dy2 * alpha2 + topology.dx3, z * thetaY - topology.dy2},
               {z * thetaX + topology.dy2 * alpha2 + topology.dx4, z * thetaY + topology.dy2},
               {z * thetaX + topology.dy2 * alpha2 - topology.dx4, z * thetaY + topology.dy2}}};
  return loft;
}

Loft MakeParaLoft(const ParaTopology& topology)
{
  const auto tanAlpha = std::tan(topology.alpha);
  const auto tanTheta = std::tan(topology.theta);
  const auto thetaX = tanTheta * std::cos(topology.phi);
  const auto thetaY = tanTheta * std::sin(topology.phi);
  Loft loft;
  loft.dz = topology.dz;
  for (std::size_t end = 0; end < 2; ++end)
  {
    const auto zSign = end == 0 ? -1. : 1.;
    Face& face = end == 0 ? loft.bottom : loft.top;
    const std::array<std::pair<G4double, G4double>, 4> signs = {
      std::pair<G4double, G4double>{-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.}};
    for (std::size_t i = 0; i < 4; ++i)
    {
      const auto x = signs[i].first * topology.dx + signs[i].second * topology.dy * tanAlpha
                     + zSign * topology.dz * thetaX;
      const auto y = signs[i].second * topology.dy + zSign * topology.dz * thetaY;
      face[i] = {x, y};
    }
  }
  return loft;
}

Loft MakeTrdLoft(const TrdTopology& topology)
{
  Loft loft;
  loft.dz = topology.dz;
  loft.bottom = {{{-topology.dx1, -topology.dy1},
                  {topology.dx1, -topology.dy1},
                  {topology.dx1, topology.dy1},
                  {-topology.dx1, topology.dy1}}};
  loft.top = {{{-topology.dx2, -topology.dy2},
               {topology.dx2, -topology.dy2},
               {topology.dx2, topology.dy2},
               {-topology.dx2, topology.dy2}}};
  return loft;
}

Loft MakeGenericTrapLoft(G4double dz, const std::vector<G4TwoVector>& vertices)
{
  Loft loft;
  loft.dz = dz;
  for (std::size_t i = 0; i < 4; ++i)
  {
    loft.bottom[i] = vertices[i];
    loft.top[i] = vertices[i + 4];
  }
  return loft;
}

std::vector<G4TwoVector> MakeVertices(const std::array<G4double, 8>& x,
                                      const std::array<G4double, 8>& y)
{
  std::vector<G4TwoVector> vertices;
  vertices.reserve(8);
  for (std::size_t i = 0; i < 8; ++i)
    vertices.emplace_back(x[i] * cm, y[i] * cm);
  return vertices;
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
    // A twisted face may not provide a supporting exit normal. Obtain the
    // local surface normal separately when the exit output is not usable.
    if (checkQueries && !validOutNormal)
      outNormal =
        solid.SurfaceNormal(data.inside[i] + distanceOutWithNormal * data.directionOut[i]);
    if (inside != kInside || outside != kOutside || surface != kSurface
        || (checkQueries
            && (normal.mag2() < .9 || normal.mag2() > 1.1 || safetyIn < 0. || safetyIn >= kInfinity
                || distanceIn < 0. || distanceIn >= kInfinity || safetyOut < 0.
                || safetyOut >= kInfinity || distanceOut < 0. || distanceOut >= kInfinity
                || std::abs(distanceOutWithNormal - distanceOut) > 1.e-9 * mm
                || outNormal.mag2() < .9 || outNormal.mag2() > 1.1
                || safetyIn > distanceIn + 1.e-9 * mm || safetyOut > distanceOut + 1.e-9 * mm)))
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
                const G4VSolid& solid, const Loft& loft, std::uint64_t calls)
{
  G4cout << "\n"
         << shape << "  implementation=" << implementation << "  topology=" << topology << G4endl;

  const auto data = MakeData(loft);
  Validate(solid, data, shape, topology);
  const auto emit = [&](const char* query, auto function) {
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

void RunTrap(const TrapTopology& topology, std::uint64_t calls)
{
  const G4Trap solid(topology.name, topology.dz, topology.theta, topology.phi, topology.dy1,
                     topology.dx1, topology.dx2, topology.alpha1, topology.dy2, topology.dx3,
                     topology.dx4, topology.alpha2);
#if defined(G4GEOM_USE_UTRAP)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Trap", topology.name, solid, MakeTrapLoft(topology), calls);
}

void RunPara(const ParaTopology& topology, std::uint64_t calls)
{
  const G4Para solid(topology.name, topology.dx, topology.dy, topology.dz, topology.alpha,
                     topology.theta, topology.phi);
#if defined(G4GEOM_USE_UPARA)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Para", topology.name, solid, MakeParaLoft(topology), calls);
}

void RunTrd(const TrdTopology& topology, std::uint64_t calls)
{
  const G4Trd solid(topology.name, topology.dx1, topology.dx2, topology.dy1, topology.dy2,
                    topology.dz);
#if defined(G4GEOM_USE_UTRD)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Trd", topology.name, solid, MakeTrdLoft(topology), calls);
}

void RunGenericTrap(const GenericTrapTopology& topology, std::uint64_t calls)
{
  constexpr G4double dz = 10. * cm;
  const G4GenericTrap solid(topology.name, dz, topology.vertices);
#if defined(G4GEOM_USE_UGENERICTRAP)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4GenericTrap", topology.name, solid,
             MakeGenericTrapLoft(dz, topology.vertices), calls);
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

  const TrapTopology trapTopologies[] = {
    {"box-like", 15. * cm, 0., 0., 10. * cm, 7.5 * cm, 7.5 * cm, 0., 10. * cm, 7.5 * cm, 7.5 * cm,
     0.},
    {"yz-rectangle", 15. * cm, 15. * deg, 0., 10. * cm, 7.5 * cm, 7.5 * cm, 0., 10. * cm, 7.5 * cm,
     7.5 * cm, 0.},
    {"trd1-like", 15. * cm, 0., 0., 10. * cm, 7.5 * cm, 7.5 * cm, 0., 10. * cm, 3.5 * cm, 3.5 * cm,
     0.},
    {"trd2-like", 15. * cm, 0., 0., 10. * cm, 7.5 * cm, 7.5 * cm, 0., 8. * cm, 3.5 * cm, 3.5 * cm,
     0.},
    {"para-like", 15. * cm, 20. * deg, 5. * deg, 10. * cm, 7.5 * cm, 7.5 * cm, 10. * deg, 10. * cm,
     7.5 * cm, 7.5 * cm, 10. * deg},
    {"xy-trapezoid", 15. * cm, 0., 0., 10. * cm, 7.5 * cm, 4. * cm, 0., 10. * cm, 7.5 * cm, 4. * cm,
     0.},
    {"general", 15. * cm, 20. * deg, 5. * deg, 10. * cm, 7.5 * cm, 10. * cm, 10. * deg, 4. * cm,
     2.5 * cm, 3.5 * cm, 10. * deg}};

  const ParaTopology paraTopologies[] = {
    {"box-like", 10. * cm, 7. * cm, 15. * cm, 0., 0., 0.},
    {"general-skew", 10. * cm, 7. * cm, 15. * cm, 20. * deg, 15. * deg, 35. * deg}};

  const TrdTopology trdTopologies[] = {
    {"box", 12. * cm, 12. * cm, 9. * cm, 9. * cm, 15. * cm},
    {"trd1", 8. * cm, 13. * cm, 9. * cm, 9. * cm, 15. * cm},
    {"trd2-increasing", 7. * cm, 13. * cm, 6. * cm, 10. * cm, 15. * cm},
    {"trd2-decreasing", 13. * cm, 7. * cm, 10. * cm, 6. * cm, 15. * cm},
    {"six-vertex-x-wedge", 12. * cm, 0., 9. * cm, 9. * cm, 15. * cm},
    {"six-vertex-y-wedge", 12. * cm, 12. * cm, 9. * cm, 0., 15. * cm},
    {"five-vertex-pyramid", 12. * cm, 0., 9. * cm, 0., 15. * cm},
    {"four-vertex-tetrahedron", 12. * cm, 0., 0., 9. * cm, 15. * cm}};

  const GenericTrapTopology genericTrapTopologies[] = {
    {"twisted",
     MakeVertices({-3., -2.5, 3., 2.5, -2., -2., 2., 2.}, {-2.5, 3., 2.5, -3., -2., 2., 2., -2.})},
    {"one-twisted-face",
     MakeVertices({-3., -3., 3., 3., -2., -2., 2., 2.}, {-3., 3., 3., -3., -2., 2., 2., -1.})},
    {"two-twisted-faces",
     MakeVertices({-3., -3., 3., 3., -2.5, -1.5, 2., 3.}, {-3., 3., 3., -3., -2., 2., 2., -2.})},
    {"planar",
     MakeVertices({-3., -3., 3., 3., -2., -2., 2., 2.}, {-3., 3., 3., -3., -2., 2., 2., -2.})},
    {"one-face-triangle",
     MakeVertices({-3., -3., 3., 2.5, -2., -2., 2., 2.}, {-2.5, -2.5, 2.5, -3., -2., 2., 2., -2.})},
    {"one-face-line", MakeVertices({-3., -3., 2.5, 2.5, -2., -2., 2., 2.},
                                   {-2.5, -2.5, -3., -3., -2., 2., 2., -2.})},
    {"one-face-point", MakeVertices({-3., -3., -3., -3., -2., -2., 2., 2.},
                                    {-2.5, -2.5, -2.5, -2.5, -2., 2., 2., -2.})},
    {"line-to-triangle", MakeVertices({-3., -3., 2.5, 2.5, -2., -2., 2., 2.},
                                      {-2.5, -2.5, -3., -3., -2., -2., 2., -2.})},
    {"point-to-triangle", MakeVertices({-3., -3., -3., -3., -2., -2., 2., 2.},
                                       {-2.5, -2.5, -2.5, -2.5, -2., -2., 2., -2.})},
    {"line-to-line",
     MakeVertices({-3., 3., 3., -3., 0., 0., 0., 0.}, {0., 0., 0., 0., 3., 3., -3., -3.})},
    {"degenerated-planar",
     MakeVertices({-.507492, -.507508, 1.522492, -.507492, -.507492, -.507508, 1.522492, -.507492},
                  {-3.634, 3.634, 3.634, -3.634, -3.634, 3.634, 3.634, -3.634})}};

  G4cout << "*********************************************************************" << G4endl;
  G4cout << "* Benchmark for G4Trap, G4Para, G4Trd and G4GenericTrap            *" << G4endl;
  G4cout << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *" << G4endl;
  G4cout << "* Configurations: regular, skewed, tapered and degenerate topologies *" << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "\n       Number of calls per method : " << calls << G4endl;
  G4cout << "       Number of sampled points   : " << kDataSize << G4endl;

  for (const auto& topology : trapTopologies)
    RunTrap(topology, calls);
  for (const auto& topology : paraTopologies)
    RunPara(topology, calls);
  for (const auto& topology : trdTopologies)
    RunTrd(topology, calls);
  for (const auto& topology : genericTrapTopologies)
    RunGenericTrap(topology, calls);
  return gSink == -1. ? 1 : 0;
}
