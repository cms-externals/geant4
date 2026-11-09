//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided under the terms and    *
// * conditions of the Geant4 Software License, included in the file  *
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
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the Geant4 collaboration.                      *
// * By using, copying, modifying or distributing the software (or    *
// * any work based on the software) you agree to acknowledge its use *
// * in resulting scientific publications, and indicate your          *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Benchmark twisted-solid navigation queries for all supported topologies.

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"
#include "G4VSolid.hh"
#include "globals.hh"

#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
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

struct FacetedGeometry
{
    G4double twist, dz, theta, phi;
    G4double dy1, dx1, dx2;
    G4double dy2, dx3, dx4;
    G4double alpha;
};

struct BoxTopology
{
    const char* name;
    G4double twist, dx, dy, dz;
};

struct TrdTopology
{
    const char* name;
    G4double dx1, dx2, dy1, dy2, dz, twist;
};

struct TrapTopology
{
    const char* name;
    G4bool simple;
    FacetedGeometry geometry;
};

enum class TubsConstructor
{
  kEndRadiiDPhi,
  kEndRadiiSegments,
  kMidRadiiDPhi,
  kMidRadiiSegments
};

struct TubsTopology
{
    const char* name;
    TubsConstructor constructor;
    G4double twist, rinner, router;
    G4double zlower, zupper;
    G4int nseg;
    G4double phi;
};

G4double Uniform(std::mt19937_64& engine, G4double low, G4double high)
{
  std::uniform_real_distribution<G4double> uniform(low, high);
  return uniform(engine);
}

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  const auto z = Uniform(engine, -1., 1.);
  const auto phi = Uniform(engine, 0., twopi);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

G4ThreeVector OutsidePoint(const G4VSolid& solid, std::mt19937_64& engine)
{
  G4ThreeVector minimum;
  G4ThreeVector maximum;
  solid.BoundingLimits(minimum, maximum);
  const auto centre = .5 * (minimum + maximum);
  const auto halfDiagonal = .5 * (maximum - minimum).mag();
  return centre + 1.5 * halfDiagonal * UnitVector(engine);
}

G4ThreeVector FacetedPoint(const FacetedGeometry& geometry, G4double fractionZ, G4double fractionX,
                           G4double fractionY)
{
  const auto interpolate = [fractionZ](G4double lower, G4double upper) {
    return (1. - fractionZ) * lower + fractionZ * upper;
  };
  const auto dy = interpolate(geometry.dy1, geometry.dy2);
  const auto dxLower = interpolate(geometry.dx1, geometry.dx3);
  const auto dxUpper = interpolate(geometry.dx2, geometry.dx4);
  const auto y = fractionY * dy;
  const auto dx = .5 * ((1. - fractionY) * dxLower + (1. + fractionY) * dxUpper);
  const auto x = y * std::tan(geometry.alpha) + fractionX * dx;
  const auto z = (2. * fractionZ - 1.) * geometry.dz;
  const auto rotation = (fractionZ - .5) * geometry.twist;
  const auto shift = z * std::tan(geometry.theta);
  const auto cosRotation = std::cos(rotation);
  const auto sinRotation = std::sin(rotation);
  return {x * cosRotation - y * sinRotation + shift * std::cos(geometry.phi),
          x * sinRotation + y * cosRotation + shift * std::sin(geometry.phi), z};
}

Data MakeFacetedData(const G4VSolid& solid, const FacetedGeometry& geometry)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto fractionZ = Uniform(engine, .12, .88);
    const auto fractionX = Uniform(engine, -.68, .68);
    const auto fractionY = Uniform(engine, -.68, .68);
    const auto pin = FacetedPoint(geometry, fractionZ, fractionX, fractionY);

    const auto surfaceZ = Uniform(engine, .1, .9);
    const auto surfaceX = Uniform(engine, -.72, .72);
    const auto surfaceY = Uniform(engine, -.72, .72);
    G4ThreeVector psurface;
    switch (i % 6)
    {
      case 0:
        psurface = FacetedPoint(geometry, 0., surfaceX, surfaceY);
        break;
      case 1:
        psurface = FacetedPoint(geometry, 1., surfaceX, surfaceY);
        break;
      case 2:
        psurface = FacetedPoint(geometry, surfaceZ, 1., surfaceY);
        break;
      case 3:
        psurface = FacetedPoint(geometry, surfaceZ, surfaceX, 1.);
        break;
      case 4:
        psurface = FacetedPoint(geometry, surfaceZ, -1., surfaceY);
        break;
      default:
        psurface = FacetedPoint(geometry, surfaceZ, surfaceX, -1.);
        break;
    }

    const auto pout = OutsidePoint(solid, engine);
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(psurface);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

G4ThreeVector TubsPoint(const G4TwistedTubs& solid, G4double fractionZ, G4double fractionRadius,
                        G4double fractionPhi)
{
  const auto z = (1. - fractionZ) * solid.GetEndZ(0) + fractionZ * solid.GetEndZ(1);
  const auto rinner =
    std::sqrt(solid.GetInnerRadius() * solid.GetInnerRadius() + z * z * solid.GetTanInnerStereo2());
  const auto router =
    std::sqrt(solid.GetOuterRadius() * solid.GetOuterRadius() + z * z * solid.GetTanOuterStereo2());
  const auto radius = (1. - fractionRadius) * rinner + fractionRadius * router;
  const auto centrePhi = std::atan2(solid.GetKappa() * z, 1.);
  const auto phi = centrePhi + .5 * fractionPhi * solid.GetDPhi();
  return {radius * std::cos(phi), radius * std::sin(phi), z};
}

Data MakeTubsData(const G4TwistedTubs& solid)
{
  std::mt19937_64 engine(0x5eed1234ULL);
  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto fractionZ = Uniform(engine, .12, .88);
    const auto fractionRadius = Uniform(engine, .15, .85);
    const auto fractionPhi = Uniform(engine, -.72, .72);
    const auto pin = TubsPoint(solid, fractionZ, fractionRadius, fractionPhi);

    const auto surfaceZ = Uniform(engine, .1, .9);
    const auto surfaceRadius = Uniform(engine, .12, .88);
    const auto surfacePhi = Uniform(engine, -.76, .76);
    G4ThreeVector psurface;
    switch (i % 6)
    {
      case 0:
        psurface = TubsPoint(solid, 0., surfaceRadius, surfacePhi);
        break;
      case 1:
        psurface = TubsPoint(solid, 1., surfaceRadius, surfacePhi);
        break;
      case 2:
        psurface = TubsPoint(solid, surfaceZ, 0., surfacePhi);
        break;
      case 3:
        psurface = TubsPoint(solid, surfaceZ, 1., surfacePhi);
        break;
      case 4:
        psurface = TubsPoint(solid, surfaceZ, surfaceRadius, -1.);
        break;
      default:
        psurface = TubsPoint(solid, surfaceZ, surfaceRadius, 1.);
        break;
    }

    const auto pout = OutsidePoint(solid, engine);
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
    const auto normal2 = normal.mag2();
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
    const auto outNormal2 = outNormal.mag2();
    EInside beforeEntry = kSurface;
    EInside afterEntry = kSurface;
    EInside beforeExit = kSurface;
    EInside afterExit = kSurface;
    if (checkQueries && distanceIn < kInfinity && distanceOut < kInfinity)
    {
      constexpr G4double probe = 1.e-4 * mm;
      const auto entry = data.outside[i] + distanceIn * data.directionIn[i];
      const auto exit = data.inside[i] + distanceOut * data.directionOut[i];
      beforeEntry = solid.Inside(entry - probe * data.directionIn[i]);
      afterEntry = solid.Inside(entry + probe * data.directionIn[i]);
      beforeExit = solid.Inside(exit - probe * data.directionOut[i]);
      afterExit = solid.Inside(exit + probe * data.directionOut[i]);
    }
    // Normals on twisted faces are local rather than globally valid, so
    // validOutNormal may legitimately be false. The magnitude is unspecified.
    (void)validOutNormal;
    const auto validDistance = [](G4double distance) {
      return std::isfinite(distance) && distance >= 0. && distance < kInfinity;
    };
    if (inside != kInside || outside != kOutside || surface != kSurface
        || (checkQueries
            && (!std::isfinite(normal2) || normal2 < .9 || normal2 > 1.1 || !validDistance(safetyIn)
                || !validDistance(distanceIn) || !validDistance(safetyOut)
                || !validDistance(distanceOut) || !std::isfinite(distanceOutWithNormal)
                || std::abs(distanceOutWithNormal - distanceOut) > 1.e-9 * mm
                || !std::isfinite(outNormal2) || outNormal2 == 0. || beforeEntry == kInside
                || afterEntry == kOutside || beforeExit == kOutside || afterExit == kInside
                || safetyIn > distanceIn + 1.e-9 * mm || safetyOut > distanceOut + 1.e-9 * mm)))
    {
      std::cerr << "invalid generated data at index " << i << " for " << shape << '/' << topology
                << ": Inside=" << inside << ", Outside=" << outside << ", Surface=" << surface
                << ", Normal2=" << normal2 << ", SafetyToIn=" << safetyIn
                << ", DistanceToIn=" << distanceIn << ", SafetyToOut=" << safetyOut
                << ", DistanceToOut=" << distanceOut
                << ", DistanceToOutWithNormal=" << distanceOutWithNormal
                << ", OutNormal2=" << outNormal2 << ", EntryTransition=" << beforeEntry << '/'
                << afterEntry << ", ExitTransition=" << beforeExit << '/' << afterExit << '\n';
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

void RunQueries(const char* shape, const std::string& topology, const G4VSolid& solid,
                const Data& data, std::uint64_t calls)
{
  G4cout << "\n" << shape << "  implementation=native  topology=" << topology << G4endl;

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

FacetedGeometry BoxGeometry(const BoxTopology& topology)
{
  return {topology.twist, topology.dz, 0.,          0.,          topology.dy, topology.dx,
          topology.dx,    topology.dy, topology.dx, topology.dx, 0.};
}

FacetedGeometry TrdGeometry(const TrdTopology& topology)
{
  return {topology.twist, topology.dz,  0.,           0., topology.dy1, topology.dx1, topology.dx1,
          topology.dy2,   topology.dx2, topology.dx2, 0.};
}

void RunBox(const BoxTopology& topology, std::uint64_t calls)
{
  const G4TwistedBox solid(topology.name, topology.twist, topology.dx, topology.dy, topology.dz);
  const auto geometry = BoxGeometry(topology);
  RunQueries("G4TwistedBox", topology.name, solid, MakeFacetedData(solid, geometry), calls);
}

void RunTrd(const TrdTopology& topology, std::uint64_t calls)
{
  const G4TwistedTrd solid(topology.name, topology.dx1, topology.dx2, topology.dy1, topology.dy2,
                           topology.dz, topology.twist);
  const auto geometry = TrdGeometry(topology);
  RunQueries("G4TwistedTrd", topology.name, solid, MakeFacetedData(solid, geometry), calls);
}

void RunTrap(const TrapTopology& topology, std::uint64_t calls)
{
  const auto& geometry = topology.geometry;
  if (topology.simple)
  {
    const G4TwistedTrap solid(topology.name, geometry.twist, geometry.dx1, geometry.dx2,
                              geometry.dy1, geometry.dz);
    RunQueries("G4TwistedTrap", topology.name, solid, MakeFacetedData(solid, geometry), calls);
  }
  else
  {
    const G4TwistedTrap solid(topology.name, geometry.twist, geometry.dz, geometry.theta,
                              geometry.phi, geometry.dy1, geometry.dx1, geometry.dx2, geometry.dy2,
                              geometry.dx3, geometry.dx4, geometry.alpha);
    RunQueries("G4TwistedTrap", topology.name, solid, MakeFacetedData(solid, geometry), calls);
  }
}

void RunTubs(const TubsTopology& topology, std::uint64_t calls)
{
  std::unique_ptr<G4TwistedTubs> solid;
  switch (topology.constructor)
  {
    case TubsConstructor::kEndRadiiDPhi:
      solid = std::make_unique<G4TwistedTubs>(topology.name, topology.twist, topology.rinner,
                                              topology.router, topology.zupper, topology.phi);
      break;
    case TubsConstructor::kEndRadiiSegments:
      solid = std::make_unique<G4TwistedTubs>(topology.name, topology.twist, topology.rinner,
                                              topology.router, topology.zupper, topology.nseg,
                                              topology.phi);
      break;
    case TubsConstructor::kMidRadiiDPhi:
      solid = std::make_unique<G4TwistedTubs>(topology.name, topology.twist, topology.rinner,
                                              topology.router, topology.zlower, topology.zupper,
                                              topology.phi);
      break;
    case TubsConstructor::kMidRadiiSegments:
      solid = std::make_unique<G4TwistedTubs>(topology.name, topology.twist, topology.rinner,
                                              topology.router, topology.zlower, topology.zupper,
                                              topology.nseg, topology.phi);
      break;
  }
  RunQueries("G4TwistedTubs", topology.name, *solid, MakeTubsData(*solid), calls);
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

  // Topologies are the discrete constructor, surface, equality and handedness
  // regimes. Continuous dimensional variants and symmetry-equivalent rotations
  // are represented once.
  const BoxTopology boxTopologies[] = {
    {"square-right-handed", 30. * deg, 10. * cm, 10. * cm, 15. * cm},
    {"rectangular-left-handed", -30. * deg, 12. * cm, 7. * cm, 15. * cm},
    {"rectangular-large-twist", 80. * deg, 12. * cm, 7. * cm, 15. * cm}};

  const TrdTopology trdTopologies[] = {
    {"constant-section", 11. * cm, 11. * cm, 8. * cm, 8. * cm, 15. * cm, 30. * deg},
    {"x-taper", 7. * cm, 13. * cm, 9. * cm, 9. * cm, 15. * cm, 30. * deg},
    {"y-taper", 11. * cm, 11. * cm, 6. * cm, 10. * cm, 15. * cm, -30. * deg},
    {"xy-expanding", 7. * cm, 13. * cm, 6. * cm, 10. * cm, 15. * cm, 30. * deg},
    {"xy-contracting", 13. * cm, 7. * cm, 10. * cm, 6. * cm, 15. * cm, -30. * deg},
    {"mixed-taper", 7. * cm, 13. * cm, 10. * cm, 6. * cm, 15. * cm, 45. * deg}};

  const TrapTopology trapTopologies[] = {
    {"regular-rectangular-caps",
     true,
     {30. * deg, 15. * cm, 0., 0., 9. * cm, 8. * cm, 8. * cm, 9. * cm, 8. * cm, 8. * cm, 0.}},
    {"regular-trapezoidal-caps",
     true,
     {35. * deg, 15. * cm, 0., 0., 9. * cm, 6. * cm, 11. * cm, 9. * cm, 6. * cm, 11. * cm, 0.}},
    {"general-rectangular-caps",
     false,
     {30. * deg, 15. * cm, 12. * deg, 35. * deg, 7. * cm, 8. * cm, 8. * cm, 10. * cm, 12. * cm,
      12. * cm, 0.}},
    {"general-trapezoidal-caps",
     false,
     {-40. * deg, 15. * cm, 10. * deg, 25. * deg, 8. * cm, 6. * cm, 10. * cm, 12. * cm, 7. * cm,
      13. * cm, -12. * deg}}};

  // Keep each effective segment width strictly below pi because navigation
  // represents the phi boundaries with two oriented half-plane tests.
  const TubsTopology tubsTopologies[] = {
    {"end-radii-dphi", TubsConstructor::kEndRadiiDPhi, 30. * deg, 4. * cm, 12. * cm, 0., 15. * cm,
     0, 75. * deg},
    {"end-radii-segments", TubsConstructor::kEndRadiiSegments, -30. * deg, 4. * cm, 12. * cm, 0.,
     15. * cm, 4, twopi},
    {"mid-radii-asymmetric-dphi", TubsConstructor::kMidRadiiDPhi, 30. * deg, 4. * cm, 12. * cm,
     -10. * cm, 18. * cm, 0, 140. * deg},
    {"mid-radii-asymmetric-segments", TubsConstructor::kMidRadiiSegments, -30. * deg, 4. * cm,
     12. * cm, -3. * cm, 18. * cm, 2, 300. * deg}};

  G4cout << "*********************************************************************" << G4endl;
  G4cout << "* Benchmark for G4TwistedBox, G4TwistedTrap, G4TwistedTrd and      *" << G4endl;
  G4cout << "* G4TwistedTubs                                                    *" << G4endl;
  G4cout << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *" << G4endl;
  G4cout << "* Configurations: handed, tapered, skewed and constructor topologies *" << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "\n       Number of calls per method : " << calls << G4endl;
  G4cout << "       Number of sampled points   : " << kDataSize << G4endl;

  for (const auto& topology : boxTopologies)
    RunBox(topology, calls);
  for (const auto& topology : trapTopologies)
    RunTrap(topology, calls);
  for (const auto& topology : trdTopologies)
    RunTrd(topology, calls);
  for (const auto& topology : tubsTopologies)
    RunTubs(topology, calls);
  return gSink == -1. ? 1 : 0;
}
