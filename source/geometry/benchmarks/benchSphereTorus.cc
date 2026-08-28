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
// Benchmark G4Sphere/G4Torus navigation queries for all topologies.

#include "G4PhysicalConstants.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Torus.hh"
#include "G4VPVParameterisation.hh"
#include "G4VSolid.hh"
#include "globals.hh"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kDataSize = 8192;
constexpr G4double kStartPhi = 0.;
volatile G4double gSink = 0.;

G4ThreeVector TorusPoint(G4double minorRadius, G4double minorPhi, G4double majorRadius,
                         G4double majorPhi);
G4ThreeVector SpherePoint(G4double radius, G4double theta, G4double phi);

class AngularParameterisation final : public G4VPVParameterisation
{
  public:

    void ComputeTransformation(G4int, G4VPhysicalVolume*) const override {}

    void ComputeDimensions(G4Sphere& solid, G4int copyNo, const G4VPhysicalVolume*) const override
    {
      const G4double startPhi = copyNo == 0 ? 30. * deg : 150. * deg;
      solid.SetStartPhiAngle(startPhi, false);
      solid.SetDeltaPhiAngle(60. * deg);
      solid.SetStartThetaAngle(60. * deg);
      solid.SetDeltaThetaAngle(60. * deg);
    }

    void ComputeDimensions(G4Torus& solid, G4int copyNo, const G4VPhysicalVolume*) const override
    {
      const G4double startPhi = copyNo == 0 ? 30. * deg : 150. * deg;
      solid.SetAllParameters(solid.GetRmin(), solid.GetRmax(), solid.GetRtor(), startPhi,
                             60. * deg);
    }
};

void ValidateParameterisationUpdates()
{
  AngularParameterisation parameterisation;
  G4Sphere sphere("parameterised-sphere", 0., 12. * cm, 0., 90. * deg, 0., pi);
  G4Torus torus("parameterised-torus", 0., 6. * cm, 24. * cm, 0., 90. * deg);

  sphere.ComputeDimensions(&parameterisation, 1, nullptr);
  torus.ComputeDimensions(&parameterisation, 1, nullptr);

  const G4ThreeVector sphereInside = G4ThreeVector(std::cos(pi), std::sin(pi), 0.) * (8. * cm);
  const G4ThreeVector sphereOldPhi =
    G4ThreeVector(std::cos(45. * deg), std::sin(45. * deg), 0.) * (8. * cm);
  const G4ThreeVector torusInside = G4ThreeVector(std::cos(pi), std::sin(pi), 0.) * (27. * cm);
  const G4ThreeVector torusOldPhi =
    G4ThreeVector(std::cos(45. * deg), std::sin(45. * deg), 0.) * (27. * cm);

  if (sphere.Inside(sphereInside) != kInside || sphere.Inside(sphereOldPhi) != kOutside)
    throw std::runtime_error("parameterised G4Sphere angular caches were not updated");
  if (torus.Inside(torusInside) != kInside || torus.Inside(torusOldPhi) != kOutside)
    throw std::runtime_error("parameterised G4Torus angular caches were not updated");
}

void ValidateTorusKnownDistances()
{
  constexpr std::array<G4double, 3> scales = {1. * mm, 1. * m, 1000. * m};
  constexpr std::array<G4double, 4> deltaPhi = {twopi, 90. * deg, pi, 270. * deg};
  for (const auto scale : scales)
  {
    const G4double rtor = scale;
    const G4double rmax = 0.25 * scale;
    const G4double gap = 0.07 * scale;
    for (const G4bool hollow : {false, true})
    {
      const G4double rmin = hollow ? 0.08 * scale : 0.;
      const G4double middle = 0.5 * (rmin + rmax);
      for (const auto dphi : deltaPhi)
      {
        G4Torus torus("known-distance-torus", rmin, rmax, rtor, kStartPhi, dphi);
        for (G4int i = 0; i < 128; ++i)
        {
          const G4double majorPhi = kStartPhi + (0.1 + 0.8 * (i + 0.5) / 128.) * dphi;
          const G4double minorPhi = twopi * (i + 0.5) / 128.;
          const G4ThreeVector normal(std::cos(minorPhi) * std::cos(majorPhi),
                                     std::cos(minorPhi) * std::sin(majorPhi), std::sin(minorPhi));
          const G4ThreeVector outer = TorusPoint(rmax, minorPhi, rtor, majorPhi);
          const G4ThreeVector inside = TorusPoint(middle, minorPhi, rtor, majorPhi);
          const G4double distanceIn = torus.DistanceToIn(outer + gap * normal, -normal);
          const G4double distanceOut = torus.DistanceToOut(inside, normal);
          const G4double tolerance = 2.e-9 * scale;
          if (std::fabs(distanceIn - gap) > tolerance
              || std::fabs(distanceOut - (rmax - middle)) > tolerance)
            throw std::runtime_error("G4Torus known exact distance validation failed");
          if (hollow)
          {
            const G4double distanceToInner = torus.DistanceToOut(inside, -normal);
            if (std::fabs(distanceToInner - (middle - rmin)) > tolerance)
              throw std::runtime_error("G4Torus inner exact distance validation failed");
          }
        }
      }
    }
  }
}

void ValidateSphereKnownDistances()
{
  constexpr std::array<G4double, 3> scales = {1. * mm, 1. * m, 1000. * m};
  constexpr std::array<G4double, 4> deltaPhi = {twopi, 90. * deg, pi, 270. * deg};
  constexpr std::array<std::array<G4double, 2>, 4> theta = {std::array<G4double, 2>{0., pi},
                                                            {0., 60. * deg},
                                                            {45. * deg, 90. * deg},
                                                            {120. * deg, 60. * deg}};
  for (const auto scale : scales)
  {
    const G4double rmax = scale;
    const G4double gap = 0.07 * scale;
    for (const G4bool hollow : {false, true})
    {
      const G4double rmin = hollow ? 0.25 * scale : 0.;
      const G4double middle = 0.5 * (rmin + rmax);
      for (const auto dphi : deltaPhi)
      {
        for (const auto& thetaRange : theta)
        {
          G4Sphere sphere("known-distance-sphere", rmin, rmax, kStartPhi, dphi, thetaRange[0],
                          thetaRange[1]);
          for (G4int i = 0; i < 32; ++i)
          {
            const G4double phi = kStartPhi + (0.02 + 0.96 * (i + 0.5) / 32.) * dphi;
            const G4double polar = thetaRange[0] + (0.02 + 0.96 * (i + 0.5) / 32.) * thetaRange[1];
            const G4ThreeVector direction = SpherePoint(1., polar, phi);
            const G4double distanceIn = sphere.DistanceToIn((rmax + gap) * direction, -direction);
            const G4double distanceOut = sphere.DistanceToOut(middle * direction, direction);
            const G4double tolerance = 2.e-9 * scale;
            if (std::fabs(distanceIn - gap) > tolerance
                || std::fabs(distanceOut - (rmax - middle)) > tolerance)
            {
              std::cerr << "G4Sphere known distance mismatch: scale=" << scale
                        << " hollow=" << hollow << " dphi=" << dphi << " stheta=" << thetaRange[0]
                        << " dtheta=" << thetaRange[1] << " distanceIn=" << distanceIn
                        << " expectedIn=" << gap << " distanceOut=" << distanceOut
                        << " expectedOut=" << rmax - middle << '\n';
              throw std::runtime_error("G4Sphere known exact distance validation failed");
            }
          }
        }
      }
    }
  }
}

struct RadialTopology
{
    const char* name;
    G4bool hollow;
};

struct PhiTopology
{
    const char* name;
    G4double dphi;
};

struct ThetaTopology
{
    const char* name;
    G4double stheta;
    G4double dtheta;
};

struct Data
{
    std::vector<G4ThreeVector> inside, outside, surface, directionIn, directionOut;
};

enum class SphereSurface
{
  kOuter,
  kInner,
  kStartPhi,
  kEndPhi,
  kStartTheta,
  kEndTheta
};

enum class TorusSurface
{
  kOuter,
  kInner,
  kStartPhi,
  kEndPhi
};

G4ThreeVector UnitVector(std::mt19937_64& engine)
{
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  const auto z = 2. * uniform(engine) - 1.;
  const auto phi = twopi * uniform(engine);
  const auto rho = std::sqrt((1. - z) * (1. + z));
  return {rho * std::cos(phi), rho * std::sin(phi), z};
}

G4ThreeVector SpherePoint(G4double radius, G4double theta, G4double phi)
{
  const auto rho = radius * std::sin(theta);
  return {rho * std::cos(phi), rho * std::sin(phi), radius * std::cos(theta)};
}

G4ThreeVector TorusPoint(G4double minorRadius, G4double minorPhi, G4double majorRadius,
                         G4double majorPhi)
{
  const auto rho = majorRadius + minorRadius * std::cos(minorPhi);
  return {rho * std::cos(majorPhi), rho * std::sin(majorPhi), minorRadius * std::sin(minorPhi)};
}

Data MakeSphereData(const RadialTopology& radialTopology, const PhiTopology& phiTopology,
                    const ThetaTopology& thetaTopology)
{
  constexpr G4double rmax = 12. * cm;
  const G4double rmin = radialTopology.hollow ? 4. * cm : 0.;
  const G4double etheta = thetaTopology.stheta + thetaTopology.dtheta;
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  std::vector<SphereSurface> surfaces = {SphereSurface::kOuter};
  if (radialTopology.hollow) surfaces.push_back(SphereSurface::kInner);
  if (phiTopology.dphi < twopi)
  {
    surfaces.push_back(SphereSurface::kStartPhi);
    surfaces.push_back(SphereSurface::kEndPhi);
  }
  if (thetaTopology.stheta > 0.) surfaces.push_back(SphereSurface::kStartTheta);
  if (etheta < pi) surfaces.push_back(SphereSurface::kEndTheta);

  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto radius = rmin + (.15 + .7 * uniform(engine)) * (rmax - rmin);
    const auto phi = kStartPhi + (.08 + .84 * uniform(engine)) * phiTopology.dphi;
    const auto theta = thetaTopology.stheta + (.08 + .84 * uniform(engine)) * thetaTopology.dtheta;
    const auto pin = SpherePoint(radius, theta, phi);
    const auto pout = SpherePoint(2.25 * rmax, theta, phi);
    G4ThreeVector psurface;

    switch (surfaces[i % surfaces.size()])
    {
      case SphereSurface::kOuter:
        psurface = SpherePoint(rmax, theta, phi);
        break;
      case SphereSurface::kInner:
        psurface = SpherePoint(rmin, theta, phi);
        break;
      case SphereSurface::kStartPhi:
        psurface = SpherePoint(radius, theta, kStartPhi);
        break;
      case SphereSurface::kEndPhi:
        psurface = SpherePoint(radius, theta, kStartPhi + phiTopology.dphi);
        break;
      case SphereSurface::kStartTheta:
        psurface = SpherePoint(radius, thetaTopology.stheta, phi);
        break;
      case SphereSurface::kEndTheta:
        psurface = SpherePoint(radius, etheta, phi);
        break;
    }
    data.inside.push_back(pin);
    data.outside.push_back(pout);
    data.surface.push_back(psurface);
    data.directionIn.push_back((pin - pout).unit());
    data.directionOut.push_back(UnitVector(engine));
  }
  return data;
}

Data MakeTorusData(const RadialTopology& radialTopology, const PhiTopology& phiTopology)
{
  constexpr G4double rmax = 6. * cm;
  constexpr G4double rtor = 24. * cm;
  const G4double rmin = radialTopology.hollow ? 2. * cm : 0.;
  std::mt19937_64 engine(0x5eed1234ULL);
  std::uniform_real_distribution<G4double> uniform(0., 1.);
  std::vector<TorusSurface> surfaces = {TorusSurface::kOuter};
  if (radialTopology.hollow) surfaces.push_back(TorusSurface::kInner);
  if (phiTopology.dphi < twopi)
  {
    surfaces.push_back(TorusSurface::kStartPhi);
    surfaces.push_back(TorusSurface::kEndPhi);
  }

  Data data;
  data.inside.reserve(kDataSize);
  data.outside.reserve(kDataSize);
  data.surface.reserve(kDataSize);
  data.directionIn.reserve(kDataSize);
  data.directionOut.reserve(kDataSize);

  for (std::size_t i = 0; i < kDataSize; ++i)
  {
    const auto minorRadius = rmin + (.15 + .7 * uniform(engine)) * (rmax - rmin);
    const auto minorPhi = twopi * uniform(engine);
    const auto majorPhi = kStartPhi + (.08 + .84 * uniform(engine)) * phiTopology.dphi;
    const auto pin = TorusPoint(minorRadius, minorPhi, rtor, majorPhi);
    const auto pout = TorusPoint(2.25 * rmax, minorPhi, rtor, majorPhi);
    G4ThreeVector psurface;

    switch (surfaces[i % surfaces.size()])
    {
      case TorusSurface::kOuter:
        psurface = TorusPoint(rmax, minorPhi, rtor, majorPhi);
        break;
      case TorusSurface::kInner:
        psurface = TorusPoint(rmin, minorPhi, rtor, majorPhi);
        break;
      case TorusSurface::kStartPhi:
        psurface = TorusPoint(minorRadius, minorPhi, rtor, kStartPhi);
        break;
      case TorusSurface::kEndPhi:
        psurface = TorusPoint(minorRadius, minorPhi, rtor, kStartPhi + phiTopology.dphi);
        break;
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
    const auto checkDistances = (i % 128) == 0;
    const auto distanceIn =
      checkDistances ? solid.DistanceToIn(data.outside[i], data.directionIn[i]) : 0.;
    const auto distanceOut =
      checkDistances ? solid.DistanceToOut(data.inside[i], data.directionOut[i]) : 0.;
    if (inside != kInside || outside != kOutside || surface != kSurface || distanceIn < 0.
        || distanceIn >= kInfinity || distanceOut < 0. || distanceOut >= kInfinity)
    {
      std::cerr << "invalid generated data at index " << i << " for " << shape << '/' << topology
                << ": Inside=" << inside << ", Outside=" << outside << ", Surface=" << surface
                << ", DistanceToIn=" << distanceIn << ", DistanceToOut=" << distanceOut << '\n';
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
  gSink += sum;
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
         << shape << "  implementation=" << implementation << "  topology=" << topology
         << G4endl;

  Validate(solid, data, shape, topology);
  const auto emit = [&](const char* query, auto function) {
    const auto ms = Time(calls, function);
    PrintTime(query, ms);
  };

  // A short untimed pass faults code/data pages in before measurement.
  for (std::size_t i = 0; i < kDataSize; ++i)
    gSink += static_cast<G4double>(solid.Inside(data.inside[i]));
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

void RunSphere(const RadialTopology& radialTopology, const PhiTopology& phiTopology,
               const ThetaTopology& thetaTopology, std::uint64_t calls)
{
  constexpr G4double rmax = 12. * cm;
  const auto rmin = radialTopology.hollow ? 4. * cm : 0.;
  const std::string topology =
    std::string(radialTopology.name) + '-' + phiTopology.name + '-' + thetaTopology.name;
  const G4Sphere solid(topology, rmin, rmax, kStartPhi, phiTopology.dphi, thetaTopology.stheta,
                       thetaTopology.dtheta);
  const auto data = MakeSphereData(radialTopology, phiTopology, thetaTopology);
#if defined(G4GEOM_USE_USPHERE)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Sphere", topology, solid, data, calls);
}

void RunTorus(const RadialTopology& radialTopology, const PhiTopology& phiTopology,
              std::uint64_t calls)
{
  constexpr G4double rmax = 6. * cm;
  constexpr G4double rtor = 24. * cm;
  const auto rmin = radialTopology.hollow ? 2. * cm : 0.;
  const std::string topology = std::string(radialTopology.name) + '-' + phiTopology.name;
  const G4Torus solid(topology, rmin, rmax, rtor, kStartPhi, phiTopology.dphi);
  const auto data = MakeTorusData(radialTopology, phiTopology);
#if defined(G4GEOM_USE_UTORUS)
  constexpr const char* implementation = "VecGeom";
#else
  constexpr const char* implementation = "native";
#endif
  RunQueries(implementation, "G4Torus", topology, solid, data, calls);
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

  ValidateParameterisationUpdates();
#if !defined(G4GEOM_USE_USPHERE)
  ValidateSphereKnownDistances();
#endif
  ValidateTorusKnownDistances();

  const RadialTopology radialTopologies[] = {{"solid", false}, {"hollow", true}};
  const PhiTopology phiTopologies[] = {{"full-phi", twopi},
                                       {"less-pi-phi", 90. * deg},
                                       {"pi-phi", pi},
                                       {"greater-pi-phi", 270. * deg}};
  const ThetaTopology thetaTopologies[] = {{"full-theta", 0., pi},
                                           {"upper-cap", 0., 60. * deg},
                                           {"band", 45. * deg, 90. * deg},
                                           {"lower-cap", 120. * deg, 60. * deg}};

  G4cout << "*********************************************************************" << G4endl;
  G4cout << "* Benchmark for G4Sphere and G4Torus                                *" << G4endl;
  G4cout << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *" << G4endl;
  G4cout << "* Configurations: solid/hollow, full/partial phi and theta          *" << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "\n       Number of calls per method : " << calls << G4endl;
  G4cout << "       Number of sampled points   : " << kDataSize << G4endl;

  for (const auto& radialTopology : radialTopologies)
  {
    for (const auto& phiTopology : phiTopologies)
    {
      for (const auto& thetaTopology : thetaTopologies)
        RunSphere(radialTopology, phiTopology, thetaTopology, calls);
      RunTorus(radialTopology, phiTopology, calls);
    }
  }
  return gSink == -1. ? 1 : 0;
}
