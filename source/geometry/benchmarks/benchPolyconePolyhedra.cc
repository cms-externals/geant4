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

#include "G4GenericPolycone.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "geomdefs.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

namespace
{
constexpr G4int kDefaultLoops = 200000;
constexpr G4int kDefaultNumPoints = 10000;

volatile G4double gDoubleSink = 0.;
volatile G4int gIntSink = 0;
G4int gNumPoints = kDefaultNumPoints;
G4bool gCsvMode = false;

enum class ProfileKind
{
  Sinusoidal,
  ZMonotonicConvex,
  ZStepped
};

const char* ProfileName(ProfileKind kind)
{
  if (kind == ProfileKind::ZMonotonicConvex)
  {
    return "z-monotonic-convex";
  }
  return (kind == ProfileKind::ZStepped) ? "z-stepped" : "sinusoidal";
}

G4bool HasProfile(ProfileKind kind, G4int zSections)
{
  return kind != ProfileKind::ZStepped || zSections >= 6;
}

struct Profile
{
    std::vector<G4double> z;
    std::vector<G4double> rmin;
    std::vector<G4double> rmax;
};

struct Sample
{
    std::vector<G4ThreeVector> inside;
    std::vector<G4ThreeVector> outside;
    std::vector<G4ThreeVector> surface;
    std::vector<G4ThreeVector> toInDirection;
    std::vector<G4ThreeVector> toOutDirection;
};

struct Config
{
    G4String solidName;
    ProfileKind profileKind;
    G4int zSections;
    G4bool fullPhi;
    G4bool hasRmin;
    std::unique_ptr<G4VSolid> solid;
};

Profile MakeProfile(G4int zSections, G4bool hasRmin, ProfileKind kind)
{
  Profile p;
  p.z.resize(zSections);
  p.rmin.resize(zSections);
  p.rmax.resize(zSections);

  for (G4int i = 0; i < zSections; ++i)
  {
    const G4double t = (zSections == 1) ? 0. : G4double(i) / G4double(zSections - 1);
    const G4double angle = twopi * t;

    if (kind == ProfileKind::ZStepped)
    {
      const G4int nLevel = (zSections + 1) / 2;
      const G4int level = std::min(i / 2, nLevel - 1);
      const G4double s = (nLevel <= 1) ? 0. : G4double(level) / G4double(nLevel - 1);
      const G4double u = 2. * s - 1.;
      p.z[i] = (-50. + 100. * s) * cm;
      p.rmax[i] = (16. + 12. * (1. - std::fabs(u))) * cm;
      p.rmin[i] = hasRmin ? (4. + 2. * std::fabs(u)) * cm : 0.;
      if ((i & 1) != 0)
      {
        p.rmax[i] += 3. * cm;
        if (hasRmin)
        {
          p.rmin[i] = std::max(0., p.rmin[i] - 1. * cm);
        }
      }
    }
    else if (kind == ProfileKind::ZMonotonicConvex)
    {
      const G4double u = 2. * t - 1.;
      p.z[i] = (-50. + 100. * t) * cm;
      p.rmax[i] = (28. - 6. * u * u) * cm;
      p.rmin[i] = hasRmin ? (4. + 2. * u * u) * cm : 0.;
    }
    else
    {
      p.z[i] = (-50. + 100. * t) * cm;
      p.rmax[i] = (24. + 5. * std::sin(angle) + 2. * std::sin(3. * angle)) * cm;
      p.rmin[i] = hasRmin ? (5. + 1.2 * std::cos(2. * angle)) * cm : 0.;
    }
    p.rmin[i] = std::min(p.rmin[i], p.rmax[i] - 5. * cm);
  }
  return p;
}

std::unique_ptr<G4VSolid> MakePolycone(const G4String& name, G4int zSections, G4double sphi,
                                       G4double dphi, G4bool hasRmin, ProfileKind kind)
{
  Profile p = MakeProfile(zSections, hasRmin, kind);
  return std::unique_ptr<G4VSolid>(
    new G4Polycone(name, sphi, dphi, zSections, &p.z[0], &p.rmin[0], &p.rmax[0]));
}

std::unique_ptr<G4VSolid> MakePolyhedra(const G4String& name, G4int zSections, G4double sphi,
                                        G4double dphi, G4bool hasRmin, ProfileKind kind)
{
  Profile p = MakeProfile(zSections, hasRmin, kind);
  constexpr G4int numSide = 12;
  return std::unique_ptr<G4VSolid>(
    new G4Polyhedra(name, sphi, dphi, numSide, zSections, &p.z[0], &p.rmin[0], &p.rmax[0]));
}

std::unique_ptr<G4VSolid> MakeGenericPolycone(const G4String& name, G4int zSections, G4double sphi,
                                              G4double dphi, G4bool hasRmin, ProfileKind kind)
{
  Profile p = MakeProfile(zSections, hasRmin, kind);
  std::vector<G4double> r;
  std::vector<G4double> z;

  r.reserve(hasRmin ? 2 * zSections : zSections + 2);
  z.reserve(hasRmin ? 2 * zSections : zSections + 2);

  for (G4int i = 0; i < zSections; ++i)
  {
    r.push_back(p.rmax[i]);
    z.push_back(p.z[i]);
  }

  if (hasRmin)
  {
    for (G4int i = zSections - 1; i >= 0; --i)
    {
      r.push_back(p.rmin[i]);
      z.push_back(p.z[i]);
    }
  }
  else
  {
    r.push_back(0.);
    z.push_back(p.z.back());
    r.push_back(0.);
    z.push_back(p.z.front());
  }

  return std::unique_ptr<G4VSolid>(
    new G4GenericPolycone(name, sphi, dphi, G4int(r.size()), &r[0], &z[0]));
}

G4double Uniform(G4double a, G4double b)
{
  return a + (b - a) * G4UniformRand();
}

G4ThreeVector RandomPoint(const G4ThreeVector& pmin, const G4ThreeVector& pmax, G4double margin)
{
  return G4ThreeVector(Uniform(pmin.x() - margin, pmax.x() + margin),
                       Uniform(pmin.y() - margin, pmax.y() + margin),
                       Uniform(pmin.z() - margin, pmax.z() + margin));
}

G4ThreeVector RandomDirection()
{
  const G4double z = Uniform(-1., 1.);
  const G4double phi = Uniform(0., twopi);
  const G4double r = std::sqrt(std::max(0., 1. - z * z));
  return G4ThreeVector(r * std::cos(phi), r * std::sin(phi), z);
}

void FillClassifiedPoints(G4VSolid& solid, Sample& sample)
{
  G4ThreeVector pmin, pmax;
  solid.BoundingLimits(pmin, pmax);
  const G4double margin = 0.25 * (pmax - pmin).mag();

  sample.inside.reserve(gNumPoints);
  sample.outside.reserve(gNumPoints);
  sample.surface.reserve(gNumPoints);
  sample.toInDirection.reserve(gNumPoints);
  sample.toOutDirection.reserve(gNumPoints);

  G4int trials = 0;
  while (G4int(sample.inside.size()) < gNumPoints || G4int(sample.outside.size()) < gNumPoints)
  {
    const G4ThreeVector p = RandomPoint(pmin, pmax, margin);
    const EInside in = solid.Inside(p);
    if (in != kOutside && G4int(sample.inside.size()) < gNumPoints)
    {
      sample.inside.push_back(p);
    }
    else if (in == kOutside && G4int(sample.outside.size()) < gNumPoints)
    {
      sample.outside.push_back(p);
    }

    if (++trials > 200 * gNumPoints)
    {
      G4Exception("benchPolyconePolyhedra::FillClassifiedPoints()", "GeomBench0001", FatalException,
                  "Could not generate enough classified benchmark points.");
    }
  }

  for (G4int i = 0; i < gNumPoints; ++i)
  {
    sample.surface.push_back(solid.GetPointOnSurface());
    sample.toInDirection.push_back((sample.inside[i] - sample.outside[i]).unit());
    sample.toOutDirection.push_back(RandomDirection());
  }
}

G4double ElapsedMs(clock_t start)
{
  return G4double(clock() - start) / CLOCKS_PER_SEC * 1000.;
}

void PrintTime(const char* label, G4double ms)
{
  if (gCsvMode) return;
  G4cout << "  " << std::left << std::setw(34) << label << std::right << std::setw(10) << ms
         << " ms" << G4endl;
}

G4double CheckInside(G4VSolid& solid, const Sample& sample, G4int nloop)
{
  G4int k = 0;
  G4int result = 0;
  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result += solid.Inside(sample.inside[k]);
    if (++k == gNumPoints) k = 0;
  }
  gIntSink += result;
  const G4double ms = ElapsedMs(start);
  PrintTime("Inside(p)", ms);
  return ms;
}

G4double CheckSurfaceNormal(G4VSolid& solid, const Sample& sample, G4int nloop)
{
  G4int k = 0;
  G4ThreeVector result;
  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result += solid.SurfaceNormal(sample.surface[k]);
    if (++k == gNumPoints) k = 0;
  }
  gDoubleSink += result.mag2();
  const G4double ms = ElapsedMs(start);
  PrintTime("SurfaceNormal(p)", ms);
  return ms;
}

G4double CheckDistanceToInSafety(G4VSolid& solid, const Sample& sample, G4int nloop)
{
  G4int k = 0;
  G4double result = 0.;
  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result += solid.DistanceToIn(sample.outside[k]);
    if (++k == gNumPoints) k = 0;
  }
  gDoubleSink += result;
  const G4double ms = ElapsedMs(start);
  PrintTime("DistanceToIn(p)", ms);
  return ms;
}

G4double CheckDistanceToInVector(G4VSolid& solid, const Sample& sample, G4int nloop)
{
  G4int k = 0;
  G4double result = 0.;
  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result += solid.DistanceToIn(sample.outside[k], sample.toInDirection[k]);
    if (++k == gNumPoints) k = 0;
  }
  gDoubleSink += result;
  const G4double ms = ElapsedMs(start);
  PrintTime("DistanceToIn(p,v)", ms);
  return ms;
}

G4double CheckDistanceToOutSafety(G4VSolid& solid, const Sample& sample, G4int nloop)
{
  G4int k = 0;
  G4double result = 0.;
  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result += solid.DistanceToOut(sample.inside[k]);
    if (++k == gNumPoints) k = 0;
  }
  gDoubleSink += result;
  const G4double ms = ElapsedMs(start);
  PrintTime("DistanceToOut(p)", ms);
  return ms;
}

G4double CheckDistanceToOutVector(G4VSolid& solid, const Sample& sample, G4bool calcNorm,
                                  G4int nloop)
{
  G4int k = 0;
  G4int trueNormCount = 0;
  G4bool validNorm = false;
  G4ThreeVector norm;
  G4double result = 0.;

  const clock_t start = clock();
  for (G4int i = 0; i < nloop; ++i)
  {
    result +=
      solid.DistanceToOut(sample.inside[k], sample.toOutDirection[k], calcNorm, &validNorm, &norm);
    if (validNorm) ++trueNormCount;
    if (++k == gNumPoints) k = 0;
  }

  gDoubleSink += result + norm.mag2();
  gIntSink += trueNormCount;
  const G4double ms = ElapsedMs(start);
  PrintTime(calcNorm ? "DistanceToOut(p,v,n)" : "DistanceToOut(p,v)", ms);
  return ms;
}

void RunCase(Config& config, G4int nloop)
{
  G4cout << "\n"
         << config.solidName << "  profile=" << ProfileName(config.profileKind)
         << "  z-planes=" << config.zSections
         << "  phi=" << (config.fullPhi ? "full" : "partial")
         << "  rmin=" << (config.hasRmin ? "non-zero" : "zero") << G4endl;

  Sample sample;
  FillClassifiedPoints(*config.solid, sample);

  CheckInside(*config.solid, sample, nloop);
  CheckSurfaceNormal(*config.solid, sample, nloop);
  CheckDistanceToInSafety(*config.solid, sample, nloop);
  CheckDistanceToInVector(*config.solid, sample, nloop);
  CheckDistanceToOutSafety(*config.solid, sample, nloop);
  CheckDistanceToOutVector(*config.solid, sample, false, nloop);
  CheckDistanceToOutVector(*config.solid, sample, true, nloop);
}

void RunCaseCsv(Config& config, G4int sections, G4int nloop)
{
  Sample sample;
  FillClassifiedPoints(*config.solid, sample);

  const G4double inside = CheckInside(*config.solid, sample, nloop);
  const G4double normal = CheckSurfaceNormal(*config.solid, sample, nloop);
  const G4double dIn = CheckDistanceToInSafety(*config.solid, sample, nloop);
  const G4double dInVec = CheckDistanceToInVector(*config.solid, sample, nloop);
  const G4double dOut = CheckDistanceToOutSafety(*config.solid, sample, nloop);
  const G4double dOutVec = CheckDistanceToOutVector(*config.solid, sample, false, nloop);
  const G4double dOutVecNorm = CheckDistanceToOutVector(*config.solid, sample, true, nloop);

  G4cout << config.solidName << "," << ProfileName(config.profileKind) << "," << sections
         << "," << config.zSections << "," << (config.fullPhi ? "full" : "partial")
         << "," << (config.hasRmin ? "nonzero" : "zero") << "," << inside << "," << normal
         << "," << dIn << "," << dInVec << "," << dOut << "," << dOutVec << ","
         << dOutVecNorm << G4endl;
}

Config MakeConfig(const G4String& solidName, ProfileKind kind, G4int zPlanes, G4bool fullPhi,
                  G4bool hasRmin)
{
  const G4double sphi = fullPhi ? 0. : 25. * deg;
  const G4double dphi = fullPhi ? twopi : 270. * deg;

  std::ostringstream suffix;
  suffix << "_" << ProfileName(kind) << "_nz" << zPlanes
         << (fullPhi ? "_full" : "_partial") << (hasRmin ? "_rmin" : "_solid");

  if (solidName == "G4Polycone")
  {
    return {"G4Polycone", kind, zPlanes, fullPhi, hasRmin,
            MakePolycone("polycone" + suffix.str(), zPlanes, sphi, dphi, hasRmin, kind)};
  }
  if (solidName == "G4Polyhedra")
  {
    return {"G4Polyhedra", kind, zPlanes, fullPhi, hasRmin,
            MakePolyhedra("polyhedra" + suffix.str(), zPlanes, sphi, dphi, hasRmin, kind)};
  }

  return {"G4GenericPolycone", kind, zPlanes, fullPhi, hasRmin,
          MakeGenericPolycone("genericpolycone" + suffix.str(), zPlanes, sphi, dphi,
                              hasRmin, kind)};
}

std::vector<Config> MakeConfigs()
{
  std::vector<Config> configs;
  const ProfileKind profileKinds[] = {ProfileKind::Sinusoidal,
                                      ProfileKind::ZMonotonicConvex,
                                      ProfileKind::ZStepped};
  const G4int zSections[] = {2, 3, 4, 8, 100};
  const G4bool phiModes[] = {true, false};
  const G4bool rminModes[] = {false, true};

  for (auto kind : profileKinds)
  {
    for (auto nz : zSections)
    {
      if (!HasProfile(kind, nz))
      {
        continue;
      }
      for (auto fullPhi : phiModes)
      {
        for (auto hasRmin : rminModes)
        {
          configs.push_back(MakeConfig("G4Polycone", kind, nz, fullPhi, hasRmin));
          configs.push_back(MakeConfig("G4Polyhedra", kind, nz, fullPhi, hasRmin));
#ifndef G4GEOM_USE_USOLIDS
          configs.push_back(MakeConfig("G4GenericPolycone", kind, nz, fullPhi, hasRmin));
#endif
        }
      }
    }
  }
#ifdef G4GEOM_USE_USOLIDS
  G4cout << G4endl
         << "==================================================================================" << G4endl
         << ">>> NOTE: Benchmark for G4GenericPolycone disabled !" << G4endl
         << ">>> VecGeom not supporting polycones with non-monotonic Z-sections configurations." << G4endl
         << "=================================================================================="
         << G4endl;
#endif
  return configs;
}

void RunScan(G4int maxSections, G4int nloop)
{
  const G4bool phiModes[] = {true, false};
  const G4bool rminModes[] = {false, true};
  const ProfileKind profileKinds[] = {ProfileKind::Sinusoidal,
                                      ProfileKind::ZMonotonicConvex,
                                      ProfileKind::ZStepped};
  const G4String solidNames[] = {"G4Polycone", "G4Polyhedra"};

  gCsvMode = true;
  G4cout << "solid,profile,sections,zplanes,phi,rmin,inside_ms,surface_normal_ms,distance_to_in_ms,"
            "distance_to_in_vec_ms,distance_to_out_ms,distance_to_out_vec_ms,"
            "distance_to_out_vec_norm_ms"
         << G4endl;

  for (auto kind : profileKinds)
  {
    for (G4int sections = 1; sections <= maxSections; ++sections)
    {
      const G4int zPlanes = sections + 1;
      if (!HasProfile(kind, zPlanes))
      {
        continue;
      }
      for (const auto& solidName : solidNames)
      {
        for (auto fullPhi : phiModes)
        {
          for (auto hasRmin : rminModes)
          {
            auto config = MakeConfig(solidName, kind, zPlanes, fullPhi, hasRmin);
            RunCaseCsv(config, sections, nloop);
          }
        }
      }
    }
  }
}
}  // namespace

int main(int argc, char** argv)
{
  G4int nloop = kDefaultLoops;
  G4bool scanMode = false;
  G4int maxSections = 100;
  if (argc > 1 && G4String(argv[1]) == "scan")
  {
    scanMode = true;
    if (argc > 2)
    {
      maxSections = std::max(1, std::atoi(argv[2]));
    }
    if (argc > 3)
    {
      nloop = std::max(1, std::atoi(argv[3]));
    }
    if (argc > 4)
    {
      gNumPoints = std::max(1, std::atoi(argv[4]));
    }
  }
  else if (argc > 1)
  {
    nloop = std::max(1, std::atoi(argv[1]));
    if (argc > 2)
    {
      gNumPoints = std::max(1, std::atoi(argv[2]));
    }
  }

  CLHEP::HepRandom::setTheSeed(1234567);

  if (scanMode)
  {
    RunScan(maxSections, nloop);
    return 0;
  }

  G4cout << "*********************************************************************" << G4endl;
  G4cout << "* Benchmark for G4Polycone and G4Polyhedra                          *" << G4endl;
  G4cout << "* Methods: Inside, SurfaceNormal, DistanceToIn and DistanceToOut     *" << G4endl;
  G4cout << "* Configurations: 2/3/4/8/100 Z planes, full/partial phi, zero/nonzero R *"
         << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "\n       Number of calls per method : " << nloop << G4endl;
  G4cout << "       Number of sampled points   : " << gNumPoints << G4endl;

  auto configs = MakeConfigs();
  for (auto& config : configs)
  {
    RunCase(config, nloop);
  }

  if (gDoubleSink < 0. || gIntSink < 0)
  {
    G4cout << "Unused result guard: " << gDoubleSink << " " << gIntSink << G4endl;
  }
  return 0;
}
