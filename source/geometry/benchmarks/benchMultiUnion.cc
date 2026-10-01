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

// Build with the neighbouring CMakeLists.txt and -DG4TARGET=benchMultiUnion.
// Usage: benchMultiUnion [calls [points [repeats]]]
//        benchMultiUnion scan [max-components [calls [points [repeats]]]]
// Default: components 1/8/64/512/1024; scan: every count from 1 to maximum.
// Query times are median wall-clock milliseconds; setup is timed once.
// Memory is measured in a separate batch: C++ new/new[] requested bytes only,
// not malloc/realloc, allocator overhead, RSS or GPU/library-private allocators.
// peak_extra_bytes is the maximum live-byte increase over the phase baseline;
// retained_bytes is the signed live-byte change at the end of that phase.
// Sum build and voxelization retained bytes for the union's owned footprint.
// The shared component box and probe buffers are excluded. The allocator shim
// still adds header/accounting overhead during timing; compare like builds.
// No worker threads are started. Do not combine with another allocator profiler.
// The two layouts deliberately use unrotated unit-half-length boxes, fixed
// spacing and analytic probes spread over the component indices. Longitudinal
// exits traverse the overlapping chain; transverse exits exercise local lookup.
// Expensive volume/area estimators and visualization queries are not included.
#include "G4AffineTransform.hh"
#include "G4Box.hh"
#include "G4MultiUnion.hh"
#include "G4VoxelLimits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

// Replace all C++17 allocation/deallocation forms, including over-aligned
// allocations. Requested bytes exclude our header and allocator overhead.
namespace Allocation
{
struct Header
{
    void* base;
    std::size_t bytes;
};
std::atomic<std::size_t> live{0}, peak{0}, calls{0}, bytes{0};
std::atomic<bool> enabled{false};

void* Allocate(std::size_t size, std::size_t alignment)
{
  alignment = std::max(alignment, alignof(Header));
  const auto overhead = sizeof(Header) + alignment - 1;
  if (size > std::numeric_limits<std::size_t>::max() - overhead) throw std::bad_alloc();
  void* base = std::malloc(size + overhead);
  if (!base) throw std::bad_alloc();
  auto address = reinterpret_cast<std::uintptr_t>(base) + sizeof(Header);
  address = (address + alignment - 1) & ~(std::uintptr_t(alignment) - 1);
  auto* header = reinterpret_cast<Header*>(address) - 1;
  new (header) Header{base, size};
  const auto current = live.fetch_add(size, std::memory_order_relaxed) + size;
  if (enabled.load(std::memory_order_relaxed))
  {
    calls.fetch_add(1, std::memory_order_relaxed);
    bytes.fetch_add(size, std::memory_order_relaxed);
    auto previous = peak.load(std::memory_order_relaxed);
    while (previous < current && !peak.compare_exchange_weak(previous, current))
    {}
  }
  return reinterpret_cast<void*>(address);
}
void Free(void* pointer) noexcept
{
  if (!pointer) return;
  auto* header = static_cast<Header*>(pointer) - 1;
  live.fetch_sub(header->bytes, std::memory_order_relaxed);
  std::free(header->base);
}
}  // namespace Allocation
void* operator new(std::size_t n)
{
  return Allocation::Allocate(n, alignof(std::max_align_t));
}
void* operator new[](std::size_t n)
{
  return ::operator new(n);
}
void operator delete(void* p) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p) noexcept
{
  Allocation::Free(p);
}
void operator delete(void* p, std::size_t) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p, std::size_t) noexcept
{
  Allocation::Free(p);
}
void* operator new(std::size_t n, std::align_val_t a)
{
  return Allocation::Allocate(n, static_cast<std::size_t>(a));
}
void* operator new[](std::size_t n, std::align_val_t a)
{
  return ::operator new(n, a);
}
void operator delete(void* p, std::align_val_t) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p, std::align_val_t) noexcept
{
  Allocation::Free(p);
}
void operator delete(void* p, std::size_t, std::align_val_t) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept
{
  Allocation::Free(p);
}
void* operator new(std::size_t n, const std::nothrow_t&) noexcept
{
  try
  {
    return ::operator new(n);
  }
  catch (...)
  {
    return nullptr;
  }
}
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept
{
  try
  {
    return ::operator new[](n);
  }
  catch (...)
  {
    return nullptr;
  }
}
void operator delete(void* p, const std::nothrow_t&) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p, const std::nothrow_t&) noexcept
{
  Allocation::Free(p);
}
void* operator new(std::size_t n, std::align_val_t a, const std::nothrow_t&) noexcept
{
  try
  {
    return ::operator new(n, a);
  }
  catch (...)
  {
    return nullptr;
  }
}
void* operator new[](std::size_t n, std::align_val_t a, const std::nothrow_t&) noexcept
{
  try
  {
    return ::operator new[](n, a);
  }
  catch (...)
  {
    return nullptr;
  }
}
void operator delete(void* p, std::align_val_t, const std::nothrow_t&) noexcept
{
  Allocation::Free(p);
}
void operator delete[](void* p, std::align_val_t, const std::nothrow_t&) noexcept
{
  Allocation::Free(p);
}

namespace
{
volatile double gSink = 0.;
struct Memory
{
    std::size_t calls, bytes, peak;
    std::int64_t retained;
};
template<class F>
Memory MeasureMemory(F operation)
{
  const auto baseline = Allocation::live.load();
  Allocation::calls = 0;
  Allocation::bytes = 0;
  Allocation::peak = baseline;
  Allocation::enabled = true;
  try
  {
    operation();
  }
  catch (...)
  {
    Allocation::enabled = false;
    throw;
  }
  Allocation::enabled = false;
  return {Allocation::calls.load(), Allocation::bytes.load(), Allocation::peak.load() - baseline,
          static_cast<std::int64_t>(Allocation::live.load()) - static_cast<std::int64_t>(baseline)};
}
template<class F>
double Time(F operation)
{
  const auto start = std::chrono::steady_clock::now();
  operation();
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
    .count();
}
struct Options
{
    int loops = 16384, points = 1024, repeats = 3, maximum = 1024;
    bool scan = false;
};
struct Probe
{
    G4ThreeVector inside, outside, surface;
    double exit;
};
G4ThreeVector Centre(int i, int n, bool chain)
{
  if (chain) return {1.5 * i, 0., 0.};
  const int side = static_cast<int>(std::ceil(std::cbrt(n)));
  return {3. * (i % side), 3. * ((i / side) % side), 3. * (i / (side * side))};
}
void Check(bool valid, const char* message)
{
  if (!valid) throw std::runtime_error(message);
}
void Near(double value, double expected)
{
  Check(std::isfinite(value) && std::abs(value - expected) < 1.e-7,
        "Probe distance/normal validation failed");
}

void Print(const Options& o, bool chain, int n, const char* label, int count, double ms,
           const Memory& memory)
{
  if (o.scan)
  {
    std::cout << "G4MultiUnion," << (chain ? "overlapping-chain" : "separated-grid") << ',' << n
              << ",\"" << label << "\"," << count << ',' << ms << ',' << memory.calls << ','
              << memory.bytes << ',' << double(memory.calls) / count << ','
              << double(memory.bytes) / count << ',' << memory.peak << ',' << memory.retained
              << '\n';
    return;
  }
  std::cout << "  " << std::left << std::setw(38) << label << std::right << std::setw(12) << ms
            << " ms  " << std::setw(12) << double(memory.calls) / count << " alloc/call  "
            << std::setw(12) << double(memory.bytes) / count << " B/call  " << std::setw(12)
            << memory.peak << " peak B  " << std::setw(12) << memory.retained << " retained B\n";
}

void RunCase(const Options& o, int n, bool chain)
{
  if (!o.scan)
    std::cout << "\nG4MultiUnion  layout=" << (chain ? "overlapping-chain" : "separated-grid")
              << "  components=" << n << '\n';
  // Reuse one box, as in a geometry containing repeated placements. Declare
  // the component first so it outlives the union (which does not own it).
  G4Box box("component", 1., 1., 1.);
  std::unique_ptr<G4MultiUnion> solid;
  auto build = [&] {
    solid = std::make_unique<G4MultiUnion>("union");
    for (int i = 0; i < n; ++i)
      solid->AddNode(box, G4Transform3D(G4RotationMatrix(), Centre(i, n, chain)));
  };
  double buildMs = Time(build);
  double voxelMs = Time([&] { solid->Voxelize(); });
  solid.reset();
  const auto buildMemory = MeasureMemory(build);
  const auto voxelMemory = MeasureMemory([&] { solid->Voxelize(); });
  Check(
    buildMemory.calls > 0 && voxelMemory.calls > 0,
    "Allocation interception did not observe Geant4 allocations; use a compatible/static build");
  Print(o, chain, n, "AddNode/build", 1, buildMs, buildMemory);
  Print(o, chain, n, "Voxelize", 1, voxelMs, voxelMemory);

  std::vector<Probe> probes;
  probes.reserve(o.points);
  for (int k = 0; k < o.points; ++k)
  {
    const int i = static_cast<int>((std::uint64_t(k) * n) / o.points);
    const auto c = Centre(i, n, chain);
    probes.push_back({c, c + G4ThreeVector(0., 1.5, 0.), c + G4ThreeVector(0., 1., 0.),
                      chain ? 1.5 * (n - 1 - i) + 1. : 1.});
  }
  const G4ThreeVector down(0., -1., 0.), along(1., 0., 0.);
  G4ThreeVector lo, hi;
  solid->BoundingLimits(lo, hi);
  const G4ThreeVector miss = hi + G4ThreeVector(2., 2., 2.);
  G4VoxelLimits limits;
  G4AffineTransform transform;
  for (const auto& p : probes)
  {
    Check(solid->Inside(p.inside) == kInside && solid->Inside(p.outside) == kOutside
            && solid->Inside(p.surface) == kSurface,
          "Probe classification failed");
    Near(solid->DistanceToIn(p.outside, down), .5);
    Near(solid->DistanceToOut(p.inside, along), p.exit);
    Check(solid->InsideNoVoxels(p.inside) == kInside, "Non-voxel classification failed");
    Near(solid->DistanceToInNoVoxels(p.outside, down), .5);
    G4ThreeVector normal;
    G4bool valid = false;
    Near(solid->DistanceToOutNoVoxels(p.inside, down, &normal), 1.);
    Near(solid->DistanceToOut(p.inside, along, true, &valid, &normal), p.exit);
    // calcNorm/validNorm are unused by G4MultiUnion; validate distance only.
    Near((solid->SurfaceNormal(p.surface) - G4ThreeVector(0., 1., 0.)).mag(), 0.);
    for (bool accurate : {false, true})
    {
      solid->SetAccurateSafety(accurate);
      const double safety = solid->DistanceToIn(p.outside);
      Check(safety >= 0. && safety <= .5 + 1.e-7, "Outside safety is not conservative");
      Near(solid->DistanceToOut(p.inside), 1.);
    }
  }
  solid->SetAccurateSafety(false);
  Check(solid->DistanceToIn(miss, along) == kInfinity, "Miss ray intersects union");

  // Sample generation, validation, warmup and output are excluded from both
  // timing and allocation scopes. Each query uses the same probe sequence.
  auto query = [&](const char* label, auto operation) {
    auto batch = [&] {
      double result = 0.;
      int k = 0;
      for (int j = 0; j < o.loops; ++j)
      {
        result += static_cast<double>(operation(probes[k]));
        if (++k == o.points) k = 0;
      }
      gSink = result;
    };
    gSink = static_cast<double>(operation(probes.front()));
    std::vector<double> times;
    times.reserve(o.repeats);
    for (int r = 0; r < o.repeats; ++r)
      times.push_back(Time(batch));
    std::sort(times.begin(), times.end());
    const auto memory = MeasureMemory(batch);
    const double median = .5 * (times[(times.size() - 1) / 2] + times[times.size() / 2]);
    Print(o, chain, n, label, o.loops, median, memory);
  };
  query("Inside(p)-inside", [&](const Probe& p) { return solid->Inside(p.inside); });
  query("Inside(p)-outside", [&](const Probe& p) { return solid->Inside(p.outside); });
  query("Inside(p)-surface", [&](const Probe& p) { return solid->Inside(p.surface); });
  query("SurfaceNormal(p)", [&](const Probe& p) { return solid->SurfaceNormal(p.surface).y(); });
  query("DistanceToIn(p)", [&](const Probe& p) { return solid->DistanceToIn(p.outside); });
  solid->SetAccurateSafety(true);
  query("DistanceToIn(p)-accurate", [&](const Probe& p) { return solid->DistanceToIn(p.outside); });
  solid->SetAccurateSafety(false);
  query("DistanceToOut(p)", [&](const Probe& p) { return solid->DistanceToOut(p.inside); });
  query("DistanceToIn(p,v)-hit",
        [&](const Probe& p) { return solid->DistanceToIn(p.outside, down); });
  query("DistanceToIn(p,v)-miss",
        [&](const Probe&) { return solid->DistanceToIn(miss, along) == kInfinity; });
  query("DistanceToOut(p,v)",
        [&](const Probe& p) { return solid->DistanceToOut(p.inside, along); });
  query("DistanceToOut(p,v,n)", [&](const Probe& p) {
    G4ThreeVector normal;
    G4bool valid = false;
    const auto d = solid->DistanceToOut(p.inside, along, true, &valid, &normal);
    return d + normal.x();  // G4MultiUnion does not set validNorm.
  });
  query("InsideNoVoxels(p)", [&](const Probe& p) { return solid->InsideNoVoxels(p.inside); });
  query("DistanceToInNoVoxels(p,v)",
        [&](const Probe& p) { return solid->DistanceToInNoVoxels(p.outside, down); });
  query("DistanceToOut(p,v)-transverse",
        [&](const Probe& p) { return solid->DistanceToOut(p.inside, down); });
  query("DistanceToOutNoVoxels(p,v)-transverse", [&](const Probe& p) {
    G4ThreeVector normal;
    return solid->DistanceToOutNoVoxels(p.inside, down, &normal);
  });
  query("BoundingLimits", [&](const Probe&) {
    G4ThreeVector a, b;
    solid->BoundingLimits(a, b);
    return a.x() + b.x();
  });
  query("CalculateExtent", [&](const Probe&) {
    double a = 0., b = 0.;
    const bool valid = solid->CalculateExtent(kXAxis, limits, transform, a, b);
    return valid ? b - a : 0.;
  });
  query("GetPointOnSurface", [&](const Probe&) { return solid->GetPointOnSurface().x(); });
}
int Positive(const char* value)
{
  std::size_t end = 0;
  const auto n = std::stoll(value, &end);
  if (end != std::string(value).size() || n < 1 || n > std::numeric_limits<int>::max())
    throw std::invalid_argument("Expected a positive integer");
  return static_cast<int>(n);
}
}  // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc > 1 && std::string(argv[1]) == "--help")
    {
      std::cout << "Usage: benchMultiUnion [calls [points [repeats]]]\n"
                   "       benchMultiUnion scan [max-components [calls [points [repeats]]]]\n";
      return 0;
    }
    Options o;
    int arg = 1;
    if (arg < argc && std::string(argv[arg]) == "scan")
    {
      o.scan = true;
      ++arg;
      if (arg < argc) o.maximum = Positive(argv[arg++]);
    }
    if (arg < argc) o.loops = Positive(argv[arg++]);
    if (arg < argc) o.points = Positive(argv[arg++]);
    if (arg < argc) o.repeats = Positive(argv[arg++]);
    if (arg != argc) throw std::invalid_argument("Too many arguments (see --help)");
    CLHEP::HepRandom::setTheSeed(1234567);
    std::cout << std::fixed << std::setprecision(3);
    if (o.scan)
      std::cout << "solid,layout,components,query,calls,time_ms,allocations,allocated_bytes,"
                   "allocations_per_call,bytes_per_call,peak_extra_bytes,retained_bytes\n";
    else
      std::cout << "*********************************************************************\n"
                   "* Benchmark for G4MultiUnion                                        *\n"
                   "* Navigation queries, construction and memory scalability           *\n"
                   "*********************************************************************\n"
                << "\n       Number of calls per method : " << o.loops
                << "\n       Number of sampled points   : " << o.points
                << "\n       Number of timing batches   : " << o.repeats
                << "\n       Query time: median; memory: separate batch, requested C++ bytes\n";
    for (bool chain : {false, true})
    {
      if (o.scan)
      {
        for (int n = 1;; ++n)
        {
          RunCase(o, n, chain);
          if (n == o.maximum) break;
        }
      }
      else
        for (int n : {1, 8, 64, 512, 1024})
          RunCase(o, n, chain);
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << "benchMultiUnion: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
