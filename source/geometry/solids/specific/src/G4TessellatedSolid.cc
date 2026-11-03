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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4TessellatedSolid implementation
//
// Closing the solid builds cached facet geometry and a hybrid median/surface-
// area hierarchy. Twelve-bin splits must improve the median cost by more than
// 10%; small nodes and depths of 24 or more use median splits. No hierarchy is
// built below eight facets or for extremely ill-conditioned facets.
//
// Nearest-surface queries visit nearer bounds first and evaluate finite facets
// at leaves. Safety is a finite-surface distance, not a bounding-box distance;
// the hierarchy returns zero within half the surface tolerance.
// SurfaceNormal first searches the surface band, then uses the nearest facet.
//
// Inside normally uses the nearest face's material side for strictly interior
// projections. On large meshes, costly nearest searches switch to an oriented
// axial ray, with finite-facet surface checks in the same traversal. Near-surface
// searches retain the nearest-face path; ambiguous crossings, tangencies and
// misses retain nearest-face/alternate-ray fallbacks. Classification supports
// concavity, cavities, handles and disconnected components.
//
// Directional traversal uses padded bounds and the existing facet intersection
// kernels. Voxel entry/exit traversal remains below 4,096/64 facets respectively,
// or when no hierarchy is available. Hierarchy query scratch uses bounded local
// stacks; the legacy fallbacks retain their original traversal and storage.
// During closure, temporary directional support bounds accelerate supporting-
// plane certification without increasing retained memory.
//
// Original author: P.R.Truscott, QinetiQ Ltd (UK) - 31.10.2004
//
// CHANGE HISTORY
// --------------
//
// 17.09.2007, P.R.Truscott, QinetiQ Ltd & Richard Holmberg
//             Random rays are predefined rather than making use of random
//             number generator at run-time.
// 12.10.2012, M.Gayer, CERN
//             Complete rewrite reducing memory requirements more than 50%
//             and speedup with voxelization of surface.
// 23.10.2016, E.Tcherniaev, CERN
//             Reimplemented CalculateExtent() to make use of G4BoundingEnvelope.
// 09.09.2026, G.Cosmo, CERN
//             Accelerated Inside, normals, finite-surface safety and exact
//             distances using a facet hierarchy and cached geometry; reduced
//             closure work and memory. Improved scaling on number of facets.
// --------------------------------------------------------------------

#include "G4TessellatedSolid.hh"

#if !defined(G4GEOM_USE_UTESSELLATEDSOLID)

#  include "G4AffineTransform.hh"
#  include "G4AutoLock.hh"
#  include "G4BoundingEnvelope.hh"
#  include "G4GeometryTolerance.hh"
#  include "G4PhysicalConstants.hh"
#  include "G4QuadrangularFacet.hh"
#  include "G4QuickRand.hh"
#  include "G4SystemOfUnits.hh"
#  include "G4TriangularFacet.hh"
#  include "G4VGraphicsScene.hh"
#  include "G4VisExtent.hh"
#  include "G4VoxelLimits.hh"

#  include "geomdefs.hh"

#  include <algorithm>
#  include <array>
#  include <iostream>
#  include <limits>
#  include <numeric>
#  include <stack>

namespace
{
  G4Mutex tessMutex = G4MUTEX_INITIALIZER;
}

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Default constructor: name "dummy", no facets.
//
G4TessellatedSolid::G4TessellatedSolid() : G4VSolid("dummy")
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Named constructor: initialises an empty solid.
//
G4TessellatedSolid::G4TessellatedSolid(const G4String& name) : G4VSolid(name)
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Constructor reserved for object persistency; initialises an empty solid.
//
G4TessellatedSolid::G4TessellatedSolid(__void__& a) : G4VSolid(a)
{
  Initialize();
  fMinExtent.set(0, 0, 0);
  fMaxExtent.set(0, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
G4TessellatedSolid::~G4TessellatedSolid()
{
  DeleteObjects();
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor.
//
G4TessellatedSolid::G4TessellatedSolid(const G4TessellatedSolid& ts) : G4VSolid(ts)
{
  Initialize();

  CopyObjects(ts);
}

///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator.
//
G4TessellatedSolid& G4TessellatedSolid::operator=(const G4TessellatedSolid& ts)
{
  if (&ts == this)
  {
    return *this;
  }

  // Copy base class data
  G4VSolid::operator=(ts);

  DeleteObjects();

  Initialize();

  CopyObjects(ts);

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::Initialize()
{
  kCarToleranceHalf = 0.5 * kCarTolerance;

  fRebuildPolyhedron = false;
  fpPolyhedron = nullptr;
  fCubicVolume = 0.;
  fSurfaceArea = 0.;

  fGeometryType = "G4TessellatedSolid";
  fSolidClosed = false;

  fMinExtent.set(kInfinity, kInfinity, kInfinity);
  fMaxExtent.set(-kInfinity, -kInfinity, -kInfinity);

  SetRandomVectors();
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DeleteObjects()
{
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    delete fFacets[i];
  }
  fFacets.clear();
  fFacetList.clear();
  fVoxels.Clear();
  fInsides.Clear();
  fVertexList.clear();
  fExtremeFacets.clear();
  fFacetPlanes.clear();
  fFacetTree.clear();
  fFacetOrder.clear();
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::CopyObjects(const G4TessellatedSolid& ts)
{
  G4ThreeVector reductionRatio;
  G4int fmaxVoxels = ts.fVoxels.GetMaxVoxels(reductionRatio);
  if (fmaxVoxels < 0)
  {
    fVoxels.SetMaxVoxels(reductionRatio);
  }
  else
  {
    fVoxels.SetMaxVoxels(fmaxVoxels);
  }

  G4int n = ts.GetNumberOfFacets();
  for (G4int i = 0; i < n; ++i)
  {
    G4VFacet* facetClone = (ts.GetFacet(i))->GetClone();
    AddFacet(facetClone);
  }
  if (ts.GetSolidClosed())
  {
    SetSolidClosed(true);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// AddFacet
//
// Add a facet to the facet list.
// Note that you can add, but you cannot delete.
//
G4bool G4TessellatedSolid::AddFacet(G4VFacet* aFacet)
{
  // Add the facet to the vector.
  //
  if (fSolidClosed)
  {
    G4Exception("G4TessellatedSolid::AddFacet()", "GeomSolids1002", JustWarning,
                "Attempt to add facets when solid is closed.");
    return false;
  }
  if (aFacet->IsDefined())
  {
    fCubicVolume = 0.;
    fSurfaceArea = 0.;
    set<G4VertexInfo, G4VertexComparator>::iterator begin = fFacetList.begin(),
                                                    end = fFacetList.end(), pos, it;
    G4ThreeVector p = aFacet->GetCircumcentre();
    G4VertexInfo value;
    value.id = (G4int)fFacetList.size();
    value.mag2 = p.x() + p.y() + p.z();

    G4bool found = false;
    if (!OutsideOfExtent(p, kCarTolerance))
    {
      G4double kCarTolerance3 = 3 * kCarTolerance;
      pos = fFacetList.lower_bound(value);

      it = pos;
      while (!found && it != end)  // Loop checking, 13.08.2015, G.Cosmo
      {
        G4int id = (*it).id;
        G4VFacet* facet = fFacets[id];
        G4ThreeVector q = facet->GetCircumcentre();
        if ((found = (facet == aFacet)))
        {
          break;
        }
        G4double dif = q.x() + q.y() + q.z() - value.mag2;
        if (dif > kCarTolerance3)
        {
          break;
        }
        it++;
      }

      if (fFacets.size() > 1)
      {
        it = pos;
        while (!found && it != begin)  // Loop checking, 13.08.2015, G.Cosmo
        {
          --it;
          G4int id = (*it).id;
          G4VFacet* facet = fFacets[id];
          G4ThreeVector q = facet->GetCircumcentre();
          found = (facet == aFacet);
          if (found)
          {
            break;
          }
          G4double dif = value.mag2 - (q.x() + q.y() + q.z());
          if (dif > kCarTolerance3)
          {
            break;
          }
        }
      }
    }

    if (!found)
    {
      fFacets.push_back(aFacet);
      fFacetList.insert(value);
    }
    return true;
  }

  G4Exception("G4TessellatedSolid::AddFacet()", "GeomSolids1002", JustWarning,
              "Attempt to add facet not properly defined.");
  aFacet->StreamInfo(G4cout);
  return false;
}

///////////////////////////////////////////////////////////////////////////////
//
// BuildFacetTree
//
// Immutable, surface-area-split hierarchy. Leaves retain the original facet
// kernels, including their treatment of tolerance bands and grazing
// intersections.
//
void G4TessellatedSolid::BuildFacetTree()
{

  fFacetPlanes.clear();
  fFacetTree.clear();
  fFacetOrder.clear();
  if (fFacets.size() < 8) return;
  std::vector<FacetTreeNode> bounds(fFacets.size());
  for (std::size_t i = 0; i < fFacets.size(); ++i)
  {
    const auto& facet = *fFacets[i];
    auto& box = bounds[i];
    box.lower = box.upper = facet.GetVertex(0);
    G4double padding = 4. * kCarTolerance;
    for (G4int j = 0; j < facet.GetNumberOfVertices(); ++j)
    {
      for (G4int axis = 0; axis < 3; ++axis)
      {
        const auto value = facet.GetVertex(j)[axis];
        box.lower[axis] = std::min(box.lower[axis], value);
        box.upper[axis] = std::max(box.upper[axis], value);
      }
    }
    // Cover the barycentric acceptance band of the triangle Intersect kernel,
    // including elongated triangles, as well as coordinate roundoff.
    //
    for (G4int j = 1; j + 1 < facet.GetNumberOfVertices(); ++j)
    {
      const auto a = facet.GetVertex(j) - facet.GetVertex(0);
      const auto b = facet.GetVertex(j + 1) - facet.GetVertex(0);
      const auto length = std::max({a.mag(), b.mag(), (a - b).mag()});
      const auto area2 = a.cross(b).mag2();

      // For extremely small/ill-conditioned facets the legacy barycentric
      // tolerance can extend far beyond the finite polygon. Keep its original
      // traversal rather than assuming that a finite box bounds accepted hits.
      //
      if (area2 <= 0. || kCarTolerance*length*length / area2 >= .0625)
      {
        return;
      }
      padding = std::max(padding, 64.*kCarTolerance*length*length*length / area2);
    }
    for (G4int axis = 0; axis < 3; ++axis)
    {
      const auto roundoff = 64. * std::numeric_limits<G4double>::epsilon()
            * std::max(std::abs(box.lower[axis]), std::abs(box.upper[axis]));
      box.padding = std::max(box.padding, padding + roundoff);
    }
  }

  // Populate every finite-face cache before making the tree available.
  // Tree queries can therefore rely on the cached local half spaces.
  //
  fFacetPlanes.resize(fFacets.size());
  for (std::size_t i = 0; i < fFacets.size(); ++i)
  {
    const auto& facet = *fFacets[i];
    auto& data = fFacetPlanes[i];
    data.anchor = facet.GetVertex(0);
    data.normal = facet.GetSurfaceNormal();
    data.count = facet.GetNumberOfVertices();
    for (G4int j = 0; j < data.count; ++j)
    {
      const auto vertex = facet.GetVertex(j);
      const auto edge = facet.GetVertex((j + 1) % data.count) - vertex;
      data.inward[j] = data.normal.cross(edge).unit();
      data.offset[j] = data.inward[j].dot(vertex - data.anchor);
    }
  }
  fFacetOrder.resize(fFacets.size());
  std::iota(fFacetOrder.begin(), fFacetOrder.end(), 0);
  fFacetTree.reserve(fFacets.size());
  BuildFacetNode(0, static_cast<G4int>(fFacets.size()), bounds);
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::BuildFacetNode(G4int first, G4int last,
                          const std::vector<FacetTreeNode>& bounds, G4int depth)
{
  const auto index = static_cast<G4int>(fFacetTree.size());
  fFacetTree.push_back(bounds[fFacetOrder[first]]);
  for (auto i = first + 1; i < last; ++i)
  {
    const auto& box = bounds[fFacetOrder[i]];
    auto& node = fFacetTree[index];
    for (G4int axis = 0; axis < 3; ++axis)
    {
      node.lower[axis] = std::min(node.lower[axis], box.lower[axis]);
      node.upper[axis] = std::max(node.upper[axis], box.upper[axis]);
    }
    node.padding = std::max(node.padding, box.padding);
  }
  if (last - first <= 4)
  {
    fFacetTree[index].first = first;
    fFacetTree[index].count = last - first;
    return index;
  }
  const auto centre = [&](G4int facet, G4int axis)
  {
    return .5 * bounds[facet].lower[axis] + .5 * bounds[facet].upper[axis];
  };
  const auto extent = fFacetTree[index].upper - fFacetTree[index].lower;
  const G4int medianAxis = extent.x() > extent.y()
      ? (extent.x() > extent.z() ? 0 : 2) : (extent.y() > extent.z() ? 1 : 2);
  G4int middle = first + (last - first) / 2;
  std::nth_element(fFacetOrder.begin() + first, fFacetOrder.begin() + middle,
                   fFacetOrder.begin() + last, [&](G4int a, G4int b)
                   {
                     const auto ca = centre(a, medianAxis),
                                cb = centre(b, medianAxis);
                     return ca < cb || (ca == cb && a < b);
                   });

  // Binned surface-area splits avoid grouping large overlapping facet bounds
  // solely by their centres. No facets are duplicated. Keep leaves of 2--4
  // facets, so the existing one-node-per-facet reservation remains sufficient.
  // After 24 levels, median splits bound the remaining depth by 31 for int
  // facet counts. Every traversal's 64-entry stack therefore remains safe.
  //
  if (last - first > 8 && depth < 24)
  {
    constexpr G4int nBins = 12;
    const auto binOf = [&](G4int facet, G4int axis, G4double lo, G4double scale)
    {
      return std::min(nBins - 1,
                      static_cast<G4int>((centre(facet, axis) - lo) * scale));
    };
    struct Bin
    {
      G4ThreeVector lower, upper;
      G4int count = 0;
      void Add(const G4ThreeVector& lo, const G4ThreeVector& hi, G4int n)
      {
        if (!n) { return; }
        if (!count)
        {
          lower = lo;
          upper = hi;
        }
        else
        {
          for (G4int axis = 0; axis < 3; ++axis)
          {
            lower[axis] = std::min(lower[axis], lo[axis]);
            upper[axis] = std::max(upper[axis], hi[axis]);
          }
        }
        count += n;
      }
      G4double Cost() const
      {
        const auto d = upper - lower;
        return count * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
      }
    };
    Bin medianLeft, medianRight;
    for (G4int i = first; i < last; ++i)
    {
      const auto facet = fFacetOrder[i];
      (i < middle ? medianLeft : medianRight).Add(bounds[facet].lower,
                                                  bounds[facet].upper, 1);
    }
    G4int bestAxis = -1, bestBin = -1;

    // Avoid disturbing a good balanced split for a marginal area improvement.
    //
    G4double bestCost = .9 * (medianLeft.Cost() + medianRight.Cost());
    G4double bestMin = 0., bestScale = 0.;
    for (G4int axis = 0; axis < 3; ++axis)
    {
      G4double lo = centre(fFacetOrder[first], axis), hi = lo;
      for (G4int i = first + 1; i < last; ++i)
      {
        const auto c = centre(fFacetOrder[i], axis);
        lo = std::min(lo, c);
        hi = std::max(hi, c);
      }
      if (!(hi > lo)) { continue; }
      const auto scale = nBins / (hi - lo);
      if (!(scale > 0.) || !std::isfinite(scale)) { continue; }
      std::array<Bin, nBins> bins, suffix;
      for (G4int i = first; i < last; ++i)
      {
        const auto facet = fFacetOrder[i];
        bins[binOf(facet, axis, lo, scale)].Add(bounds[facet].lower,
                                                bounds[facet].upper, 1);
      }
      Bin rightBox;
      for (G4int i = nBins - 1; i >= 0; --i)
      {
        rightBox.Add(bins[i].lower, bins[i].upper, bins[i].count);
        suffix[i] = rightBox;
      }
      Bin leftBox;
      for (G4int i = 0; i < nBins - 1; ++i)
      {
        leftBox.Add(bins[i].lower, bins[i].upper, bins[i].count);
        if (leftBox.count < 2 || suffix[i + 1].count < 2) { continue; }
        const auto cost = leftBox.Cost() + suffix[i + 1].Cost();
        if (cost < bestCost)
        {
          bestCost = cost;
          bestAxis = axis;
          bestBin = i;
          bestMin = lo;
          bestScale = scale;
        }
      }
    }
    if (bestAxis >= 0)
    {
      const auto split =
        std::partition(fFacetOrder.begin()+first, fFacetOrder.begin()+last, [&](G4int facet)
        {
          return binOf(facet, bestAxis, bestMin, bestScale) <= bestBin;
        });
      middle = static_cast<G4int>(split - fFacetOrder.begin());
    }
  }
  BuildFacetNode(first, middle, bounds, depth + 1);  // Left child immediately follows parent.
  const auto right = BuildFacetNode(middle, last, bounds, depth + 1);
  fFacetTree[index].right = right;

  return index;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::FacetBoxDistance(const FacetTreeNode& node,
                                              const G4ThreeVector& p)
{
  G4double squared = 0.;
  for (G4int axis = 0; axis < 3; ++axis)
  {
    const auto d = std::max({node.lower[axis] - p[axis],
                             p[axis] - node.upper[axis], 0.});
    squared += d * d;
  }
  return squared;
}

///////////////////////////////////////////////////////////////////////////////
//
// FacetBoxEntry
//
// In reciprocal mode v contains inverse components; zero still denotes a
// parallel slab. Callers retain division mode if an inverse is not finite.
//
template<G4bool reciprocal>
G4double G4TessellatedSolid::FacetBoxEntry(const FacetTreeNode& node, const G4ThreeVector& p,
                                           const G4ThreeVector& v, G4double limit)
{
  G4double enter = 0., leave = limit;
  for (G4int axis = 0; axis < 3; ++axis)
  {
    // Explicit parallel-slab handling, also for negative zero.
    if (v[axis] == 0.)
    {
      if (p[axis] < node.lower[axis] - node.padding
       || p[axis] > node.upper[axis] + node.padding)
      {
        return kInfinity;
      }
      continue;
    }
    auto a = node.lower[axis] - node.padding - p[axis];
    auto b = node.upper[axis] + node.padding - p[axis];
    if constexpr (reciprocal)
    {
      a *= v[axis];
      b *= v[axis];
    }
    else
    {
      a /= v[axis];
      b /= v[axis];
    }
    if (a > b) { std::swap(a, b); }
    enter = std::max(enter, a);
    leave = std::min(leave, b);
    if (enter > leave) { return kInfinity; }
  }
  return enter;
}

///////////////////////////////////////////////////////////////////////////////
//
// InteriorSurfaceFacet
//
// Return a surface facet only when its projection is away from the edges.
//
G4int G4TessellatedSolid::InteriorSurfaceFacet(const G4ThreeVector& p) const
{
  if (!HasFacetTree()) { return -1; }
  const auto index = SurfaceFacet(p);
  return index >= 0 && FacetProjectionInside(index, p) ? index : -1;
}

///////////////////////////////////////////////////////////////////////////////
//
// SurfaceFacet
//
// SurfaceNormal only needs a facet in the surface band, not the global
// nearest distance. A point-overlap traversal avoids distance ordering and
// stops at the first verified finite facet. Off-surface normals still use
// the nearest-facet search below.
//
G4int G4TessellatedSolid::SurfaceFacet(const G4ThreeVector& p) const
{

  std::array<G4int, 64> stack;
  G4int size = 0;
  stack[size++] = 0;
  while (size)
  {
    const auto current = stack[--size];
    const auto& node = fFacetTree[current];
    if (p.x() < node.lower.x() - node.padding
     || p.x() > node.upper.x() + node.padding
     || p.y() < node.lower.y() - node.padding
     || p.y() > node.upper.y() + node.padding
     || p.z() < node.lower.z() - node.padding
     || p.z() > node.upper.z() + node.padding)
    {
      continue;
    }
    if (node.count)
    {
      for (auto i = node.first; i < node.first + node.count; ++i)
      {
        const auto candidate = fFacetOrder[i];
        if (FacetDistance(candidate, p, kCarTolerance) <= kCarToleranceHalf)
        {
          return candidate;
        }
      }
    }
    else
    {
      stack[size++] = current + 1;
      stack[size++] = node.right;
    }
  }
  return -1;
}

///////////////////////////////////////////////////////////////////////////////
//
// NearestFacet
//
// Construction bounds the depth below 64, including median fallback levels.
// One pending sibling per level fits without allocation or shared scratch.
//
template<G4bool bounded>
G4double G4TessellatedSolid::NearestFacet(const G4ThreeVector& p, G4int& index) const
{
  struct Pending
  {
    G4int index;
    G4double bound;
  };
  std::array<Pending, 64> stack;
  G4int size = 0;
  stack[size++] = {0, 0.};
  G4double distance = kInfinity;
  index = -1;
  G4int visited = 0, tested = 0;
  G4bool budgetActive = bounded;
  while (size != 0)
  {
    if constexpr (bounded)
    {
      if (budgetActive && ++visited > 128)
      {
        index = -2;  // Classification may switch to a ray; distance is incomplete.
        return kInfinity;
      }
    }
    const auto pending = stack[--size];
    if (pending.bound >= distance * distance) { continue; }
    const auto current = pending.index;
    const auto& node = fFacetTree[current];
    if constexpr (bounded)
    {
      // A well-separated subtree already justifies ray classification;
      // descending to its first leaf would add no useful surface information.
      if (budgetActive && pending.bound > 0.
          && pending.bound > (node.upper - node.lower).mag2())
      {
        index = -2;
        return kInfinity;
      }
    }
    if (node.count != 0)
    {
      if constexpr (bounded)
      {
        if (tested == 0 && budgetActive)
        {
          // Near-surface probes favour the existing nearest-face search.
          // Switch algorithms only when the point is well separated from the
          // first leaf relative to that leaf's size.
          const auto extent = node.upper - node.lower;
          const auto diagonalSquared = extent.mag2();
          budgetActive = pending.bound > 1.e-4 * diagonalSquared;
        }
      }
      for (auto i = node.first; i < node.first + node.count; ++i)
      {
        if constexpr (bounded)
        {
          if (budgetActive && ++tested > 8)
          {
            index = -2;
            return kInfinity;
          }
        }
        const auto candidate = fFacetOrder[i];
        const auto d = FacetDistance(candidate, p, distance);
        if (d < distance || (d == distance && candidate < index))
        {
          index = candidate;
          distance = d;
          if (d <= kCarToleranceHalf) { return 0.; }
        }
      }
    }
    else
    {
      const auto left = current + 1, right = node.right;
      const auto dl = FacetBoxDistance(fFacetTree[left], p);
      const auto dr = FacetBoxDistance(fFacetTree[right], p);
      if (dl < dr)
      {
        if (dr < distance * distance) stack[size++] = {right, dr};
        if (dl < distance * distance) stack[size++] = {left, dl};
      }
      else
      {
        if (dl < distance * distance) stack[size++] = {left, dl};
        if (dr < distance * distance) stack[size++] = {right, dr};
      }
    }
  }
  return distance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::FacetDistance(G4int index, const G4ThreeVector& p,
                                           G4double limit) const
{
  const auto& data = fFacetPlanes[index];
  const auto distance = std::abs(data.normal.dot(p - data.anchor));
  if (distance > limit + 4. * kCarTolerance) { return kInfinity; }

  // Preserve the original kernels around the classification tolerance.
  // Strictly interior projections elsewhere need only a plane distance.
  //
  if ((distance < .25 * kCarTolerance || distance > 4. * kCarTolerance)
      && FacetProjectionInside(index, p))
  {
    return distance;
  }
  return fFacets[index]->Distance(p, limit);
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4TessellatedSolid::FacetProjectionInside(G4int index,
                                                 const G4ThreeVector& p) const
{
  const auto& data = fFacetPlanes[index];
  const auto relative = p - data.anchor;
  for (G4int i = 0; i < data.count; ++i)
  {
    if (data.inward[i].dot(relative) - data.offset[i] <= 4. * kCarTolerance)
    {
      return false;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//
template<G4bool classify>
G4double G4TessellatedSolid::IntersectFacetTree(const G4ThreeVector& p, const G4ThreeVector& v,
                                                G4bool outgoing, G4int& index,
                                                G4ThreeVector& normal) const
{
  std::array<G4int, 64> stack;

  // Classification visits both face orientations. Reuse child slab entries
  // when popping siblings and compute reciprocal ray components only once.
  //
  std::array<G4double, classify ? 64 : 1> entries;
  G4ThreeVector inverse;
  if constexpr (classify)
  {
    for (G4int axis = 0; axis < 3; ++axis)
    {
      inverse[axis] = v[axis] == 0. ? 0. : 1. / v[axis];
    }
  }
  const auto boxEntry = [&](const FacetTreeNode& node, G4double limit)
  {
    if constexpr (classify)
    {
      if (std::isfinite(inverse.x()) && std::isfinite(inverse.y()) && std::isfinite(inverse.z()))
      {
        return FacetBoxEntry<true>(node, p, inverse, limit);
      }
    }
    return FacetBoxEntry(node, p, v, limit);
  };
  if constexpr (classify)
  {
    entries[0] = boxEntry(fFacetTree[0], kInfinity);
  }
  G4int size = 0;
  stack[size++] = 0;
  G4double minimum = kInfinity;
  index = -1;
  G4bool ambiguous = false;
  while (size != 0)
  {
    const auto current = stack[--size];
    const auto& node = fFacetTree[current];
    const auto limit = classify ? minimum + kCarToleranceHalf : minimum;
    if constexpr (classify)
    {
      if (entries[size] == kInfinity || entries[size] > limit) { continue; }
    }
    else if (FacetBoxEntry(node, p, v, limit) == kInfinity) { continue; }
    if (node.count != 0)
    {
      for (auto i = node.first; i < node.first + node.count; ++i)
      {
        const auto candidate = fFacetOrder[i];
        auto& facet = *fFacets[candidate];
        if constexpr (classify)
        {
          // A surface-band point lies within the padded facet box, hence
          // its ray entry is zero. Check finite distance in the same traversal.
          //
          if (entries[size] == 0.
              && FacetDistance(candidate, p, kCarTolerance) <= kCarToleranceHalf)
          {
            index = -2;  // Verified surface point, independently of ray orientation.
            return 0.;
          }
          outgoing = fFacetPlanes[candidate].normal.dot(v) > 0.;
        }
        G4double distance, fromSurface;
        G4ThreeVector n;
        if (!facet.Intersect(p, v, outgoing, distance, fromSurface, n))
        {
          continue;
        }
        if (outgoing)
        {
          if (fromSurface > 0. && fromSurface <= kCarToleranceHalf
           && facet.Distance(p, kCarTolerance) <= kCarToleranceHalf)
          {
            distance = 0.;
          }
        }
        else if (std::abs(distance) <= kCarToleranceHalf)
        {
          distance = 0.;
        }
        if constexpr (classify)
        {
          if (distance >= 0. && std::abs(distance-minimum) <= kCarToleranceHalf
              && normal.dot(v) * n.dot(v) < 0.)
          {
            ambiguous = true;
          }
          if (distance >= 0. && distance < minimum - kCarToleranceHalf)
          {
            ambiguous = false;
          }
        }
        if (distance >= 0. && distance < minimum)
        {
          minimum = distance;
          index = candidate;
          normal = n;
        }
      }
    }
    else
    {
      const auto left = current + 1, right = node.right;
      const auto dl = boxEntry(fFacetTree[left], limit);
      const auto dr = boxEntry(fFacetTree[right], limit);
      const auto push = [&](G4int child, G4double entry)
      {
        if (entry == kInfinity) { return; }
        if constexpr (classify)
        {
          entries[size] = entry;
        }
        stack[size++] = child;
      };
      if (dl < dr)
      {
        push(right, dr);
        push(left, dl);
      }
      else
      {
        push(left, dl);
        push(right, dr);
      }
    }
  }
  if constexpr (classify)
  {
    if (ambiguous) { normal.set(0., 0., 0.); }
  }
  return minimum;
}

///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::InsideFacetTree(const G4ThreeVector& p) const
{
  if (OutsideOfExtent(p, kCarTolerance)) { return kOutside; }
  G4int index;
  auto distance = fFacets.size() >= 4096
                ? NearestFacet<true>(p, index) : NearestFacet(p, index);
  if (index == -2)
  {
    // A costly nearest search can be replaced by a surface-band test and the first crossing,
    // not an exact global nearest distance. Prefer the nearby outward surface;
    // if that ray misses, reverse it before falling back to nearest distance.
    // An axial ray makes two slab tests parallel. Prefer the largest centre
    // displacement; at the centre retain an oblique ray to avoid mesh seams.
    //
    const auto relative = p - (.5 * fMinExtent + .5 * fMaxExtent);
    G4int axis = 0;
    if (std::abs(relative[1]) > std::abs(relative[axis])) axis = 1;
    if (std::abs(relative[2]) > std::abs(relative[axis])) axis = 2;
    G4ThreeVector direction;
    direction[axis] = relative[axis] < 0. ? -1. : 1.;
    if (relative.mag2() == 0.) direction = fRandir.front();
    G4ThreeVector normal;
    G4int hit;
    auto crossingDistance = IntersectFacetTree<true>(p, direction, false, hit, normal);
    if (crossingDistance == kInfinity)
    {
      direction = -direction;
      crossingDistance = IntersectFacetTree<true>(p, direction, false, hit, normal);
    }
    if (hit == -2) { return kSurface; }
    if (crossingDistance > 4. * kCarTolerance && crossingDistance < kInfinity)
    {
      const auto crossing = normal.dot(direction);
      if (crossing < -1.e-14) { return kOutside; }
      if (crossing > 1.e-14) { return kInside; }
    }
    // Tangencies, tied crossings and misses retain nearest-face classification
    // and the existing alternate-ray fallback below.
    //
    distance = NearestFacet(p, index);
  }
  if (distance <= kCarToleranceHalf) { return kSurface; }
  if (index >= 0 && FacetProjectionInside(index, p))
  {
    const auto& data = fFacetPlanes[index];
    const auto signedDistance = (p - data.anchor).dot(data.normal);

    // A nearest point strictly within a face determines the local material
    // side, for concave surfaces, cavities and disconnected components too.
    //
    if (std::abs(signedDistance) > kCarToleranceHalf)
    {
      return signedDistance < 0. ? kInside : kOutside;
    }
  }

  // Edges/vertices do not have a unique face normal. Classify by the nearest
  // oriented ray crossing; try another direction if two crossings tie.
  //
  for (const auto& direction : fRandir)
  {
    G4ThreeVector n;
    G4int inIndex, outIndex;
    const auto entry = IntersectFacetTree(p, direction, false, inIndex, n);
    const auto exit = IntersectFacetTree(p, direction, true, outIndex, n);
    if (entry == kInfinity && exit == kInfinity) { return kOutside; }
    if (entry < exit - kCarToleranceHalf) { return kOutside; }
    if (exit < entry - kCarToleranceHalf) { return kInside; }
  }
  return InsideNoVoxels(p);
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::SetAllUsingStack(const std::vector<G4int>& voxel,
                                           const std::vector<G4int>& max, G4bool status,
                                           G4SurfBits& checked)
{
  vector<G4int> xyz = voxel;
  stack<vector<G4int>> pos;
  pos.push(xyz);
  G4int filled = 0;

  vector<G4int> candidates;

  while (!pos.empty())  // Loop checking, 13.08.2015, G.Cosmo
  {
    xyz = pos.top();
    pos.pop();
    G4int index = fVoxels.GetVoxelsIndex(xyz);
    if (!checked[index])
    {
      checked.SetBitNumber(index, true);

      if (fVoxels.IsEmpty(index))
      {
        ++filled;

        fInsides.SetBitNumber(index, status);

        for (auto i = 0; i <= 2; ++i)
        {
          if (xyz[i] < max[i] - 1)
          {
            xyz[i]++;
            pos.push(xyz);
            xyz[i]--;
          }

          if (xyz[i] > 0)
          {
            xyz[i]--;
            pos.push(xyz);
            xyz[i]++;
          }
        }
      }
    }
  }
  return filled;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::PrecalculateInsides()
{
  vector<G4int> voxel(3), maxVoxels(3);
  for (auto i = 0; i <= 2; ++i)
  {
    maxVoxels[i] = (G4int)fVoxels.GetBoundary(i).size();
  }
  G4int size = maxVoxels[0] * maxVoxels[1] * maxVoxels[2];

  G4SurfBits checked(size - 1);
  fInsides.Clear();
  fInsides.ResetBitNumber(size - 1);

  G4ThreeVector point;
  for (voxel[2] = 0; voxel[2] < maxVoxels[2] - 1; ++voxel[2])
  {
    for (voxel[1] = 0; voxel[1] < maxVoxels[1] - 1; ++voxel[1])
    {
      for (voxel[0] = 0; voxel[0] < maxVoxels[0] - 1; ++voxel[0])
      {
        G4int index = fVoxels.GetVoxelsIndex(voxel);
        if (!checked[index] && fVoxels.IsEmpty(index))
        {
          for (auto i = 0; i <= 2; ++i)
          {
            point[i] = fVoxels.GetBoundary(i)[voxel[i]];
          }
          auto inside = (G4bool)(InsideNoVoxels(point) == kInside);
          SetAllUsingStack(voxel, maxVoxels, inside, checked);
        }
        else
        {
          checked.SetBitNumber(index);
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::Voxelize()
{
#  ifdef G4SPECSDEBUG
  G4cout << "Voxelizing..." << G4endl;
#  endif

  // Large meshes use the facet tree for every query. Retain voxel entry
  // traversal on smaller meshes, where it is faster, but not unused safety
  // lists or precomputed inside classifications.
  //
  if (UsesFacetTreeForEntry())
  {
    fVoxels.Clear();
    fInsides.Clear();
    return;
  }
  fVoxels.Voxelize(fFacets, !HasFacetTree());

  if (!HasFacetTree() && fVoxels.Empty().GetNbits() != 0u)
  {
#  ifdef G4SPECSDEBUG
    G4cout << "Precalculating Insides..." << G4endl;
#  endif
    PrecalculateInsides();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// SetExtremeFacets
//
// Find facets whose supporting planes contain every vertex on the material
// side. Preserve the exact IsInside predicate, including its equality case.
// The facet hierarchy rejects whole groups of vertices away from the plane.
//
void G4TessellatedSolid::SetExtremeFacets()
{
  fExtremeFacets.clear();
  if (fVertexList.empty()) return;
  std::array<G4ThreeVector, 6> points;
  points.fill(fVertexList.front());
  for (const auto& vertex : fVertexList)
  {
    for (G4int axis = 0; axis < 3; ++axis)
    {
      if (vertex[axis] < points[2 * axis][axis]) points[2 * axis] = vertex;
      if (vertex[axis] > points[2 * axis + 1][axis]) points[2 * axis + 1] = vertex;
    }
  }
  const auto eligible = [&](const G4VFacet* facet)
  {
    // Strict IsInside can reject even a face's own vertices through roundoff.
    //
    for (G4int i = 0; i < facet->GetNumberOfVertices(); ++i)
    {
      if (!facet->IsInside(facet->GetVertex(i))) { return false; }
    }
    return std::all_of(points.begin(), points.end(),
                       [facet](const G4ThreeVector& p) { return facet->IsInside(p); });
  };
  // A box alone is a loose support bound for a curved patch. Temporarily
  // bound each large-tree node along its mean facet normal as well. Removing
  // that component of the query normal leaves only a small tangential box
  // bound. These construction-only data do not increase retained memory.
  struct Support
  {
    G4ThreeVector normal;
    G4double upper = 0.;
    G4int first = 0, last = 0;
  };
  std::vector<Support> support;
  constexpr auto roundoff = 64. * std::numeric_limits<G4double>::epsilon();
  if (HasFacetTree() && fFacets.size() >= 4096
      && static_cast<std::size_t>(std::count_if(fFacets.begin(), fFacets.end(), eligible))
           > fFacets.size() / 16)
  {
    support.resize(fFacetTree.size());
    for (G4int i = static_cast<G4int>(support.size()) - 1; i >= 0; --i)
    {
      const auto& node = fFacetTree[i];
      auto& bound = support[i];
      if (node.count)
      {
        bound.first = node.first;
        bound.last = node.first + node.count;
        for (auto j = bound.first; j < bound.last; ++j)
        {
          bound.normal += fFacets[fFacetOrder[j]]->GetSurfaceNormal();
        }
      }
      else
      {
        const auto& left = support[i + 1];
        const auto& right = support[node.right];
        bound.first = left.first;
        bound.last = right.last;
        bound.normal = (left.last - left.first) * left.normal
                     + (right.last - right.first) * right.normal;
      }
      bound.normal = bound.normal.unit();
      const auto centre = .5 * node.lower + .5 * node.upper;
      bound.upper = -kInfinity;
      for (auto j = bound.first; j < bound.last; ++j)
      {
        const auto& facet = *fFacets[fFacetOrder[j]];
        for (G4int k = 0; k < facet.GetNumberOfVertices(); ++k)
        {
          const auto vertex = facet.GetVertex(k);
          G4double error = 0.;
          for (G4int axis = 0; axis < 3; ++axis)
          {
            error += (std::abs(vertex[axis]) + std::abs(centre[axis]))
                   * std::abs(bound.normal[axis]);
          }
          bound.upper = std::max(bound.upper, (vertex - centre).dot(bound.normal)
                              + roundoff * error + std::numeric_limits<G4double>::min());
        }
      }
    }
  }
  for (auto* facet : fFacets)
  {
    const auto inside = [facet](const G4ThreeVector& vertex)
    {
      return facet->IsInside(vertex);
    };
    if (!eligible(facet)) { continue; }
    if (!HasFacetTree())
    {
      if (std::all_of(fVertexList.begin(), fVertexList.end(), inside))
      {
        fExtremeFacets.insert(facet);
      }
      continue;
    }
    const auto anchor = facet->GetVertex(0);
    const auto normal = facet->GetSurfaceNormal();
    std::array<G4int, 64> stack;
    G4int pending = 0;
    stack[pending++] = 0;
    G4bool extreme = true;
    while (pending && extreme)
    {
      const auto index = stack[--pending];
      const auto& node = fFacetTree[index];
      G4ThreeVector corner;
      G4double error = 0.;
      G4bool nonpositive = true;
      for (G4int axis = 0; axis < 3; ++axis)
      {
        corner[axis] = normal[axis] >= 0. ? node.upper[axis] : node.lower[axis];

        // Compare signs without multiplying: a positive subnormal product
        // must not round to zero and certify a box as behind the plane.
        //
        nonpositive = nonpositive && (normal[axis] == 0.
              || (normal[axis] > 0. ? corner[axis] <= anchor[axis]
                                    : corner[axis] >= anchor[axis]));
        error += (std::abs(corner[axis]) + std::abs(anchor[axis]))
                * std::abs(normal[axis]);
      }
      error = error * roundoff + std::numeric_limits<G4double>::min();

      // Only prune boxes strictly behind the plane. At the plane, use
      // precisely the original IsInside predicate, including roundoff.
      // If every contribution is nonpositive, their sum cannot round
      // positive. This also prunes exactly coplanar axis-aligned faces.
      //
      if (nonpositive || (corner - anchor).dot(normal) < -error) { continue; }
      if (!support.empty())
      {
        const auto& bound = support[index];
        const auto lambda = normal.dot(bound.normal);
        if (lambda > 0.)
        {
          const auto centre = .5 * node.lower + .5 * node.upper;
          const auto tangent = normal - lambda * bound.normal;
          G4double upper = (centre - anchor).dot(normal) + lambda * bound.upper;
          G4double uncertainty = std::abs(lambda * bound.upper);
          for (G4int axis = 0; axis < 3; ++axis)
          {
            const auto edge = tangent[axis] >= 0. ? node.upper[axis] : node.lower[axis];
            upper += tangent[axis] * (edge - centre[axis]);
            uncertainty += (std::abs(normal[axis]) + std::abs(lambda * bound.normal[axis]))
                         * (std::abs(edge) + std::abs(anchor[axis]) + 2. * std::abs(centre[axis]));
          }
          if (upper < -roundoff * uncertainty - std::numeric_limits<G4double>::min())
          {
            continue;
          }
        }
      }
      if (node.count)
      {
        for (G4int i = node.first; i < node.first + node.count && extreme; ++i)
        {
          const auto& other = *fFacets[fFacetOrder[i]];
          for (G4int j = 0; j < other.GetNumberOfVertices(); ++j)
          {
            if (!facet->IsInside(other.GetVertex(j)))
            {
              extreme = false;
              break;
            }
          }
        }
      }
      else
      {
        stack[pending++] = node.right;
        stack[pending++] = index + 1;
      }
    }
    if (extreme) { fExtremeFacets.insert(facet); }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::CreateVertexList()
{
  // The algorithm:
  // we will have additional vertexListSorted, where all the items will be
  // sorted by magnitude of vertice vector.
  // New candidate for fVertexList - we will determine the position fo first
  // item which would be within its magnitude - 0.5*kCarTolerance.
  // We will go trough until we will reach > +0.5 kCarTolerance.
  // Comparison (q-p).mag() < 0.5*kCarTolerance will be made.
  // They can be just stored in std::vector, with custom insertion based
  // on binary search.

  set<G4VertexInfo, G4VertexComparator> vertexListSorted;
  set<G4VertexInfo, G4VertexComparator>::iterator pos, it;
  auto begin = vertexListSorted.begin(), end = vertexListSorted.end();
  G4ThreeVector p;
  G4VertexInfo value;

  fVertexList.clear();
  std::size_t size = fFacets.size();

  G4double kCarTolerance24 = kCarTolerance * kCarTolerance / 4.0;
  G4double kCarTolerance3 = 3 * kCarTolerance;
  vector<G4int> newIndex(100);

  for (std::size_t k = 0; k < size; ++k)
  {
    G4VFacet& facet = *fFacets[k];
    G4int max = facet.GetNumberOfVertices();

    for (G4int i = 0; i < max; ++i)
    {
      p = facet.GetVertex(i);
      value.id = (G4int)fVertexList.size();
      value.mag2 = p.x() + p.y() + p.z();

      G4bool found = false;
      G4int id = 0;
      if (!OutsideOfExtent(p, kCarTolerance))
      {
        pos = vertexListSorted.lower_bound(value);
        it = pos;
        while (it != end)  // Loop checking, 13.08.2015, G.Cosmo
        {
          id = (*it).id;
          G4ThreeVector q = fVertexList[id];
          G4double dif = (q - p).mag2();
          found = (dif < kCarTolerance24);
          if (found)
          {
            break;
          }
          dif = q.x() + q.y() + q.z() - value.mag2;
          if (dif > kCarTolerance3)
          {
            break;
          }
          ++it;
        }

        if (!found && (fVertexList.size() > 1))
        {
          it = pos;
          while (it != begin)  // Loop checking, 13.08.2015, G.Cosmo
          {
            --it;
            id = (*it).id;
            G4ThreeVector q = fVertexList[id];
            G4double dif = (q - p).mag2();
            found = (dif < kCarTolerance24);
            if (found)
            {
              break;
            }
            dif = value.mag2 - (q.x() + q.y() + q.z());
            if (dif > kCarTolerance3)
            {
              break;
            }
          }
        }
      }

      if (!found)
      {
#  ifdef G4SPECSDEBUG
        G4cout << p.x() << ":" << p.y() << ":" << p.z() << G4endl;
        G4cout << "Adding new vertex #" << i << " of facet " << k << " id " << value.id << G4endl;
        G4cout << "===" << G4endl;
#  endif
        fVertexList.push_back(p);
        vertexListSorted.insert(value);
        begin = vertexListSorted.begin();
        end = vertexListSorted.end();
        newIndex[i] = value.id;
        //
        // Now update the maximum x, y and z limits of the volume.
        //
        if (value.id == 0)
        {
          fMinExtent = fMaxExtent = p;
        }
        else
        {
          if (p.x() > fMaxExtent.x())
          {
            fMaxExtent.setX(p.x());
          }
          else if (p.x() < fMinExtent.x())
          {
            fMinExtent.setX(p.x());
          }
          if (p.y() > fMaxExtent.y())
          {
            fMaxExtent.setY(p.y());
          }
          else if (p.y() < fMinExtent.y())
          {
            fMinExtent.setY(p.y());
          }
          if (p.z() > fMaxExtent.z())
          {
            fMaxExtent.setZ(p.z());
          }
          else if (p.z() < fMinExtent.z())
          {
            fMinExtent.setZ(p.z());
          }
        }
      }
      else
      {
#  ifdef G4SPECSDEBUG
        G4cout << p.x() << ":" << p.y() << ":" << p.z() << G4endl;
        G4cout << "Vertex #" << i << " of facet " << k << " found, redirecting to " << id << G4endl;
        G4cout << "===" << G4endl;
#  endif
        newIndex[i] = id;
      }
    }
    // only now it is possible to change vertices pointer
    //
    facet.SetVertices(&fVertexList);
    for (G4int i = 0; i < max; ++i)
    {
      facet.SetVertexIndex(i, newIndex[i]);
    }
  }
  vector<G4ThreeVector>(fVertexList).swap(fVertexList);

#  ifdef G4SPECSDEBUG
  G4double previousValue = 0.;
  for (auto res = vertexListSorted.cbegin();
            res != vertexListSorted.cend(); ++res)
  {
    G4int id = (*res).id;
    G4ThreeVector vec = fVertexList[id];
    G4double mvalue = vec.x() + vec.y() + vec.z();
    if (previousValue && (previousValue - 1e-9 > mvalue))
      G4cout << "Error in CreateVertexList: previousValue " << previousValue
             << " is smaller than mvalue " << mvalue << G4endl;
    previousValue = mvalue;
  }
#  endif
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DisplayAllocatedMemory()
{
  G4int without = AllocatedMemoryWithoutVoxels();
  G4int with = AllocatedMemory();
  G4double ratio = (G4double)with / without;
  G4cout << "G4TessellatedSolid - Allocated memory without voxel overhead "
         << without << "; with "
         << with << "; ratio: " << ratio << G4endl;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::SetSolidClosed(const G4bool t)
{
  if (t && fSolidClosed) { return; }
  if (!t && fSolidClosed)
  {
    // Detach indexed vertices before a subsequent closure rebuilds their pool.
    for (auto& facet : fFacets)
    {
      auto* copy = facet->GetClone();
      delete facet;
      facet = copy;
      const auto centre = facet->GetCircumcentre();
      fFacetList.insert(
        {static_cast<G4int>(fFacetList.size()), centre.x() + centre.y() + centre.z()});
    }
    fVertexList.clear();
    fExtremeFacets.clear();
  }
  if (t)
  {
#  ifdef G4SPECSDEBUG
    G4cout << "Creating vertex list..." << G4endl;
#  endif
    // The duplicate-facet index is needed only while accepting new facets.
    // Recreate it on reopen, rather than retaining it throughout navigation.
    //
    fFacetList.clear();
    CreateVertexList();

#  ifdef G4SPECSDEBUG
    G4cout << "Setting extreme facets..." << G4endl;
#  endif
    BuildFacetTree();
    SetExtremeFacets();

#  ifdef G4SPECSDEBUG
    G4cout << "Voxelizing..." << G4endl;
#  endif
    Voxelize();

#  ifdef G4SPECSDEBUG
    DisplayAllocatedMemory();
#  endif

#  ifdef G4SPECSDEBUG
    G4cout << "Checking Structure..." << G4endl;
#  endif
    G4int irep = CheckStructure();
    if (irep != 0)
    {
      if ((irep & 1) != 0)
      {
        std::ostringstream message;
        message << "Defects in solid: " << GetName()
                << " - negative cubic volume, please check orientation of facets!";
        G4Exception("G4TessellatedSolid::SetSolidClosed()",
                    "GeomSolids1001", JustWarning, message);
      }
      if ((irep & 2) != 0)
      {
        std::ostringstream message;
        message << "Defects in solid: " << GetName()
                << " - some facets have wrong orientation!";
        G4Exception("G4TessellatedSolid::SetSolidClosed()",
                    "GeomSolids1001", JustWarning, message);
      }
      if ((irep & 4) != 0)
      {
        std::ostringstream message;
        message << "Defects in solid: " << GetName()
                << " - there are holes in the surface!";
        G4Exception("G4TessellatedSolid::SetSolidClosed()",
                    "GeomSolids1001", JustWarning, message);
      }
    }
  }
  if (!t)
  {
    fVoxels.Clear();
    fInsides.Clear();
    fFacetPlanes.clear();
    fFacetTree.clear();
    fFacetOrder.clear();
  }
  fSolidClosed = t;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetSolidClosed
//
// Used to determine whether the solid is closed to adding further fFacets.
//
G4bool G4TessellatedSolid::GetSolidClosed() const
{
  return fSolidClosed;
}

///////////////////////////////////////////////////////////////////////////////
//
// CheckStructure
//
// Checks structure of the solid. Return value is a sum of the following
// defect indicators, if any (0 means no defects):
//   1 - cubic volume is negative, wrong orientation of facets
//   2 - some facets have wrong orientation
//   4 - holes in the surface
//
G4int G4TessellatedSolid::CheckStructure() const
{
  G4int nedge = 0;
  std::size_t nface = fFacets.size();

  // Calculate volume
  //
  G4double volume = 0.;
  for (std::size_t i = 0; i < nface; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    nedge += facet.GetNumberOfVertices();
    volume += facet.GetArea() * (facet.GetVertex(0).dot(facet.GetSurfaceNormal()));
  }
  auto ivolume = static_cast<G4int>(volume <= 0.);

  // Create sorted vector of edges
  //
  std::vector<int64_t> iedge(nedge);
  G4int kk = 0;
  for (std::size_t i = 0; i < nface; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    G4int nnode = facet.GetNumberOfVertices();
    for (G4int k = 0; k < nnode; ++k)
    {
      int64_t i1 = facet.GetVertexIndex((k == 0) ? nnode - 1 : k - 1);
      int64_t i2 = facet.GetVertexIndex(k);
      auto inverse = static_cast<int64_t>(i2 > i1);
      if (inverse != 0)
      {
        std::swap(i1, i2);
      }
      iedge[kk++] = i1 * 1000000000 + i2 * 2 + inverse;
    }
  }
  std::sort(iedge.begin(), iedge.end());

  // Check edges, correct structure should consist of paired edges
  // with different orientation
  //
  G4int iorder = 0;
  G4int ihole = 0;
  G4int i = 0;
  while (i < nedge - 1)
  {
    if (iedge[i + 1] - iedge[i] == 1)  // paired edges with different orientation
    {
      i += 2;
    }
    else if (iedge[i + 1] == iedge[i])  // paired edges with the same orientation
    {
      iorder = 2;
      i += 2;
    }
    else  // unpaired edge
    {
      ihole = 4;
      i++;
    }
  }
  return ivolume + iorder + ihole;
}

///////////////////////////////////////////////////////////////////////////////
//
// operator+=
//
// This operator allows the user to add two tessellated solids together, so
// that the solid on the left then includes all of the facets in the solid
// on the right.  Note that copies of the facets are generated, rather than
// using the original facet set of the solid on the right.
//
G4TessellatedSolid& G4TessellatedSolid::operator+=(const G4TessellatedSolid& right)
{
  G4int size = right.GetNumberOfFacets();
  for (G4int i = 0; i < size; ++i)
  {
    AddFacet(right.GetFacet(i)->GetClone());
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetNumberOfFacets
//
G4int G4TessellatedSolid::GetNumberOfFacets() const
{
  return (G4int)fFacets.size();
}

///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::InsideVoxels(const G4ThreeVector& p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, kCarTolerance))
  {
    return kOutside;
  }

  vector<G4int> startingVoxel(3);
  fVoxels.GetVoxel(startingVoxel, p);

  const G4double dirTolerance = 1.0E-14;

  const vector<G4int>& startingCandidates = fVoxels.GetCandidates(startingVoxel);
  std::size_t limit = startingCandidates.size();
  if (limit == 0 && (fInsides.GetNbits() != 0u))
  {
    G4int index = fVoxels.GetPointIndex(p);
    EInside location = fInsides[index] ? kInside : kOutside;
    return location;
  }

  G4double minDist = kInfinity;

  for (std::size_t i = 0; i < limit; ++i)
  {
    G4int candidate = startingCandidates[i];
    G4VFacet& facet = *fFacets[candidate];
    G4double dist = facet.Distance(p, minDist);
    if (dist < minDist)
    {
      minDist = dist;
    }
    if (dist <= kCarToleranceHalf)
    {
      return kSurface;
    }
  }

  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003. In essence,
  // we're trying to determine whether we're inside the volume by projecting
  // a few rays and determining if the first surface crossed is has a normal
  // vector between 0 to pi/2 (out-going) or pi/2 to pi (in-going).
  // We should also avoid rays which are nearly within the plane of the
  // tessellated surface, and therefore produce rays randomly.
  // For the moment, this is a bit over-engineered (belt-braces-and-ducttape).
  //
  G4double distOut = kInfinity;
  G4double distIn = kInfinity;
  G4double distO = 0.0;
  G4double distI = 0.0;
  G4double distFromSurfaceO = 0.0;
  G4double distFromSurfaceI = 0.0;
  G4ThreeVector normalO, normalI;
  G4bool crossingO = false;
  G4bool crossingI = false;
  EInside location = kOutside;
  G4int sm = 0;

  G4bool nearParallel = false;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    // We loop until we find direction where the vector is not nearly parallel
    // to the surface of any facet since this causes ambiguities.  The usual
    // case is that the angles should be sufficiently different, but there
    // are 20 random directions to select from - hopefully sufficient.
    //
    distOut = distIn = kInfinity;
    const G4ThreeVector& v = fRandir[sm];
    ++sm;
    //
    // This code could be voxelized by the same algorithm, which is used for
    // DistanceToOut(). We will traverse through fVoxels. we will call
    // intersect only for those, which would be candidates and was not
    // checked before.
    //
    G4ThreeVector currentPoint = p;
    G4ThreeVector direction = v.unit();
    vector<G4int> curVoxel(3);
    curVoxel = startingVoxel;
    G4double shiftBonus = kCarTolerance;

    G4bool crossed = false;
    G4bool started = true;

    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int>& candidates =
        started ? startingCandidates : fVoxels.GetCandidates(curVoxel);
      started = false;
      if (auto candidatesCount = (G4int)candidates.size())
      {
        for (G4int i = 0; i < candidatesCount; ++i)
        {
          G4int candidate = candidates[i];
          // bits.SetBitNumber(candidate);
          G4VFacet& facet = *fFacets[candidate];

          crossingO = facet.Intersect(p, v, true, distO, distFromSurfaceO, normalO);
          crossingI = facet.Intersect(p, v, false, distI, distFromSurfaceI, normalI);

          if (crossingO || crossingI)
          {
            crossed = true;

            nearParallel = (crossingO && std::fabs(normalO.dot(v)) < dirTolerance)
                           || (crossingI && std::fabs(normalI.dot(v)) < dirTolerance);
            if (!nearParallel)
            {
              if (crossingO && distO > 0.0 && distO < distOut)
              {
                distOut = distO;
              }
              if (crossingI && distI > 0.0 && distI < distIn)
              {
                distIn = distI;
              }
            }
            else
            {
              break;
            }
          }
        }
        if (nearParallel)
        {
          break;
        }
      }
      else
      {
        if (!crossed)
        {
          G4int index = fVoxels.GetVoxelsIndex(curVoxel);
          G4bool inside = fInsides[index];
          location = inside ? kInside : kOutside;
          return location;
        }
      }

      G4double shift = fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity)
      {
        break;
      }

      currentPoint += direction * (shift + shiftBonus);
    } while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));

  } while (nearParallel && sm != fMaxTries);
  //
  // Here we loop through the facets to find out if there is an intersection
  // between the ray and that facet.  The test if performed separately whether
  // the ray is entering the facet or exiting.
  //
#  ifdef G4VERBOSE
  if (sm == fMaxTries)
  {
    //
    // We've run out of random vector directions. If nTries is set sufficiently
    // low (nTries <= 0.5*maxTries) then this would indicate that there is
    // something wrong with geometry.
    //
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Cannot determine whether point is inside or outside volume!" << G4endl
            << "Solid name       = " << GetName() << G4endl
            << "Geometry Type    = " << fGeometryType << G4endl
            << "Number of facets = " << fFacets.size() << G4endl << "Position:" << G4endl << G4endl
            << "p.x() = " << p.x() / mm << " mm" << G4endl << "p.y() = " << p.y() / mm << " mm"
            << G4endl << "p.z() = " << p.z() / mm << " mm";
    message.precision(oldprc);
    G4Exception("G4TessellatedSolid::Inside()", "GeomSolids1002", JustWarning, message);
  }
#  endif

  // In the next if-then-elseif G4String the logic is as follows:
  // (1) You don't hit anything so cannot be inside volume, provided volume
  //     constructed correctly!
  // (2) Distance to inside (ie. nearest facet such that you enter facet) is
  //     shorter than distance to outside (nearest facet such that you exit
  //     facet) - on condition of safety distance - therefore we're outside.
  // (3) Distance to outside is shorter than distance to inside therefore
  //     we're inside.
  //
  if (distIn == kInfinity && distOut == kInfinity)
  {
    location = kOutside;
  }
  else if (distIn <= distOut - kCarToleranceHalf)
  {
    location = kOutside;
  }
  else if (distOut <= distIn - kCarToleranceHalf)
  {
    location = kInside;
  }

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::InsideNoVoxels(const G4ThreeVector& p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, kCarTolerance))
  {
    return kOutside;
  }

  const G4double dirTolerance = 1.0E-14;

  G4double minDist = kInfinity;
  //
  // Check if we are close to a surface
  //
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    G4double dist = facet.Distance(p, minDist);
    if (dist < minDist)
    {
      minDist = dist;
    }
    if (dist <= kCarToleranceHalf)
    {
      return kSurface;
    }
  }
  //
  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003. In essence, we're
  // trying to determine whether we're inside the volume by projecting a few
  // rays and determining if the first surface crossed is has a normal vector
  // between 0 to pi/2 (out-going) or pi/2 to pi (in-going). We should also
  // avoid rays which are nearly within the plane of the tessellated surface,
  // and therefore produce rays randomly. For the moment, this is a bit
  // over-engineered (belt-braces-and-ducttape).
  //
#  if G4SPECSDEBUG
  G4int nTry = 7;
#  else
  G4int nTry = 3;
#  endif
  G4double distOut = kInfinity;
  G4double distIn = kInfinity;
  G4double distO = 0.0;
  G4double distI = 0.0;
  G4double distFromSurfaceO = 0.0;
  G4double distFromSurfaceI = 0.0;
  G4ThreeVector normalO(0.0, 0.0, 0.0);
  G4ThreeVector normalI(0.0, 0.0, 0.0);
  G4bool crossingO = false;
  G4bool crossingI = false;
  EInside location = kOutside;
  EInside locationprime = kOutside;
  G4int sm = 0;

  for (G4int i = 0; i < nTry; ++i)
  {
    G4bool nearParallel = false;
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      //
      // We loop until we find direction where the vector is not nearly parallel
      // to the surface of any facet since this causes ambiguities.  The usual
      // case is that the angles should be sufficiently different, but there
      // are 20 random directions to select from - hopefully sufficient.
      //
      distOut = distIn = kInfinity;
      G4ThreeVector v = fRandir[sm];
      sm++;
      auto f = fFacets.cbegin();

      do  // Loop checking, 13.08.2015, G.Cosmo
      {
        //
        // Here we loop through the facets to find out if there is an
        // intersection between the ray and that facet. The test if performed
        // separately whether the ray is entering the facet or exiting.
        //
        crossingO = ((*f)->Intersect(p, v, true, distO, distFromSurfaceO, normalO));
        crossingI = ((*f)->Intersect(p, v, false, distI, distFromSurfaceI, normalI));
        if (crossingO || crossingI)
        {
          nearParallel = (crossingO && std::fabs(normalO.dot(v)) < dirTolerance)
                         || (crossingI && std::fabs(normalI.dot(v)) < dirTolerance);
          if (!nearParallel)
          {
            if (crossingO && distO > 0.0 && distO < distOut)
            {
              distOut = distO;
            }
            if (crossingI && distI > 0.0 && distI < distIn)
            {
              distIn = distI;
            }
          }
        }
      } while (!nearParallel && ++f != fFacets.cend());
    } while (nearParallel && sm != fMaxTries);

#  ifdef G4VERBOSE
    if (sm == fMaxTries)
    {
      //
      // We've run out of random vector directions. If nTries is set
      // sufficiently low (nTries <= 0.5*maxTries) then this would indicate
      // that there is something wrong with geometry.
      //
      std::ostringstream message;
      G4long oldprc = message.precision(16);
      message << "Cannot determine whether point is inside or outside volume!" << G4endl
              << "Solid name       = " << GetName() << G4endl
              << "Geometry Type    = " << fGeometryType << G4endl
              << "Number of facets = " << fFacets.size() << G4endl << "Position:" << G4endl
              << G4endl << "p.x() = " << p.x() / mm << " mm" << G4endl << "p.y() = " << p.y() / mm
              << " mm" << G4endl << "p.z() = " << p.z() / mm << " mm";
      message.precision(oldprc);
      G4Exception("G4TessellatedSolid::Inside()", "GeomSolids1002", JustWarning, message);
    }
#  endif
    //
    // In the next if-then-elseif G4String the logic is as follows:
    // (1) You don't hit anything so cannot be inside volume, provided volume
    //     constructed correctly!
    // (2) Distance to inside (ie. nearest facet such that you enter facet) is
    //     shorter than distance to outside (nearest facet such that you exit
    //     facet) - on condition of safety distance - therefore we're outside.
    // (3) Distance to outside is shorter than distance to inside therefore
    // we're inside.
    //
    if (distIn == kInfinity && distOut == kInfinity)
    {
      locationprime = kOutside;
    }
    else if (distIn <= distOut - kCarToleranceHalf)
    {
      locationprime = kOutside;
    }
    else if (distOut <= distIn - kCarToleranceHalf)
    {
      locationprime = kInside;
    }

    if (i == 0)
    {
      location = locationprime;
    }
  }

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetFacetIndex
//
// Return index of the facet closest to the point p, normally the point should
// be located on the surface. Return -1 if no facet selected.
//
G4int G4TessellatedSolid::GetFacetIndex(const G4ThreeVector& p) const
{
  if (HasFacetTree())
  {
    G4int index;
    NearestFacet(p, index);
    return index;
  }
  G4int index = -1;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    vector<G4int> curVoxel(3);
    fVoxels.GetVoxel(curVoxel, p);
    const vector<G4int>& candidates = fVoxels.GetCandidates(curVoxel);
    if (auto limit = (G4int)candidates.size())
    {
      G4double minDist = kInfinity;
      for (G4int i = 0; i < limit; ++i)
      {
        G4int candidate = candidates[i];
        G4VFacet& facet = *fFacets[candidate];
        G4double dist = facet.Distance(p, minDist);
        if (dist <= kCarToleranceHalf)
        {
          return index = candidate;
        }
        if (dist < minDist)
        {
          minDist = dist;
          index = candidate;
        }
      }
    }
  }
  else
  {
    G4double minDist = kInfinity;
    std::size_t size = fFacets.size();
    for (std::size_t i = 0; i < size; ++i)
    {
      G4VFacet& facet = *fFacets[i];
      G4double dist = facet.Distance(p, minDist);
      if (dist < minDist)
      {
        minDist = dist;
        index = (G4int)i;
      }
    }
  }
  return index;
}

///////////////////////////////////////////////////////////////////////////////
//
// Normal
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.
//
G4bool G4TessellatedSolid::Normal(const G4ThreeVector& p, G4ThreeVector& aNormal) const
{
  if (HasFacetTree())
  {
    const auto surface = SurfaceFacet(p);
    if (surface >= 0)
    {
      aNormal = fFacets[surface]->GetSurfaceNormal();
      return true;
    }

    G4int index;
    const auto distance = NearestFacet(p, index);
    if (index >= 0)
    {
      aNormal = fFacets[index]->GetSurfaceNormal();
      return distance <= kCarToleranceHalf;
    }
  }
  G4double minDist;
  G4VFacet* facet = nullptr;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    vector<G4int> curVoxel(3);
    fVoxels.GetVoxel(curVoxel, p);
    const vector<G4int>& candidates = fVoxels.GetCandidates(curVoxel);

    if (auto limit = (G4int)candidates.size())
    {
      minDist = kInfinity;
      for (G4int i = 0; i < limit; ++i)
      {
        G4int candidate = candidates[i];
        G4VFacet& fct = *fFacets[candidate];
        G4double dist = fct.Distance(p, minDist);
        if (dist < minDist)
        {
          minDist = dist;
        }
        if (dist <= kCarToleranceHalf)
        {
          aNormal = fct.GetSurfaceNormal();
          return true;
        }
      }
    }
  }
  minDist = MinDistanceFacet(p, facet);

  if (minDist != kInfinity)
  {
    if (facet != nullptr)
    {
      aNormal = facet->GetSurfaceNormal();
    }
    return minDist <= kCarToleranceHalf;
  }

#  ifdef G4VERBOSE
  std::ostringstream message;
  message << "Point p is not on surface !?" << G4endl
          << "          No facets found for point: " << p << " !" << G4endl
          << "          Returning approximated value for normal.";

  G4Exception("G4TessellatedSolid::SurfaceNormal(p)", "GeomSolids1002", JustWarning, message);
#  endif

  aNormal = (p.z() > 0 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(0, 0, -1));
  return false;
}

///////////////////////////////////////////////////////////////////////////////
//
// DistanceToInNoVoxels
//
// Return the distance along the normalised vector v to the shape,
// from the point at offset p. If there is no intersection, return
// kInfinity. The first intersection resulting from 'leaving' a
// surface/volume is discarded. Hence, this is tolerant of points on
// surface of shape.
//
G4double G4TessellatedSolid::DistanceToInNoVoxels(const G4ThreeVector& p, const G4ThreeVector& v,
                                                  G4double /*aPstep*/) const
{
  G4double minDist = kInfinity;
  G4double dist = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

#  if G4SPECSDEBUG
  if (Inside(p) == kInside)
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is already inside!?" << G4endl << "Position:" << G4endl << G4endl
            << "   p.x() = " << p.x() / mm << " mm" << G4endl << "   p.y() = " << p.y() / mm
            << " mm" << G4endl << "   p.z() = " << p.z() / mm << " mm" << G4endl
            << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc);
    G4Exception("G4TriangularFacet::DistanceToIn(p,v)", "GeomSolids1002", JustWarning, message);
  }
#  endif

  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    if (facet.Intersect(p, v, false, dist, distFromSurface, normal))
    {
      //
      // set minDist to the new distance to current facet if distFromSurface
      // is in positive direction and point is not at surface. If the point is
      // within 0.5*kCarTolerance of the surface, then force distance to be
      // zero and leave member function immediately (for efficiency), as
      // proposed by & credit to Akira Okumura.
      //
      if (distFromSurface > kCarToleranceHalf && dist >= 0.0 && dist < minDist)
      {
        minDist = dist;
      }
      else
      {
        if (-kCarToleranceHalf <= dist && dist <= kCarToleranceHalf)
        {
          return 0.0;
        }

        if (distFromSurface > -kCarToleranceHalf && distFromSurface < kCarToleranceHalf)
        {
          minDist = dist;
        }
      }
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToOutNoVoxels(const G4ThreeVector& p, const G4ThreeVector& v,
                                                   G4ThreeVector& aNormalVector, G4bool& aConvex,
                                                   G4double /*aPstep*/) const
{
  G4double minDist = kInfinity;
  G4double dist = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal, minNormal;

#  if G4SPECSDEBUG
  if (Inside(p) == kOutside)
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is already outside!?" << G4endl << "Position:" << G4endl << G4endl
            << "   p.x() = " << p.x() / mm << " mm" << G4endl << "   p.y() = " << p.y() / mm
            << " mm" << G4endl << "   p.z() = " << p.z() / mm << " mm" << G4endl
            << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc);
    G4Exception("G4TriangularFacet::DistanceToOut(p)", "GeomSolids1002", JustWarning, message);
  }
#  endif

  G4bool isExtreme = false;
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    if (facet.Intersect(p, v, true, dist, distFromSurface, normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= kCarToleranceHalf
          && facet.Distance(p, kCarTolerance) <= kCarToleranceHalf)
      {
        // We are on a surface. Return zero.
        aConvex = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
        // Normal(p, aNormalVector);
        // aNormalVector = facet.GetSurfaceNormal();
        aNormalVector = normal;
        return 0.0;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist = dist;
        minNormal = normal;
        isExtreme = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
      }
    }
  }
  if (minDist < kInfinity)
  {
    aNormalVector = minNormal;
    aConvex = isExtreme;
    return minDist;
  }

  // No intersection found
  aConvex = false;
  Normal(p, aNormalVector);
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DistanceToOutCandidates(const std::vector<G4int>& candidates,
                                                 const G4ThreeVector& aPoint,
                                                 const G4ThreeVector& direction, G4double& minDist,
                                                 G4ThreeVector& minNormal,
                                                 G4int& minCandidate) const
{
  auto candidatesCount = (G4int)candidates.size();
  G4double dist = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

  for (G4int i = 0; i < candidatesCount; ++i)
  {
    G4int candidate = candidates[i];
    G4VFacet& facet = *fFacets[candidate];
    if (facet.Intersect(aPoint, direction, true, dist, distFromSurface, normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= kCarToleranceHalf
          && facet.Distance(aPoint, kCarTolerance) <= kCarToleranceHalf)
      {
        // We are on a surface
        //
        minDist = 0.0;
        minNormal = normal;
        minCandidate = candidate;
        break;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist = dist;
        minNormal = normal;
        minCandidate = candidate;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToOutCore(const G4ThreeVector& aPoint,
                                               const G4ThreeVector& aDirection,
                                               G4ThreeVector& aNormalVector, G4bool& aConvex,
                                               G4double aPstep) const
{
  if (HasFacetTree() && fFacets.size() >= 64)
  {
    G4int index;
    const auto distance = IntersectFacetTree(aPoint, aDirection, true, index, aNormalVector);
    aConvex = index >= 0
           && fExtremeFacets.find(fFacets[index]) != fExtremeFacets.end();
    if (index >= 0) { return distance; }
    Normal(aPoint, aNormalVector);
    return 0.;
  }

  G4double minDistance;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = kInfinity;

    G4ThreeVector currentPoint = aPoint;
    G4ThreeVector direction = aDirection.unit();
    G4double totalShift = 0.;
    vector<G4int> curVoxel(3);
    if (!fVoxels.Contains(aPoint))
    {
      return 0.;
    }

    fVoxels.GetVoxel(curVoxel, currentPoint);

    G4double shiftBonus = kCarTolerance;

    const vector<G4int>* old = nullptr;

    G4int minCandidate = -1;
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int>& candidates = fVoxels.GetCandidates(curVoxel);
      if (old == &candidates)
      {
        ++old;
      }
      if (old != &candidates && !candidates.empty())
      {
        DistanceToOutCandidates(candidates, aPoint, direction, minDistance, aNormalVector,
                                minCandidate);
        if (minDistance <= totalShift)
        {
          break;
        }
      }

      G4double shift = fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity)
      {
        break;
      }

      totalShift += shift;
      if (minDistance <= totalShift)
      {
        break;
      }

      currentPoint += direction * (shift + shiftBonus);

      old = &candidates;
    } while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));

    if (minCandidate < 0)
    {
      // No intersection found
      minDistance = 0.;
      aConvex = false;
      Normal(aPoint, aNormalVector);
    }
    else
    {
      aConvex = (fExtremeFacets.find(fFacets[minCandidate]) != fExtremeFacets.end());
    }
  }
  else
  {
    minDistance = DistanceToOutNoVoxels(aPoint, aDirection, aNormalVector, aConvex, aPstep);
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToInCandidates(const std::vector<G4int>& candidates,
                                                    const G4ThreeVector& aPoint,
                                                    const G4ThreeVector& direction) const
{
  auto candidatesCount = (G4int)candidates.size();
  G4double dist = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

  G4double minDistance = kInfinity;
  for (G4int i = 0; i < candidatesCount; ++i)
  {
    G4int candidate = candidates[i];
    G4VFacet& facet = *fFacets[candidate];
    if (facet.Intersect(aPoint, direction, false, dist, distFromSurface, normal))
    {
      //
      // Set minDist to the new distance to current facet if distFromSurface is
      // in positive direction and point is not at surface. If the point is
      // within 0.5*kCarTolerance of the surface, then force distance to be
      // zero and leave member function immediately (for efficiency), as
      // proposed by & credit to Akira Okumura.
      //
      if ((distFromSurface > kCarToleranceHalf) && (dist >= 0.0) && (dist < minDistance))
      {
        minDistance = dist;
      }
      else
      {
        if (-kCarToleranceHalf <= dist && dist <= kCarToleranceHalf)
        {
          return 0.0;
        }
        if (distFromSurface > -kCarToleranceHalf && distFromSurface < kCarToleranceHalf)
        {
          minDistance = dist;
        }
      }
    }
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToInCore(const G4ThreeVector& aPoint,
                                              const G4ThreeVector& aDirection,
                                              G4double aPstep) const
{
  if (UsesFacetTreeForEntry())
  {
    G4int index;
    G4ThreeVector normal;
    return IntersectFacetTree(aPoint, aDirection, false, index, normal);
  }

  G4double minDistance;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = kInfinity;
    G4ThreeVector currentPoint = aPoint;
    G4ThreeVector direction = aDirection.unit();
    G4double shift = fVoxels.DistanceToFirst(currentPoint, direction);
    if (shift == kInfinity)
    {
      return shift;
    }
    G4double shiftBonus = kCarTolerance;
    if (shift != 0.0)
    {
      currentPoint += direction * (shift + shiftBonus);
    }
    G4double totalShift = shift;

    vector<G4int> curVoxel(3);

    fVoxels.GetVoxel(curVoxel, currentPoint);
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int>& candidates = fVoxels.GetCandidates(curVoxel);
      if (!candidates.empty())
      {
        G4double distance = DistanceToInCandidates(candidates, aPoint, direction);
        if (minDistance > distance)
        {
          minDistance = distance;
        }
        if (distance < totalShift)
        {
          break;
        }
      }

      shift = fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity)
      {
        break;
      }

      totalShift += shift;
      if (minDistance < totalShift)
      {
        break;
      }

      currentPoint += direction * (shift + shiftBonus);
    } while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));
  }
  else
  {
    minDistance = DistanceToInNoVoxels(aPoint, aDirection, aPstep);
  }

  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::MinDistanceFacet(const G4ThreeVector& p,
                                              G4VFacet*& minFacet) const
{
  if (HasFacetTree())
  {
    G4int index;
    const auto distance = NearestFacet(p, index);
    minFacet = index >= 0 ? fFacets[index] : nullptr;
    return distance;
  }

  G4double minDist = kInfinity;
  if (fVoxels.GetCountOfVoxels() <= 1)
  {
    for (auto* facet : fFacets)
    {
      const auto distance = facet->Distance(p, minDist);
      if (distance < minDist)
      {
        minDist = distance;
        minFacet = facet;
      }
    }
    return minDist;
  }

  G4int size = fVoxels.GetVoxelBoxesSize();
  vector<pair<G4int, G4double>> voxelsSorted(size);

  for (G4int i = 0; i < size; ++i)
  {
    const G4VoxelBox& voxelBox = fVoxels.GetVoxelBox(i);

    G4ThreeVector pointShifted = p - voxelBox.pos;
    G4double safety = fVoxels.MinDistanceToBox(pointShifted, voxelBox.hlen);
    voxelsSorted[i] = {i, safety};
  }

  std::sort(voxelsSorted.begin(), voxelsSorted.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

  for (G4int i = 0; i < size; ++i)
  {
    const pair<G4int, G4double>& inf = voxelsSorted[i];
    G4double dist = inf.second;
    if (dist > minDist)
    {
      break;
    }

    const vector<G4int>& candidates = fVoxels.GetVoxelBoxCandidates(inf.first);
    auto csize = (G4int)candidates.size();
    for (G4int j = 0; j < csize; ++j)
    {
      G4int candidate = candidates[j];
      G4VFacet& facet = *fFacets[candidate];
      dist = facet.Distance(p, minDist);
      if (dist < minDist)
      {
        minDist = dist;
        minFacet = &facet;
      }
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::SafetyFromOutside(const G4ThreeVector& p, G4bool) const
{
  if (HasFacetTree())
  {
    G4int index;
    return NearestFacet(p, index);
  }

#  if G4SPECSDEBUG
  if (Inside(p) == kInside)
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is already inside!?" << G4endl << "Position:" << G4endl << G4endl
            << "p.x() = " << p.x() / mm << " mm" << G4endl << "p.y() = " << p.y() / mm << " mm"
            << G4endl << "p.z() = " << p.z() / mm << " mm" << G4endl
            << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc);
    G4Exception("G4TriangularFacet::DistanceToIn(p)", "GeomSolids1002", JustWarning, message);
  }
#  endif

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    if (!OutsideOfExtent(p, kCarTolerance))
    {
      vector<G4int> startingVoxel(3);
      fVoxels.GetVoxel(startingVoxel, p);
      const vector<G4int>& candidates = fVoxels.GetCandidates(startingVoxel);
      if (candidates.empty() && (fInsides.GetNbits() != 0u))
      {
        G4int index = fVoxels.GetPointIndex(p);
        if (fInsides[index])
        {
          return 0.;
        }
      }
    }
  }
  G4VFacet* facet = nullptr;
  return MinDistanceFacet(p, facet);
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::SafetyFromInside(const G4ThreeVector& p, G4bool) const
{
  if (HasFacetTree())
  {
    if (OutsideOfExtent(p, kCarTolerance)) { return 0.; }
    G4int index;
    return NearestFacet(p, index);
  }

#if G4SPECSDEBUG
  if (Inside(p) == kOutside)
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is already outside!?" << G4endl << "Position:" << G4endl << G4endl
            << "p.x() = " << p.x() / mm << " mm" << G4endl << "p.y() = " << p.y() / mm << " mm"
            << G4endl << "p.z() = " << p.z() / mm << " mm" << G4endl
            << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc);
    G4Exception("G4TriangularFacet::DistanceToOut(p)", "GeomSolids1002", JustWarning, message);
  }
#endif

  if (OutsideOfExtent(p, kCarTolerance)) { return 0.; }
  G4VFacet* facet = nullptr;
  return MinDistanceFacet(p, facet);
}

///////////////////////////////////////////////////////////////////////////////
//
// GetEntityType
//
// Provide identification of the class of an object
//
G4GeometryType G4TessellatedSolid::GetEntityType() const
{
  return fGeometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
// IsFaceted
//
G4bool G4TessellatedSolid::IsFaceted() const
{
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream& G4TessellatedSolid::StreamInfo(std::ostream& os) const
{
  os << G4endl;
  os << "Solid name       = " << GetName() << G4endl;
  os << "Geometry Type    = " << fGeometryType << G4endl;
  os << "Number of facets = " << fFacets.size() << G4endl;

  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    os << "FACET #          = " << i + 1 << G4endl;
    G4VFacet& facet = *fFacets[i];
    facet.StreamInfo(os);
  }
  os << G4endl;

  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4TessellatedSolid::Clone() const
{
  return new G4TessellatedSolid(*this);
}

///////////////////////////////////////////////////////////////////////////////
//
// Inside
//
// This method must return:
//    * kOutside if the point at offset p is outside the shape
//      boundaries plus kCarTolerance/2,
//    * kSurface if the point is <= kCarTolerance/2 from a surface, or
//    * kInside otherwise.
//
EInside G4TessellatedSolid::Inside(const G4ThreeVector& aPoint) const
{
  if (HasFacetTree())
  {
    return InsideFacetTree(aPoint);
  }
  EInside location;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    location = InsideVoxels(aPoint);
  }
  else
  {
    location = InsideNoVoxels(aPoint);
  }
  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::SurfaceNormal(const G4ThreeVector& p) const
{
  G4ThreeVector n;
  Normal(p, n);
  return n;
}

///////////////////////////////////////////////////////////////////////////////
//
// DistanceToIn
//
// Calculate distance to nearest surface of shape from an outside point p. The
// distance can be an underestimate.
//
G4double G4TessellatedSolid::DistanceToIn(const G4ThreeVector& p) const
{
  return SafetyFromOutside(p, false);
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  G4double dist = DistanceToInCore(p, v, kInfinity);
#  ifdef G4SPECSDEBUG
  if (dist < kInfinity)
  {
    if (Inside(p + dist * v) != kSurface)
    {
      std::ostringstream message;
      message << "Invalid response from facet in solid '" << GetName() << "'," << G4endl
              << "at point: " << p << "and direction: " << v;
      G4Exception("G4TessellatedSolid::DistanceToIn(p,v)", "GeomSolids1002", JustWarning, message);
    }
  }
#  endif
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// DistanceToOut
//
// Calculate distance to nearest surface of shape from an inside
// point. The distance can be an underestimate.
//
G4double G4TessellatedSolid::DistanceToOut(const G4ThreeVector& p) const
{
  return SafetyFromInside(p, false);
}

///////////////////////////////////////////////////////////////////////////////
//
// DistanceToOut
//
// Return distance along the normalised vector v to the shape, from a
// point at an offset p inside or on the surface of the
// shape. Intersections with surfaces, when the point is not greater
// than kCarTolerance/2 from a surface, must be ignored.
//     If calcNorm is true, then it must also set validNorm to either
//     * true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//     * false, if the solid does not lie entirely behind or on the
//       exiting surface.
// If calcNorm is false, then validNorm and n are unused.
//
G4double G4TessellatedSolid::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                           const G4bool calcNorm, G4bool* validNorm,
                                           G4ThreeVector* norm) const
{
  G4ThreeVector n;
  G4bool valid;

  G4double dist = DistanceToOutCore(p, v, n, valid);
  if (calcNorm)
  {
    *norm = n;
    *validNorm = valid;
  }
#  ifdef G4SPECSDEBUG
  if (dist < kInfinity)
  {
    if (Inside(p + dist * v) != kSurface)
    {
      std::ostringstream message;
      message << "Invalid response from facet in solid '" << GetName() << "'," << G4endl
              << "at point: " << p << "and direction: " << v;
      G4Exception("G4TessellatedSolid::DistanceToOut(p,v,..)", "GeomSolids1002", JustWarning,
                  message);
    }
  }
#  endif
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

///////////////////////////////////////////////////////////////////////////////
//
G4Polyhedron* G4TessellatedSolid::CreatePolyhedron() const
{
  auto nVertices = (G4int)fVertexList.size();
  auto nFacets = (G4int)fFacets.size();
  auto polyhedron = new G4Polyhedron(nVertices, nFacets);
  for (auto i = 0; i < nVertices; ++i)
  {
    polyhedron->SetVertex(i + 1, fVertexList[i]);
  }

  for (auto i = 0; i < nFacets; ++i)
  {
    G4VFacet* facet = fFacets[i];
    G4int v[4] = {0};
    G4int n = facet->GetNumberOfVertices();
    if (n > 4)
    {
      n = 4;
    }
    for (auto j = 0; j < n; ++j)
    {
      v[j] = facet->GetVertexIndex(j) + 1;
    }
    polyhedron->SetFacet(i + 1, v[0], v[1], v[2], v[3]);
  }
  polyhedron->SetReferences();

  return polyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4TessellatedSolid::GetPolyhedron() const
{
  G4RecursiveAutoLock l(&GetPolyhedronMutex());
  if (fpPolyhedron == nullptr || fRebuildPolyhedron
      || fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation()
           != fpPolyhedron->GetNumberOfRotationSteps())
  {
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
  }
  return fpPolyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
// BoundingLimits
//
// Get bounding box
//
void G4TessellatedSolid::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin = fMinExtent;
  pMax = fMaxExtent;

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " << GetName() << " !"
            << "\npMin = " << pMin << "\npMax = " << pMax;
    G4Exception("G4TessellatedSolid::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// CalculateExtent
//
// Calculate extent under transform and specified limit
//
G4bool G4TessellatedSolid::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                           const G4AffineTransform& pTransform, G4double& pMin,
                                           G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin, bmax);
  G4BoundingEnvelope bbox(bmin, bmax);

  // Use simple bounding-box to help in the case of complex meshes
  //
  return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);

#  if 0
  // Precise extent computation (disabled by default for this shape)
  //
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return (pMin < pMax) ? true : false;
  }

  // The extent is calculated as cumulative extent of the pyramids
  // formed by facets and the center of the bounding box.
  //
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  G4ThreeVectorList base;
  G4ThreeVectorList apex(1);
  std::vector<const G4ThreeVectorList *> pyramid(2);
  pyramid[0] = &base;
  pyramid[1] = &apex;
  apex[0] = (bmin+bmax)*0.5;

  // main loop along facets
  pMin =  kInfinity;
  pMax = -kInfinity;
  for (G4int i=0; i<GetNumberOfFacets(); ++i)
  {
    G4VFacet* facet = GetFacet(i);
    if (std::abs((facet->GetSurfaceNormal()).dot(facet->GetVertex(0)-apex[0]))
        < kCarToleranceHalf) { continue; }

    G4int nv = facet->GetNumberOfVertices();
    base.resize(nv);
    for (G4int k=0; k<nv; ++k) { base[k] = facet->GetVertex(k); }

    G4double emin,emax;
    G4BoundingEnvelope benv(pyramid);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax))
    {
      continue;
    }
    if (emin < pMin) { pMin = emin; }
    if (emax > pMax) { pMax = emax; }
    if (eminlim > pMin && emaxlim < pMax) { break; } // max possible extent
  }
  return (pMin < pMax);
#  endif
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinXExtent() const
{
  return fMinExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxXExtent() const
{
  return fMaxExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinYExtent() const
{
  return fMinExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxYExtent() const
{
  return fMaxExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinZExtent() const
{
  return fMinExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxZExtent() const
{
  return fMaxExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
G4VisExtent G4TessellatedSolid::GetExtent() const
{
  return {fMinExtent.x(), fMaxExtent.x(), fMinExtent.y(),
          fMaxExtent.y(), fMinExtent.z(), fMaxExtent.z()};
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetCubicVolume()
{
  if (fCubicVolume != 0.)
  {
    return fCubicVolume;
  }

  // For explanation of the following algorithm see:
  // https://en.wikipedia.org/wiki/Polyhedron#Volume
  // http://wwwf.imperial.ac.uk/~rn/centroid.pdf

  G4AutoLock l(&tessMutex);
  if (fCubicVolume != 0.)
  {
    l.unlock();
    return fCubicVolume;
  }
  G4double volume = 0.;
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    G4double area = facet.GetArea();
    G4ThreeVector unit_normal = facet.GetSurfaceNormal();
    volume += area * (facet.GetVertex(0).dot(unit_normal));
  }
  fCubicVolume = volume / 3.;
  l.unlock();

  return fCubicVolume;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetSurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    return fSurfaceArea;
  }

  G4AutoLock l(&tessMutex);
  if (fSurfaceArea != 0.)
  {
    l.unlock();
    return fSurfaceArea;
  }
  G4double area = 0.;
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    area += facet.GetArea();
  }
  fSurfaceArea = area;
  l.unlock();

  return fSurfaceArea;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::GetPointOnSurface() const
{
  // Select randomly a facet proportionally to its area and return a
  // random point on it.
  G4double area = 0.;
  for (const auto& facet : fFacets)
  {
    area += facet->GetArea();
  }

  if (fFacets.empty() || area <= 0.)
  {
    return {};
  }

  G4double select = area * G4QuickRand();
  G4double cumulative = 0.;
  for (const auto& facet : fFacets)
  {
    cumulative += facet->GetArea();
    if (select < cumulative)
    {
      return facet->GetPointOnFace();
    }
  }
  return fFacets.back()->GetPointOnFace();
}

///////////////////////////////////////////////////////////////////////////////
//
// SetRandomVectorSet
//
// This is a set of predefined random vectors (if that isn't a contradition
// in terms!) used to generate rays from a user-defined point.  The member
// function Inside uses these to determine whether the point is inside or
// outside of the tessellated solid.  All vectors should be unit vectors.
//
void G4TessellatedSolid::SetRandomVectors()
{
  fRandir.resize(20);
  fRandir[0] = G4ThreeVector(-0.9577428892113370, 0.2732676269591740, 0.0897405271949221);
  fRandir[1] = G4ThreeVector(-0.8331264504940770, -0.5162067214954600, -0.1985722492445700);
  fRandir[2] = G4ThreeVector(-0.1516671651108820, 0.9666292616127460, 0.2064580868390110);
  fRandir[3] = G4ThreeVector(0.6570250350323190, -0.6944539025883300, 0.2933460081893360);
  fRandir[4] = G4ThreeVector(-0.4820456281280320, -0.6331060000098690, -0.6056474264406270);
  fRandir[5] = G4ThreeVector(0.7629032554236800, 0.1016854697539910, -0.6384658864065180);
  fRandir[6] = G4ThreeVector(0.7689540409061150, 0.5034929891988220, 0.3939600142169160);
  fRandir[7] = G4ThreeVector(0.5765188359255740, 0.5997271636278330, -0.5549354566343150);
  fRandir[8] = G4ThreeVector(0.6660632777862070, -0.6362809868288380, 0.3892379937580790);
  fRandir[9] = G4ThreeVector(0.3824415020414780, 0.6541792713761380, -0.6525243125110690);
  fRandir[10] = G4ThreeVector(-0.5107726564526760, 0.6020905056811610, 0.6136760679616570);
  fRandir[11] = G4ThreeVector(0.7459135439578050, 0.6618796061649330, 0.0743530220183488);
  fRandir[12] = G4ThreeVector(0.1536405855311580, 0.8117477913978260, -0.5634359711967240);
  fRandir[13] = G4ThreeVector(0.0744395301705579, -0.8707110101772920, -0.4861286795736560);
  fRandir[14] = G4ThreeVector(-0.1665874645185400, 0.6018553940549240, -0.7810369397872780);
  fRandir[15] = G4ThreeVector(0.7766902003633100, 0.6014617505959970, -0.1870724331097450);
  fRandir[16] = G4ThreeVector(-0.8710128685847430, -0.1434320216603030, -0.4698551243971010);
  fRandir[17] = G4ThreeVector(0.8901082092766820, -0.4388411398893870, 0.1229871120030100);
  fRandir[18] = G4ThreeVector(-0.6430417431544370, -0.3295938228697690, 0.6912779675984150);
  fRandir[19] = G4ThreeVector(0.6331124368380410, 0.6306211461665000, 0.4488714875425340);

  fMaxTries = 20;
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::AllocatedMemoryWithoutVoxels()
{
  G4int base = sizeof(*this);
  base += fVertexList.capacity() * sizeof(G4ThreeVector);
  base += fRandir.capacity() * sizeof(G4ThreeVector);
  base += fFacetTree.capacity() * sizeof(FacetTreeNode);
  base += fFacetOrder.capacity() * sizeof(G4int);
  base += fFacetPlanes.capacity() * sizeof(FacetPlaneData);

  std::size_t limit = fFacets.size();
  for (std::size_t i = 0; i < limit; ++i)
  {
    G4VFacet& facet = *fFacets[i];
    base += facet.AllocatedMemory();
  }

  base += fFacets.capacity() * sizeof(G4VFacet*);
  base += fFacetList.size() * sizeof(G4VertexInfo);
  base += fExtremeFacets.size() * sizeof(G4VFacet*);
  return base;
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::AllocatedMemory()
{
  G4int size = AllocatedMemoryWithoutVoxels();
  G4int sizeInsides = fInsides.GetNbytes();
  G4int sizeVoxels = fVoxels.AllocatedMemory();
  size += sizeInsides + sizeVoxels;
  return size;
}

#endif
