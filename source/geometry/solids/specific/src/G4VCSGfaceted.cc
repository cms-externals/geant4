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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4VCSGfaceted implementation; a virtual class of a CSG type shape
// that is built entirely out of G4VCSGface faces.
//
// Author: David C. Williams (UCSC), 1998
// --------------------------------------------------------------------

#include "G4VCSGfaceted.hh"

#include "G4AffineTransform.hh"
#include "G4AutoLock.hh"
#include "G4Polyhedron.hh"
#include "G4QuickRand.hh"
#include "G4SolidExtentList.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"
#include "G4VoxelLimits.hh"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
  G4Mutex vcsgMutex = G4MUTEX_INITIALIZER;
}

//
// Constructor
//
G4VCSGfaceted::G4VCSGfaceted(const G4String& name)
  : G4VSolid(name), fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.)
{}

//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VCSGfaceted::G4VCSGfaceted(__void__& a)
  : G4VSolid(a), fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.)
{}

//
// Destructor
//
G4VCSGfaceted::~G4VCSGfaceted()
{
  DeleteStuff();
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

//
// Copy constructor
//
G4VCSGfaceted::G4VCSGfaceted(const G4VCSGfaceted& source) : G4VSolid(source)
{
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  CopyStuff(source);
}

//
// Assignment operator
//
G4VCSGfaceted& G4VCSGfaceted::operator=(const G4VCSGfaceted& source)
{
  if (&source == this)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(source);

  // Copy data
  //
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  DeleteStuff();
  CopyStuff(source);

  return *this;
}

//
// CopyStuff (protected)
//
// Copy the contents of source
//
void G4VCSGfaceted::CopyStuff(const G4VCSGfaceted& source)
{
  numFace = source.numFace;
  if (numFace == 0)
  {
    return;
  }  // odd, but permissable?

  faces = new G4VCSGface*[numFace];

  G4VCSGface **face = faces, **sourceFace = source.faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    *face = (*sourceFace)->Clone();
  } while (++sourceFace, ++face < faces + numFace);
  fCubicVolume = source.fCubicVolume;
  fSurfaceArea = source.fSurfaceArea;
  fRebuildPolyhedron = false;
  fpPolyhedron = nullptr;
  fZSectionMin = source.fZSectionMin;
  fZSectionMax = source.fZSectionMax;
  fRSectionMax = source.fRSectionMax;
  fZSlabBoundary = source.fZSlabBoundary;
  fZSlabCandidate = source.fZSlabCandidate;
  fZGlobalCandidate = source.fZGlobalCandidate;
}

//
// DeleteStuff (protected)
//
// Delete all allocated objects
//
void G4VCSGfaceted::DeleteStuff()
{
  if (numFace != 0)
  {
    G4VCSGface** face = faces;
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      delete *face;
    } while (++face < faces + numFace);

    delete[] faces;
  }
  numFace = 0;
  faces = nullptr;
  fZSectionMin.clear();
  fZSectionMax.clear();
  fRSectionMax.clear();
  fZSlabBoundary.clear();
  fZSlabCandidate.clear();
  fZGlobalCandidate.clear();
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

//
// BuildZSectionIndex
//
void G4VCSGfaceted::BuildZSectionIndex()
{
  fZSectionMin.clear();
  fZSectionMax.clear();
  fRSectionMax.clear();
  fZSlabBoundary.clear();
  fZSlabCandidate.clear();
  fZGlobalCandidate.clear();

  fZSectionMin.reserve(numFace);
  fZSectionMax.reserve(numFace);
  fRSectionMax.reserve(numFace);

  const G4ThreeVector zAxis(0., 0., 1.);
  const G4ThreeVector minusZAxis(0., 0., -1.);
  const G4ThreeVector xAxis(1., 0., 0.);
  const G4ThreeVector minusXAxis(-1., 0., 0.);
  const G4ThreeVector yAxis(0., 1., 0.);
  const G4ThreeVector minusYAxis(0., -1., 0.);

  G4VCSGface** face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    const G4double zMax = (*face)->Extent(zAxis);
    const G4double zMin = -(*face)->Extent(minusZAxis);
    fZSectionMin.push_back(std::min(zMin, zMax));
    fZSectionMax.push_back(std::max(zMin, zMax));
    const G4double xAbs = std::max((*face)->Extent(xAxis),
                                  (*face)->Extent(minusXAxis));
    const G4double yAbs = std::max((*face)->Extent(yAxis),
                                  (*face)->Extent(minusYAxis));
    fRSectionMax.push_back(std::sqrt(xAbs * xAbs + yAbs * yAbs));
    fZSlabBoundary.push_back(fZSectionMin.back());
    fZSlabBoundary.push_back(fZSectionMax.back());
  } while (++face < faces + numFace);

  if (numFace < kMinZSectionIndexFaces)
  {
    // Retain only face extents; a slab index would cost more than a short scan.
    fZSlabBoundary.clear();
    return;
  }

  std::sort(fZSlabBoundary.begin(), fZSlabBoundary.end());
  fZSlabBoundary.erase(std::unique(fZSlabBoundary.begin(), fZSlabBoundary.end()),
                       fZSlabBoundary.end());

  if (fZSlabBoundary.size() < 2)
  {
    fZSlabBoundary.clear();
    return;
  }

  fZSlabCandidate.resize(fZSlabBoundary.size() - 1);
  for (G4int i = 0; i < numFace; ++i)
  {
    std::vector<std::size_t> overlappingSlabs;
    overlappingSlabs.reserve(4);
    for (std::size_t slab = 0; slab < fZSlabCandidate.size(); ++slab)
    {
      const G4double zMin = fZSlabBoundary[slab];
      const G4double zMax = fZSlabBoundary[slab + 1];
      if (fZSectionMax[i] >= zMin && fZSectionMin[i] <= zMax)
      {
        overlappingSlabs.push_back(slab);
      }
    }
    if (overlappingSlabs.size() > 8)
    {
      // Broad faces are cheaper to visit once than to duplicate into many slabs.
      fZGlobalCandidate.push_back(i);
      continue;
    }
    for (auto slab : overlappingSlabs)
    {
      fZSlabCandidate[slab].push_back(i);
    }
  }
}

//
// FindZSlabIndex
//
std::size_t G4VCSGfaceted::FindZSlabIndex(const std::vector<G4double>& boundaries,
                                          G4double z)
{
  auto upper = std::upper_bound(boundaries.begin(), boundaries.end(), z);
  if (upper == boundaries.begin())
  {
    return 0;
  }
  if (upper == boundaries.end())
  {
    return boundaries.size() - 2;
  }
  return std::size_t((upper - boundaries.begin()) - 1);
}

//
// FindZSlabIndex
//
std::size_t G4VCSGfaceted::FindZSlabIndex(const std::vector<G4double>& boundaries,
                                          G4double z, G4double edgeTolerance)
{
  if (boundaries.size() < 2)
  {
    return 0;
  }
  if (z <= boundaries.front())
  {
    z = boundaries.front() + edgeTolerance;
  }
  else if (z >= boundaries.back())
  {
    z = boundaries.back() - edgeTolerance;
  }
  return FindZSlabIndex(boundaries, z);
}

//
// BuildZVoxelizer
//
void G4VCSGfaceted::BuildZVoxelizer(const std::vector<G4double>& sectionZMin,
                                    const std::vector<G4double>& sectionZMax,
                                    G4double minHalfZ, std::vector<G4double>& slabBoundary,
                                    std::vector<std::vector<G4int>>& slabCandidate)
{
  slabBoundary.clear();
  slabCandidate.clear();
  if (sectionZMin.empty() || sectionZMin.size() != sectionZMax.size())
  {
    return;
  }

  std::vector<G4double> zMin(sectionZMin.size());
  std::vector<G4double> zMax(sectionZMin.size());
  slabBoundary.reserve(2 * sectionZMin.size());
  for (std::size_t i = 0; i < sectionZMin.size(); ++i)
  {
    const G4double halfZ = std::max(0.5 * (sectionZMax[i] - sectionZMin[i]), minHalfZ);
    const G4double zc = 0.5 * (sectionZMin[i] + sectionZMax[i]);
    zMin[i] = zc - halfZ;
    zMax[i] = zc + halfZ;
    slabBoundary.push_back(zMin[i]);
    slabBoundary.push_back(zMax[i]);
  }

  std::sort(slabBoundary.begin(), slabBoundary.end());
  slabBoundary.erase(std::unique(slabBoundary.begin(), slabBoundary.end()),
                     slabBoundary.end());
  if (slabBoundary.size() < 2)
  {
    slabBoundary.clear();
    return;
  }

  slabCandidate.resize(slabBoundary.size() - 1);
  for (std::size_t slab = 0; slab < slabCandidate.size(); ++slab)
  {
    const G4double lo = slabBoundary[slab];
    const G4double hi = slabBoundary[slab + 1];
    for (std::size_t i = 0; i < zMin.size(); ++i)
    {
      if (zMax[i] >= lo && zMin[i] <= hi)
      {
        slabCandidate[slab].push_back(G4int(i));
      }
    }
  }
}

//
// FindZCandidates
//
const std::vector<G4int>& G4VCSGfaceted::FindZCandidates(
  const std::vector<G4double>& boundaries,
  const std::vector<std::vector<G4int>>& candidates,
  G4double z) const
{
  static const std::vector<G4int> empty;
  if (candidates.empty() || boundaries.size() < 2)
  {
    return empty;
  }
  return candidates[FindZSlabIndex(boundaries, z, 0.25 * kCarTolerance)];
}

//
// GetZCandidates
//
const std::vector<G4int>& G4VCSGfaceted::GetZCandidates(G4double z) const
{
  if (fZSlabCandidate.empty())
  {
    static const std::vector<G4int> empty;
    return empty;
  }

  return fZSlabCandidate[GetZSlabIndex(z)];
}

//
// CalculateExtent
//
G4bool G4VCSGfaceted::CalculateExtent(const EAxis axis, const G4VoxelLimits& voxelLimit,
                                      const G4AffineTransform& transform, G4double& min,
                                      G4double& max) const
{
  G4SolidExtentList extentList(axis, voxelLimit);

  //
  // Loop over all faces, checking min/max extent as we go.
  //
  G4VCSGface** face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    (*face)->CalculateExtent(axis, voxelLimit, transform, extentList);
  } while (++face < faces + numFace);

  //
  // Return min/max value
  //
  return extentList.GetExtent(min, max);
}

//
// Inside
//
// It could be a good idea to override this virtual
// member to add first a simple test (such as spherical
// test or whatnot) and to call this version only if
// the simplier test fails.
//
EInside G4VCSGfaceted::Inside(const G4ThreeVector& p) const
{
  return HasZSectionIndex() ? InsideZ(p) : InsideNoZ(p);
}

//
// SurfaceNormal
//
G4ThreeVector G4VCSGfaceted::SurfaceNormal(const G4ThreeVector& p) const
{
  G4ThreeVector answer;
  G4VCSGface** face = faces;
  G4double best = kInfinity;

  if (numFace < kMinZSectionIndexFaces)
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double distance = kInfinity;
      G4ThreeVector normal = (*face)->Normal(p, &distance);
      if (distance < best)
      {
        best = distance;
        answer = normal;
      }
    } while (++face < faces + numFace);

    return answer;
  }

  const G4bool useZSections = HasZSectionIndex();

  if (useZSections)
  {
    const G4double z = p.z();
    const auto& candidates = GetZCandidates(z);
    if (candidates.empty())
    {
      return answer;
    }
    for (auto i : fZGlobalCandidate)
    {
      G4double distance = kInfinity;
      G4ThreeVector normal = faces[i]->Normal(p, &distance);
      if (distance < best)
      {
        best = distance;
        answer = normal;
      }
    }
    for (auto i : candidates)
    {
      G4double distance = kInfinity;
      G4ThreeVector normal = faces[i]->Normal(p, &distance);
      if (distance < best)
      {
        best = distance;
        answer = normal;
      }
    }
    return answer;
  }

  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double distance = kInfinity;
    G4ThreeVector normal = (*face)->Normal(p, &distance);
    if (distance < best)
    {
      best = distance;
      answer = normal;
    }
  } while (++face < faces + numFace);

  return answer;
}

//
// DistanceToIn(p,v)
//
G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector& p,
                                     const G4ThreeVector& v) const
{
  G4double distance = kInfinity;
  G4double distFromSurface = kInfinity;
  G4VCSGface** face = faces;
  G4VCSGface* bestFace = *face;
  const G4double tolerance = kCarTolerance / 2;

  if (numFace < kMinZSectionIndexFaces)
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if ((*face)->Intersect(p, v, false, tolerance, faceDistance, faceDistFromSurface, faceNormal,
                             faceAllBehind))
      {
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          bestFace = *face;
          if (distFromSurface <= 0)
          {
            return 0;
          }
        }
      }
    } while (++face < faces + numFace);

    if (distance < kInfinity && distFromSurface < tolerance)
    {
      if (bestFace->Distance(p, false) < tolerance)
      {
        distance = 0;
      }
    }

    return distance;
  }

  if (HasZSectionIndex())
  {
    // Walk slabs in ray order and stop once they cannot improve the best hit.
    G4StackVisited visited(numFace);
    auto checkFace = [&](G4int i) -> G4bool
    {
      if (visited.TestAndSet(i))
      {
        return false;
      }
      if (!RayIntersectsFaceZ(i, p, v, distance, tolerance))
      {
        return false;
      }

      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if (faces[i]->Intersect(p, v, false, tolerance, faceDistance,
                              faceDistFromSurface, faceNormal, faceAllBehind))
      {
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          bestFace = faces[i];
          return distFromSurface <= 0;
        }
      }
      return false;
    };

    for (auto i : fZGlobalCandidate)
    {
      if (checkFace(i))
      {
        return 0;
      }
    }

    const G4double z = p.z();
    const G4double vz = v.z();
    if (std::fabs(vz) < DBL_MIN)
    {
      for (auto i : GetZCandidates(z))
      {
        if (checkFace(i))
        {
          return 0;
        }
      }
    }
    else if (!fZSlabCandidate.empty())
    {
      const G4double zFirst = fZSlabBoundary.front();
      const G4double zLast = fZSlabBoundary.back();
      if (!((vz > 0. && z > zLast + tolerance) || (vz < 0. && z < zFirst - tolerance)))
      {
        std::size_t slab = (vz > 0.)
                             ? (z < zFirst ? 0 : GetZSlabIndex(z))
                             : (z > zLast ? fZSlabCandidate.size() - 1 : GetZSlabIndex(z));
        while (true)
        {
          const G4double zEnter = (vz > 0.) ? fZSlabBoundary[slab] : fZSlabBoundary[slab + 1];
          const G4double rayEnter = (vz > 0.)
                                      ? (z < zEnter ? (zEnter - z) / vz : 0.)
                                      : (z > zEnter ? (zEnter - z) / vz : 0.);
          if (rayEnter > distance + tolerance)
          {
            break;
          }
          for (auto i : fZSlabCandidate[slab])
          {
            if (checkFace(i))
            {
              return 0;
            }
          }
          if (vz > 0.)
          {
            if (++slab >= fZSlabCandidate.size())
            {
              break;
            }
          }
          else if (slab-- == 0)
          {
            break;
          }
        }
      }
    }
  }
  else
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if ((*face)->Intersect(p, v, false, tolerance, faceDistance,
                             faceDistFromSurface, faceNormal, faceAllBehind))
      {
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          bestFace = *face;
          if (distFromSurface <= 0)
          {
            return 0;
          }
        }
      }
    } while (++face < faces + numFace);
  }

  if (distance < kInfinity && distFromSurface < tolerance)
  {
    if (bestFace->Distance(p, false) < tolerance)
    {
      distance = 0;
    }
  }

  return distance;
}

//
// DistanceToIn(p)
//
G4double G4VCSGfaceted::DistanceToIn(const G4ThreeVector& p) const
{
  return DistanceTo(p, false);
}

//
// DistanceToOut(p,v)
//
G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                      const G4bool calcNorm, G4bool* validNorm,
                                      G4ThreeVector* n) const
{
  G4bool allBehind = true;
  G4double distance = kInfinity;
  G4double distFromSurface = kInfinity;
  G4ThreeVector normal;

  G4VCSGface** face = faces;
  G4VCSGface* bestFace = *face;
  const G4double tolerance = kCarTolerance / 2;

  if (numFace < kMinZSectionIndexFaces)
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if ((*face)->Intersect(p, v, true, tolerance, faceDistance,
                             faceDistFromSurface, faceNormal, faceAllBehind))
      {
        if ((distance < kInfinity) || (!faceAllBehind))
        {
          allBehind = false;
        }
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          normal = faceNormal;
          bestFace = *face;
          if (distFromSurface <= 0.)
          {
            break;
          }
        }
      }
    } while (++face < faces + numFace);

    if (distance < kInfinity)
    {
      if (distFromSurface <= 0.)
      {
        distance = 0.;
      }
      else if (distFromSurface < tolerance)
      {
        if (bestFace->Distance(p, true) < tolerance)
        {
          distance = 0.;
        }
      }

      if (calcNorm)
      {
        *validNorm = allBehind;
        *n = normal;
      }
    }
    else
    {
      if (Inside(p) == kSurface)
      {
        distance = 0.;
      }
      if (calcNorm)
      {
        *validNorm = false;
      }
    }

    return distance;
  }

  if (HasZSectionIndex())
  {
    if (calcNorm)
    {
      allBehind = false;
    }

    // Walk slabs in ray order while suppressing faces shared by adjacent slabs.
    G4StackVisited visited(numFace);
    auto checkFace = [&](G4int i) -> G4bool
    {
      if (visited.TestAndSet(i))
      {
        return false;
      }
      if (!RayIntersectsFaceZ(i, p, v, distance, tolerance, false))
      {
        return false;
      }

      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if (faces[i]->Intersect(p, v, true, tolerance, faceDistance,
                              faceDistFromSurface, faceNormal, faceAllBehind))
      {
        if ((distance < kInfinity) || (!faceAllBehind))
        {
          allBehind = false;
        }
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          normal = faceNormal;
          bestFace = faces[i];
          return distFromSurface <= 0.;
        }
      }
      return false;
    };

    for (auto i : fZGlobalCandidate)
    {
      if (checkFace(i))
      {
        break;
      }
    }

    const G4double z = p.z();
    const G4double vz = v.z();
    if (distFromSurface > 0.)
    {
      if (std::fabs(vz) < DBL_MIN)
      {
        for (auto i : GetZCandidates(z))
        {
          if (checkFace(i))
          {
            break;
          }
        }
      }
      else if (!fZSlabCandidate.empty())
      {
        const G4double zFirst = fZSlabBoundary.front();
        const G4double zLast = fZSlabBoundary.back();
        if (!((vz > 0. && z > zLast + tolerance)
           || (vz < 0. && z < zFirst - tolerance)))
        {
          std::size_t slab = (vz > 0.)
                               ? (z < zFirst ? 0 : GetZSlabIndex(z))
                               : (z > zLast ? fZSlabCandidate.size() - 1 : GetZSlabIndex(z));
          while (true)
          {
            const G4double zEnter = (vz > 0.) ? fZSlabBoundary[slab] : fZSlabBoundary[slab + 1];
            const G4double rayEnter = (vz > 0.)
                                        ? (z < zEnter ? (zEnter - z) / vz : 0.)
                                        : (z > zEnter ? (zEnter - z) / vz : 0.);
            if (rayEnter > distance + tolerance)
            {
              break;
            }
            for (auto i : fZSlabCandidate[slab])
            {
              if (checkFace(i))
              {
                break;
              }
            }
            if (distFromSurface <= 0.)
            {
              break;
            }
            if (vz > 0.)
            {
              if (++slab >= fZSlabCandidate.size())
              {
                break;
              }
            }
            else if (slab-- == 0)
            {
              break;
            }
          }
        }
      }
    }
  }
  else
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double faceDistance, faceDistFromSurface;
      G4ThreeVector faceNormal;
      G4bool faceAllBehind;
      if ((*face)->Intersect(p, v, true, tolerance, faceDistance,
                             faceDistFromSurface, faceNormal, faceAllBehind))
      {
        if ((distance < kInfinity) || (!faceAllBehind))
        {
          allBehind = false;
        }
        if (faceDistance < distance)
        {
          distance = faceDistance;
          distFromSurface = faceDistFromSurface;
          normal = faceNormal;
          bestFace = *face;
          if (distFromSurface <= 0.)
          {
            break;
          }
        }
      }
    } while (++face < faces + numFace);
  }

  if (distance < kInfinity)
  {
    if (distFromSurface <= 0.)
    {
      distance = 0.;
    }
    else if (distFromSurface < tolerance)
    {
      if (bestFace->Distance(p, true) < tolerance)
      {
        distance = 0.;
      }
    }

    if (calcNorm)
    {
      *validNorm = allBehind;
      *n = normal;
    }
  }
  else
  {
    if (Inside(p) == kSurface)
    {
      distance = 0.;
    }
    if (calcNorm)
    {
      *validNorm = false;
    }
  }

  return distance;
}

//
// DistanceToOut(p)
//
G4double G4VCSGfaceted::DistanceToOut(const G4ThreeVector& p) const
{
  return DistanceTo(p, true);
}

//
// DistanceTo
//
// Protected routine called by DistanceToIn and DistanceToOut
//
G4double G4VCSGfaceted::DistanceTo(const G4ThreeVector& p, const G4bool outgoing) const
{
  G4VCSGface** face = faces;
  G4double best = kInfinity;
  const G4double tolerance = kCarTolerance / 2;

  if (numFace < kMinZSectionIndexFaces)
  {
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double distance = (*face)->Distance(p, outgoing);
      if (distance < best)
      {
        best = distance;
      }
    } while (++face < faces + numFace);

    return (best < 0.5 * kCarTolerance) ? 0. : best;
  }

  const G4bool useZSections = HasZSectionIndex();

  if (useZSections)
  {
    const G4double z = p.z();

    // Visit slabs by increasing Z separation until none can improve the safety.
    G4StackVisited visited(numFace);
    auto checkFace = [&](G4int i)
    {
      if (visited.TestAndSet(i))
      {
        return;
      }
      if (!FaceCanImprovePoint(i, z, best, tolerance))
      {
        return;
      }
      const G4double distance = faces[i]->Distance(p, outgoing);
      if (distance < best)
      {
        best = distance;
      }
    };

    for (auto i : fZGlobalCandidate)
    {
      checkFace(i);
      if (best <= 0.)
      {
        return 0.;
      }
    }

    VisitZSlabsBySafety(fZSlabBoundary, z, best, tolerance,
      [&](std::size_t slab)
      {
        for (auto i : fZSlabCandidate[slab])
        {
          checkFace(i);
        }
        return best <= 0.;
      });
    if (best <= 0.)
    {
      return 0.;
    }

    return (best < 0.5 * kCarTolerance) ? 0. : best;
  }

  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double distance = (*face)->Distance(p, outgoing);
    if (distance < best)
    {
      best = distance;
    }
  } while (++face < faces + numFace);

  return (best < 0.5 * kCarTolerance) ? 0. : best;
}

//
// DescribeYourselfTo
//
void G4VCSGfaceted::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

//
// GetExtent
//
// Define the sides of the box into which our solid instance would fit.
//
G4VisExtent G4VCSGfaceted::GetExtent() const
{
  static const G4ThreeVector xMax(1, 0, 0), xMin(-1, 0, 0),
                             yMax(0, 1, 0), yMin(0, -1, 0),
                             zMax(0, 0, 1), zMin(0, 0, -1);
  static const G4ThreeVector* axes[6] = {&xMin, &xMax, &yMin,
                                         &yMax, &zMin, &zMax};

  G4double answers[6] = {-kInfinity, -kInfinity, -kInfinity,
                         -kInfinity, -kInfinity, -kInfinity};

  G4VCSGface** face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    const G4ThreeVector** axis = axes + 5;
    G4double* answer = answers + 5;
    do  // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double testFace = (*face)->Extent(**axis);
      if (testFace > *answer)
      {
        *answer = testFace;
      }
    } while (--axis, --answer >= answers);

  } while (++face < faces + numFace);

  return {-answers[0], answers[1],
          -answers[2], answers[3],
          -answers[4], answers[5]};
}

//
// GetEntityType
//
G4GeometryType G4VCSGfaceted::GetEntityType() const
{
  return {"G4CSGfaceted"};
}

//
// Stream object contents to an output stream
//
std::ostream& G4VCSGfaceted::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4VCSGfaceted\n"
     << " Parameters: \n"
     << "    number of faces: " << numFace << "\n"
     << "-----------------------------------------------------------\n";

  return os;
}

//
// GetCubVolStatistics
//
G4int G4VCSGfaceted::GetCubVolStatistics() const
{
  return fStatistics;
}

//
// GetCubVolEpsilon
//
G4double G4VCSGfaceted::GetCubVolEpsilon() const
{
  return fCubVolEpsilon;
}

//
// SetCubVolStatistics
//
void G4VCSGfaceted::SetCubVolStatistics(G4int st)
{
  fCubicVolume = 0.;
  fStatistics = st;
}

//
// SetCubVolEpsilon
//
void G4VCSGfaceted::SetCubVolEpsilon(G4double ep)
{
  fCubicVolume = 0.;
  fCubVolEpsilon = ep;
}

//
// GetAreaStatistics
//
G4int G4VCSGfaceted::GetAreaStatistics() const
{
  return fStatistics;
}

//
// GetAreaAccuracy
//
G4double G4VCSGfaceted::GetAreaAccuracy() const
{
  return fAreaAccuracy;
}

//
// SetAreaStatistics
//
void G4VCSGfaceted::SetAreaStatistics(G4int st)
{
  fSurfaceArea = 0.;
  fStatistics = st;
}

//
// SetAreaAccuracy
//
void G4VCSGfaceted::SetAreaAccuracy(G4double ep)
{
  fSurfaceArea = 0.;
  fAreaAccuracy = ep;
}

//
// GetCubicVolume
//
G4double G4VCSGfaceted::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&vcsgMutex);
    if (fCubicVolume == 0)
    {
      fCubicVolume = EstimateCubicVolume(fStatistics, fCubVolEpsilon);
    }
    l.unlock();
  }
  return fCubicVolume;
}

//
// GetSurfaceArea
//
G4double G4VCSGfaceted::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&vcsgMutex);
    if (fSurfaceArea == 0)
    {
      fSurfaceArea = EstimateSurfaceArea(fStatistics, fAreaAccuracy);
    }
    l.unlock();
  }
  return fSurfaceArea;
}

//
// GetPolyhedron
//
G4Polyhedron* G4VCSGfaceted::GetPolyhedron() const
{
  if (fpPolyhedron == nullptr || fRebuildPolyhedron
      || fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation()
           != fpPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&polyhedronMutex);
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fpPolyhedron;
}

//
// GetPointOnSurfaceGeneric proportional to Areas of faces
// in case of GenericPolycone or GenericPolyhedra
//
G4ThreeVector G4VCSGfaceted::GetPointOnSurfaceGeneric() const
{
  G4double area = 0.;
  G4VCSGface** face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    area += (*face)->SurfaceArea();
  } while (++face < faces + numFace);

  if (area <= 0.)
  {
    return {};
  }

  G4double chose = area * G4QuickRand();
  G4double cumulative = 0.;
  face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    cumulative += (*face)->SurfaceArea();
    if (chose < cumulative)
    {
      return (*face)->GetPointOnFace();
    }
  } while (++face < faces + numFace);

  return faces[numFace - 1]->GetPointOnFace();
}
