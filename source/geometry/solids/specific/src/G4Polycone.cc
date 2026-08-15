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
// Implementation of G4Polycone, a CSG polycone
//
// Author: David C. Williams (UCSC), 1998
// --------------------------------------------------------------------

#include "G4Polycone.hh"

#if !defined(G4GEOM_USE_UPOLYCONE)

#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4EnclosingCylinder.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeomTools.hh"
#include "G4PolyPhiFace.hh"
#include "G4PolyconeSide.hh"
#include "G4QuickRand.hh"
#include "G4ReduciblePolygon.hh"
#include "G4VPVParameterisation.hh"
#include "G4VoxelLimits.hh"

#include <algorithm>

namespace
{
  G4Mutex surface_elementsMutex = G4MUTEX_INITIALIZER;
  G4Mutex polyconeMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

// Constructor (GEANT3 style parameters)
//
G4Polycone::G4Polycone(const G4String& name, G4double phiStart,
                       G4double phiTotal, G4int numZPlanes,
                       const G4double zPlane[], const G4double rInner[],
                       const G4double rOuter[])
  : G4VCSGfaceted(name)
{
  //
  // Some historical ugliness
  //
  original_parameters = new G4PolyconeHistorical();

  original_parameters->Start_angle = phiStart;
  original_parameters->Opening_angle = phiTotal;
  original_parameters->Num_z_planes = numZPlanes;
  original_parameters->Z_values = new G4double[numZPlanes];
  original_parameters->Rmin = new G4double[numZPlanes];
  original_parameters->Rmax = new G4double[numZPlanes];

  for (G4int i = 0; i < numZPlanes; ++i)
  {
    if (rInner[i] > rOuter[i])
    {
      DumpInfo();
      std::ostringstream message;
      message << "Cannot create a Polycone with rInner > rOuter for the same Z" << G4endl
              << "        rInner > rOuter for the same Z !" << G4endl << "        rMin[" << i
              << "] = " << rInner[i] << " -- rMax[" << i << "] = " << rOuter[i];
      G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002", FatalErrorInArgument, message);
    }
    if ((i < numZPlanes - 1) && (zPlane[i] == zPlane[i + 1]))
    {
      if ((rInner[i] > rOuter[i + 1]) || (rInner[i + 1] > rOuter[i]))
      {
        DumpInfo();
        std::ostringstream message;
        message << "Cannot create a Polycone with no contiguous segments." << G4endl
                << "        Segments are not contiguous !" << G4endl << "        rMin[" << i
                << "] = " << rInner[i] << " -- rMax[" << i + 1 << "] = " << rOuter[i + 1] << G4endl
                << "        rMin[" << i + 1 << "] = " << rInner[i + 1] << " -- rMax[" << i
                << "] = " << rOuter[i];
        G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002", FatalErrorInArgument, message);
      }
    }
    original_parameters->Z_values[i] = zPlane[i];
    original_parameters->Rmin[i] = rInner[i];
    original_parameters->Rmax[i] = rOuter[i];
  }

  //
  // Build RZ polygon using special PCON/PGON GEANT3 constructor
  //
  auto rz = new G4ReduciblePolygon(rInner, rOuter, zPlane, numZPlanes);

  //
  // Do the real work
  //
  Create(phiStart, phiTotal, rz);

  delete rz;
}

// Constructor (generic parameters)
//
G4Polycone::G4Polycone(const G4String& name, G4double phiStart,
                       G4double phiTotal, G4int numRZ,
                       const G4double r[], const G4double z[])
  : G4VCSGfaceted(name)
{
  auto rz = new G4ReduciblePolygon(r, z, numRZ);

  Create(phiStart, phiTotal, rz);

  // Set original_parameters struct for consistency
  //

  G4bool convertible = SetOriginalParameters(rz);

  if (!convertible)
  {
    std::ostringstream message;
    message << "Polycone " << GetName() << "cannot be converted" << G4endl
            << "to Polycone with (Rmin,Rmaz,Z) parameters!";
    G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002", FatalException, message,
                "Use G4GenericPolycone instead!");
  }
  else
  {
    G4cout << "INFO: Converting polycone " << GetName() << G4endl
           << "to optimized polycone with (Rmin,Rmaz,Z) parameters !" << G4endl;
  }
  delete rz;
}

// Create
//
// Generic create routine, called by each constructor after
// conversion of arguments
//
void G4Polycone::Create(G4double phiStart, G4double phiTotal,
                        G4ReduciblePolygon* rz)
{
  //
  // Perform checks of rz values
  //
  if (rz->Amin() < 0.0)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        All R values must be >= 0 !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002", FatalErrorInArgument, message);
  }

  G4double rzArea = rz->Area();
  if (rzArea < -kCarTolerance)
  {
    rz->ReverseOrder();
  }
  else if (rzArea < kCarTolerance)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        R/Z cross section is zero or near zero: " << rzArea;
    G4Exception("G4Polycone::Create()", "GeomSolids0002", FatalErrorInArgument, message);
  }

  if ((!rz->RemoveDuplicateVertices(kCarTolerance))
      || (!rz->RemoveRedundantVertices(kCarTolerance)))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        Too few unique R/Z values !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002", FatalErrorInArgument, message);
  }

  if (rz->CrossesItself(1 / kInfinity))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        R/Z segments cross !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002", FatalErrorInArgument, message);
  }

  numCorner = rz->NumVertices();

  startPhi = phiStart;
  while (startPhi < 0.)  // Loop checking, 13.08.2015, G.Cosmo
  {
    startPhi += twopi;
  }
  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if ((phiTotal <= 0) || (phiTotal > twopi * (1 - DBL_EPSILON)))
  {
    phiIsOpen = false;
    startPhi = 0.;
    endPhi = twopi;
  }
  else
  {
    phiIsOpen = true;
    endPhi = startPhi + phiTotal;
  }

  //
  // Allocate corner array.
  //
  corners = new G4PolyconeSideRZ[numCorner];

  //
  // Copy corners
  //
  G4ReduciblePolygonIterator iterRZ(rz);

  G4PolyconeSideRZ* next = corners;
  iterRZ.Begin();
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    next->r = iterRZ.GetA();
    next->z = iterRZ.GetB();
  } while (++next, iterRZ.Next());

  //
  // Allocate face pointer array
  //
  numFace = phiIsOpen ? numCorner + 2 : numCorner;
  faces = new G4VCSGface*[numFace];

  //
  // Construct conical faces
  //
  // But! Don't construct a face if both points are at zero radius!
  //
  G4PolyconeSideRZ *corner = corners, *prev = corners + numCorner - 1, *nextNext;
  G4VCSGface** face = faces;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    next = corner + 1;
    if (next >= corners + numCorner)
    {
      next = corners;
    }
    nextNext = next + 1;
    if (nextNext >= corners + numCorner)
    {
      nextNext = corners;
    }

    if (corner->r < 1 / kInfinity && next->r < 1 / kInfinity)
    {
      continue;
    }

    //
    // We must decide here if we can dare declare one of our faces
    // as having a "valid" normal (i.e. allBehind = true). This
    // is never possible if the face faces "inward" in r.
    //
    G4bool allBehind;
    if (corner->z > next->z)
    {
      allBehind = false;
    }
    else
    {
      //
      // Otherwise, it is only true if the line passing
      // through the two points of the segment do not
      // split the r/z cross section
      //
      allBehind = !rz->BisectedBy(corner->r, corner->z, next->r, next->z, kCarTolerance);
    }

    *face++ = new G4PolyconeSide(prev, corner, next, nextNext, startPhi, endPhi - startPhi,
                                 phiIsOpen, allBehind);
  } while (prev = corner, corner = next, corner > corners);

  if (phiIsOpen)
  {
    //
    // Construct phi open edges
    //
    *face++ = new G4PolyPhiFace(rz, startPhi, 0, endPhi);
    *face++ = new G4PolyPhiFace(rz, endPhi, 0, startPhi);
  }

  //
  // We might have dropped a face or two: recalculate numFace
  //
  numFace = (G4int)(face - faces);

  //
  // Make enclosingCylinder
  //
  enclosingCylinder = new G4EnclosingCylinder(rz, phiIsOpen, phiStart, phiTotal);

  BuildZSectionIndex();
  BuildRZIndex();
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Polycone::G4Polycone(__void__& a)
  : G4VCSGfaceted(a), startPhi(0.), endPhi(0.), numCorner(0)
{
}

// Destructor
//
G4Polycone::~G4Polycone()
{
  delete[] corners;
  delete original_parameters;
  delete enclosingCylinder;
  delete fElements;
  delete fpPolyhedron;
  corners = nullptr;
  original_parameters = nullptr;
  enclosingCylinder = nullptr;
  fElements = nullptr;
  fpPolyhedron = nullptr;
}

// Copy constructor
//
G4Polycone::G4Polycone(const G4Polycone& source) : G4VCSGfaceted(source)
{
  CopyStuff(source);
}

// Assignment operator
//
G4Polycone& G4Polycone::operator=(const G4Polycone& source)
{
  if (this == &source)
  {
    return *this;
  }

  G4VCSGfaceted::operator=(source);

  delete[] corners;
  delete original_parameters;

  delete enclosingCylinder;

  CopyStuff(source);

  return *this;
}

// CopyStuff
//
void G4Polycone::CopyStuff(const G4Polycone& source)
{
  //
  // Simple stuff
  //
  startPhi = source.startPhi;
  endPhi = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  fZPlanesMonotonic = source.fZPlanesMonotonic;
  fSinStartPhi = source.fSinStartPhi;
  fCosStartPhi = source.fCosStartPhi;
  fSinEndPhi = source.fSinEndPhi;
  fCosEndPhi = source.fCosEndPhi;
  numCorner = source.numCorner;

  //
  // The corner array
  //
  corners = new G4PolyconeSideRZ[numCorner];

  G4PolyconeSideRZ *corn = corners, *sourceCorn = source.corners;
  do  // Loop checking, 13.08.2015, G.Cosmo
  {
    *corn = *sourceCorn;
  } while (++sourceCorn, ++corn < corners + numCorner);

  //
  // Original parameters
  //
  if (source.original_parameters != nullptr)
  {
    original_parameters = new G4PolyconeHistorical(*source.original_parameters);
  }

  //
  // Enclosing cylinder
  //
  enclosingCylinder = new G4EnclosingCylinder(*source.enclosingCylinder);

  //
  // Surface elements
  //
  delete fElements;
  fElements = nullptr;

  //
  // Polyhedron
  //
  fRebuildPolyhedron = false;
  delete fpPolyhedron;
  fpPolyhedron = nullptr;

  fRZSection = source.fRZSection;
  fRZCone = source.fRZCone;
  fRZVoxel = CopyZCandidateData(source.fRZVoxel);
}

// Reset
//
G4bool G4Polycone::Reset()
{
  //
  // Clear old setup
  //
  G4VCSGfaceted::DeleteStuff();
  delete[] corners;
  delete enclosingCylinder;
  delete fElements;
  corners = nullptr;
  fElements = nullptr;
  enclosingCylinder = nullptr;
  fRZSection.clear();
  fRZCone.clear();
  fRZVoxel.reset();

  //
  // Rebuild polycone
  //
  auto rz =
    new G4ReduciblePolygon(original_parameters->Rmin, original_parameters->Rmax,
                           original_parameters->Z_values, original_parameters->Num_z_planes);
  Create(original_parameters->Start_angle, original_parameters->Opening_angle, rz);
  delete rz;

  return false;
}

// BuildRZIndex
//
void G4Polycone::BuildRZIndex()
{
  fRZSection.clear();
  fRZCone.clear();
  fRZVoxel.reset();
  fSinStartPhi = phiIsOpen ? std::sin(startPhi) : 0.;
  fCosStartPhi = phiIsOpen ? std::cos(startPhi) : 1.;
  fSinEndPhi = phiIsOpen ? std::sin(endPhi) : 0.;
  fCosEndPhi = phiIsOpen ? std::cos(endPhi) : 1.;

  // Standard ordered Z planes permit direct segment lookup and exact kernels.
  fZPlanesMonotonic = original_parameters != nullptr
                    && original_parameters->Num_z_planes >= 2;
  if (fZPlanesMonotonic)
  {
    const auto* z = original_parameters->Z_values;
    for (G4int i = 1; i < original_parameters->Num_z_planes; ++i)
    {
      if (z[i] <= z[i - 1])
      {
        fZPlanesMonotonic = false;
        break;
      }
    }
  }

  if (numCorner < 3)
  {
    return;
  }

  // Small solids are faster without an allocated candidate table.
  fRZVoxel = MakeZCandidateData(numFace);

  fRZSection.reserve(numCorner);
  for (G4int i = 0; i < numCorner; ++i)
  {
    const auto& a = corners[i];
    const auto& b = corners[(i + 1) % numCorner];
    if (!HasRZEdge(a, b))
    {
      continue;
    }
    RZSection sec;
    sec.r1 = b.r;
    sec.z1 = b.z;
    if (!InitRZSection(sec, a, b))
    {
      continue;
    }
    if (sec.z0 <= sec.z1)
    {
      // Record whether this edge alone defines a supporting boundary.
      G4bool hasPositive = false;
      G4bool hasNegative = false;
      for (G4int j = 0; j < numCorner; ++j)
      {
        const auto& c = corners[j];
        const G4double side = (c.r - sec.r0)
                            * sec.rNorm + (c.z - sec.z0) * sec.zNorm;
        if (side > kCarTolerance)
        {
          hasPositive = true;
        }
        else if (side < -kCarTolerance)
        {
          hasNegative = true;
        }
        if (hasPositive && hasNegative)
        {
          break;
        }
      }
      sec.allBehind = !(hasPositive && hasNegative);
    }
    const G4double rr[2] = {sec.r0, sec.r1};
    const G4double zz[2] = {sec.z0, sec.z1};
    AddRZSectionCandidate(sec, fRZSection, fRZVoxel.get());
    fRZCone.emplace_back(rr, zz);
  }

  BuildRZVoxel(fRZSection, fRZVoxel);
}

// InsideZPlane
//
EInside G4Polycone::InsideZPlane(const G4ThreeVector& p) const
{
  if (!HasZPlaneInsideFastPath())
  {
    return InsideRZOnly(p);
  }

  const auto* pars = original_parameters;
  const G4int n = pars->Num_z_planes;
  const auto* zPlane = pars->Z_values;
  const auto* rMin = pars->Rmin;
  const auto* rMax = pars->Rmax;
  const G4double z = p.z();
  const G4double tol = 0.5 * kCarTolerance;

  if (z < zPlane[0] - tol || z > zPlane[n - 1] + tol)
  {
    return kOutside;
  }

  G4int iz = 0;

  // Avoid binary-search setup for the usual small number of Z planes.
  if (n <= kSmallZPlaneFastPath)
  {
    iz = n - 2;
    for (G4int i = 0; i < n - 1; ++i)
    {
      if (z <= zPlane[i + 1] + tol)
      {
        iz = i;
        break;
      }
    }
  }
  else
  {
    const auto* upper = std::upper_bound(zPlane + 1, zPlane + n, z + tol);
    iz = std::min(n - 2, G4int(upper - zPlane) - 1);
  }

  const G4double dz = zPlane[iz + 1] - zPlane[iz];
  const G4double t = (std::fabs(dz) > 0.) ? (z - zPlane[iz]) / dz : 0.;
  const G4double rIn = rMin[iz] + t * (rMin[iz + 1] - rMin[iz]);
  const G4double rOut = rMax[iz] + t * (rMax[iz + 1] - rMax[iz]);
  const G4double rho2 = p.x() * p.x() + p.y() * p.y();

  const G4double rOutPlus = rOut + tol;
  if (rho2 > rOutPlus * rOutPlus)
  {
    return kOutside;
  }
  if (rIn > tol)
  {
    const G4double rInMinus = std::max(0., rIn - tol);
    if (rho2 < rInMinus * rInMinus)
    {
      return kOutside;
    }
  }

  const G4bool onZ = std::fabs(z - zPlane[0]) <= tol || std::fabs(z - zPlane[n - 1]) <= tol;
  const G4double rOutMinus = std::max(0., rOut - tol);
  const G4bool onOuter = rho2 >= rOutMinus * rOutMinus;
  G4bool onInner = false;
  if (rIn > tol)
  {
    const G4double rInPlus = rIn + tol;
    onInner = rho2 <= rInPlus * rInPlus;
  }
  return (onZ || onOuter || onInner) ? kSurface : kInside;
}

// InsideRZOnly
//
EInside G4Polycone::InsideRZOnly(const G4ThreeVector& p) const
{
  const G4double rho = p.perp();
  const G4double z = p.z();
  const G4double tol = 0.5 * kCarTolerance;
  EInside answer = kInside;
  G4double best = kInfinity;

  auto check = [&](G4int i)
  {
    const auto& sec = fRZSection[i];
    const G4double dr = rho - sec.r0;
    const G4double dz = z - sec.z0;
    const G4double q = dr * sec.rS + dz * sec.zS;
    G4double norm = dr * sec.rNorm + dz * sec.zNorm;
    G4double distOut2 = 0.;
    if (q < 0.)
    {
      distOut2 = q * q;
    }
    else if (q > sec.length)
    {
      distOut2 = (q - sec.length) * (q - sec.length);
    }
    const G4double distance = std::sqrt(norm * norm + distOut2);
    EInside result = kInside;
    if (std::fabs(norm) < tol && distOut2 < tol * tol)
    {
      result = kSurface;
    }
    else if (norm > tol)
    {
      result = kOutside;
    }
    if (result == kSurface)
    {
      answer = kSurface;
      best = 0.;
      return;
    }
    if (distance < best)
    {
      best = distance;
      answer = result;
    }
  };

  if (fRZVoxel == nullptr || fRZVoxel->candidate.empty())
  {
    for (G4int i = 0; i < G4int(fRZSection.size()); ++i)
    {
      check(i);
    }
  }
  else
  {
    for (auto i : fRZVoxel->global)
    {
      check(i);
    }
    for (auto i : FindZCandidates(fRZVoxel->boundary, fRZVoxel->candidate, z))
    {
      check(i);
    }
  }

  return answer;
}

// IntersectRZSection
//
G4bool G4Polycone::IntersectRZSection(G4int i, const G4ThreeVector& p, const G4ThreeVector& v,
                                      G4bool outgoing, G4double surfTolerance,
                                      G4double& distance, G4double& distFromSurface,
                                      G4ThreeVector* normal) const
{
  const auto& sec = fRZSection[i];
  const auto& cone = fRZCone[i];
  G4double s1 = 0., s2 = 0.;
  const G4int nside = cone.LineHitsCone(p, v, &s1, &s2);
  if (nside == 0)
  {
    return false;
  }

  const G4double normSign = outgoing ? +1. : -1.;
  auto distanceAway = [&](G4double& distOutside2)
  {
    const G4double rx = p.perp();
    const G4double deltaR = rx - sec.r0;
    const G4double deltaZ = p.z() - sec.z0;
    const G4double answer = deltaR * sec.rNorm + deltaZ * sec.zNorm;
    const G4double q = deltaR * sec.rS + deltaZ * sec.zS;
    if (q < 0.)
    {
      distOutside2 = q * q;
    }
    else if (q > sec.length)
    {
      distOutside2 = sqr(q - sec.length);
    }
    else
    {
      distOutside2 = 0.;
    }
    return answer;
  };

  auto tryHit = [&](G4double root)
  {
    const G4ThreeVector hit = p + root * v;
    const G4double rx = hit.perp();
    if (!cone.HitOn(rx, hit.z()))
    {
      return false;
    }
    if (phiIsOpen && PhiState(hit) == kOutside)
    {
      return false;
    }
    if (rx < DBL_MIN)
    {
      const G4double nz = sec.zNorm < 0. ? -1. : 1.;
      if (normSign * v.z() * nz <= 0.)
      {
        return false;
      }
      if (normal != nullptr)
      {
        *normal = G4ThreeVector(0., 0., nz);
      }
    }
    else
    {
      const G4double rvdot = v.x() * hit.x() + v.y() * hit.y();
      if (normSign * (sec.rNorm * rvdot / rx + sec.zNorm * v.z()) <= 0.)
      {
        return false;
      }
      if (normal != nullptr)
      {
        *normal = G4ThreeVector(sec.rNorm * hit.x() / rx,
                                sec.rNorm * hit.y() / rx, sec.zNorm);
      }
    }

    G4double pr = p.perp();
    if (pr < DBL_MIN)
    {
      pr = DBL_MIN;
    }
    const G4double prvdot = v.x() * p.x() + v.y() * p.y();
    if (normSign * (sec.rNorm * prvdot / pr + sec.zNorm * v.z()) > 0.)
    {
      G4double distOutside2 = 0.;
      distFromSurface = -normSign * distanceAway(distOutside2);
      if (distOutside2 < surfTolerance * surfTolerance
       && distFromSurface > -surfTolerance)
      {
        distance = root;
        return true;
      }
    }
    else
    {
      distFromSurface = root;
    }

    if (root > 0.)
    {
      distance = root;
      return true;
    }
    return false;
  };

  if (tryHit(s1))
  {
    return true;
  }
  return nside > 1 && tryHit(s2);
}

// IntersectPhiPlane
//
G4bool G4Polycone::IntersectPhiPlane(const G4ThreeVector& p, const G4ThreeVector& v,
                                     G4bool outgoing, G4double surfTolerance,
                                     G4double& distance, G4double& distFromSurface,
                                     G4ThreeVector* normal) const
{
  if (!phiIsOpen)
  {
    return false;
  }

  const G4double sStart = fSinStartPhi, cStart = fCosStartPhi;
  const G4double sEnd = fSinEndPhi, cEnd = fCosEndPhi;
  const G4ThreeVector nStart(sStart, -cStart, 0.);
  const G4ThreeVector nEnd(-sEnd, cEnd, 0.);
  const G4double normSign = outgoing ? 1. : -1.;

  auto tryPlane = [&](const G4ThreeVector& planeNormal)
  {
    const G4double denom = v.dot(planeNormal);
    if (std::fabs(denom) < DBL_MIN)
    {
      return false;
    }
    const G4double root = -p.dot(planeNormal) / denom;
    if (root < -surfTolerance || root > distance + surfTolerance)
    {
      return false;
    }
    if (normSign * denom <= 0.)
    {
      return false;
    }
    const G4ThreeVector hit = p + root * v;
    if (InsideZPlane(hit) == kOutside)
    {
      return false;
    }
    const G4ThreeVector probe = hit + (2. * surfTolerance) * v;
    const EInside probeState = (InsideZPlane(probe) == kOutside
                            || PhiState(probe) == kOutside) ? kOutside : kInside;
    if ((outgoing && probeState != kOutside) || (!outgoing && probeState == kOutside))
    {
      return false;
    }
    distance = std::max(0., root);
    distFromSurface = std::fabs(p.dot(planeNormal));
    if (normal != nullptr)
    {
      *normal = planeNormal;
    }
    return true;
  };

  G4bool hit = tryPlane(nStart);
  hit = tryPlane(nEnd) || hit;
  return hit;
}

// DistanceToSmallProfile
//
G4double G4Polycone::DistanceToFullPhiTwoPlaneNoNormal(const G4ThreeVector& p,
                                                       const G4ThreeVector& v,
                                                       G4bool outgoing) const
{
  const G4double tolerance = kCarTolerance / 2.;
  const G4double* zPlane = original_parameters->Z_values;
  const G4double* rMin = original_parameters->Rmin;
  const G4double* rMax = original_parameters->Rmax;
  const G4double px = p.x(), py = p.y(), pz = p.z();
  const G4double vx = v.x(), vy = v.y(), vz = v.z();

  const G4double z0 = zPlane[0];
  const G4double z1 = zPlane[1];
  const G4double invDz = 1. / (z1 - z0);
  const G4bool hasInner = rMin[0] > tolerance || rMin[1] > tolerance;
  const G4double zMid = 0.5 * (z0 + z1);
  const G4double dzHalf = 0.5 * (z1 - z0);
  const G4double rminMinus = rMin[0], rminPlus = rMin[1];
  const G4double rmaxMinus = rMax[0], rmaxPlus = rMax[1];
  const G4double zc = pz - zMid;
  const G4double halfCarTolerance = 0.5 * kCarTolerance;
  const G4double kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  const G4double halfRadTolerance = 0.5 * kRadTolerance;

  // The common solid-cone entry needs only the outer conical surface and caps.
  if (!outgoing && !hasInner)
  {
    const G4double tanRMax = (rmaxPlus - rmaxMinus) * 0.5 / dzHalf;
    const G4double secRMax = std::sqrt(1. + tanRMax * tanRMax);
    const G4double rMaxAv = 0.5 * (rmaxMinus + rmaxPlus);
    const G4double tolIDz = dzHalf - halfCarTolerance;
    const G4double tolODz = dzHalf + halfCarTolerance;

    if (std::fabs(zc) >= tolIDz)
    {
      if (zc * vz < 0.)
      {
        G4double sd = (std::fabs(zc) - dzHalf) / std::fabs(vz);
        if (sd < 0.)
        {
          sd = 0.;
        }
        const G4double x = px + sd * vx;
        const G4double y = py + sd * vy;
        const G4double rho2 = x * x + y * y;
        const G4double rAtZ = (vz > 0.) ? rmaxMinus : rmaxPlus;
        const G4double tolIRMax = rAtZ - halfRadTolerance * secRMax;
        const G4double tolIRMax2 = (tolIRMax > 0.) ? tolIRMax * tolIRMax : 0.;
        if (rho2 <= tolIRMax2)
        {
          return sd;
        }
      }
      else
      {
        return kInfinity;
      }
    }

    const G4double t1 = vx * vx + vy * vy;
    const G4double t2 = px * vx + py * vy;
    const G4double t3 = px * px + py * py;
    const G4double rout = tanRMax * zc + rMaxAv;
    const G4double tanVzMax = tanRMax * vz;
    const G4double nt1 = t1 - tanVzMax * tanVzMax;
    const G4double nt2 = t2 - tanVzMax * rout;
    const G4double nt3 = t3 - rout * rout;
    if (std::fabs(nt1) > kRadTolerance)
    {
      const G4double bq = nt2 / nt1;
      const G4double cq = nt3 / nt1;
      const G4double disc = bq * bq - cq;
      if ((nt3 > rout * rout * kRadTolerance * kRadTolerance * secRMax * secRMax)
       || (rout < 0.))
      {
        if (disc >= 0.)
        {
          const G4double sqrtDisc = std::sqrt(disc);
          G4double sd = kInfinity;
          if ((rout < 0.) && (nt3 <= 0.))
          {
            sd = (bq > 0.) ? cq / (-bq - sqrtDisc) : -bq + sqrtDisc;
          }
          else if ((bq <= 0.) && (cq >= 0.))
          {
            sd = cq / (-bq + sqrtDisc);
          }
          else if (cq <= 0.)
          {
            sd = -bq + sqrtDisc;
            if ((sd < 0.) && (sd > -halfRadTolerance))
            {
              sd = 0.;
            }
          }
          else
          {
            return kInfinity;
          }
          if (sd >= 0.)
          {
            const G4double zi = zc + sd * vz;
            if (std::fabs(zi) <= tolODz)
            {
              return sd;
            }
          }
        }
      }
      else if ((t3 > halfRadTolerance * halfRadTolerance) && (nt2 < 0.)
               && (disc >= 0.) && (std::fabs(zc) <= tolIDz))
      {
        return 0.;
      }
    }
    else if (std::fabs(nt2) > kRadTolerance)
    {
      const G4double sd = -0.5 * nt3 / nt2;
      if (sd < 0.)
      {
        return kInfinity;
      }
      const G4double zi = zc + sd * vz;
      if ((std::fabs(zi) <= tolODz) && (nt2 < 0.))
      {
        return sd;
      }
    }
    return kInfinity;
  }

  // Starting inside fixes the first exit to an inner/outer wall or a Z cap.
  if (outgoing)
  {
    G4double snxt = kInfinity;
    if (vz > 0.)
    {
      const G4double pdist = dzHalf - zc;
      if (pdist > halfCarTolerance)
      {
        snxt = pdist / vz;
      }
      else
      {
        return 0.;
      }
    }
    else if (vz < 0.)
    {
      const G4double pdist = dzHalf + zc;
      if (pdist > halfCarTolerance)
      {
        snxt = -pdist / vz;
      }
      else
      {
        return 0.;
      }
    }

    const G4double t1 = vx * vx + vy * vy;
    const G4double t2 = px * vx + py * vy;
    const G4double t3 = px * px + py * py;
    const G4double tanRMax = (rmaxPlus - rmaxMinus) * 0.5 / dzHalf;
    const G4double secRMax = std::sqrt(1. + tanRMax * tanRMax);
    const G4double rMaxAv = 0.5 * (rmaxMinus + rmaxPlus);
    const G4double rout = tanRMax * zc + rMaxAv;
    const G4double tanVzMax = tanRMax * vz;
    const G4double nt1 = t1 - tanVzMax * tanVzMax;
    const G4double nt2 = t2 - tanVzMax * rout;
    const G4double nt3 = t3 - rout * rout;
    G4double srd = kInfinity;

    G4double deltaRoi2 = 1.;
    if (vz > 0.)
    {
      deltaRoi2 =
        snxt * snxt * t1 + 2. * snxt * t2 + t3 - rmaxPlus * (rmaxPlus + kRadTolerance * secRMax);
    }
    else if (vz < 0.)
    {
      deltaRoi2 =
        snxt * snxt * t1 + 2. * snxt * t2 + t3 - rmaxMinus * (rmaxMinus + kRadTolerance * secRMax);
    }

    if ((nt1 != 0.) && (deltaRoi2 > 0.))
    {
      const G4double bq = nt2 / nt1;
      const G4double cq = nt3 / nt1;
      const G4double disc = bq * bq - cq;
      if (disc >= 0.)
      {
        const G4double sqrtDisc = std::sqrt(disc);
        if ((nt3 > -halfRadTolerance) && (nt2 >= 0.))
        {
          return 0.;
        }
        srd = (bq > 0.) ? -bq - sqrtDisc : cq / (-bq + sqrtDisc);
        G4double zi = zc + srd * vz;
        G4double ri = tanRMax * zi + rMaxAv;
        if ((ri < 0.) || (srd < halfRadTolerance))
        {
          const G4double sr2 = (bq > 0.) ? cq / (-bq - sqrtDisc) : -bq + sqrtDisc;
          zi = zc + sr2 * vz;
          ri = tanRMax * zi + rMaxAv;
          srd = ((ri >= 0.) && (sr2 > halfRadTolerance)) ? sr2 : kInfinity;
        }
      }
      else
      {
        return 0.;
      }
    }
    else if ((nt2 != 0.) && (deltaRoi2 > 0.))
    {
      return 0.;
    }

    if (hasInner)
    {
      const G4double tanRMin = (rminPlus - rminMinus) * 0.5 / dzHalf;
      const G4double rMinAv = 0.5 * (rminMinus + rminPlus);
      const G4double rin = tanRMin * zc + rMinAv;
      const G4double tanVzMin = tanRMin * vz;
      const G4double inNt1 = t1 - tanVzMin * tanVzMin;
      if (inNt1 != 0.)
      {
        const G4double inNt2 = t2 - tanVzMin * rin;
        const G4double inNt3 = t3 - rin * rin;
        const G4double bq = inNt2 / inNt1;
        const G4double cq = inNt3 / inNt1;
        const G4double disc = bq * bq - cq;
        if (disc >= 0.)
        {
          const G4double sqrtDisc = std::sqrt(disc);
          if (inNt3 < kRadTolerance * (rin + 0.25 * kRadTolerance))
          {
            if (inNt2 < 0.)
            {
              return 0.;
            }
          }
          else
          {
            G4double sr2 = (bq > 0.) ? -bq - sqrtDisc : cq / (-bq + sqrtDisc);
            G4double zi = zc + sr2 * vz;
            G4double ri = tanRMin * zi + rMinAv;
            if ((ri < 0.) || (sr2 < halfRadTolerance))
            {
              const G4double sr3 = (bq > 0.) ? cq / (-bq - sqrtDisc) : -bq + sqrtDisc;
              if (sr3 > halfRadTolerance)
              {
                zi = zc + sr3 * vz;
                ri = tanRMin * zi + rMinAv;
                if ((ri >= 0.) && (sr3 < srd))
                {
                  srd = sr3;
                }
              }
            }
            else if ((sr2 < srd) && (sr2 > halfCarTolerance))
            {
              srd = sr2;
            }
          }
        }
      }
    }

    if (srd < snxt)
    {
      snxt = srd;
    }
    return (snxt < halfCarTolerance) ? 0. : snxt;
  }

  const G4double outerSlope = (rMax[1] - rMax[0]) * invDz;
  const G4double innerSlope = (rMin[1] - rMin[0]) * invDz;
  const G4double zMin = z0 - tolerance, zMax = z1 + tolerance;
  G4double best = kInfinity;

  // The remaining entry case is hollow; test all four analytic boundaries.
  auto radiusAt = [=](G4double z, G4double r0, G4double slope)
  {
    return r0 + (z - z0) * slope;
  };
  auto radialValid = [&](G4double z, G4double rho2)
  {
    const G4double rOuter = radiusAt(z, rMax[0], outerSlope) + tolerance;
    if (rho2 > rOuter * rOuter)
    {
      return false;
    }
    const G4double rInner = std::max(0., radiusAt(z, rMin[0], innerSlope) - tolerance);
    if (rho2 < rInner * rInner)
    {
      return false;
    }
    return true;
  };

  auto tryZPlane = [&](G4double z, G4double normalZ)
  {
    if (std::fabs(vz) <= DBL_MIN)
    {
      return;
    }
    const G4double root = (z - pz) / vz;
    const G4double normalDotV = normalZ * vz;
    if (root < -tolerance || root > best + tolerance
     || normalDotV >= 0.)
    {
      return;
    }
    const G4double x = px + root * vx;
    const G4double y = py + root * vy;
    if (radialValid(z, x * x + y * y))
    {
      best = std::max(0., root);
    }
  };

  auto tryCone = [&](const G4double* radius, G4double slope, G4bool inner)
  {
    const G4double rAtP = radiusAt(pz, radius[0], slope);
    const G4double svz = slope * vz;
    const G4double a = vx * vx + vy * vy - svz * svz;
    const G4double b = 2. * (px * vx + py * vy - svz * rAtP);
    const G4double c = px * px + py * py - rAtP * rAtP;
    G4double roots[2];
    G4int nRoots = 0;
    if (std::fabs(a) < DBL_EPSILON)
    {
      if (std::fabs(b) > DBL_EPSILON)
      {
        roots[nRoots++] = -c / b;
      }
    }
    else
    {
      G4double radical = b * b - 4. * a * c;
      const G4double eps = DBL_EPSILON * std::fabs(b);
      if (radical >= -eps)
      {
        if (radical <= eps)
        {
          roots[nRoots++] = -0.5 * b / a;
        }
        else
        {
          radical = std::sqrt(radical);
          const G4double q = -0.5 * (b + (b < 0. ? -radical : radical));
          roots[nRoots++] = q / a;
          roots[nRoots++] = c / q;
          if (roots[1] < roots[0])
          {
            std::swap(roots[0], roots[1]);
          }
        }
      }
    }

    const G4double rSign = inner ? -1. : 1.;
    for (G4int k = 0; k < nRoots; ++k)
    {
      const G4double root = roots[k];
      if (root < -tolerance || root > best + tolerance)
      {
        continue;
      }
      const G4double z = pz + root * vz;
      if (z < zMin || z > zMax)
      {
        continue;
      }
      const G4double normalDotV = rSign * (a * root + 0.5 * b);
      if (normalDotV >= 0.)
      {
        continue;
      }
      const G4double rHit = radiusAt(z, radius[0], slope);
      if (rHit < DBL_MIN)
      {
        continue;
      }
      const G4double x = px + root * vx;
      const G4double y = py + root * vy;
      const G4double rho2 = x * x + y * y;
      if (!inner)
      {
        const G4double rInner = std::max(0., radiusAt(z, rMin[0], innerSlope) - tolerance);
        if (rho2 < rInner * rInner)
        {
          continue;
        }
      }
      else
      {
        const G4double rOuter = radiusAt(z, rMax[0], outerSlope) + tolerance;
        if (rho2 > rOuter * rOuter)
        {
          continue;
        }
      }
      best = std::max(0., root);
    }
  };

  tryZPlane(z0, -1.);
  tryZPlane(z1, 1.);
  tryCone(rMax, outerSlope, false);
  tryCone(rMin, innerSlope, true);
  return best;
}

G4double G4Polycone::DistanceToRZ(const G4ThreeVector& p, const G4ThreeVector& v,
                                  G4bool outgoing, G4bool calcNorm, G4ThreeVector* normal,
                                  G4bool* validNorm) const
{
  if (fRZSection.empty() || fRZCone.size() != fRZSection.size())
  {
    return kInfinity;
  }

  const G4double tolerance = kCarTolerance / 2;
  G4double distance = kInfinity;
  G4double distFromSurface = kInfinity;
  G4ThreeVector bestNormal;
  G4bool allBehind = true;
  const G4bool useVisited = fRZVoxel != nullptr && !fRZVoxel->candidate.empty();
  G4StackVisited checked(useVisited ? G4int(fRZSection.size()) : 0);

  auto checkSection = [&](G4int i)
  {
    if (useVisited)
    {
      if (checked.TestAndSet(i))
      {
        return;
      }
    }
    const auto& sec = fRZSection[i];
    if (std::fabs(v.z()) < DBL_MIN)
    {
      if (p.z() < sec.zMin - tolerance || p.z() > sec.zMax + tolerance)
      {
        return;
      }
    }
    else
    {
      G4double z1 = (sec.zMin - tolerance - p.z()) / v.z();
      G4double z2 = (sec.zMax + tolerance - p.z()) / v.z();
      if (z1 > z2)
      {
        std::swap(z1, z2);
      }
      if (z2 < -tolerance || z1 > distance + tolerance)
      {
        return;
      }
    }
    G4double faceDistance = kInfinity;
    G4double faceDistFromSurface = kInfinity;
    G4ThreeVector faceNormal;
    if (IntersectRZSection(i, p, v, outgoing, tolerance, faceDistance, faceDistFromSurface,
                           calcNorm ? &faceNormal : nullptr))
    {
      if ((distance < kInfinity) || (!sec.allBehind))
      {
        allBehind = false;
      }
      if (faceDistance < distance)
      {
        distance = faceDistance;
        distFromSurface = faceDistFromSurface;
        if (calcNorm)
        {
          bestNormal = faceNormal;
        }
      }
    }
  };

  auto checkPhi = [&]()
  {
    G4double phiDistance = distance;
    G4double phiDistFromSurface = kInfinity;
    G4ThreeVector phiNormal;
    if (IntersectPhiPlane(p, v, outgoing, tolerance, phiDistance, phiDistFromSurface,
                          calcNorm ? &phiNormal : nullptr)
        && phiDistance < distance)
    {
      distance = phiDistance;
      distFromSurface = phiDistFromSurface;
      if (calcNorm)
      {
        bestNormal = phiNormal;
      }
      allBehind = false;
    }
  };

  if (fRZVoxel != nullptr && !fRZVoxel->candidate.empty() && fRZVoxel->boundary.size() >= 2)
  {
    checkPhi();
    if (!calcNorm && distFromSurface <= 0.)
    {
      return 0.;
    }
    for (auto i : fRZVoxel->global)
    {
      checkSection(i);
    }
    if (!calcNorm && distFromSurface <= 0.)
    {
      return 0.;
    }

    if (std::fabs(v.z()) < DBL_MIN)
    {
      for (auto i : FindZCandidates(fRZVoxel->boundary, fRZVoxel->candidate, p.z()))
      {
        checkSection(i);
      }
      if (!calcNorm && distFromSurface <= 0.)
      {
        return 0.;
      }
    }
    else
    {
      const auto upper = std::upper_bound(fRZVoxel->boundary.begin(),
                                          fRZVoxel->boundary.end(), p.z());
      G4int slab = 0;
      if (upper == fRZVoxel->boundary.begin())
      {
        slab = 0;
      }
      else if (upper == fRZVoxel->boundary.end())
      {
        slab = G4int(fRZVoxel->boundary.size()) - 2;
      }
      else
      {
        slab = G4int((upper - fRZVoxel->boundary.begin()) - 1);
      }

      const G4int nSlab = G4int(fRZVoxel->boundary.size()) - 1;
      const G4int step = (v.z() > 0.) ? 1 : -1;
      for (; slab >= 0 && slab < nSlab; slab += step)
      {
        const G4double enterZ = (step > 0)
                              ? fRZVoxel->boundary[slab]
                              : fRZVoxel->boundary[slab + 1];
        const G4double enterDistance = (enterZ - p.z()) / v.z();
        if (enterDistance > distance + tolerance)
        {
          break;
        }

        for (auto i : fRZVoxel->candidate[slab])
        {
          checkSection(i);
        }
        if (!calcNorm && distFromSurface <= 0.)
        {
          return 0.;
        }

        const G4double exitZ = (step > 0)
                             ? fRZVoxel->boundary[slab + 1]
                             : fRZVoxel->boundary[slab];
        const G4double exitDistance = (exitZ - p.z()) / v.z();
        if (exitDistance > distance + tolerance)
        {
          break;
        }
      }
    }
  }
  else
  {
    checkPhi();
    if (!calcNorm && distFromSurface <= 0.)
    {
      return 0.;
    }
    for (G4int i = 0; i < G4int(fRZSection.size()); ++i)
    {
      checkSection(i);
      if (!calcNorm && distFromSurface <= 0.)
      {
        return 0.;
      }
    }
  }

  if (distance < kInfinity)
  {
    if (distFromSurface <= 0.)
    {
      distance = 0.;
    }
    else if (distFromSurface < tolerance
          && DistanceToRZ(p, outgoing) < tolerance)
    {
      distance = 0.;
    }
    if (normal != nullptr)
    {
      *normal = bestNormal;
    }
    if (validNorm != nullptr)
    {
      *validNorm = allBehind;
    }
  }
  return distance;
}

// SurfaceNormalRZ
//
G4ThreeVector G4Polycone::SurfaceNormalRZ(const G4ThreeVector& p) const
{
  const G4double rho = p.perp();
  const G4double z = p.z();
  G4double best = kInfinity;
  const RZSection* bestSection = nullptr;
  auto check = [&](G4int i)
  {
    const auto& sec = fRZSection[i];
    G4double dzBound = 0.;
    if (z < sec.zMin)
    {
      dzBound = sec.zMin - z;
    }
    else if (z > sec.zMax)
    {
      dzBound = z - sec.zMax;
    }
    if (dzBound * dzBound > best + kCarTolerance)
    {
      return;
    }
    const G4double dr = rho - sec.r0;
    const G4double dz = z - sec.z0;
    const G4double q = std::max(0., std::min(sec.length,
                                             dr * sec.rS + dz * sec.zS));
    const G4double rr = dr - q * sec.rS;
    const G4double zz = dz - q * sec.zS;
    const G4double dist2 = rr * rr + zz * zz;
    if (dist2 < best)
    {
      best = dist2;
      bestSection = &sec;
    }
  };
  if (fRZVoxel != nullptr && !fRZVoxel->candidate.empty())
  {
    for (auto i : fRZVoxel->global)
    {
      check(i);
    }
    for (auto i : FindZCandidates(fRZVoxel->boundary, fRZVoxel->candidate, z))
    {
      check(i);
    }
  }
  if (bestSection == nullptr)
  {
    for (G4int i = 0; i < G4int(fRZSection.size()); ++i)
    {
      check(i);
    }
  }
  if (bestSection == nullptr)
  {
    return G4VCSGfaceted::SurfaceNormal(p);
  }
  G4ThreeVector phiNormal;
  const G4double phiDistance = DistanceToPhiBoundary(p, &phiNormal);
  if (phiDistance * phiDistance < best)
  {
    return phiNormal;
  }
  if (rho > 0.)
  {
    return {bestSection->rNorm * p.x() / rho, bestSection->rNorm * p.y() / rho,
            bestSection->zNorm};
  }
  return G4ThreeVector(0., 0., bestSection->zNorm).unit();
}

// Inside
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
EInside G4Polycone::Inside(const G4ThreeVector& p) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->MustBeOutside(p))
  {
    return kOutside;
  }

  //
  // Long answer
  //
  return fRZSection.empty() ? G4VCSGfaceted::Inside(p) : InsideRZ(p);
}

// SurfaceNormal
//
G4ThreeVector G4Polycone::SurfaceNormal(const G4ThreeVector& p) const
{
  return fRZSection.empty() ? G4VCSGfaceted::SurfaceNormal(p)
                            : SurfaceNormalRZ(p);
}

// DistanceToIn
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
G4double G4Polycone::DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->ShouldMiss(p, v))
  {
    return kInfinity;
  }

  //
  // Long answer
  //
  if (!fRZSection.empty())
  {
    if (!phiIsOpen && fZPlanesMonotonic && original_parameters != nullptr
        && original_parameters->Num_z_planes == 2)
    {
      if (const G4double d = DistanceToFullPhiTwoPlaneNoNormal(p, v, false);
          d < kInfinity)
      {
        return d;
      }
    }
    if (const G4double d = DistanceToRZ(p, v, false); d < kInfinity)
    {
      return d;
    }
  }
  return G4VCSGfaceted::DistanceToIn(p, v);
}

// DistanceToIn
//
G4double G4Polycone::DistanceToIn(const G4ThreeVector& p) const
{
  if (fRZSection.empty())
  {
    return G4VCSGfaceted::DistanceToIn(p);
  }
  if (!phiIsOpen)
  {
    return DistanceToRZ(p, false);
  }
  const EInside phi = PhiState(p);
  if (phi != kOutside)
  {
    return DistanceToRZ(p, false);
  }
  const G4double phiSafety = DistanceToPhiBoundary(p);
  const EInside rz = InsideZPlane(p);
  if (rz != kOutside)
  {
    return phiSafety;
  }
  return std::max(DistanceToRZ(p, false), phiSafety);
}

// DistanceToOut
//
G4double G4Polycone::DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm, G4bool* validNorm,
                                   G4ThreeVector* n) const
{
  if (!fRZSection.empty())
  {
    G4ThreeVector normal;
    G4bool rzValidNorm = false;
    G4double d = kInfinity;
    if (!calcNorm && !phiIsOpen && fZPlanesMonotonic && original_parameters != nullptr
        && original_parameters->Num_z_planes == 2)
    {
      d = DistanceToFullPhiTwoPlaneNoNormal(p, v, true);
    }
    if (d == kInfinity)
    {
      d = DistanceToRZ(p, v, true, calcNorm,
                       calcNorm ? &normal : nullptr, &rzValidNorm);
    }
    if (d < kInfinity)
    {
      if (calcNorm)
      {
        if (validNorm != nullptr)
        {
          *validNorm = rzValidNorm;
        }
        if (n != nullptr)
        {
          *n = normal;
        }
      }
      return d;
    }
  }
  return G4VCSGfaceted::DistanceToOut(p, v, calcNorm, validNorm, n);
}

// DistanceToOut
//
G4double G4Polycone::DistanceToOut(const G4ThreeVector& p) const
{
  if (fRZSection.empty())
  {
    return G4VCSGfaceted::DistanceToOut(p);
  }
  if (!phiIsOpen)
  {
    return DistanceToRZ(p, true);
  }
  return std::min(DistanceToRZ(p, true), DistanceToPhiBoundary(p));
}

// Get bounding box
//
void G4Polycone::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double rmin = kInfinity, rmax = -kInfinity;
  G4double zmin = kInfinity, zmax = -kInfinity;

  for (G4int i = 0; i < GetNumRZCorner(); ++i)
  {
    G4PolyconeSideRZ corner = GetCorner(i);
    if (corner.r < rmin)
    {
      rmin = corner.r;
    }
    if (corner.r > rmax)
    {
      rmax = corner.r;
    }
    if (corner.z < zmin)
    {
      zmin = corner.z;
    }
    if (corner.z > zmax)
    {
      zmax = corner.z;
    }
  }

  if (IsOpen())
  {
    G4TwoVector vmin, vmax;
    G4GeomTools::DiskExtent(rmin, rmax, GetSinStartPhi(), GetCosStartPhi(),
                            GetSinEndPhi(), GetCosEndPhi(), vmin, vmax);
    pMin.set(vmin.x(), vmin.y(), zmin);
    pMax.set(vmax.x(), vmax.y(), zmax);
  }
  else
  {
    pMin.set(-rmax, -rmax, zmin);
    pMax.set(rmax, rmax, zmax);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " << GetName() << " !"
            << "\npMin = " << pMin << "\npMax = " << pMax;
    G4Exception("G4Polycone::BoundingLimits()",
                "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

// Calculate extent under transform and specified limit
//
G4bool G4Polycone::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform, G4double& pMin,
                                   G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin, bmax);
  G4BoundingEnvelope bbox(bmin, bmax);
#  ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
#  endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis, pVoxelLimit, pTransform, pMin, pMax))
  {
    return exist = pMin < pMax;
  }

  // To find the extent, RZ contour of the polycone is subdivided
  // in triangles. The extent is calculated as cumulative extent of
  // all sub-polycones formed by rotation of triangles around Z
  //
  G4TwoVectorList contourRZ;
  G4TwoVectorList triangles;
  std::vector<G4int> iout;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  // get RZ contour, ensure anticlockwise order of corners
  for (G4int i = 0; i < GetNumRZCorner(); ++i)
  {
    G4PolyconeSideRZ corner = GetCorner(i);
    contourRZ.emplace_back(corner.r, corner.z);
  }
  G4GeomTools::RemoveRedundantVertices(contourRZ, iout, 2 * kCarTolerance);
  G4double area = G4GeomTools::PolygonArea(contourRZ);
  if (area < 0.)
  {
    std::reverse(contourRZ.begin(), contourRZ.end());
  }

  // triangulate RZ countour
  if (!G4GeomTools::TriangulatePolygon(contourRZ, triangles))
  {
    std::ostringstream message;
    message << "Triangulation of RZ contour has failed for solid: " << GetName() << " !"
            << "\nExtent has been calculated using boundary box";
    G4Exception("G4Polycone::CalculateExtent()", "GeomMgt1002", JustWarning, message);
    return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
  }

  // set trigonometric values
  const G4int NSTEPS = 24;  // number of steps for whole circle
  G4double astep = twopi / NSTEPS;  // max angle for one step

  G4double sphi = GetStartPhi();
  G4double ephi = GetEndPhi();
  G4double dphi = IsOpen() ? ephi - sphi : twopi;
  G4int ksteps = (dphi <= astep) ? 1 : (G4int)((dphi - deg) / astep) + 1;
  G4double ang = dphi / ksteps;

  G4double sinHalf = std::sin(0.5 * ang);
  G4double cosHalf = std::cos(0.5 * ang);
  G4double sinStep = 2. * sinHalf * cosHalf;
  G4double cosStep = 1. - 2. * sinHalf * sinHalf;

  G4double sinStart = GetSinStartPhi();
  G4double cosStart = GetCosStartPhi();
  G4double sinEnd = GetSinEndPhi();
  G4double cosEnd = GetCosEndPhi();

  // define vectors and arrays
  std::vector<const G4ThreeVectorList*> polygons;
  polygons.resize(ksteps + 2);
  G4ThreeVectorList pols[NSTEPS + 2];
  for (G4int k = 0; k < ksteps + 2; ++k)
  {
    pols[k].resize(6);
  }
  for (G4int k = 0; k < ksteps + 2; ++k)
  {
    polygons[k] = &pols[k];
  }
  G4double r0[6], z0[6];  // contour with original edges of triangle
  G4double r1[6];  // shifted radii of external edges of triangle

  // main loop along triangles
  pMin = kInfinity;
  pMax = -kInfinity;
  G4int ntria = (G4int)triangles.size() / 3;
  for (G4int i = 0; i < ntria; ++i)
  {
    G4int i3 = i * 3;
    for (G4int k = 0; k < 3; ++k)
    {
      G4int e0 = i3 + k, e1 = (k < 2) ? e0 + 1 : i3;
      G4int k2 = k * 2;
      // set contour with original edges of triangle
      r0[k2 + 0] = triangles[e0].x();
      z0[k2 + 0] = triangles[e0].y();
      r0[k2 + 1] = triangles[e1].x();
      z0[k2 + 1] = triangles[e1].y();
      // set shifted radii
      r1[k2 + 0] = r0[k2 + 0];
      r1[k2 + 1] = r0[k2 + 1];
      if (z0[k2 + 1] - z0[k2 + 0] <= 0)
      {
        continue;
      }
      r1[k2 + 0] /= cosHalf;
      r1[k2 + 1] /= cosHalf;
    }

    // rotate countour, set sequence of 6-sided polygons
    G4double sinCur = sinStart * cosHalf + cosStart * sinHalf;
    G4double cosCur = cosStart * cosHalf - sinStart * sinHalf;
    for (G4int j = 0; j < 6; ++j)
    {
      pols[0][j].set(r0[j] * cosStart, r0[j] * sinStart, z0[j]);
    }
    for (G4int k = 1; k < ksteps + 1; ++k)
    {
      for (G4int j = 0; j < 6; ++j)
      {
        pols[k][j].set(r1[j] * cosCur, r1[j] * sinCur, z0[j]);
      }
      G4double sinTmp = sinCur;
      sinCur = sinCur * cosStep + cosCur * sinStep;
      cosCur = cosCur * cosStep - sinTmp * sinStep;
    }
    for (G4int j = 0; j < 6; ++j)
    {
      pols[ksteps + 1][j].set(r0[j] * cosEnd, r0[j] * sinEnd, z0[j]);
    }

    // set sub-envelope and adjust extent
    G4double emin, emax;
    G4BoundingEnvelope benv(polygons);
    if (!benv.CalculateExtent(pAxis, pVoxelLimit, pTransform, emin, emax))
    {
      continue;
    }
    if (emin < pMin)
    {
      pMin = emin;
    }
    if (emax > pMax)
    {
      pMax = emax;
    }
    if (eminlim > pMin && emaxlim < pMax)
    {
      return true;
    }  // max possible extent
  }
  return (pMin < pMax);
}

// ComputeDimensions
//
void G4Polycone::ComputeDimensions(G4VPVParameterisation* p, const G4int n,
                                   const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this, n, pRep);
}

// GetEntityType
//
G4GeometryType G4Polycone::GetEntityType() const
{
  return {"G4Polycone"};
}

// Make a clone of the object
//
G4VSolid* G4Polycone::Clone() const
{
  return new G4Polycone(*this);
}

//
// Stream object contents to an output stream
//
std::ostream& G4Polycone::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Polycone\n"
     << " Parameters: \n"
     << "    starting phi angle : " << startPhi / degree << " degrees \n"
     << "    ending phi angle   : " << endPhi / degree << " degrees \n";
  G4int i = 0;

  G4int numPlanes = original_parameters->Num_z_planes;
  os << "    number of Z planes: " << numPlanes << "\n"
     << "              Z values: \n";
  for (i = 0; i < numPlanes; ++i)
  {
    os << "              Z plane " << i << ": " << original_parameters->Z_values[i] << "\n";
  }
  os << "              Tangent distances to inner surface (Rmin): \n";
  for (i = 0; i < numPlanes; ++i)
  {
    os << "              Z plane " << i << ": " << original_parameters->Rmin[i] << "\n";
  }
  os << "              Tangent distances to outer surface (Rmax): \n";
  for (i = 0; i < numPlanes; ++i)
  {
    os << "              Z plane " << i << ": " << original_parameters->Rmax[i] << "\n";
  }

  os << "    number of RZ points: " << numCorner << "\n"
     << "              RZ values (corners): \n";
  for (i = 0; i < numCorner; ++i)
  {
    os << "                         " << corners[i].r << ", " << corners[i].z << "\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

// Return volume
//
G4double G4Polycone::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&polyconeMutex);
    G4double total = 0.;
    G4int nrz = GetNumRZCorner();
    G4PolyconeSideRZ a = GetCorner(nrz - 1);
    for (G4int i = 0; i < nrz; ++i)
    {
      G4PolyconeSideRZ b = GetCorner(i);
      total += (b.r * b.r + b.r * a.r + a.r * a.r) * (b.z - a.z);
      a = b;
    }
    fCubicVolume = std::abs(total) * (GetEndPhi() - GetStartPhi()) / 6.;
    l.unlock();
  }
  return fCubicVolume;
}

// Return surface area
//
G4double G4Polycone::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&polyconeMutex);
    // phi cut area
    G4int nrz = GetNumRZCorner();
    G4double scut = 0.;
    if (IsOpen())
    {
      G4PolyconeSideRZ a = GetCorner(nrz - 1);
      for (G4int i = 0; i < nrz; ++i)
      {
        G4PolyconeSideRZ b = GetCorner(i);
        scut += a.r * b.z - a.z * b.r;
        a = b;
      }
      scut = std::abs(scut);
    }
    // lateral surface area
    G4double slat = 0;
    G4PolyconeSideRZ a = GetCorner(nrz - 1);
    for (G4int i = 0; i < nrz; ++i)
    {
      G4PolyconeSideRZ b = GetCorner(i);
      G4double h = std::sqrt((b.r - a.r) * (b.r - a.r) + (b.z - a.z) * (b.z - a.z));
      slat += (b.r + a.r) * h;
      a = b;
    }
    slat *= (GetEndPhi() - GetStartPhi()) / 2.;
    fSurfaceArea = scut + slat;
    l.unlock();
  }
  return fSurfaceArea;
}

// Set vector of surface elements, auxiliary method for sampling
// random points on surface
//
void G4Polycone::SetSurfaceElements() const
{
  fElements = new std::vector<G4Polycone::surface_element>;
  G4double total = 0.;
  G4int nrz = GetNumRZCorner();

  // set lateral surface elements
  G4double dphi = GetEndPhi() - GetStartPhi();
  G4int ia = nrz - 1;
  for (G4int ib = 0; ib < nrz; ++ib)
  {
    G4PolyconeSideRZ a = GetCorner(ia);
    G4PolyconeSideRZ b = GetCorner(ib);
    G4Polycone::surface_element selem;
    selem.i0 = ia;
    selem.i1 = ib;
    selem.i2 = -1;
    ia = ib;
    if (a.r == 0. && b.r == 0.)
    {
      continue;
    }
    G4double h = std::sqrt((b.r - a.r) * (b.r - a.r) + (b.z - a.z) * (b.z - a.z));
    total += 0.5 * dphi * (b.r + a.r) * h;
    selem.area = total;
    fElements->push_back(selem);
  }

  // set elements for phi cuts
  if (IsOpen())
  {
    G4TwoVectorList contourRZ;
    std::vector<G4int> triangles;
    for (G4int i = 0; i < nrz; ++i)
    {
      G4PolyconeSideRZ corner = GetCorner(i);
      contourRZ.emplace_back(corner.r, corner.z);
    }
    G4GeomTools::TriangulatePolygon(contourRZ, triangles);
    auto ntria = (G4int)triangles.size();
    for (G4int i = 0; i < ntria; i += 3)
    {
      G4Polycone::surface_element selem;
      selem.i0 = triangles[i];
      selem.i1 = triangles[i + 1];
      selem.i2 = triangles[i + 2];
      G4PolyconeSideRZ a = GetCorner(selem.i0);
      G4PolyconeSideRZ b = GetCorner(selem.i1);
      G4PolyconeSideRZ c = GetCorner(selem.i2);
      G4double stria = std::abs(G4GeomTools::TriangleArea(a.r, a.z, b.r, b.z, c.r, c.z));
      total += stria;
      selem.area = total;
      fElements->push_back(selem);  // start phi
      total += stria;
      selem.area = total;
      selem.i0 += nrz;
      fElements->push_back(selem);  // end phi
    }
  }
}

// Generate random point on surface
//
G4ThreeVector G4Polycone::GetPointOnSurface() const
{
  // Set surface elements
  if (fElements == nullptr)
  {
    G4AutoLock l(&surface_elementsMutex);
    if (fElements == nullptr)
    {
      SetSurfaceElements();
    }
    l.unlock();
  }

  // Select surface element
  G4Polycone::surface_element selem;
  selem = fElements->back();
  G4double select = selem.area * G4QuickRand();
  auto it = std::lower_bound(fElements->begin(), fElements->end(), select,
    [](const G4Polycone::surface_element& x, G4double val) -> G4bool
  { return x.area < val; });

  // Generate random point
  G4double r = 0, z = 0, phi = 0;
  G4double u = G4QuickRand();
  G4double v = G4QuickRand();
  G4int i0 = (*it).i0;
  G4int i1 = (*it).i1;
  G4int i2 = (*it).i2;
  if (i2 < 0)  // lateral surface
  {
    G4PolyconeSideRZ p0 = GetCorner(i0);
    G4PolyconeSideRZ p1 = GetCorner(i1);
    if (p1.r < p0.r)
    {
      p0 = GetCorner(i1);
      p1 = GetCorner(i0);
    }
    if (p1.r - p0.r < kCarTolerance)  // cylindrical surface
    {
      r = (p1.r - p0.r) * u + p0.r;
      z = (p1.z - p0.z) * u + p0.z;
    }
    else  // conical surface
    {
      r = std::sqrt(p1.r * p1.r * u + p0.r * p0.r * (1. - u));
      z = p0.z + (p1.z - p0.z) * (r - p0.r) / (p1.r - p0.r);
    }
    phi = (GetEndPhi() - GetStartPhi()) * v + GetStartPhi();
  }
  else  // phi cut
  {
    G4int nrz = GetNumRZCorner();
    phi = (i0 < nrz) ? GetStartPhi() : GetEndPhi();
    if (i0 >= nrz)
    {
      i0 -= nrz;
    }
    G4PolyconeSideRZ p0 = GetCorner(i0);
    G4PolyconeSideRZ p1 = GetCorner(i1);
    G4PolyconeSideRZ p2 = GetCorner(i2);
    if (u + v > 1.)
    {
      u = 1. - u;
      v = 1. - v;
    }
    r = (p1.r - p0.r) * u + (p2.r - p0.r) * v + p0.r;
    z = (p1.z - p0.z) * u + (p2.z - p0.z) * v + p0.z;
  }
  return {r * std::cos(phi), r * std::sin(phi), z};
}

// CreatePolyhedron
//
G4Polyhedron* G4Polycone::CreatePolyhedron() const
{
  std::vector<G4TwoVector> rz(numCorner);
  for (G4int i = 0; i < numCorner; ++i)
  {
    rz[i].set(corners[i].r, corners[i].z);
  }
  return new G4PolyhedronPcon(startPhi, endPhi - startPhi, rz);
}

// SetOriginalParameters
//
G4bool G4Polycone::SetOriginalParameters(G4ReduciblePolygon* rz)
{
  G4int numPlanes = numCorner;
  G4bool isConvertible = true;
  G4double Zmax = rz->Bmax();
  rz->StartWithZMin();

  // Prepare vectors for storage
  //
  std::vector<G4double> Z;
  std::vector<G4double> Rmin;
  std::vector<G4double> Rmax;

  G4int countPlanes = 1;
  G4int icurr = 0;
  G4int icurl = 0;

  // first plane Z=Z[0]
  //
  Z.push_back(corners[0].z);
  G4double Zprev = Z[0];
  if (Zprev == corners[1].z)
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back(corners[1].r);
    icurr = 1;
  }
  else if (Zprev == corners[numPlanes - 1].z)
  {
    Rmin.push_back(corners[numPlanes - 1].r);
    Rmax.push_back(corners[0].r);
    icurl = numPlanes - 1;
  }
  else
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back(corners[0].r);
  }

  // next planes until last
  //
  G4int inextr = 0, inextl = 0;
  for (G4int i = 0; i < numPlanes - 2; ++i)
  {
    inextr = 1 + icurr;
    inextl = (icurl <= 0) ? numPlanes - 1 : icurl - 1;

    if ((static_cast<G4int>(corners[inextr].z >= Zmax)
       & static_cast<G4int>(corners[inextl].z >= Zmax)) != 0)
    {
      break;
    }

    G4double Zleft = corners[inextl].z;
    G4double Zright = corners[inextr].z;
    if (Zright > Zleft)  // Next plane will be Zleft
    {
      Z.push_back(Zleft);
      countPlanes++;
      G4double difZr = corners[inextr].z - corners[icurr].z;
      G4double difZl = corners[inextl].z - corners[icurl].z;

      if (std::fabs(difZl) < kCarTolerance)
      {
        if (std::fabs(difZr) < kCarTolerance)
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r
                         + (Zleft - corners[icurr].z) / difZr
                             * (corners[inextr].r - corners[icurr].r));
        }
      }
      else if (difZl >= kCarTolerance)
      {
        if (std::fabs(difZr) < kCarTolerance)
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r
                         + (Zleft - corners[icurr].z) / difZr
                             * (corners[inextr].r - corners[icurr].r));
        }
      }
      else
      {
        isConvertible = false;
        break;
      }
      icurl = (icurl == 0) ? numPlanes - 1 : icurl - 1;
    }
    else if (std::fabs(Zright - Zleft) < kCarTolerance)  // Zright=Zleft
    {
      Z.push_back(Zleft);
      ++countPlanes;
      ++icurr;

      icurl = (icurl == 0) ? numPlanes - 1 : icurl - 1;

      Rmin.push_back(corners[inextl].r);
      Rmax.push_back(corners[inextr].r);
    }
    else  // Zright<Zleft
    {
      Z.push_back(Zright);
      ++countPlanes;

      G4double difZr = corners[inextr].z - corners[icurr].z;
      G4double difZl = corners[inextl].z - corners[icurl].z;
      if (std::fabs(difZr) < kCarTolerance)
      {
        if (std::fabs(difZl) < kCarTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[icurl].r
                         + (Zright - corners[icurl].z) / difZl
                             * (corners[inextl].r - corners[icurl].r));
          Rmax.push_back(corners[inextr].r);
        }
        ++icurr;
      }  // plate
      else if (difZr >= kCarTolerance)
      {
        if (std::fabs(difZl) < kCarTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurr].r);
        }
        else
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurl].r
                         + (Zright - corners[icurl].z) / difZl
                             * (corners[inextl].r - corners[icurl].r));
        }
        ++icurr;
      }
      else
      {
        isConvertible = false;
        break;
      }
    }
  }  // end for loop

  // last plane Z=Zmax
  //
  Z.push_back(Zmax);
  ++countPlanes;
  inextr = 1 + icurr;
  inextl = (icurl <= 0) ? numPlanes - 1 : icurl - 1;

  Rmax.push_back(corners[inextr].r);
  Rmin.push_back(corners[inextl].r);

  // Set original parameters Rmin,Rmax,Z
  //
  if (isConvertible)
  {
    original_parameters = new G4PolyconeHistorical;
    original_parameters->Z_values = new G4double[countPlanes];
    original_parameters->Rmin = new G4double[countPlanes];
    original_parameters->Rmax = new G4double[countPlanes];

    for (G4int j = 0; j < countPlanes; ++j)
    {
      original_parameters->Z_values[j] = Z[j];
      original_parameters->Rmax[j] = Rmax[j];
      original_parameters->Rmin[j] = Rmin[j];
    }
    original_parameters->Start_angle = startPhi;
    original_parameters->Opening_angle = endPhi - startPhi;
    original_parameters->Num_z_planes = countPlanes;
  }
  else  // Set parameters(r,z) with Rmin==0 as convention
  {
#  ifdef G4SPECSDEBUG
    std::ostringstream message;
    message << "Polycone " << GetName() << G4endl
            << "cannot be converted to Polycone with (Rmin,Rmaz,Z) parameters!";
    G4Exception("G4Polycone::SetOriginalParameters()", "GeomSolids0002", JustWarning, message);
#  endif
    original_parameters = new G4PolyconeHistorical;
    original_parameters->Z_values = new G4double[numPlanes];
    original_parameters->Rmin = new G4double[numPlanes];
    original_parameters->Rmax = new G4double[numPlanes];

    for (G4int j = 0; j < numPlanes; ++j)
    {
      original_parameters->Z_values[j] = corners[j].z;
      original_parameters->Rmax[j] = corners[j].r;
      original_parameters->Rmin[j] = 0.0;
    }
    original_parameters->Start_angle = startPhi;
    original_parameters->Opening_angle = endPhi - startPhi;
    original_parameters->Num_z_planes = numPlanes;
  }
  return isConvertible;
}

#endif
