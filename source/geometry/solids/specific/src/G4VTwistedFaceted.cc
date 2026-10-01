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
// G4VTwistedFaceted implementation
//
// Author: Oliver Link (CERN), 27.10.2004 - Created
// --------------------------------------------------------------------

#include "G4VTwistedFaceted.hh"

#include "G4AffineTransform.hh"
#include "G4AutoLock.hh"
#include "G4BoundingEnvelope.hh"
#include "G4ClippablePolygon.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polyhedron.hh"
#include "G4QuickRand.hh"
#include "G4SolidExtentList.hh"
#include "G4SystemOfUnits.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VTwistSurfaceWorkspace.hh"
#include "G4VisExtent.hh"
#include "G4VoxelLimits.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace
{
  G4Mutex vtwMutex = G4MUTEX_INITIALIZER;
}

//=====================================================================
//* constructors ------------------------------------------------------

G4VTwistedFaceted::G4VTwistedFaceted(const G4String& pname,  // Name of instance
                                     G4double PhiTwist,  // twist angle
                                     G4double pDz,  // half z length
                                     G4double pTheta,  // direction between end planes
                                     G4double pPhi,  // defined by polar and azim. angles
                                     G4double pDy1,  // half y length at -pDz
                                     G4double pDx1,  // half x length at -pDz,-pDy
                                     G4double pDx2,  // half x length at -pDz,+pDy
                                     G4double pDy2,  // half y length at +pDz
                                     G4double pDx3,  // half x length at +pDz,-pDy
                                     G4double pDx4,  // half x length at +pDz,+pDy
                                     G4double pAlph)  // tilt angle
  : G4VSolid(pname)
{
  G4double pDytmp;
  G4double fDxUp;
  G4double fDxDown;

  fDx1 = pDx1;
  fDx2 = pDx2;
  fDx3 = pDx3;
  fDx4 = pDx4;
  fDy1 = pDy1;
  fDy2 = pDy2;
  fDz = pDz;

  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  // maximum values
  //
  fDxDown = (fDx1 > fDx2 ? fDx1 : fDx2);
  fDxUp = (fDx3 > fDx4 ? fDx3 : fDx4);
  fDx = (fDxUp > fDxDown ? fDxUp : fDxDown);
  fDy = (fDy1 > fDy2 ? fDy1 : fDy2);

  // planarity check
  //
  if (fDx1 != fDx2 && fDx3 != fDx4)
  {
    pDytmp = fDy1 * (fDx3 - fDx4) / (fDx1 - fDx2);
    if (std::fabs(pDytmp - fDy2) > kCarTolerance)
    {
      std::ostringstream message;
      message << "Not planar surface in untwisted Trapezoid: " << GetName() << G4endl << "fDy2 is "
              << fDy2 << " but should be " << pDytmp << ".";
      G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "GeomSolids0002", FatalErrorInArgument,
                  message);
    }
  }

#ifdef G4TWISTDEBUG
  if (fDx1 == fDx2 && fDx3 == fDx4)
  {
    G4cout << "Trapezoid is a box" << G4endl;
  }

#endif

  if ((fDx1 == fDx2 && fDx3 != fDx4) || (fDx1 != fDx2 && fDx3 == fDx4))
  {
    std::ostringstream message;
    message << "Not planar surface in untwisted Trapezoid: " << GetName() << G4endl
            << "One endcap is rectangular, the other is a trapezoid." << G4endl
            << "For planarity reasons they have to be rectangles or trapezoids "
            << "on both sides.";
    G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "GeomSolids0002", FatalErrorInArgument,
                message);
  }

  // twist angle
  //
  fPhiTwist = PhiTwist;

  // tilt angle
  //
  fAlph = pAlph;
  fTAlph = std::tan(fAlph);

  fTheta = pTheta;
  fPhi = pPhi;

  // dx in surface equation
  //
  fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi);

  // dy in surface equation
  //
  fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi);

  if ((fDx1 <= 2 * kCarTolerance) || (fDx2 <= 2 * kCarTolerance) || (fDx3 <= 2 * kCarTolerance)
      || (fDx4 <= 2 * kCarTolerance) || (fDy1 <= 2 * kCarTolerance) || (fDy2 <= 2 * kCarTolerance)
      || (fDz <= 2 * kCarTolerance) || (std::fabs(fPhiTwist) <= 2 * kAngTolerance)
      || (std::fabs(fPhiTwist) >= pi / 2) || (std::fabs(fAlph) >= pi / 2) || fTheta >= pi / 2
      || fTheta < 0)
  {
    std::ostringstream message;
    message << "Invalid dimensions. Too small, or twist angle too big: " << GetName() << G4endl
            << "fDx 1-4 = " << fDx1 / cm << ", " << fDx2 / cm << ", " << fDx3 / cm << ", "
            << fDx4 / cm << " cm" << G4endl << "fDy 1-2 = " << fDy1 / cm << ", " << fDy2 / cm
            << ", "
            << " cm" << G4endl << "fDz = " << fDz / cm << " cm" << G4endl << " twistangle "
            << fPhiTwist / deg << " deg" << G4endl << " phi,theta = " << fPhi / deg << ", "
            << fTheta / deg << " deg";
    G4Exception("G4TwistedTrap::G4VTwistedFaceted()", "GeomSolids0002", FatalErrorInArgument,
                message);
  }
  CreateSurfaces();
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4VTwistedFaceted::G4VTwistedFaceted(__void__& a)
  : G4VSolid(a)
{
}

//=====================================================================
//* destructor --------------------------------------------------------

G4VTwistedFaceted::~G4VTwistedFaceted()
{
  DeleteSurfaces();
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

//=====================================================================
//* Copy constructor --------------------------------------------------

G4VTwistedFaceted::G4VTwistedFaceted(const G4VTwistedFaceted& rhs)
  : G4VSolid(rhs),
    fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea),
    fTheta(rhs.fTheta),
    fPhi(rhs.fPhi),
    fDy1(rhs.fDy1),
    fDx1(rhs.fDx1),
    fDx2(rhs.fDx2),
    fDy2(rhs.fDy2),
    fDx3(rhs.fDx3),
    fDx4(rhs.fDx4),
    fDz(rhs.fDz),
    fDx(rhs.fDx),
    fDy(rhs.fDy),
    fAlph(rhs.fAlph),
    fTAlph(rhs.fTAlph),
    fdeltaX(rhs.fdeltaX),
    fdeltaY(rhs.fdeltaY),
    fPhiTwist(rhs.fPhiTwist)
{
  CreateSurfaces();
}

//=====================================================================
//* Assignment operator -----------------------------------------------

G4VTwistedFaceted& G4VTwistedFaceted::operator=(const G4VTwistedFaceted& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  DeleteSurfaces();

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fTheta = rhs.fTheta;
  fPhi = rhs.fPhi;
  fDy1 = rhs.fDy1;
  fDx1 = rhs.fDx1;
  fDx2 = rhs.fDx2;
  fDy2 = rhs.fDy2;
  fDx3 = rhs.fDx3;
  fDx4 = rhs.fDx4;
  fDz = rhs.fDz;
  fDx = rhs.fDx;
  fDy = rhs.fDy;
  fAlph = rhs.fAlph;
  fTAlph = rhs.fTAlph;
  fdeltaX = rhs.fdeltaX;
  fdeltaY = rhs.fdeltaY;
  fPhiTwist = rhs.fPhiTwist;
  fCubicVolume = rhs.fCubicVolume;
  fSurfaceArea = rhs.fSurfaceArea;
  fRebuildPolyhedron = false;
  delete fpPolyhedron;
  fpPolyhedron = nullptr;

  CreateSurfaces();

  return *this;
}

//=====================================================================
//* Extent ------------------------------------------------------------

void G4VTwistedFaceted::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double cosPhi = std::cos(fPhi);
  G4double sinPhi = std::sin(fPhi);
  G4double tanTheta = std::tan(fTheta);
  G4double tanAlpha = fTAlph;

  G4double xmid1 = fDy1 * tanAlpha;
  G4double x1 = std::abs(xmid1 + fDx1);
  G4double x2 = std::abs(xmid1 - fDx1);
  G4double x3 = std::abs(xmid1 + fDx2);
  G4double x4 = std::abs(xmid1 - fDx2);
  G4double xmax1 = std::max(std::max(std::max(x1, x2), x3), x4);
  G4double rmax1 = std::sqrt(xmax1 * xmax1 + fDy1 * fDy1);

  G4double xmid2 = fDy2 * tanAlpha;
  G4double x5 = std::abs(xmid2 + fDx3);
  G4double x6 = std::abs(xmid2 - fDx3);
  G4double x7 = std::abs(xmid2 + fDx4);
  G4double x8 = std::abs(xmid2 - fDx4);
  G4double xmax2 = std::max(std::max(std::max(x5, x6), x7), x8);
  G4double rmax2 = std::sqrt(xmax2 * xmax2 + fDy2 * fDy2);

  G4double x0 = fDz * tanTheta * cosPhi;
  G4double y0 = fDz * tanTheta * sinPhi;
  G4double xmin = std::min(-x0 - rmax1, x0 - rmax2);
  G4double ymin = std::min(-y0 - rmax1, y0 - rmax2);
  G4double xmax = std::max(-x0 + rmax1, x0 + rmax2);
  G4double ymax = std::max(-y0 + rmax1, y0 + rmax2);
  pMin.set(xmin, ymin, -fDz);
  pMax.set(xmax, ymax, fDz);
}

//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool G4VTwistedFaceted::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                          const G4AffineTransform& pTransform, G4double& pMin,
                                          G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin, bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin, bmax);
  return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
}

//=====================================================================
//* Inside ------------------------------------------------------------

EInside G4VTwistedFaceted::Inside(const G4ThreeVector& p) const
{
  EInside tmpin = kOutside;
  const G4double halfTolerance = 0.5 * kCarTolerance;
  const G4double absZ = std::fabs(p.z());
  const G4bool withinZTolerance = absZ <= fDz + halfTolerance;

  // No point beyond the endcap tolerance can be inside the solid. Check this
  // before performing the untwisting transformation and its trigonometry.
  if (!withinZTolerance)
  {
    return kOutside;
  }
  const G4bool withinZInterior = absZ <= fDz - halfTolerance;

  G4double phi = p.z() / (2 * fDz) * fPhiTwist;  // rotate the point to z=0
  G4double cphi = std::cos(-phi);
  G4double sphi = std::sin(-phi);

  G4double px = p.x() + fdeltaX * (-phi / fPhiTwist);  // shift
  G4double py = p.y() + fdeltaY * (-phi / fPhiTwist);

  G4double posx = px * cphi - py * sphi;  // rotation
  G4double posy = px * sphi + py * cphi;

  G4double xMin = Xcoef(posy, phi, fTAlph) - 2 * Xcoef(posy, phi, 0.);
  G4double xMax = Xcoef(posy, phi, fTAlph);

  G4double yMax = GetValueB(phi) / 2.;  // b(phi)/2 is limit
  G4double yMin = -yMax;

#ifdef G4TWISTDEBUG

  G4cout << "inside called: p = " << p << G4endl;
  G4cout << "fDx1 = " << fDx1 << G4endl;
  G4cout << "fDx2 = " << fDx2 << G4endl;
  G4cout << "fDx3 = " << fDx3 << G4endl;
  G4cout << "fDx4 = " << fDx4 << G4endl;

  G4cout << "fDy1 = " << fDy1 << G4endl;
  G4cout << "fDy2 = " << fDy2 << G4endl;

  G4cout << "fDz  = " << fDz << G4endl;

  G4cout << "Tilt angle alpha = " << fAlph << G4endl;
  G4cout << "phi,theta = " << fPhi << " , " << fTheta << G4endl;

  G4cout << "Twist angle = " << fPhiTwist << G4endl;

  G4cout << "posx = " << posx << G4endl;
  G4cout << "posy = " << posy << G4endl;
  G4cout << "xMin = " << xMin << G4endl;
  G4cout << "xMax = " << xMax << G4endl;
  G4cout << "yMin = " << yMin << G4endl;
  G4cout << "yMax = " << yMax << G4endl;

#endif

  if (posx <= xMax - halfTolerance && posx >= xMin + halfTolerance)
  {
    if (posy <= yMax - halfTolerance && posy >= yMin + halfTolerance)
    {
      tmpin = withinZInterior ? kInside : kSurface;
    }
    else if (posy <= yMax + halfTolerance && posy >= yMin - halfTolerance)
    {
      tmpin = kSurface;
    }
  }
  else if (posx <= xMax + halfTolerance && posx >= xMin - halfTolerance)
  {
    if (posy <= yMax + halfTolerance && posy >= yMin - halfTolerance)
    {
      tmpin = kSurface;
    }
  }

#ifdef G4TWISTDEBUG
  G4cout << "inside = " << tmpin << G4endl;
#endif

  return tmpin;
}

//=====================================================================
//* FindBoundarySurface -----------------------------------------------

G4int G4VTwistedFaceted::FindBoundarySurface(const G4ThreeVector& p) const
{
  const G4double halfTolerance = 0.5 * kCarTolerance;
  if (std::fabs(p.z()) > fDz + halfTolerance)
  {
    return -1;
  }

  const G4double phi = p.z() / (2 * fDz) * fPhiTwist;
  const G4double cphi = std::cos(-phi);
  const G4double sphi = std::sin(-phi);
  const G4double px = p.x() + fdeltaX * (-phi / fPhiTwist);
  const G4double py = p.y() + fdeltaY * (-phi / fPhiTwist);
  const G4double posx = px * cphi - py * sphi;
  const G4double posy = px * sphi + py * cphi;
  const G4double xMin = Xcoef(posy, phi, fTAlph) - 2 * Xcoef(posy, phi, 0.);
  const G4double xMax = Xcoef(posy, phi, fTAlph);
  const G4double yMax = GetValueB(phi) / 2.;
  const G4double yMin = -yMax;

  if (posx > xMax + halfTolerance || posx < xMin - halfTolerance || posy > yMax + halfTolerance
      || posy < yMin - halfTolerance)
  {
    return -1;
  }

  // A direct normal is used only for a single face whose equation is
  // satisfied to floating-point roundoff. A point that is merely in the
  // geometry-tolerance shell, or is close to an edge/corner, retains the
  // general closest-surface calculation and its ordered tie semantics.
  const G4double separations[6] = {posx - xMax, posy - yMax, posx - xMin,
                                   posy - yMin, p.z() + fDz, p.z() - fDz};
  constexpr G4double roundoffFactor = 128.;
  const G4double epsilon = std::numeric_limits<G4double>::epsilon();
  G4int toleranceBoundary = -1;
  G4int toleranceBoundaryCount = 0;
  for (G4int i = 0; i < 6; ++i)
  {
    if (std::fabs(separations[i]) <= halfTolerance)
    {
      toleranceBoundary = i;
      ++toleranceBoundaryCount;
    }
  }
  if (toleranceBoundaryCount != 1)
  {
    return -1;
  }

  G4double coordinate = p.z();
  G4double boundary = fDz;
  if (toleranceBoundary == 0 || toleranceBoundary == 2)
  {
    coordinate = posx;
    boundary = toleranceBoundary == 0 ? xMax : xMin;
  }
  else if (toleranceBoundary == 1 || toleranceBoundary == 3)
  {
    coordinate = posy;
    boundary = toleranceBoundary == 1 ? yMax : yMin;
  }
  const G4double scale = std::max({1., std::fabs(coordinate), std::fabs(boundary)});
  return std::fabs(separations[toleranceBoundary]) <= roundoffFactor * epsilon * scale
           ? toleranceBoundary
           : -1;
}

//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4VTwistedFaceted::SurfaceNormal(const G4ThreeVector& p) const
{
  const G4int surfaceIndex = FindBoundarySurface(p);
  if (surfaceIndex >= 0)
  {
    return GetSurface(surfaceIndex)->GetNormal(p, true);
  }
  return LegacySurfaceNormal(p);
}

//=====================================================================
//* LegacySurfaceNormal -----------------------------------------------

G4ThreeVector G4VTwistedFaceted::LegacySurfaceNormal(const G4ThreeVector& p) const
{
  //
  // return the normal unit vector to the Hyperbolical Surface at a point
  // p on (or nearly on) the surface
  //
  // Which of the three or four surfaces are we closest to?
  //

  G4double distance = kInfinity;

  G4ThreeVector xx;
  G4ThreeVector bestxx;
  G4int i;
  G4int besti = -1;
  for (i = 0; i < 6; i++)
  {
    auto* surface = GetSurface(i);
    G4double tmpdistance = surface->DistanceTo(p, xx);
    if (tmpdistance < distance)
    {
      // The first exact hit is also the ordered tie winner, so no later
      // surface can alter the result.
      if (tmpdistance == 0.)
      {
        return surface->GetNormal(xx, true);
      }
      distance = tmpdistance;
      bestxx = xx;
      besti = i;
    }
  }

  return GetSurface(besti)->GetNormal(bestxx, true);
}

//=====================================================================
//* DistanceToBoundaryFallback ----------------------------------------

G4double G4VTwistedFaceted::DistanceToBoundaryFallback(const G4ThreeVector& p,
                                                       const G4ThreeVector& v, G4bool entering,
                                                       G4bool finitePoint, G4int selectedFace,
                                                       G4double selectedLower,
                                                       G4double selectedUpper, G4bool calcNorm,
                                                       G4bool* validNorm, G4ThreeVector* norm) const
{
  G4VTwistSurface* surfaces[] = {fSide0, fSide90, fSide180, fSide270, fLowerEndcap, fUpperEndcap};
  G4VTwistSurface::ExactQueryWorkspace workspace;
  G4ThreeVector xx;

  const auto querySurface = [&](G4VTwistSurface* surface)
  {
    return entering ? surface->DistanceToIn(p, v, xx, workspace)
                    : surface->DistanceToOut(p, v, xx, workspace);
  };

  if (selectedFace >= 0)
  {
    const G4double selectedDistance = querySurface(surfaces[selectedFace]);
    if (selectedDistance >= std::max(0., selectedLower) && selectedDistance <= selectedUpper)
    {
      if (!entering && calcNorm)
      {
        *norm = surfaces[selectedFace]->GetNormal(xx, true);
        *validNorm = surfaces[selectedFace]->IsValidNorm();
      }
      return selectedDistance;
    }
  }

  G4double distance = kInfinity;
  G4int besti = -1;
  G4ThreeVector bestxx;
  for (G4int i = 0; i < 6; ++i)
  {
#ifdef G4TWISTDEBUG
    if (entering)
    {
      G4cout << G4endl << "surface " << i << ": " << G4endl << G4endl;
    }
#endif
    if (finitePoint && i < 4 && !MayCrossLateralFace(i, v, entering))
    {
      continue;
    }
    if (finitePoint && std::isfinite(v.z()))
    {
      // Flat endcaps have outward normals -Z and +Z, respectively.
      const G4bool pointsAway = entering ? ((i == 4 && v.z() <= 0.) || (i == 5 && v.z() >= 0.))
                                         : ((i == 4 && v.z() >= 0.) || (i == 5 && v.z() <= 0.));
      if (pointsAway)
      {
        continue;
      }
    }

    const G4double candidate = querySurface(surfaces[i]);
#ifdef G4TWISTDEBUG
    if (entering)
    {
      G4cout << "Solid DistanceToIn : distance = " << candidate << G4endl;
      G4cout << "intersection point = " << xx << G4endl;
    }
#endif
    if (candidate < distance)
    {
      distance = candidate;
      if (!entering)
      {
        bestxx = xx;
        besti = i;
      }
    }
  }

  if (!entering && calcNorm && besti != -1)
  {
    *norm = surfaces[besti]->GetNormal(bestxx, true);
    *validNorm = surfaces[besti]->IsValidNorm();
  }

#ifdef G4TWISTDEBUG
  if (entering)
  {
    G4cout << "best distance = " << distance << G4endl;
  }
#endif
  return distance;
}

//=====================================================================
//* DistanceToIn (p, v) -----------------------------------------------

G4double G4VTwistedFaceted::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  // DistanceToIn (p, v):
  // Calculate distance to surface of shape from `outside'
  // along with the v, allowing for tolerance.
  // The function returns kInfinity if no intersection or
  // just grazing within tolerance.

  //
  // Calculate DistanceToIn(p,v)
  //

  const EInside currentside = Inside(p);

  if (currentside == kSurface)
  {
    // particle is just on a boundary.
    // if the particle is entering to the volume, return 0
    //
    const G4int surfaceIndex = FindBoundarySurface(p);
    G4ThreeVector normal =
      surfaceIndex >= 0 ? GetSurface(surfaceIndex)->GetNormal(p, true) : LegacySurfaceNormal(p);
    if (normal * v < 0)
    {
      return 0;
    }
  }
  else if (currentside == kOutside)
  {
    G4ThreeVector boxMinimum;
    G4ThreeVector boxMaximum;
    if (fQueryCacheReady)
    {
      boxMinimum = fQueryBoxMin;
      boxMaximum = fQueryBoxMax;
    }
    else
    {
      BoundingLimits(boxMinimum, boxMaximum);
      const G4ThreeVector margin(0.5 * kCarTolerance, 0.5 * kCarTolerance, 0.5 * kCarTolerance);
      boxMinimum -= margin;
      boxMaximum += margin;
    }
    // The point is outside, and a ray missing the tolerance-expanded
    // conservative box cannot reach any constituent surface.
    if (!RayIntersectsBox(p, v, boxMinimum, boxMaximum))
    {
      return kInfinity;
    }
  }

  G4int selectedFace = -1;
  G4double selectedLower = 0.;
  G4double selectedUpper = 0.;
  G4double selectedSlopeLower = 0.;
  const G4bool finitePoint = IsFinite(p);
  if (currentside == kOutside && finitePoint)
  {
    selectedFace =
      SelectFirstTransition(p, v, true, selectedLower, selectedUpper, selectedSlopeLower);
    if (selectedFace >= 0 && selectedFace < 4)
    {
      G4double refinedDistance = kInfinity;
      G4ThreeVector refinedNormal;
      if (RefineFirstTransition(p, v, selectedFace, true, selectedLower, selectedUpper,
                                selectedSlopeLower, refinedDistance, refinedNormal))
      {
        return refinedDistance;
      }
    }
  }

  return DistanceToBoundaryFallback(p, v, true, finitePoint, selectedFace, selectedLower,
                                    selectedUpper, false, nullptr, nullptr);
}

//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4VTwistedFaceted::DistanceToIn(const G4ThreeVector& p) const
{
  // DistanceToIn(p):
  // Calculate distance to surface of shape from `outside',
  // allowing for tolerance
  //

  //
  // Calculate DistanceToIn(p)
  //

  const EInside currentside = Inside(p);

  switch (currentside)
  {
    case (kInside):
    {
    }
    case (kSurface):
    {
      return 0;
    }
    case (kOutside):
    {
      G4ThreeVector boxMinimum;
      G4ThreeVector boxMaximum;
      if (fQueryCacheReady)
      {
        boxMinimum = fQueryBoxMin;
        boxMaximum = fQueryBoxMax;
      }
      else
      {
        BoundingLimits(boxMinimum, boxMaximum);
        const G4ThreeVector margin(0.5 * kCarTolerance, 0.5 * kCarTolerance, 0.5 * kCarTolerance);
        boxMinimum -= margin;
        boxMaximum += margin;
      }
      const G4double boxDistance = DistanceToBox(p, boxMinimum, boxMaximum);
      // Outside the expanded enclosing box, its distance is already a safe
      // lower bound on the distance to the solid.
      if (boxDistance > 0.)
      {
        return boxDistance;
      }

      // Initialize
      //
      G4double distance = kInfinity;

      // Find intersections and choose nearest one
      //
      G4VTwistSurface* surfaces[] = { fSide0,   fSide90,      fSide180,
                                      fSide270, fLowerEndcap, fUpperEndcap };

      G4ThreeVector xx;
      for (const auto& surface : surfaces)
      {
        G4double tmpdistance = surface->DistanceTo(p, xx);
        if (tmpdistance < distance)
        {
          distance = tmpdistance;
        }
      }
      return distance;
    }
    default:
    {
      G4Exception("G4VTwistedFaceted::DistanceToIn(p)", "GeomSolids0003",
                  FatalException, "Unknown point location!");
    }
  }  // switch end

  return 0.;
}

//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double G4VTwistedFaceted::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                          const G4bool calcNorm, G4bool* validNorm,
                                          G4ThreeVector* norm) const
{
  // DistanceToOut (p, v):
  // Calculate distance to surface of shape from `inside'
  // along with the v, allowing for tolerance.
  // The function returns kInfinity if no intersection or
  // just grazing within tolerance.

  //
  // Calculate DistanceToOut(p,v)
  //

  const EInside currentside = Inside(p);

  if (currentside == kSurface)
  {
    // particle is just on a boundary.
    // if the particle is exiting from the volume, return 0
    //
    const G4int surfaceIndex = FindBoundarySurface(p);
    G4ThreeVector normal =
      surfaceIndex >= 0 ? GetSurface(surfaceIndex)->GetNormal(p, true) : LegacySurfaceNormal(p);
    if (normal * v > 0)
    {
      if (calcNorm)
      {
        *norm = normal;
        *validNorm = true;
      }
      // timer.Stop();
      return 0;
    }
  }

  G4int selectedFace = -1;
  G4double selectedLower = 0.;
  G4double selectedUpper = 0.;
  G4double selectedSlopeLower = 0.;
  const G4bool finitePoint = IsFinite(p);
  if (currentside == kInside && finitePoint)
  {
    selectedFace =
      SelectFirstTransition(p, v, false, selectedLower, selectedUpper, selectedSlopeLower);
    if (selectedFace >= 0 && selectedFace < 4)
    {
      G4double refinedDistance = kInfinity;
      G4ThreeVector refinedNormal;
      if (RefineFirstTransition(p, v, selectedFace, false, selectedLower, selectedUpper,
                                selectedSlopeLower, refinedDistance, refinedNormal))
      {
        if (calcNorm)
        {
          *norm = refinedNormal;
          *validNorm = GetSurface(selectedFace)->IsValidNorm();
        }
        return refinedDistance;
      }
    }
  }

  return DistanceToBoundaryFallback(p, v, false, finitePoint, selectedFace, selectedLower,
                                    selectedUpper, calcNorm, validNorm, norm);
}

//=====================================================================
//* DistanceToOut (p) -------------------------------------------------

G4double G4VTwistedFaceted::DistanceToOut(const G4ThreeVector& p) const
{
  // DistanceToOut(p):
  // Calculate distance to surface of shape from `inside',
  // allowing for tolerance

  //
  // Calculate DistanceToOut(p)
  //

  EInside currentside = Inside(p);
  G4double retval = kInfinity;

  switch (currentside)
  {
    case (kOutside):
    {
#ifdef G4SPECSDEBUG
      G4int oldprc = G4cout.precision(16);
      G4cout << G4endl;
      DumpInfo();
      G4cout << "Position:" << G4endl << G4endl;
      G4cout << "p.x() = " << p.x() / mm << " mm" << G4endl;
      G4cout << "p.y() = " << p.y() / mm << " mm" << G4endl;
      G4cout << "p.z() = " << p.z() / mm << " mm" << G4endl << G4endl;
      G4cout.precision(oldprc);
      G4Exception("G4VTwistedFaceted::DistanceToOut(p)", "GeomSolids1002",
                  JustWarning, "Point p is outside !?");
#endif
      break;
    }
    case (kSurface):
    {
      retval = 0;
      break;
    }
    case (kInside):
    {
      if (fQueryCacheReady && IsFinite(p))
      {
        const G4double untwist = -fUntwistRate * p.z();
        const G4double cosUntwist = std::cos(untwist);
        const G4double sinUntwist = std::sin(untwist);
        const G4double shiftedX = p.x() - fUntwistShiftX * p.z();
        const G4double shiftedY = p.y() - fUntwistShiftY * p.z();
        const G4ThreeVector q(shiftedX * cosUntwist - shiftedY * sinUntwist,
                              shiftedX * sinUntwist + shiftedY * cosUntwist, p.z());

        G4double planeMargin = fDz - std::fabs(p.z());
        for (const auto& face : fUntwistedFaces)
        {
          planeMargin = std::min(planeMargin, face.offset - face.normal.dot(q));
        }
        // The untwisting map is Lipschitz on every straight segment between
        // points of the solid. Thus this normalized plane margin is a proven
        // lower bound on the corresponding distance in world coordinates.
        retval = std::max(0., planeMargin - fUntwistSafetyRoundoff) * fUntwistInvLipschitz;
        break;
      }

      // Initialize
      //
      G4double distance = kInfinity;

      // find intersections and choose nearest one
      //
      G4VTwistSurface* surfaces[] = { fSide0,   fSide90,      fSide180,
                                      fSide270, fLowerEndcap, fUpperEndcap };

      G4ThreeVector xx;
      for (const auto& surface : surfaces)
      {
        G4double tmpdistance = surface->DistanceTo(p, xx);
        if (tmpdistance < distance)
        {
          distance = tmpdistance;
        }
      }
      retval = distance;
      break;
    }
    default:
    {
      G4Exception("G4VTwistedFaceted::DistanceToOut(p)", "GeomSolids0003",
                  FatalException, "Unknown point location!");
      break;
    }
  }  // switch end

  return retval;
}

//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4VTwistedFaceted::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4VTwistedFaceted\n"
     << " Parameters: \n"
     << "  polar angle theta = " << fTheta / degree << " deg" << G4endl
     << "  azimuthal angle phi = " << fPhi / degree << " deg" << G4endl
     << "  tilt angle  alpha = " << fAlph / degree << " deg" << G4endl
     << "  TWIST angle = " << fPhiTwist / degree << " deg" << G4endl
     << "  Half length along y (lower endcap) = " << fDy1 / cm << " cm" << G4endl
     << "  Half length along x (lower endcap, bottom) = " << fDx1 / cm << " cm" << G4endl
     << "  Half length along x (lower endcap, top) = " << fDx2 / cm << " cm" << G4endl
     << "  Half length along y (upper endcap) = " << fDy2 / cm << " cm" << G4endl
     << "  Half length along x (upper endcap, bottom) = " << fDx3 / cm << " cm" << G4endl
     << "  Half length along x (upper endcap, top) = " << fDx4 / cm << " cm" << G4endl
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4VTwistedFaceted::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4VTwistedFaceted::GetExtent() const
{
  G4double maxRad = std::sqrt(fDx * fDx + fDy * fDy);

  return {-maxRad, maxRad, -maxRad, maxRad, -fDz, fDz};
}

//=====================================================================
//* CreateSurfaces ----------------------------------------------------

void G4VTwistedFaceted::CreateSurfaces()
{
  // create 6 surfaces of TwistedTub.

  if (fDx1 == fDx2 && fDx3 == fDx4)  // special case : Box
  {
    fSide0 = new G4TwistBoxSide("0deg", fPhiTwist, fDz, fTheta, fPhi, fDy1, fDx1, fDx1, fDy2, fDx3,
                                fDx3, fAlph, 0. * deg);
    fSide180 = new G4TwistBoxSide("180deg", fPhiTwist, fDz, fTheta, fPhi + pi, fDy1, fDx1, fDx1,
                                  fDy2, fDx3, fDx3, fAlph, 180. * deg);
  }
  else  // default general case
  {
    fSide0 = new G4TwistTrapAlphaSide("0deg", fPhiTwist, fDz, fTheta, fPhi, fDy1, fDx1, fDx2, fDy2,
                                      fDx3, fDx4, fAlph, 0. * deg);
    fSide180 = new G4TwistTrapAlphaSide("180deg", fPhiTwist, fDz, fTheta, fPhi + pi, fDy1, fDx2,
                                        fDx1, fDy2, fDx4, fDx3, fAlph, 180. * deg);
  }

  // create parallel sides
  //
  fSide90 = new G4TwistTrapParallelSide("90deg", fPhiTwist, fDz, fTheta, fPhi, fDy1, fDx1, fDx2,
                                        fDy2, fDx3, fDx4, fAlph, 0. * deg);
  fSide270 = new G4TwistTrapParallelSide("270deg", fPhiTwist, fDz, fTheta, fPhi + pi, fDy1, fDx2,
                                         fDx1, fDy2, fDx4, fDx3, fAlph, 180. * deg);

  // create endcaps
  //
  fUpperEndcap =
    new G4TwistTrapFlatSide("UpperCap", fPhiTwist, fDx3, fDx4, fDy2, fDz, fAlph, fPhi, fTheta, 1);
  fLowerEndcap =
    new G4TwistTrapFlatSide("LowerCap", fPhiTwist, fDx1, fDx2, fDy1, fDz, fAlph, fPhi, fTheta, -1);

  // Set neighbour surfaces

  fSide0->SetNeighbours(fSide270, fLowerEndcap, fSide90, fUpperEndcap);
  fSide90->SetNeighbours(fSide0, fLowerEndcap, fSide180, fUpperEndcap);
  fSide180->SetNeighbours(fSide90, fLowerEndcap, fSide270, fUpperEndcap);
  fSide270->SetNeighbours(fSide180, fLowerEndcap, fSide0, fUpperEndcap);
  fUpperEndcap->SetNeighbours(fSide180, fSide270, fSide0, fSide90);
  fLowerEndcap->SetNeighbours(fSide180, fSide270, fSide0, fSide90);

  CreateQueryCache();
}

//=====================================================================
//* DeleteSurfaces ----------------------------------------------------

void G4VTwistedFaceted::DeleteSurfaces()
{
  delete fLowerEndcap;
  delete fUpperEndcap;
  delete fSide0;
  delete fSide90;
  delete fSide180;
  delete fSide270;

  fLowerEndcap = nullptr;
  fUpperEndcap = nullptr;
  fSide0 = nullptr;
  fSide90 = nullptr;
  fSide180 = nullptr;
  fSide270 = nullptr;
  fQueryCacheReady = false;
}

//=====================================================================
//* CreateQueryCache --------------------------------------------------

void G4VTwistedFaceted::CreateQueryCache()
{
  fQueryCacheReady = false;

  const G4double ta = fTAlph;
  const std::array<G4ThreeVector, 8> vertices = {
    G4ThreeVector(-fDx1 - fDy1 * ta, -fDy1, -fDz), G4ThreeVector(fDx1 - fDy1 * ta, -fDy1, -fDz),
    G4ThreeVector(-fDx2 + fDy1 * ta, fDy1, -fDz),  G4ThreeVector(fDx2 + fDy1 * ta, fDy1, -fDz),
    G4ThreeVector(-fDx3 - fDy2 * ta, -fDy2, fDz),  G4ThreeVector(fDx3 - fDy2 * ta, -fDy2, fDz),
    G4ThreeVector(-fDx4 + fDy2 * ta, fDy2, fDz),   G4ThreeVector(fDx4 + fDy2 * ta, fDy2, fDz)};

  // Vertex order follows the constituent surface order +X,+Y,-X,-Y.
  constexpr G4int faceVertices[4][4] = {{1, 3, 7, 5}, {2, 6, 7, 3}, {0, 4, 6, 2}, {0, 1, 5, 4}};
  const G4double epsilon = std::numeric_limits<G4double>::epsilon();
  G4double qRadiusMax = 0.;
  G4double planeResidualMax = 0.;
  for (const auto& vertex : vertices)
  {
    qRadiusMax = std::max(qRadiusMax, std::hypot(vertex.x(), vertex.y()));
  }

  for (G4int i = 0; i < 4; ++i)
  {
    const auto& q0 = vertices[faceVertices[i][0]];
    const auto& q1 = vertices[faceVertices[i][1]];
    const auto& q3 = vertices[faceVertices[i][3]];
    G4ThreeVector normal = (q1 - q0).cross(q3 - q0);
    const G4double normalMag = normal.mag();
    if (!(normalMag > 0.) || !std::isfinite(normalMag))
    {
      return;
    }
    normal /= normalMag;
    G4double offset = normal.dot(q0);
    if (offset < 0.)
    {
      normal = -normal;
      offset = -offset;
    }
    if (!(offset > 0.) || !std::isfinite(offset))
    {
      return;
    }

    G4double qJMin = kInfinity;
    G4double qJMax = -kInfinity;
    for (G4int j = 0; j < 4; ++j)
    {
      const auto& q = vertices[faceVertices[i][j]];
      const G4double planeResidual = normal.dot(q) - offset;
      const G4double scale = std::max({1., std::fabs(offset), q.mag()});
      if (std::fabs(planeResidual) > 4096. * epsilon * scale)
      {
        // Near-planar input accepted within geometry tolerance retains the
        // established navigation rather than using a plane-based bound.
        return;
      }
      planeResidualMax = std::max(planeResidualMax, std::fabs(planeResidual));
      const G4double qJ = -normal.x() * q.y() + normal.y() * q.x();
      qJMin = std::min(qJMin, qJ);
      qJMax = std::max(qJMax, qJ);
    }
    fUntwistedFaces[i] = {normal, offset, qJMin, qJMax};
  }

  G4double safetyScale = std::max({1., std::fabs(fDz), qRadiusMax});
  for (const auto& face : fUntwistedFaces)
  {
    safetyScale = std::max(safetyScale, std::fabs(face.offset));
  }
  fUntwistSafetyRoundoff = std::nextafter(
    std::max(512. * epsilon * safetyScale,
             planeResidualMax + 512. * epsilon * (planeResidualMax + safetyScale + 1.)),
    kInfinity);

  fUntwistRate = fPhiTwist / (2. * fDz);
  fUntwistShiftX = fdeltaX / (2. * fDz);
  fUntwistShiftY = fdeltaY / (2. * fDz);
  const G4double halfTwist = 0.5 * std::fabs(fPhiTwist);
  fSinHalfTwist = std::sin(halfTwist);
  fCosHalfTwist = std::cos(halfTwist);

  const G4double shiftRate = std::hypot(fUntwistShiftX, fUntwistShiftY);
  const G4double shearBound = std::fabs(fUntwistRate) * qRadiusMax + shiftRate;
  const G4double lipschitzSquared =
    1. + 0.5 * shearBound * shearBound + 0.5 * shearBound * std::sqrt(shearBound * shearBound + 4.);
  if (!(lipschitzSquared >= 1.) || !std::isfinite(lipschitzSquared))
  {
    return;
  }
  fUntwistInvLipschitz = std::nextafter((1. - 512. * epsilon) / std::sqrt(lipschitzSquared), 0.);
  BoundingLimits(fQueryBoxMin, fQueryBoxMax);
  const G4ThreeVector boxMargin(0.5 * kCarTolerance, 0.5 * kCarTolerance, 0.5 * kCarTolerance);
  fQueryBoxMin -= boxMargin;
  fQueryBoxMax += boxMargin;

  // Certify that the analytic plane orientation agrees with each concrete
  // side implementation before it is allowed to suppress a root solve.
  for (G4int i = 0; i < 4; ++i)
  {
    G4ThreeVector q(0., 0., 0.);
    for (G4int j = 0; j < 4; ++j)
    {
      q += vertices[faceVertices[i][j]];
    }
    q *= 0.25;
    const G4double angle = fUntwistRate * q.z();
    const G4double ca = std::cos(angle);
    const G4double sa = std::sin(angle);
    const G4ThreeVector point(ca * q.x() - sa * q.y() + fUntwistShiftX * q.z(),
                              sa * q.x() + ca * q.y() + fUntwistShiftY * q.z(), q.z());
    const auto& face = fUntwistedFaces[i];
    const G4double worldNx = ca * face.normal.x() - sa * face.normal.y();
    const G4double worldNy = sa * face.normal.x() + ca * face.normal.y();
    const G4double qJ = -face.normal.x() * q.y() + face.normal.y() * q.x();
    const G4ThreeVector expected(worldNx, worldNy,
                                 face.normal.z() - worldNx * fUntwistShiftX
                                   - worldNy * fUntwistShiftY - fUntwistRate * qJ);
    const G4ThreeVector actual = GetSurface(i)->GetNormal(point, true);
    if (!IsFinite(actual) || expected.unit().dot(actual) < 1. - 1024. * epsilon)
    {
      return;
    }
  }

  fQueryCacheReady = true;
}

//=====================================================================
//* MayCrossLateralFace -----------------------------------------------

G4bool G4VTwistedFaceted::MayCrossLateralFace(G4int surfaceIndex, const G4ThreeVector& v,
                                              G4bool entering) const
{
  if (!fQueryCacheReady || surfaceIndex < 0 || surfaceIndex >= 4 || !IsFinite(v))
  {
    return true;
  }

  const auto& face = fUntwistedFaces[surfaceIndex];
  const G4double wx = v.x() - fUntwistShiftX * v.z();
  const G4double wy = v.y() - fUntwistShiftY * v.z();
  const G4double a = face.normal.x() * wx + face.normal.y() * wy;
  const G4double b = face.normal.x() * wy - face.normal.y() * wx;

  // phi is confined to [-|twist|/2,+|twist|/2], where sine is monotonic
  // and cosine lies in [cos(|twist|/2),1]. Independent interval terms make
  // this deliberately conservative; uncertain faces are always evaluated.
  const G4double trigLower = (a >= 0. ? a * fCosHalfTwist : a) - std::fabs(b) * fSinHalfTwist;
  const G4double trigUpper = (a >= 0. ? a : a * fCosHalfTwist) + std::fabs(b) * fSinHalfTwist;
  const G4double qJCoefficient = -fUntwistRate * v.z();
  const G4double qJ0 = qJCoefficient * face.qJMin;
  const G4double qJ1 = qJCoefficient * face.qJMax;
  const G4double zTerm = face.normal.z() * v.z();
  G4double lower = trigLower + std::min(qJ0, qJ1) + zTerm;
  G4double upper = trigUpper + std::max(qJ0, qJ1) + zTerm;

  const G4double magnitude = std::fabs(trigLower) + std::fabs(trigUpper)
                             + std::max(std::fabs(qJ0), std::fabs(qJ1)) + std::fabs(zTerm) + 1.;
  const G4double guard = 256. * std::numeric_limits<G4double>::epsilon() * magnitude;
  lower = std::nextafter(lower - guard, -kInfinity);
  upper = std::nextafter(upper + guard, kInfinity);
  if (!std::isfinite(lower) || !std::isfinite(upper))
  {
    return true;
  }
  return entering ? lower < 0. : upper > 0.;
}

//=====================================================================
//* SelectFirstTransition ---------------------------------------------

G4int G4VTwistedFaceted::SelectFirstTransition(const G4ThreeVector& p, const G4ThreeVector& v,
                                               G4bool entering, G4double& lower, G4double& upper,
                                               G4double& slopeLower) const
{
  const G4double speed2 = v.mag2();
  const G4double epsilon = std::numeric_limits<G4double>::epsilon();
  if (!fQueryCacheReady || !(fUntwistInvLipschitz > 0.) || !IsFinite(p) || !IsFinite(v)
      || !(speed2 > 0.) || !std::isfinite(speed2) || std::fabs(speed2 - 1.) > 256. * epsilon)
  {
    return -1;
  }

  // Intersect an outward-rounded, tolerance-expanded cached box.  This is
  // used only to bound a proof; failure or uncertainty retains the general
  // surface loop and can therefore never suppress a physical intersection.
  const G4ThreeVector boxMinimum(std::nextafter(fQueryBoxMin.x(), -kInfinity),
                                 std::nextafter(fQueryBoxMin.y(), -kInfinity),
                                 std::nextafter(fQueryBoxMin.z(), -kInfinity));
  const G4ThreeVector boxMaximum(std::nextafter(fQueryBoxMax.x(), kInfinity),
                                 std::nextafter(fQueryBoxMax.y(), kInfinity),
                                 std::nextafter(fQueryBoxMax.z(), kInfinity));
  G4double boxEnter = -kInfinity;
  G4double boxLeave = kInfinity;
  const auto updateSlab = [&](G4double coordinate, G4double direction, G4double minimum,
                              G4double maximum)
  {
    if (direction == 0.)
    {
      return coordinate >= minimum && coordinate <= maximum;
    }
    G4double first = (minimum - coordinate) / direction;
    G4double last = (maximum - coordinate) / direction;
    if (first > last)
    {
      std::swap(first, last);
    }
    first = std::nextafter(first, -kInfinity);
    last = std::nextafter(last, kInfinity);
    if (!std::isfinite(first) || !std::isfinite(last))
    {
      return false;
    }
    boxEnter = std::max(boxEnter, first);
    boxLeave = std::min(boxLeave, last);
    return boxEnter <= boxLeave;
  };
  if (!updateSlab(p.x(), v.x(), boxMinimum.x(), boxMaximum.x())
      || !updateSlab(p.y(), v.y(), boxMinimum.y(), boxMaximum.y())
      || !updateSlab(p.z(), v.z(), boxMinimum.z(), boxMaximum.z()))
  {
    return -1;
  }

  G4double intervalLower = entering ? std::max(0., boxEnter) : 0.;
  G4double intervalUpper = boxLeave;
  G4bool endsAtCap = false;
  G4int capIndex = -1;
  G4double capDistance = kInfinity;
  if (!entering && v.z() > 0.)
  {
    capDistance = (fDz - p.z()) / v.z();
    capIndex = 5;
  }
  else if (!entering && v.z() < 0.)
  {
    capDistance = (-fDz - p.z()) / v.z();
    capIndex = 4;
  }
  if (!entering && capDistance >= 0. && capDistance <= intervalUpper && std::isfinite(capDistance))
  {
    intervalUpper = capDistance;
    endsAtCap = true;
  }
  if (!(intervalUpper > intervalLower) || !std::isfinite(intervalLower)
      || !std::isfinite(intervalUpper))
  {
    return -1;
  }

  // Entry rays may originate arbitrarily far from the solid.  Re-anchor the
  // analytic functions at the beginning of the bounded box interval so that
  // coefficient formation and subsequent a0+a1*t evaluation cannot amplify
  // a large cancellation involving the original point.
  const G4double anchor = entering ? intervalLower : 0.;
  const G4double intervalWidth = intervalUpper - anchor;
  const G4ThreeVector origin(std::fma(anchor, v.x(), p.x()), std::fma(anchor, v.y(), p.y()),
                             std::fma(anchor, v.z(), p.z()));
  if (!(intervalWidth > 0.) || !std::isfinite(intervalWidth) || !IsFinite(origin))
  {
    return -1;
  }
  G4double reanchorUncertainty = 0.;
  if (entering && anchor > 0.)
  {
    const G4double cancellationScale = std::max({1., std::fabs(p.x()) + std::fabs(anchor * v.x()),
                                                 std::fabs(p.y()) + std::fabs(anchor * v.y()),
                                                 std::fabs(p.z()) + std::fabs(anchor * v.z())});
    reanchorUncertainty = std::nextafter(64. * epsilon * cancellationScale, kInfinity);
    // The slab endpoint itself was obtained by floating-point division.
    // FMA removes the final multiply-add rounding but not that endpoint's
    // uncertainty.  Refuse launches for which it is material geometrically.
    if (!std::isfinite(reanchorUncertainty) || reanchorUncertainty > 0.25 * kCarTolerance)
    {
      return -1;
    }
  }

  const G4double rx0 = origin.x() - fUntwistShiftX * origin.z();
  const G4double ry0 = origin.y() - fUntwistShiftY * origin.z();
  const G4double wx = v.x() - fUntwistShiftX * v.z();
  const G4double wy = v.y() - fUntwistShiftY * v.z();
  const G4double phi0 = fUntwistRate * origin.z();
  const G4double omega = fUntwistRate * v.z();
  const G4double shiftedPositionScale =
    std::fabs(origin.x()) + std::fabs(fUntwistShiftX) * std::fabs(origin.z())
    + std::fabs(origin.y()) + std::fabs(fUntwistShiftY) * std::fabs(origin.z())
    + std::fabs(origin.z());
  const G4double shiftedDirectionScale =
    std::fabs(v.x()) + std::fabs(fUntwistShiftX) * std::fabs(v.z()) + std::fabs(v.y())
    + std::fabs(fUntwistShiftY) * std::fabs(v.z()) + std::fabs(v.z());
  const G4double coefficientScale =
    shiftedPositionScale + intervalWidth * shiftedDirectionScale + 1.;
  const G4double coefficientUncertainty =
    std::nextafter(2048. * epsilon * coefficientScale, kInfinity);
  if (!std::isfinite(coefficientUncertainty))
  {
    return -1;
  }
  std::array<RayFaceFunction, 6> faces{};
  for (G4int i = 0; i < 4; ++i)
  {
    const auto& face = fUntwistedFaces[i];
    faces[i] = { face.normal.x() * rx0 + face.normal.y() * ry0,
                 face.normal.x() * wx + face.normal.y() * wy,
                 face.normal.x() * ry0 - face.normal.y() * rx0,
                 face.normal.x() * wy - face.normal.y() * wx,
                 face.normal.z() * origin.z() - face.offset,
                 face.normal.z() * v.z(),
                 phi0,
                 omega };
  }
  // Endcaps use the same <=0 interior convention as the lateral planes.
  faces[4] = { 0., 0., 0., 0., -origin.z() - fDz, -v.z(), phi0, omega };
  faces[5] = { 0., 0., 0., 0., origin.z() - fDz, v.z(), phi0, omega };

  const G4double speed = std::sqrt(speed2);
  G4double lipschitzUpper = 1. / fUntwistInvLipschitz;
  lipschitzUpper =
    std::nextafter(lipschitzUpper + 512. * epsilon * (std::fabs(lipschitzUpper) + 1.), kInfinity);
  G4double evaluationGuard =
    fUntwistSafetyRoundoff + coefficientUncertainty + reanchorUncertainty * lipschitzUpper;
  evaluationGuard =
    std::nextafter(evaluationGuard + 512. * epsilon * (std::fabs(evaluationGuard) + 1.), kInfinity);
  G4double faceGuard = 4. * kCarTolerance * lipschitzUpper + evaluationGuard;
  faceGuard = std::nextafter(faceGuard + 512. * epsilon * (std::fabs(faceGuard) + 1.), kInfinity);
  G4double parallelGuard = std::sin(0.001) * lipschitzUpper * speed;
  parallelGuard =
    std::nextafter(parallelGuard + 512. * epsilon * (std::fabs(parallelGuard) + 1.), kInfinity);
  if (!std::isfinite(lipschitzUpper) || !std::isfinite(evaluationGuard) || !std::isfinite(faceGuard)
      || !std::isfinite(parallelGuard))
  {
    return -1;
  }

  G4int nodeCount = 0;
  const G4int faceCount = entering ? 6 : 4;
  const G4double lateralUpperLimit =
    !entering && endsAtCap ? std::nextafter(capDistance - 8. * kCarTolerance / speed, -kInfinity)
                           : kInfinity;
  std::array<FaceRange, 6> ranges{};
  SelectedTransition selected;
  const TransitionSearchStatus status = SearchFirstTransition(
    faces, faceCount, entering, 0., intervalWidth, faceGuard, parallelGuard, evaluationGuard,
    0.25 * kCarTolerance / speed, lateralUpperLimit, 0, nodeCount, ranges, selected);
  if (status == TransitionSearchStatus::kSelected)
  {
    lower = std::nextafter(anchor + selected.lower, -kInfinity);
    upper = std::nextafter(anchor + selected.upper, kInfinity);
    slopeLower = selected.slopeLower;
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(slopeLower > 0.)
        || !std::isfinite(slopeLower))
    {
      return -1;
    }
    return selected.face;
  }
  if (!entering && status == TransitionSearchStatus::kNoTransition && endsAtCap)
  {
    lower = std::nextafter(capDistance, -kInfinity);
    upper = std::nextafter(capDistance, kInfinity);
    return capIndex;
  }
  return -1;
}

//=====================================================================
//* RefineFirstTransition ---------------------------------------------

G4bool G4VTwistedFaceted::RefineFirstTransition(const G4ThreeVector& p, const G4ThreeVector& v,
                                                G4int surfaceIndex, G4bool entering, G4double lower,
                                                G4double upper, G4double selectedSlopeLower,
                                                G4double& distance, G4ThreeVector& normal) const
{
  // Flat caps intentionally retain their established exact implementation.
  const G4double epsilon = std::numeric_limits<G4double>::epsilon();
  const G4double speed2 = v.mag2();
  if (!fQueryCacheReady || surfaceIndex < 0 || surfaceIndex >= 4 || !(upper > lower)
      || !std::isfinite(lower) || !std::isfinite(upper) || !(selectedSlopeLower > 0.)
      || !std::isfinite(selectedSlopeLower) || !(speed2 > 0.) || !std::isfinite(speed2)
      || std::fabs(speed2 - 1.) > 256. * epsilon)
  {
    return false;
  }

  const G4double anchor = std::max(0., lower);
  const G4double width = upper - anchor;
  const G4ThreeVector origin(std::fma(anchor, v.x(), p.x()), std::fma(anchor, v.y(), p.y()),
                             std::fma(anchor, v.z(), p.z()));
  if (!(width > 0.) || !std::isfinite(width) || !IsFinite(origin))
  {
    return false;
  }

  const G4double cancellationScale = std::max({1., std::fabs(p.x()) + std::fabs(anchor * v.x()),
                                               std::fabs(p.y()) + std::fabs(anchor * v.y()),
                                               std::fabs(p.z()) + std::fabs(anchor * v.z())});
  const G4double reanchorUncertainty = std::nextafter(64. * epsilon * cancellationScale, kInfinity);
  if (!std::isfinite(reanchorUncertainty) || reanchorUncertainty > 0.25 * kCarTolerance)
  {
    return false;
  }

  const G4double rx0 = std::fma(-fUntwistShiftX, origin.z(), origin.x());
  const G4double ry0 = std::fma(-fUntwistShiftY, origin.z(), origin.y());
  const G4double wx = std::fma(-fUntwistShiftX, v.z(), v.x());
  const G4double wy = std::fma(-fUntwistShiftY, v.z(), v.y());
  const auto& untwistedFace = fUntwistedFaces[surfaceIndex];
  const RayFaceFunction face = {
    std::fma(untwistedFace.normal.x(), rx0, untwistedFace.normal.y() * ry0),
    std::fma(untwistedFace.normal.x(), wx, untwistedFace.normal.y() * wy),
    std::fma(untwistedFace.normal.x(), ry0, -untwistedFace.normal.y() * rx0),
    std::fma(untwistedFace.normal.x(), wy, -untwistedFace.normal.y() * wx),
    std::fma(untwistedFace.normal.z(), origin.z(), -untwistedFace.offset),
    untwistedFace.normal.z() * v.z(),
    fUntwistRate * origin.z(),
    fUntwistRate * v.z() };

  const G4double shiftedPositionScale =
    std::fabs(origin.x()) + std::fabs(fUntwistShiftX) * std::fabs(origin.z())
    + std::fabs(origin.y()) + std::fabs(fUntwistShiftY) * std::fabs(origin.z())
    + std::fabs(origin.z());
  const G4double shiftedDirectionScale =
    std::fabs(v.x()) + std::fabs(fUntwistShiftX) * std::fabs(v.z()) + std::fabs(v.y())
    + std::fabs(fUntwistShiftY) * std::fabs(v.z()) + std::fabs(v.z());
  const G4double coefficientScale = shiftedPositionScale + width * shiftedDirectionScale + 1.;
  // These re-anchored coefficients each use only a short, fused operation
  // chain.  Keep a separate bound from the broader selector envelope so the
  // final root-error proof is not dominated by unrelated faces.
  const G4double coefficientUncertainty =
    std::nextafter(128. * epsilon * coefficientScale, kInfinity);
  G4double lipschitzUpper = 1. / fUntwistInvLipschitz;
  lipschitzUpper =
    std::nextafter(lipschitzUpper + 512. * epsilon * (std::fabs(lipschitzUpper) + 1.), kInfinity);
  G4double evaluationGuard =
    fUntwistSafetyRoundoff + coefficientUncertainty + reanchorUncertainty * lipschitzUpper;
  evaluationGuard =
    std::nextafter(evaluationGuard + 512. * epsilon * (std::fabs(evaluationGuard) + 1.), kInfinity);
  const G4double speed = std::sqrt(speed2);
  G4double parallelGuard = std::sin(0.001) * lipschitzUpper * speed;
  parallelGuard =
    std::nextafter(parallelGuard + 512. * epsilon * (std::fabs(parallelGuard) + 1.), kInfinity);
  if (!std::isfinite(coefficientUncertainty) || !std::isfinite(lipschitzUpper)
      || !std::isfinite(evaluationGuard) || !std::isfinite(parallelGuard))
  {
    return false;
  }

  const auto orient = [entering](const PointRange& range)
  {
    return entering ? PointRange{-range.upper, -range.lower} : range;
  };

  // Re-certify monotonicity after re-anchoring and outward widening of the
  // selected interval.  The smaller of the two independent lower bounds is
  // carried through the residual and transition checks below.
  const G4double midpoint = 0.5 * width;
  const G4double phase = std::fma(face.omega, midpoint, face.phi0);
  FaceRange fullRange{};
  if (!EvaluateInterval(face, midpoint, 0.5 * width, std::sin(phase), std::cos(phase),
                        evaluationGuard, fullRange))
  {
    return false;
  }
  const G4double reconstructedSlopeLower =
    entering ? -fullRange.derivativeUpper : fullRange.derivativeLower;
  const G4double slopeLower = std::min(selectedSlopeLower, reconstructedSlopeLower);
  if (!(slopeLower > parallelGuard) || !std::isfinite(slopeLower))
  {
    return false;
  }
  PointRange lowerFaceRange{};
  PointRange upperFaceRange{};
  G4double lowerValue = 0.;
  G4double upperValue = 0.;
  G4double endpointDerivative = 0.;
  if (!EvaluatePointWithCentral(face, 0., evaluationGuard, lowerFaceRange, lowerValue,
                                endpointDerivative)
      || !EvaluatePointWithCentral(face, width, evaluationGuard, upperFaceRange, upperValue,
                                   endpointDerivative))
  {
    return false;
  }
  const PointRange atLower = orient(lowerFaceRange);
  const PointRange atUpper = orient(upperFaceRange);
  if (!(atLower.upper <= 0.) || !(atUpper.lower >= 0.))
  {
    return false;
  }

  // H is oriented to be monotonically increasing for both entry and exit.
  // Raw values choose Newton proposals only; a bracket endpoint is changed
  // exclusively when its outward-expanded PointRange proves the sign.
  G4double bracketLower = 0.;
  G4double bracketUpper = width;
  const G4double orientedLowerValue = entering ? -lowerValue : lowerValue;
  const G4double orientedUpperValue = entering ? -upperValue : upperValue;
  G4double current = midpoint;
  const G4double secantDenominator = orientedUpperValue - orientedLowerValue;
  if (orientedLowerValue <= 0. && orientedUpperValue >= 0. && secantDenominator > 0.
      && std::isfinite(secantDenominator))
  {
    const G4double secant = (-orientedLowerValue / secantDenominator) * width;
    if (secant > bracketLower && secant < bracketUpper && std::isfinite(secant))
    {
      current = secant;
    }
  }
  PointRange currentRange{};
  G4bool stationary = false;
  constexpr G4int maxIterations = 12;
  for (G4int iteration = 0; iteration < maxIterations; ++iteration)
  {
    G4double rawValue = 0.;
    G4double rawDerivative = 0.;
    PointRange currentFaceRange{};
    if (!EvaluatePointWithCentral(face, current, evaluationGuard, currentFaceRange, rawValue,
                                  rawDerivative))
    {
      return false;
    }
    const G4double orientedValue = entering ? -rawValue : rawValue;
    const G4double orientedDerivative = entering ? -rawDerivative : rawDerivative;
    currentRange = orient(currentFaceRange);

    if (currentRange.upper <= 0.)
    {
      bracketLower = current;
    }
    else if (currentRange.lower >= 0.)
    {
      bracketUpper = current;
    }

    G4double next = bracketLower + 0.5 * (bracketUpper - bracketLower);
    if (orientedDerivative > 0. && std::isfinite(orientedDerivative))
    {
      const G4double newtonStep = current - orientedValue / orientedDerivative;
      if (newtonStep > bracketLower && newtonStep < bracketUpper && std::isfinite(newtonStep))
      {
        const G4double stationaryGuard =
          8. * epsilon
          * std::max({1., std::fabs(current), std::fabs(bracketLower), std::fabs(bracketUpper)});
        if (std::fabs(newtonStep - current) <= stationaryGuard)
        {
          stationary = true;
          break;
        }
        next = newtonStep;
      }
    }
    if (!(next > bracketLower && next < bracketUpper))
    {
      return false;
    }
    current = next;
  }
  if (!stationary)
  {
    return false;
  }

  const G4double residualMagnitude =
    std::max(std::fabs(currentRange.lower), std::fabs(currentRange.upper));
  const G4double residualParameterError = std::nextafter(residualMagnitude / slopeLower, kInfinity);
  distance = anchor + current;
  const G4double anchorSumParameterError =
    std::nextafter(32. * epsilon * (std::fabs(anchor) + std::fabs(current) + 1.), kInfinity);
  const G4double parameterError =
    std::nextafter(residualParameterError + anchorSumParameterError, kInfinity);
  const G4double spatialError = std::nextafter(parameterError * speed, kInfinity);
  // The exact-face and orientation checks below independently validate the
  // central root. Accept its conservative error envelope up to the established
  // half Cartesian surface tolerance.
  if (!std::isfinite(distance) || distance < lower || distance > upper
      || !std::isfinite(spatialError) || spatialError > 0.5 * kCarTolerance)
  {
    return false;
  }

  // The legacy lateral solver snaps a root within half a Cartesian tolerance
  // to zero.  Do not make that threshold decision from an approximate root.
  if (distance - parameterError <= 0.5 * kCarTolerance)
  {
    return false;
  }

  const G4ThreeVector intersection(std::fma(distance, v.x(), p.x()),
                                   std::fma(distance, v.y(), p.y()),
                                   std::fma(distance, v.z(), p.z()));
  if (!IsFinite(intersection) || FindBoundarySurface(intersection) != surfaceIndex)
  {
    return false;
  }

  normal = GetSurface(surfaceIndex)->GetNormal(intersection, true);
  const G4double normalDot = normal.dot(v);
  if (!IsFinite(normal) || !std::isfinite(normalDot)
      || (entering ? normalDot >= 0. : normalDot <= 0.))
  {
    return false;
  }

  return true;
}

//=====================================================================
//* GetLateralFaceArea ------------------------------------------------

G4double G4VTwistedFaceted::GetLateralFaceArea(const G4TwoVector& p1, const G4TwoVector& p2,
                                               const G4TwoVector& p3, const G4TwoVector& p4) const
{
  constexpr G4int NSTEP = 100;
  constexpr G4double dt = 1. / NSTEP;

  G4double h = 2. * fDz;
  G4double hh = h * h;
  G4double hTanTheta = h * std::tan(fTheta);
  G4double x1 = p1.x();
  G4double y1 = p1.y();
  G4double x21 = p2.x() - p1.x();
  G4double y21 = p2.y() - p1.y();
  G4double x31 = p3.x() - p1.x();
  G4double y31 = p3.y() - p1.y();
  G4double x42 = p4.x() - p2.x();
  G4double y42 = p4.y() - p2.y();
  G4double x43 = p4.x() - p3.x();
  G4double y43 = p4.y() - p3.y();

  // check if face is plane (just in case)
  G4double lmax =
    std::max(std::max(std::abs(x21), std::abs(y21)), std::max(std::abs(x43), std::abs(y43)));
  G4double eps = lmax * kCarTolerance;
  if (std::abs(fPhiTwist) < kCarTolerance && std::abs(x21 * y43 - y21 * x43) < eps)
  {
    G4double x0 = hTanTheta * std::cos(fPhi);
    G4double y0 = hTanTheta * std::sin(fPhi);
    G4ThreeVector A(p4.x() - p1.x() + x0, p4.y() - p1.y() + y0, h);
    G4ThreeVector B(p3.x() - p2.x() + x0, p3.y() - p2.y() + y0, h);
    return (A.cross(B)).mag() * 0.5;
  }

  // twisted face
  G4double area = 0.;
  for (G4int i = 0; i < NSTEP; ++i)
  {
    G4double t = (i + 0.5) * dt;
    G4double I = x21 + (x42 - x31) * t;
    G4double J = y21 + (y42 - y31) * t;
    G4double II = I * I;
    G4double JJ = J * J;
    G4double IIJJ = hh * (I * I + J * J);

    G4double ang = fPhi + fPhiTwist * (0.5 - t);
    G4double A = fPhiTwist * (II + JJ) + x21 * y43 - x43 * y21;
    G4double B = fPhiTwist * (I * (x1 + x31 * t) + J * (y1 + y31 * t))
                 + hTanTheta * (I * std::sin(ang) - J * std::cos(ang)) + (I * y31 - J * x31);

    G4double invAA = 1. / (A * A);
    G4double sqrtAA = 2. * std::abs(A);
    G4double invSqrtAA = 1. / sqrtAA;

    G4double aa = A * A;
    G4double bb = 2. * A * B;
    G4double cc = IIJJ + B * B;

    G4double R1 = std::sqrt(aa + bb + cc);
    G4double R0 = std::sqrt(cc);
    G4double log1 = std::log(std::abs(sqrtAA * R1 + 2. * aa + bb));
    G4double log0 = std::log(std::abs(sqrtAA * R0 + bb));

    area += 0.5 * R1 + 0.25 * bb * invAA * (R1 - R0) + IIJJ * invSqrtAA * (log1 - log0);
  }
  return area * dt;
}

//=====================================================================
//* GetCubicVolume ----------------------------------------------------

G4double G4VTwistedFaceted::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&vtwMutex);
    fCubicVolume = ((fDx1 + fDx2 + fDx3 + fDx4) * (fDy1 + fDy2)
                    + (fDx4 + fDx3 - fDx2 - fDx1) * (fDy2 - fDy1) / 3)
                   * fDz;
    l.unlock();
  }
  return fCubicVolume;
}

//=====================================================================
//* GetSurfaceArea ----------------------------------------------------

G4double G4VTwistedFaceted::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&vtwMutex);
    G4TwoVector vv[8];
    vv[0] = G4TwoVector(-fDx1 - fDy1 * fTAlph, -fDy1);
    vv[1] = G4TwoVector(fDx1 - fDy1 * fTAlph, -fDy1);
    vv[2] = G4TwoVector(-fDx2 + fDy1 * fTAlph, fDy1);
    vv[3] = G4TwoVector(fDx2 + fDy1 * fTAlph, fDy1);
    vv[4] = G4TwoVector(-fDx3 - fDy2 * fTAlph, -fDy2);
    vv[5] = G4TwoVector(fDx3 - fDy2 * fTAlph, -fDy2);
    vv[6] = G4TwoVector(-fDx4 + fDy2 * fTAlph, fDy2);
    vv[7] = G4TwoVector(fDx4 + fDy2 * fTAlph, fDy2);
    fSurfaceArea = 2. * (fDy1 * (fDx1 + fDx2) + fDy2 * (fDx3 + fDx4))
                   + GetLateralFaceArea(vv[0], vv[1], vv[4], vv[5])
                   + GetLateralFaceArea(vv[1], vv[3], vv[5], vv[7])
                   + GetLateralFaceArea(vv[3], vv[2], vv[7], vv[6])
                   + GetLateralFaceArea(vv[2], vv[0], vv[6], vv[4]);
    l.unlock();
  }
  return fSurfaceArea;
}

//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4VTwistedFaceted::GetEntityType() const
{
  return {"G4VTwistedFaceted"};
}

//=====================================================================
//* GetPolyhedron -----------------------------------------------------

G4Polyhedron* G4VTwistedFaceted::GetPolyhedron() const
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

//=====================================================================
//* GetPointOnSurface -------------------------------------------------

G4ThreeVector G4VTwistedFaceted::GetPointOnSurface() const
{
  G4double phi = fPhiTwist * (G4QuickRand() - 0.5);
  G4double u, umin, umax;  //  variable for twisted surfaces
  G4double y;  //  variable for flat surface (top and bottom)

  // Compute the areas. Attention: Only correct for trapezoids
  // where the twisting is done along the z-axis. In the general case
  // the computed surface area is more difficult. However this simplification
  // does not affect the tracking through the solid.

  G4double a1 = fSide0->GetSurfaceArea();
  G4double a2 = fSide90->GetSurfaceArea();
  G4double a3 = fSide180->GetSurfaceArea();
  G4double a4 = fSide270->GetSurfaceArea();
  G4double a5 = fLowerEndcap->GetSurfaceArea();
  G4double a6 = fUpperEndcap->GetSurfaceArea();

#ifdef G4TWISTDEBUG
  G4cout << "Surface 0   deg = " << a1 << G4endl;
  G4cout << "Surface 90  deg = " << a2 << G4endl;
  G4cout << "Surface 180 deg = " << a3 << G4endl;
  G4cout << "Surface 270 deg = " << a4 << G4endl;
  G4cout << "Surface Lower   = " << a5 << G4endl;
  G4cout << "Surface Upper   = " << a6 << G4endl;
#endif

  G4double chose = (a1 + a2 + a3 + a4 + a5 + a6) * G4QuickRand();

  if (chose < a1)
  {
    umin = fSide0->GetBoundaryMin(phi);
    umax = fSide0->GetBoundaryMax(phi);
    u = umin + (umax - umin) * G4QuickRand();
    return fSide0->SurfacePoint(phi, u, true);  // point on 0 deg surface
  }
  if ((chose >= a1) && (chose < a1 + a2))
  {
    umin = fSide90->GetBoundaryMin(phi);
    umax = fSide90->GetBoundaryMax(phi);
    u = umin + (umax - umin) * G4QuickRand();
    return fSide90->SurfacePoint(phi, u, true);  // point on 90 deg surface
  }
  if ((chose >= a1 + a2) && (chose < a1 + a2 + a3))
  {
    umin = fSide180->GetBoundaryMin(phi);
    umax = fSide180->GetBoundaryMax(phi);
    u = umin + (umax - umin) * G4QuickRand();
    return fSide180->SurfacePoint(phi, u, true);  // point on 180 deg surface
  }
  if ((chose >= a1 + a2 + a3) && (chose < a1 + a2 + a3 + a4))
  {
    umin = fSide270->GetBoundaryMin(phi);
    umax = fSide270->GetBoundaryMax(phi);
    u = umin + (umax - umin) * G4QuickRand();
    return fSide270->SurfacePoint(phi, u, true);  // point on 270 deg surface
  }
  if ((chose >= a1 + a2 + a3 + a4) && (chose < a1 + a2 + a3 + a4 + a5))
  {
    y = fDy1 * (2. * G4QuickRand() - 1.);
    umin = fLowerEndcap->GetBoundaryMin(y);
    umax = fLowerEndcap->GetBoundaryMax(y);
    u = umin + (umax - umin) * G4QuickRand();
    return fLowerEndcap->SurfacePoint(u, y, true);  // point on lower endcap
  }
  y = fDy2 * (2. * G4QuickRand() - 1.);
  umin = fUpperEndcap->GetBoundaryMin(y);
  umax = fUpperEndcap->GetBoundaryMax(y);
  u = umin + (umax - umin) * G4QuickRand();
  return fUpperEndcap->SurfacePoint(u, y, true);  // point on upper endcap
}

//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4VTwistedFaceted::CreatePolyhedron() const
{
  // number of meshes
  const G4int k = G4int(G4Polyhedron::GetNumberOfRotationSteps() * std::abs(fPhiTwist) / twopi) + 2;
  const G4int n = k;

  const G4int nnodes = 4 * (k - 1) * (n - 2) + 2 * k * k;
  const G4int nfaces = 4 * (k - 1) * (n - 1) + 2 * (k - 1) * (k - 1);

  auto ph = new G4Polyhedron;
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];
  auto xyz = new G4double3[nnodes];  // number of nodes
  auto faces = new G4int4[nfaces];  // number of faces

  fLowerEndcap->GetFacets(k, k, xyz, faces, 0);
  fUpperEndcap->GetFacets(k, k, xyz, faces, 1);
  fSide270->GetFacets(k, n, xyz, faces, 2);
  fSide0->GetFacets(k, n, xyz, faces, 3);
  fSide90->GetFacets(k, n, xyz, faces, 4);
  fSide180->GetFacets(k, n, xyz, faces, 5);

  ph->createPolyhedron(nnodes, nfaces, xyz, faces);

  return ph;
}

//=====================================================================
//* RayIntersectsBox --------------------------------------------------

G4bool G4VTwistedFaceted::RayIntersectsBox(const G4ThreeVector& p, const G4ThreeVector& v,
                                           const G4ThreeVector& minimum,
                                           const G4ThreeVector& maximum) const
{
  // Preserve the exact surface implementation's handling of non-numeric
  // inputs instead of making a rejection decision from unordered values.
  if (std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z()) || std::isnan(v.x())
      || std::isnan(v.y()) || std::isnan(v.z()))
  {
    return true;
  }

  G4double rayMinimum = 0.;
  G4double rayMaximum = kInfinity;

  const auto intersectsSlab = [&rayMinimum, &rayMaximum](G4double position, G4double direction,
                                                         G4double lower, G4double upper)
  {
    if (direction == 0.)
    {
      return position >= lower && position <= upper;
    }

    G4double slabMinimum;
    G4double slabMaximum;
    if (direction > 0.)
    {
      slabMinimum = (lower - position) / direction;
      slabMaximum = (upper - position) / direction;
    }
    else
    {
      slabMinimum = (upper - position) / direction;
      slabMaximum = (lower - position) / direction;
    }
    rayMinimum = std::max(rayMinimum, slabMinimum);
    rayMaximum = std::min(rayMaximum, slabMaximum);
    return rayMinimum <= rayMaximum;
  };

  return intersectsSlab(p.x(), v.x(), minimum.x(), maximum.x())
         && intersectsSlab(p.y(), v.y(), minimum.y(), maximum.y())
         && intersectsSlab(p.z(), v.z(), minimum.z(), maximum.z());
}

//=====================================================================
//* DistanceToBox -----------------------------------------------------

G4double G4VTwistedFaceted::DistanceToBox(const G4ThreeVector& p, const G4ThreeVector& minimum,
                                          const G4ThreeVector& maximum) const
{
  // Leave unordered inputs to the established constituent-surface queries.
  if (std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z()))
  {
    return 0.;
  }

  const auto separation = [](G4double position, G4double lower, G4double upper)
  {
    if (position < lower)
    {
      return lower - position;
    }
    if (position > upper)
    {
      return position - upper;
    }
    return 0.;
  };

  const G4double dx = separation(p.x(), minimum.x(), maximum.x());
  const G4double dy = separation(p.y(), minimum.y(), maximum.y());
  const G4double dz = separation(p.z(), minimum.z(), maximum.z());

  return std::min(std::hypot(dx, dy, dz), kInfinity);
}

//=====================================================================
//* GetSurface --------------------------------------------------------

G4VTwistSurface* G4VTwistedFaceted::GetSurface(G4int surfaceIndex) const
{
  switch (surfaceIndex)
  {
    case 0:
      return fSide0;
    case 1:
      return fSide90;
    case 2:
      return fSide180;
    case 3:
      return fSide270;
    case 4:
      return fLowerEndcap;
    case 5:
      return fUpperEndcap;
    default:
      return nullptr;
  }
}

//=====================================================================
//* EvaluatePoint -----------------------------------------------------

G4VTwistedFaceted::PointRange G4VTwistedFaceted::
EvaluatePoint(const G4VTwistedFaceted::RayFaceFunction& face,
              G4double t, G4double sine, G4double cosine,
              G4double evaluationGuard) const
{
  const G4double a = std::fma(face.a1, t, face.a0);
  const G4double b = std::fma(face.b1, t, face.b0);
  const G4double c = std::fma(face.c1, t, face.c0);
  const G4double value = std::fma(a, cosine, std::fma(b, sine, c));
  const G4double magnitude = std::fabs(a) + std::fabs(b) + std::fabs(c) + std::fabs(value) + 1.;
  const G4double roundoff =
    evaluationGuard + 512. * std::numeric_limits<G4double>::epsilon() * magnitude;
  return { std::nextafter(value - roundoff, -kInfinity),
           std::nextafter(value + roundoff, kInfinity) };
}

//=====================================================================
//* EvaluatePointWithCentral ------------------------------------------

G4bool G4VTwistedFaceted::
EvaluatePointWithCentral(const G4VTwistedFaceted::RayFaceFunction& face,
                         G4double t, G4double evaluationGuard,
                         G4VTwistedFaceted::PointRange& range, G4double& value,
                         G4double& derivative) const
{
  const G4double phase = std::fma(face.omega, t, face.phi0);
  const G4double sine = std::sin(phase);
  const G4double cosine = std::cos(phase);
  const G4double a = std::fma(face.a1, t, face.a0);
  const G4double b = std::fma(face.b1, t, face.b0);
  const G4double c = std::fma(face.c1, t, face.c0);
  value = std::fma(a, cosine, std::fma(b, sine, c));
  derivative =
    std::fma(face.a1 + face.omega * b, cosine, std::fma(face.b1 - face.omega * a, sine, face.c1));
  const G4double magnitude = std::fabs(a) + std::fabs(b) + std::fabs(c) + std::fabs(value) + 1.;
  const G4double roundoff =
    evaluationGuard + 512. * std::numeric_limits<G4double>::epsilon() * magnitude;
  range = { std::nextafter(value - roundoff, -kInfinity),
            std::nextafter(value + roundoff, kInfinity) };
  return std::isfinite(value) && std::isfinite(derivative) && std::isfinite(range.lower)
         && std::isfinite(range.upper);
}

//=====================================================================
//* EvaluateInterval --------------------------------------------------

G4bool G4VTwistedFaceted::
EvaluateInterval(const G4VTwistedFaceted::RayFaceFunction& face, G4double midpoint,
                 G4double radius, G4double sine, G4double cosine,
                 G4double evaluationGuard, G4VTwistedFaceted::FaceRange& result) const
{
  const G4double epsilon = std::numeric_limits<G4double>::epsilon();
  const G4double a = std::fma(face.a1, midpoint, face.a0);
  const G4double b = std::fma(face.b1, midpoint, face.b0);
  const G4double c = std::fma(face.c1, midpoint, face.c0);
  const G4double value = std::fma(a, cosine, std::fma(b, sine, c));
  const G4double derivative =
    std::fma(face.a1 + face.omega * b, cosine, std::fma(face.b1 - face.omega * a, sine, face.c1));

  // For F=(a0+a1*t)cos(phi)+(b0+b1*t)sin(phi)+(c0+c1*t),
  // F''=omega*((2*b1-omega*A)cos(phi)
  //            -(2*a1+omega*B)sin(phi)).  Bound it independently of
  // trigonometric phase over this interval, expanding every bound outwards.
  G4double maxA = std::fabs(a) + std::fabs(face.a1) * radius;
  G4double maxB = std::fabs(b) + std::fabs(face.b1) * radius;
  const G4double abMagnitude = maxA + maxB + std::fabs(a) + std::fabs(b) + 1.;
  const G4double abRoundoff = 512. * epsilon * abMagnitude;
  maxA = std::nextafter(maxA + abRoundoff, kInfinity);
  maxB = std::nextafter(maxB + abRoundoff, kInfinity);
  const G4double curvatureTerms = 2. * std::fabs(face.b1) + std::fabs(face.omega) * maxA
                                  + 2. * std::fabs(face.a1) + std::fabs(face.omega) * maxB;
  G4double curvature = std::fabs(face.omega) * curvatureTerms;
  curvature = std::nextafter(
    curvature + 512. * epsilon * (std::fabs(curvature) + curvatureTerms + 1.), kInfinity);

  const G4double derivativeRadius = curvature * radius;
  const G4double derivativeMagnitude = std::fabs(derivative) + derivativeRadius + std::fabs(face.a1)
                                       + std::fabs(face.b1) + std::fabs(face.c1)
                                       + std::fabs(face.omega) * (std::fabs(a) + std::fabs(b)) + 1.;
  const G4double derivativeRoundoff = 1024. * epsilon * derivativeMagnitude;
  const G4double linearRadius =
    std::nextafter((std::fabs(derivative) + derivativeRoundoff) * radius, kInfinity);
  const G4double remainder = std::nextafter(0.5 * curvature * radius * radius, kInfinity);
  const G4double valueRadius = std::nextafter(linearRadius + remainder, kInfinity);
  const G4double valueMagnitude =
    std::fabs(value) + valueRadius + std::fabs(a) + std::fabs(b) + std::fabs(c) + 1.;
  const G4double valueRoundoff = evaluationGuard + 1024. * epsilon * valueMagnitude;
  result.lower = std::nextafter(value - valueRadius - valueRoundoff, -kInfinity);
  result.upper = std::nextafter(value + valueRadius + valueRoundoff, kInfinity);

  result.derivativeLower =
    std::nextafter(derivative - derivativeRadius - derivativeRoundoff, -kInfinity);
  result.derivativeUpper =
    std::nextafter(derivative + derivativeRadius + derivativeRoundoff, kInfinity);
  return std::isfinite(result.lower) && std::isfinite(result.upper)
         && std::isfinite(result.derivativeLower) && std::isfinite(result.derivativeUpper);
}

//=====================================================================
//* SearchFirstTransition ---------------------------------------------

G4VTwistedFaceted::TransitionSearchStatus G4VTwistedFaceted::
SearchFirstTransition(const std::array<G4VTwistedFaceted::RayFaceFunction, 6>& faces,
                      G4int faceCount, G4bool entering, G4double lower,
                      G4double upper, G4double faceGuard, G4double lateralParallelGuard,
                      G4double evaluationGuard, G4double minimumWidth,
                      G4double lateralUpperLimit, G4int depth, G4int& nodeCount,
                      std::array<G4VTwistedFaceted::FaceRange, 6>& ranges,
                      G4VTwistedFaceted::SelectedTransition& selected) const
{
  constexpr G4int maxDepth = 32;
  constexpr G4int maxNodes = 96;
  if (++nodeCount > maxNodes || !(upper > lower) || !std::isfinite(lower) || !std::isfinite(upper))
  {
    return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
  }

  const G4double midpoint = lower + 0.5 * (upper - lower);
  const G4double radius = 0.5 * (upper - lower);
  if (!(midpoint > lower && midpoint < upper))
  {
    return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
  }
  const G4double phase = std::fma(faces[0].omega, midpoint, faces[0].phi0);
  const G4double sine = std::sin(phase);
  const G4double cosine = std::cos(phase);
  if (!std::isfinite(sine) || !std::isfinite(cosine))
  {
    return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
  }

  for (G4int i = 0; i < faceCount; ++i)
  {
    if (!EvaluateInterval(faces[i], midpoint, radius, sine, cosine, evaluationGuard, ranges[i]))
    {
      return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
    }
  }

  if (entering)
  {
    // If one constraint is strictly positive throughout, the ray cannot be
    // inside the solid anywhere in this interval and hence cannot enter it.
    for (G4int i = 0; i < faceCount; ++i)
    {
      if (ranges[i].lower > faceGuard)
      {
        return G4VTwistedFaceted::TransitionSearchStatus::kNoTransition;
      }
    }
  }
  else
  {
    // From an interior starting point an exit is impossible while every
    // lateral constraint is strictly negative over the complete interval.
    G4bool allInterior = true;
    for (G4int i = 0; i < faceCount; ++i)
    {
      allInterior = allInterior && ranges[i].upper < -faceGuard;
    }
    if (allInterior)
    {
      return G4VTwistedFaceted::TransitionSearchStatus::kNoTransition;
    }
  }

  const G4double lowerPhase = std::fma(faces[0].omega, lower, faces[0].phi0);
  const G4double upperPhase = std::fma(faces[0].omega, upper, faces[0].phi0);
  const G4double lowerSine = std::sin(lowerPhase);
  const G4double lowerCosine = std::cos(lowerPhase);
  const G4double upperSine = std::sin(upperPhase);
  const G4double upperCosine = std::cos(upperPhase);
  if (!std::isfinite(lowerSine) || !std::isfinite(lowerCosine) || !std::isfinite(upperSine)
      || !std::isfinite(upperCosine))
  {
    return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
  }

  for (G4int candidate = 0; candidate < faceCount; ++candidate)
  {
    const PointRange atLower =
      EvaluatePoint(faces[candidate], lower, lowerSine, lowerCosine, evaluationGuard);
    const PointRange atUpper =
      EvaluatePoint(faces[candidate], upper, upperSine, upperCosine, evaluationGuard);
    const G4double parallelGuard = candidate < 4 ? lateralParallelGuard : 0.;
    const G4bool orientation = entering ? (atLower.lower >= 0. && atUpper.upper <= 0.
                                           && ranges[candidate].derivativeUpper < -parallelGuard)
                                        : (atLower.upper <= 0. && atUpper.lower >= 0.
                                           && ranges[candidate].derivativeLower > parallelGuard);
    if (!orientation)
    {
      continue;
    }

    G4bool separatedFromEdges = true;
    for (G4int other = 0; other < faceCount; ++other)
    {
      if (other != candidate && ranges[other].upper >= -faceGuard)
      {
        separatedFromEdges = false;
        break;
      }
    }
    const G4bool separatedFromEndcap = entering || upper <= lateralUpperLimit;
    if (separatedFromEdges && separatedFromEndcap)
    {
      const G4double slopeLower =
        entering ? -ranges[candidate].derivativeUpper : ranges[candidate].derivativeLower;
      selected = {candidate, lower, upper, slopeLower};
      return G4VTwistedFaceted::TransitionSearchStatus::kSelected;
    }
  }

  if (depth >= maxDepth || upper - lower <= minimumWidth)
  {
    return G4VTwistedFaceted::TransitionSearchStatus::kFallback;
  }

  const TransitionSearchStatus first = SearchFirstTransition(
    faces, faceCount, entering, lower, midpoint, faceGuard, lateralParallelGuard, evaluationGuard,
    minimumWidth, lateralUpperLimit, depth + 1, nodeCount, ranges, selected);
  if (first != G4VTwistedFaceted::TransitionSearchStatus::kNoTransition)
  {
    return first;
  }
  return SearchFirstTransition(faces, faceCount, entering, midpoint, upper, faceGuard,
                               lateralParallelGuard, evaluationGuard, minimumWidth,
                               lateralUpperLimit, depth + 1, nodeCount, ranges, selected);
}
