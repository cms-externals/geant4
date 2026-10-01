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
// Implementation for G4Para class
//
// 21.03.95 P.Kent: Modified for `tolerant' geom
// 31.10.96 V.Grichine: Modifications according G4Box/Tubs before to commit
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal
// 29.05.17 E.Tcherniaev: complete revision, speed-up
////////////////////////////////////////////////////////////////////////////

#include "G4Para.hh"

#if !defined(G4GEOM_USE_UPARA)

#  include "G4AffineTransform.hh"
#  include "G4AutoLock.hh"
#  include "G4BoundingEnvelope.hh"
#  include "G4QuickRand.hh"
#  include "G4VGraphicsScene.hh"
#  include "G4VPVParameterisation.hh"
#  include "G4VoxelLimits.hh"

namespace
{
G4Mutex paraMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
//  Constructor - set & check half widths

G4Para::G4Para(const G4String& pName, G4double pDx, G4double pDy, G4double pDz, G4double pAlpha,
               G4double pTheta, G4double pPhi)
  : G4CSGSolid(pName), halfCarTolerance(0.5 * kCarTolerance)
{
  SetAllParameters(pDx, pDy, pDz, pAlpha, pTheta, pPhi);
  fRebuildPolyhedron = false;  // default value for G4CSGSolid
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor - design of trapezoid based on 8 vertices

G4Para::G4Para(const G4String& pName, const G4ThreeVector pt[8])
  : G4CSGSolid(pName), halfCarTolerance(0.5 * kCarTolerance)
{
  // Find dimensions and trigonometric values
  //
  fDx = (pt[3].x() - pt[2].x()) * 0.5;
  fDy = (pt[2].y() - pt[1].y()) * 0.5;
  fDz = pt[7].z();
  CheckParameters();  // check dimensions

  fTalpha = (pt[2].x() + pt[3].x() - pt[1].x() - pt[0].x()) * 0.25 / fDy;
  fTthetaCphi = (pt[4].x() + fDy * fTalpha + fDx) / fDz;
  fTthetaSphi = (pt[4].y() + fDy) / fDz;
  MakePlanes();

  // Recompute vertices
  //
  G4ThreeVector v[8];
  G4double DyTalpha = fDy * fTalpha;
  G4double DzTthetaSphi = fDz * fTthetaSphi;
  G4double DzTthetaCphi = fDz * fTthetaCphi;
  v[0].set(-DzTthetaCphi - DyTalpha - fDx, -DzTthetaSphi - fDy, -fDz);
  v[1].set(-DzTthetaCphi - DyTalpha + fDx, -DzTthetaSphi - fDy, -fDz);
  v[2].set(-DzTthetaCphi + DyTalpha - fDx, -DzTthetaSphi + fDy, -fDz);
  v[3].set(-DzTthetaCphi + DyTalpha + fDx, -DzTthetaSphi + fDy, -fDz);
  v[4].set(DzTthetaCphi - DyTalpha - fDx, DzTthetaSphi - fDy, fDz);
  v[5].set(DzTthetaCphi - DyTalpha + fDx, DzTthetaSphi - fDy, fDz);
  v[6].set(DzTthetaCphi + DyTalpha - fDx, DzTthetaSphi + fDy, fDz);
  v[7].set(DzTthetaCphi + DyTalpha + fDx, DzTthetaSphi + fDy, fDz);

  // Compare with original vertices
  //
  for (G4int i = 0; i < 8; ++i)
  {
    G4double delx = std::abs(pt[i].x() - v[i].x());
    G4double dely = std::abs(pt[i].y() - v[i].y());
    G4double delz = std::abs(pt[i].z() - v[i].z());
    G4double discrepancy = std::max(std::max(delx, dely), delz);
    if (discrepancy > 0.1 * kCarTolerance)
    {
      std::ostringstream message;
      G4long oldprc = message.precision(16);
      message << "Invalid vertice coordinates for Solid: " << GetName() << "\nVertix #" << i
              << ", discrepancy = " << discrepancy << "\n  original   : " << pt[i]
              << "\n  recomputed : " << v[i];
      G4cout.precision(oldprc);
      G4Exception("G4Para::G4Para()", "GeomSolids0002", FatalException, message);
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Para::G4Para(__void__& a) : G4CSGSolid(a), halfCarTolerance(0.5 * kCarTolerance)
{
  SetAllParameters(1., 1., 1., 0., 0., 0.);
  fRebuildPolyhedron = false;  // default value for G4CSGSolid
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Para::G4Para(const G4Para& rhs)
  : G4CSGSolid(rhs),
    halfCarTolerance(rhs.halfCarTolerance),
    fDx(rhs.fDx),
    fDy(rhs.fDy),
    fDz(rhs.fDz),
    fTalpha(rhs.fTalpha),
    fTthetaCphi(rhs.fTthetaCphi),
    fTthetaSphi(rhs.fTthetaSphi)
{
  for (G4int i = 0; i < 4; ++i)
  {
    fPlanes[i] = rhs.fPlanes[i];
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Para& G4Para::operator=(const G4Para& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4CSGSolid::operator=(rhs);

  // Copy data
  //
  halfCarTolerance = rhs.halfCarTolerance;
  fDx = rhs.fDx;
  fDy = rhs.fDy;
  fDz = rhs.fDz;
  fTalpha = rhs.fTalpha;
  fTthetaCphi = rhs.fTthetaCphi;
  fTthetaSphi = rhs.fTthetaSphi;
  for (G4int i = 0; i < 4; ++i)
  {
    fPlanes[i] = rhs.fPlanes[i];
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - set and check half-widths

void G4Para::SetAllParameters(G4double pDx, G4double pDy, G4double pDz, G4double pAlpha,
                              G4double pTheta, G4double pPhi)
{
  // Reset data of the base class
  fCubicVolume = 0;
  fSurfaceArea = 0;
  fRebuildPolyhedron = true;

  // Set parameters
  fDx = pDx;
  fDy = pDy;
  fDz = pDz;
  fTalpha = std::tan(pAlpha);
  fTthetaCphi = std::tan(pTheta) * std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta) * std::sin(pPhi);

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Check dimensions

void G4Para::CheckParameters()
{
  if (fDx < 2 * kCarTolerance || fDy < 2 * kCarTolerance || fDz < 2 * kCarTolerance)
  {
    std::ostringstream message;
    message << "Invalid (too small or negative) dimensions for Solid: " << GetName() << "\n  X - "
            << fDx << "\n  Y - " << fDy << "\n  Z - " << fDz;
    G4Exception("G4Para::CheckParameters()", "GeomSolids0002", FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Set side planes

void G4Para::MakePlanes()
{
  G4ThreeVector vx(1, 0, 0);
  G4ThreeVector vy(fTalpha, 1, 0);
  G4ThreeVector vz(fTthetaCphi, fTthetaSphi, 1);

  // Set -Y & +Y planes
  //
  G4ThreeVector ynorm = (vx.cross(vz)).unit();

  fPlanes[0].a = 0.;
  fPlanes[0].b = ynorm.y();
  fPlanes[0].c = ynorm.z();
  fPlanes[0].d = fPlanes[0].b * fDy;  // point (0,fDy,0) is on plane

  fPlanes[1].a = 0.;
  fPlanes[1].b = -fPlanes[0].b;
  fPlanes[1].c = -fPlanes[0].c;
  fPlanes[1].d = fPlanes[0].d;

  // Set -X & +X planes
  //
  G4ThreeVector xnorm = (vz.cross(vy)).unit();

  fPlanes[2].a = xnorm.x();
  fPlanes[2].b = xnorm.y();
  fPlanes[2].c = xnorm.z();
  fPlanes[2].d = fPlanes[2].a * fDx;  // point (fDx,0,0) is on plane

  fPlanes[3].a = -fPlanes[2].a;
  fPlanes[3].b = -fPlanes[2].b;
  fPlanes[3].c = -fPlanes[2].c;
  fPlanes[3].d = fPlanes[2].d;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification

void G4Para::ComputeDimensions(G4VPVParameterisation* p, const G4int n,
                               const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this, n, pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Para::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dz = GetZHalfLength();
  G4double dx = GetXHalfLength();
  G4double dy = GetYHalfLength();

  G4double x0 = dz * fTthetaCphi;
  G4double x1 = dy * GetTanAlpha();
  G4double xmin =
    std::min(std::min(std::min(-x0 - x1 - dx, -x0 + x1 - dx), x0 - x1 - dx), x0 + x1 - dx);
  G4double xmax =
    std::max(std::max(std::max(-x0 - x1 + dx, -x0 + x1 + dx), x0 - x1 + dx), x0 + x1 + dx);

  G4double y0 = dz * fTthetaSphi;
  G4double ymin = std::min(-y0 - dy, y0 - dy);
  G4double ymax = std::max(-y0 + dy, y0 + dy);

  pMin.set(xmin, ymin, -dz);
  pMax.set(xmax, ymax, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " << GetName() << " !"
            << "\npMin = " << pMin << "\npMax = " << pMax;
    G4Exception("G4Para::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Para::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
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

  // Set bounding envelope (benv) and calculate extent
  //
  G4double dz = GetZHalfLength();
  G4double dx = GetXHalfLength();
  G4double dy = GetYHalfLength();

  G4double x0 = dz * fTthetaCphi;
  G4double x1 = dy * GetTanAlpha();
  G4double y0 = dz * fTthetaSphi;

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0].set(-x0 - x1 - dx, -y0 - dy, -dz);
  baseA[1].set(-x0 - x1 + dx, -y0 - dy, -dz);
  baseA[2].set(-x0 + x1 + dx, -y0 + dy, -dz);
  baseA[3].set(-x0 + x1 - dx, -y0 + dy, -dz);

  baseB[0].set(+x0 - x1 - dx, y0 - dy, dz);
  baseB[1].set(+x0 - x1 + dx, y0 - dy, dz);
  baseB[2].set(+x0 + x1 + dx, y0 + dy, dz);
  baseB[3].set(+x0 + x1 - dx, y0 + dy, dz);

  std::vector<const G4ThreeVectorList*> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin, bmax, polygons);
  exist = benv.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Determine where is point p, inside/on_surface/outside
//

EInside G4Para::Inside(const G4ThreeVector& p) const
{
  G4double xx = fPlanes[2].a * p.x() + fPlanes[2].b * p.y() + fPlanes[2].c * p.z();
  G4double dx = std::abs(xx) + fPlanes[2].d;

  G4double yy = fPlanes[0].b * p.y() + fPlanes[0].c * p.z();
  G4double dy = std::abs(yy) + fPlanes[0].d;
  G4double dxy = std::max(dx, dy);

  G4double dz = std::abs(p.z()) - fDz;
  G4double dist = std::max(dxy, dz);

  if (dist > halfCarTolerance)
  {
    return kOutside;
  }

  return (dist > -halfCarTolerance) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Determine side where point is, and return corresponding normal

G4ThreeVector G4Para::SurfaceNormal(const G4ThreeVector& p) const
{
  G4int nsurf = 0;  // number of surfaces where p is placed

  // Check Z faces
  //
  G4double nz = 0;
  G4double dz = std::abs(p.z()) - fDz;
  if (std::abs(dz) <= halfCarTolerance)
  {
    nz = (p.z() < 0) ? -1 : 1;
    ++nsurf;
  }

  // Check Y faces
  //
  G4double ny = 0;
  G4double yy = fPlanes[0].b * p.y() + fPlanes[0].c * p.z();
  if (std::abs(fPlanes[0].d + yy) <= halfCarTolerance)
  {
    ny = fPlanes[0].b;
    nz += fPlanes[0].c;
    ++nsurf;
  }
  else if (std::abs(fPlanes[1].d - yy) <= halfCarTolerance)
  {
    ny = fPlanes[1].b;
    nz += fPlanes[1].c;
    ++nsurf;
  }

  // Check X faces
  //
  G4double nx = 0;
  G4double xx = fPlanes[2].a * p.x() + fPlanes[2].b * p.y() + fPlanes[2].c * p.z();
  if (std::abs(fPlanes[2].d + xx) <= halfCarTolerance)
  {
    nx = fPlanes[2].a;
    ny += fPlanes[2].b;
    nz += fPlanes[2].c;
    ++nsurf;
  }
  else if (std::abs(fPlanes[3].d - xx) <= halfCarTolerance)
  {
    nx = fPlanes[3].a;
    ny += fPlanes[3].b;
    nz += fPlanes[3].c;
    ++nsurf;
  }

  // Return normal
  //
  if (nsurf == 1)
  {
    return {nx, ny, nz};
  }
  if (nsurf != 0)
  {
    return G4ThreeVector(nx, ny, nz).unit();  // edge or corner
  }

  // Point is not on the surface
  //
#  ifdef G4CSGDEBUG
  std::ostringstream message;
  G4int oldprc = message.precision(16);
  message << "Point p is not on surface (!?) of solid: " << GetName() << G4endl;
  message << "Position:\n";
  message << "   p.x() = " << p.x() / mm << " mm\n";
  message << "   p.y() = " << p.y() / mm << " mm\n";
  message << "   p.z() = " << p.z() / mm << " mm";
  G4cout.precision(oldprc);
  G4Exception("G4Para::SurfaceNormal(p)", "GeomSolids1002", JustWarning, message);
  DumpInfo();
#  endif

  return ApproxSurfaceNormal(p);
}

//////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Para::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  G4double dist = -DBL_MAX;
  G4int iside = 0;
  for (G4int i = 0; i < 4; ++i)
  {
    G4double d = fPlanes[i].a * p.x() + fPlanes[i].b * p.y() + fPlanes[i].c * p.z() + fPlanes[i].d;
    if (d > dist)
    {
      dist = d;
      iside = i;
    }
  }

  G4double distz = std::abs(p.z()) - fDz;
  if (dist > distz)
  {
    return {fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c};
  }

  return {0, 0, (G4double)((p.z() < 0) ? -1 : 1)};
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside
//  - return kInfinity if no intersection

G4double G4Para::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  // Transform to the oblique box frame, where the solid is bounded by
  // |x| <= fDx, |y| <= fDy and |z| <= fDz.
  //
  G4double pz = p.z();
  G4double py = p.y() - fTthetaSphi * pz;
  G4double px = p.x() - fTthetaCphi * pz - fTalpha * py;
  G4double vz = v.z();
  G4double vy = v.y() - fTthetaSphi * vz;
  G4double vx = v.x() - fTthetaCphi * vz - fTalpha * vy;

  // Check if the point is on a surface and travelling away. The X and Y
  // distances are scaled to the corresponding physical plane normals.
  //
  if ((std::abs(pz) - fDz) >= -halfCarTolerance && pz * vz >= 0)
  {
    return kInfinity;
  }
  if ((std::abs(py) - fDy) * -fPlanes[0].b >= -halfCarTolerance && py * vy >= 0)
  {
    return kInfinity;
  }
  if ((std::abs(px) - fDx) * -fPlanes[2].a >= -halfCarTolerance && px * vx >= 0)
  {
    return kInfinity;
  }

  // Intersect the ray with the three pairs of box-frame slabs.
  //
  G4double invz = (vz == 0) ? DBL_MAX : -1. / vz;
  G4double invy = (vy == 0) ? DBL_MAX : -1. / vy;
  G4double invx = (vx == 0) ? DBL_MAX : -1. / vx;

  G4double dz = std::copysign(fDz, invz);
  G4double dy = std::copysign(fDy, invy);
  G4double dx = std::copysign(fDx, invx);

  G4double tzmin = (pz - dz) * invz;
  G4double tzmax = (pz + dz) * invz;
  G4double tymin = std::max(tzmin, (py - dy) * invy);
  G4double tymax = std::min(tzmax, (py + dy) * invy);
  G4double tmin = std::max(tymin, (px - dx) * invx);
  G4double tmax = std::min(tymax, (px + dx) * invx);

  // Find distance
  //
  if (tmax <= tmin + halfCarTolerance)
  {
    return kInfinity;
  }  // touch or no hit

  return (tmin < halfCarTolerance) ? 0. : tmin;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// - returns 0 is point inside

G4double G4Para::DistanceToIn(const G4ThreeVector& p) const
{
  G4double xx = fPlanes[2].a * p.x() + fPlanes[2].b * p.y() + fPlanes[2].c * p.z();
  G4double dx = std::abs(xx) + fPlanes[2].d;

  G4double yy = fPlanes[0].b * p.y() + fPlanes[0].c * p.z();
  G4double dy = std::abs(yy) + fPlanes[0].d;
  G4double dxy = std::max(dx, dy);

  G4double dz = std::abs(p.z()) - fDz;
  G4double dist = std::max(dxy, dz);

  return (dist > 0) ? dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside and, if required,
// find normal at exit point
// - when leaving the surface, return 0

G4double G4Para::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                               const G4bool calcNorm, G4bool* validNorm, G4ThreeVector* n) const
{
  // Transform to the oblique box frame.
  //
  G4double pz = p.z();
  G4double py = p.y() - fTthetaSphi * pz;
  G4double px = p.x() - fTthetaCphi * pz - fTalpha * py;
  G4double vz = v.z();
  G4double vy = v.y() - fTthetaSphi * vz;
  G4double vx = v.x() - fTthetaCphi * vz - fTalpha * vy;

  // Check if the point is on a surface and travelling away.
  //
  if ((std::abs(pz) - fDz) >= -halfCarTolerance && pz * vz > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0, 0, (pz < 0) ? -1 : 1);
    }
    return 0.;
  }
  if ((std::abs(py) - fDy) * -fPlanes[0].b >= -halfCarTolerance && py * vy > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      const auto& plane = (py < 0) ? fPlanes[0] : fPlanes[1];
      n->set(plane.a, plane.b, plane.c);
    }
    return 0.;
  }
  if ((std::abs(px) - fDx) * -fPlanes[2].a >= -halfCarTolerance && px * vx > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      const auto& plane = (px < 0) ? fPlanes[2] : fPlanes[3];
      n->set(plane.a, plane.b, plane.c);
    }
    return 0.;
  }

  // Find the nearest exit slab.
  //
  G4double tz = (vz == 0) ? DBL_MAX : (std::copysign(fDz, vz) - pz) / vz;
  G4double ty = (vy == 0) ? DBL_MAX : (std::copysign(fDy, vy) - py) / vy;
  G4double tx = (vx == 0) ? DBL_MAX : (std::copysign(fDx, vx) - px) / vx;
  G4double tmax = std::min(std::min(tz, ty), tx);

  if (calcNorm)
  {
    *validNorm = true;
    // Preserve the historical Z, Y, X priority for coincident exits.
    //
    if (tmax == tz)
    {
      n->set(0, 0, (vz < 0) ? -1 : 1);
    }
    else
    {
      G4int iside = (tmax == ty) ? ((vy < 0) ? 0 : 1) : ((vx < 0) ? 2 : 3);
      n->set(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c);
    }
  }
  return tmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - returns 0 is point outside

G4double G4Para::DistanceToOut(const G4ThreeVector& p) const
{
#  ifdef G4CSGDEBUG
  if (Inside(p) == kOutside)
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: " << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x() / mm << " mm\n";
    message << "   p.y() = " << p.y() / mm << " mm\n";
    message << "   p.z() = " << p.z() / mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4Para::DistanceToOut(p)", "GeomSolids1002", JustWarning, message);
    DumpInfo();
  }
#  endif
  G4double xx = fPlanes[2].a * p.x() + fPlanes[2].b * p.y() + fPlanes[2].c * p.z();
  G4double dx = std::abs(xx) + fPlanes[2].d;

  G4double yy = fPlanes[0].b * p.y() + fPlanes[0].c * p.z();
  G4double dy = std::abs(yy) + fPlanes[0].d;
  G4double dxy = std::max(dx, dy);

  G4double dz = std::abs(p.z()) - fDz;
  G4double dist = std::max(dxy, dz);

  return (dist < 0) ? -dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Para::GetEntityType() const
{
  return {"G4Para"};
}

//////////////////////////////////////////////////////////////////////////
//
// IsFaceted

G4bool G4Para::IsFaceted() const
{
  return true;
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Para::Clone() const
{
  return new G4Para(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Para::StreamInfo(std::ostream& os) const
{
  G4double alpha = std::atan(fTalpha);
  G4double theta = std::atan(std::sqrt(fTthetaCphi * fTthetaCphi + fTthetaSphi * fTthetaSphi));
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);

  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Para\n"
     << " Parameters:\n"
     << "    half length X: " << fDx / mm << " mm\n"
     << "    half length Y: " << fDy / mm << " mm\n"
     << "    half length Z: " << fDz / mm << " mm\n"
     << "    alpha: " << alpha / degree << "degrees\n"
     << "    theta: " << theta / degree << "degrees\n"
     << "    phi: " << phi / degree << "degrees\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Return a point randomly and uniformly selected on the solid surface

G4ThreeVector G4Para::GetPointOnSurface() const
{
  G4double sxy = fDx * fDy;
  G4double sxz = fDx * fDz * std::sqrt(1. + sqr(fTthetaSphi));
  G4double syz =
    fDy * fDz * std::sqrt(1. + sqr(fTalpha) + sqr(fTalpha * fTthetaSphi - fTthetaCphi));

  G4double select = (sxy + sxz + syz) * G4QuickRand();
  G4double u = 2. * G4QuickRand() - 1.;
  G4double v = 2. * G4QuickRand() - 1.;

  G4double x, y, z;
  if (select < sxy)
  {
    x = u * fDx;
    y = v * fDy;
    z = (select < 0.5 * sxy) ? -fDz : fDz;
  }
  else if (select < sxy + sxz)
  {
    x = u * fDx;
    y = (select < sxy + 0.5 * sxz) ? -fDy : fDy;
    z = v * fDz;
  }
  else
  {
    x = (select < sxy + sxz + 0.5 * syz) ? -fDx : fDx;
    y = u * fDy;
    z = v * fDz;
  }
  return {x + y * fTalpha + z * fTthetaCphi, y + z * fTthetaSphi, z};
}

//////////////////////////////////////////////////////////////////////////
//
// Get volume

G4double G4Para::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&paraMutex);
    fCubicVolume = 8 * fDx * fDy * fDz;
    l.unlock();
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Get surface area

G4double G4Para::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&paraMutex);
    G4double sxy = fDx * fDy;
    G4double sxz = fDx * fDz * std::sqrt(1. + sqr(fTthetaSphi));
    G4double syz =
      fDy * fDz * std::sqrt(1. + sqr(fTalpha) + sqr(fTalpha * fTthetaSphi - fTthetaCphi));
    fSurfaceArea = 8 * (sxy + sxz + syz);
    l.unlock();
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Para::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4Polyhedron* G4Para::CreatePolyhedron() const
{
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha = std::atan(fTalpha);
  G4double theta = std::atan(std::sqrt(fTthetaCphi * fTthetaCphi + fTthetaSphi * fTthetaSphi));

  return new G4PolyhedronPara(fDx, fDy, fDz, alpha, theta, phi);
}
#endif
