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
// G4TwistedTubs implementation
//
// Author: Kotoyo Hoshina (Chiba University), 01.08.2002 - created.
//         Oliver Link (CERN), 13.11.2003 - Integration in Geant4
//                from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include "G4TwistedTubs.hh"

#include "G4AffineTransform.hh"
#include "G4AutoLock.hh"
#include "G4BoundingEnvelope.hh"
#include "G4ClippablePolygon.hh"
#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polyhedron.hh"
#include "G4QuickRand.hh"
#include "G4SystemOfUnits.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisExtent.hh"
#include "G4VoxelLimits.hh"

namespace
{
G4Mutex twtubsMutex = G4MUTEX_INITIALIZER;
}  // namespace

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistedTubs::G4TwistedTubs(const G4String& pname, G4double twistedangle, G4double endinnerrad,
                             G4double endouterrad, G4double halfzlen, G4double dphi)
  : G4VSolid(pname), fDPhi(dphi)
{
  if (endinnerrad < DBL_MIN)
  {
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument,
                "Invalid end-inner-radius!");
  }

  G4double sinhalftwist = std::sin(0.5 * twistedangle);

  G4double endinnerradX = endinnerrad * sinhalftwist;
  G4double innerrad = std::sqrt(endinnerrad * endinnerrad - endinnerradX * endinnerradX);

  G4double endouterradX = endouterrad * sinhalftwist;
  G4double outerrad = std::sqrt(endouterrad * endouterrad - endouterradX * endouterradX);

  // temporary treatment!!
  SetFields(twistedangle, innerrad, outerrad, -halfzlen, halfzlen);
  CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname, G4double twistedangle, G4double endinnerrad,
                             G4double endouterrad, G4double halfzlen, G4int nseg, G4double totphi)
  : G4VSolid(pname)
{
  if (nseg == 0)
  {
    std::ostringstream message;
    message << "Invalid number of segments." << G4endl << "        nseg = " << nseg;
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument, message);
  }
  if (totphi == DBL_MIN || endinnerrad < DBL_MIN)
  {
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument,
                "Invalid total-phi or end-inner-radius!");
  }

  G4double sinhalftwist = std::sin(0.5 * twistedangle);

  G4double endinnerradX = endinnerrad * sinhalftwist;
  G4double innerrad = std::sqrt(endinnerrad * endinnerrad - endinnerradX * endinnerradX);

  G4double endouterradX = endouterrad * sinhalftwist;
  G4double outerrad = std::sqrt(endouterrad * endouterrad - endouterradX * endouterradX);

  // temporary treatment!!
  fDPhi = totphi / nseg;
  SetFields(twistedangle, innerrad, outerrad, -halfzlen, halfzlen);
  CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname, G4double twistedangle, G4double innerrad,
                             G4double outerrad, G4double negativeEndz, G4double positiveEndz,
                             G4double dphi)
  : G4VSolid(pname), fDPhi(dphi)
{
  if (innerrad < DBL_MIN)
  {
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument,
                "Invalid end-inner-radius!");
  }

  SetFields(twistedangle, innerrad, outerrad, negativeEndz, positiveEndz);
  CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname, G4double twistedangle, G4double innerrad,
                             G4double outerrad, G4double negativeEndz, G4double positiveEndz,
                             G4int nseg, G4double totphi)
  : G4VSolid(pname)
{
  if (nseg == 0)
  {
    std::ostringstream message;
    message << "Invalid number of segments." << G4endl << "        nseg = " << nseg;
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument, message);
  }
  if (totphi == DBL_MIN || innerrad < DBL_MIN)
  {
    G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002", FatalErrorInArgument,
                "Invalid total-phi or end-inner-radius!");
  }

  fDPhi = totphi / nseg;
  SetFields(twistedangle, innerrad, outerrad, negativeEndz, positiveEndz);
  CreateSurfaces();
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistedTubs::G4TwistedTubs(__void__& a)
  : G4VSolid(a)
{
}

//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedTubs::~G4TwistedTubs()
{
  DeleteSurfaces();
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

//=====================================================================
//* Copy constructor --------------------------------------------------

G4TwistedTubs::G4TwistedTubs(const G4TwistedTubs& rhs)
  : G4VSolid(rhs),
    fPhiTwist(rhs.fPhiTwist),
    fInnerRadius(rhs.fInnerRadius),
    fOuterRadius(rhs.fOuterRadius),
    fDPhi(rhs.fDPhi),
    fZHalfLength(rhs.fZHalfLength),
    fInnerStereo(rhs.fInnerStereo),
    fOuterStereo(rhs.fOuterStereo),
    fTanInnerStereo(rhs.fTanInnerStereo),
    fTanOuterStereo(rhs.fTanOuterStereo),
    fKappa(rhs.fKappa),
    fInnerRadius2(rhs.fInnerRadius2),
    fOuterRadius2(rhs.fOuterRadius2),
    fTanInnerStereo2(rhs.fTanInnerStereo2),
    fTanOuterStereo2(rhs.fTanOuterStereo2),
    fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea)
{
  for (auto i = 0; i < 2; ++i)
  {
    fEndZ[i] = rhs.fEndZ[i];
    fEndInnerRadius[i] = rhs.fEndInnerRadius[i];
    fEndOuterRadius[i] = rhs.fEndOuterRadius[i];
    fEndPhi[i] = rhs.fEndPhi[i];
    fEndZ2[i] = rhs.fEndZ2[i];
  }
  CreateSurfaces();
}

//=====================================================================
//* Assignment operator -----------------------------------------------

G4TwistedTubs& G4TwistedTubs::operator=(const G4TwistedTubs& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Release the old constituent surfaces before replacing the geometry.
  // CreateSurfaces() requires null pointers and creates a fresh, mutually
  // linked set of surfaces from the copied dimensions.
  DeleteSurfaces();

  // Copy data
  //
  fPhiTwist = rhs.fPhiTwist;
  fInnerRadius = rhs.fInnerRadius;
  fOuterRadius = rhs.fOuterRadius;
  fDPhi = rhs.fDPhi;
  fZHalfLength = rhs.fZHalfLength;
  fInnerStereo = rhs.fInnerStereo;
  fOuterStereo = rhs.fOuterStereo;
  fTanInnerStereo = rhs.fTanInnerStereo;
  fTanOuterStereo = rhs.fTanOuterStereo;
  fKappa = rhs.fKappa;
  fInnerRadius2 = rhs.fInnerRadius2;
  fOuterRadius2 = rhs.fOuterRadius2;
  fTanInnerStereo2 = rhs.fTanInnerStereo2;
  fTanOuterStereo2 = rhs.fTanOuterStereo2;
  fCubicVolume = rhs.fCubicVolume;
  fSurfaceArea = rhs.fSurfaceArea;

  for (auto i = 0; i < 2; ++i)
  {
    fEndZ[i] = rhs.fEndZ[i];
    fEndInnerRadius[i] = rhs.fEndInnerRadius[i];
    fEndOuterRadius[i] = rhs.fEndOuterRadius[i];
    fEndPhi[i] = rhs.fEndPhi[i];
    fEndZ2[i] = rhs.fEndZ2[i];
  }

  CreateSurfaces();
  fRebuildPolyhedron = false;
  delete fpPolyhedron;
  fpPolyhedron = nullptr;

  return *this;
}

//=====================================================================
//* DeleteSurfaces ----------------------------------------------------

void G4TwistedTubs::DeleteSurfaces()
{
  delete fLowerEndcap;
  delete fUpperEndcap;
  delete fLatterTwisted;
  delete fFormerTwisted;
  delete fInnerHype;
  delete fOuterHype;

  fLowerEndcap = nullptr;
  fUpperEndcap = nullptr;
  fLatterTwisted = nullptr;
  fFormerTwisted = nullptr;
  fInnerHype = nullptr;
  fOuterHype = nullptr;
}

//=====================================================================
//* BoundingLimits ----------------------------------------------------

void G4TwistedTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  // Find bounding tube
  G4double rmin = GetInnerRadius();
  G4double rmax = GetEndOuterRadius();

  G4double zmin = std::min(GetEndZ(0), GetEndZ(1));
  G4double zmax = std::max(GetEndZ(0), GetEndZ(1));

  G4double dphi = 0.5 * GetDPhi();
  G4double sphi = std::min(GetEndPhi(0), GetEndPhi(1)) - dphi;
  G4double ephi = std::max(GetEndPhi(0), GetEndPhi(1)) + dphi;
  G4double totalphi = ephi - sphi;

  // Find bounding box
  if (dphi <= 0 || totalphi >= CLHEP::twopi)
  {
    pMin.set(-rmax, -rmax, zmin);
    pMax.set(rmax, rmax, zmax);
  }
  else
  {
    G4TwoVector vmin, vmax;
    G4GeomTools::DiskExtent(rmin, rmax, sphi, totalphi, vmin, vmax);
    pMin.set(vmin.x(), vmin.y(), zmin);
    pMax.set(vmax.x(), vmax.y(), zmax);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " << GetName() << " !"
            << "\npMin = " << pMin << "\npMax = " << pMax;
    G4Exception("G4TwistedTubs::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool G4TwistedTubs::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
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

EInside G4TwistedTubs::Inside(const G4ThreeVector& p) const
{
  const G4double halftol = 0.5 * G4GeometryTolerance::GetInstance()->GetRadialTolerance();

  // Follow the same arithmetic and tolerance comparisons as the outer
  // hyperbolic surface, avoiding virtual dispatch and generic area-code
  // construction. The boundary lines and tolerance rotations are immutable
  // properties cached when the constituent surfaces are created.
  if (p.mag() < DBL_MIN)
  {
    return kOutside;
  }

  auto* outerHype = static_cast<G4TwistTubsHypeSide*>(fOuterHype);
  const G4double outerRho = outerHype->GetRhoAtPZ(p);
  const G4double pointRho = p.getRho();
  const G4double outerDistance = outerRho - pointRho;
  if (outerDistance < -halftol)
  {
    return kOutside;
  }

  EInside outerState = kInside;
  // A direct-object persistency reader can restore a legacy object without
  // calling CreateSurfaces(), leaving the newly introduced derived caches at
  // their in-class defaults. Use the established surface implementation in
  // that case. This fallback is deliberately read-only so shared queries stay
  // thread-safe.
  const G4bool cacheReady = fInsideHalfCarTolerance != 0.
                          && fPhiBoundaryDirection[0].z() != 0.
                          && fPhiBoundaryDirection[1].z() != 0.;
  if (!cacheReady)
  {
    outerState = outerHype->Inside(p);
  }
  else
  {
    const G4ThreeVector pointXY = G4ThreeVector(p.x(), p.y(), 0.).unit();
    const auto boundaryAtZ = [&p](const G4ThreeVector& direction,
                                  const G4ThreeVector& origin)
    {
      return ((p.z() - origin.z()) / direction.z()) * direction + origin;
    };
    const auto sideOfBoundary = [this, &pointXY](const G4ThreeVector& boundary)
    {
      const G4ThreeVector boundaryXY = G4ThreeVector(boundary.x(), boundary.y(), 0.).unit();
      const G4ThreeVector inverseBoundary = fPhiToleranceRotation[0] * boundaryXY;
      const G4ThreeVector rotatedBoundary = fPhiToleranceRotation[1] * boundaryXY;
      const G4double cross = pointXY.x() * boundaryXY.y() - pointXY.y() * boundaryXY.x();
      if (pointXY.x() * inverseBoundary.y() - pointXY.y() * inverseBoundary.x() > 0. && cross >= 0.)
      {
        return 1;
      }
      if (pointXY.x() * rotatedBoundary.y() - pointXY.y() * rotatedBoundary.x() < 0. && cross <= 0.)
      {
        return -1;
      }
      return 0;
    };

    G4bool areaBoundary = false;
    G4bool areaOutside = false;
    const G4ThreeVector lowerBoundary =
      boundaryAtZ(fPhiBoundaryDirection[0], fPhiBoundaryOrigin[0]);
    const G4int lowerSide = sideOfBoundary(lowerBoundary);
    if (lowerSide >= 0)
    {
      areaBoundary = true;
      if (lowerSide > 0)
      {
        areaOutside = true;
      }
    }
    else
    {
      const G4ThreeVector upperBoundary =
        boundaryAtZ(fPhiBoundaryDirection[1], fPhiBoundaryOrigin[1]);
      const G4int upperSide = sideOfBoundary(upperBoundary);
      if (upperSide <= 0)
      {
        areaBoundary = true;
        if (upperSide < 0)
        {
          areaOutside = true;
        }
      }
    }

    if (p.z() < fEndZ[0] + fInsideHalfCarTolerance)
    {
      areaBoundary = true;
      if (p.z() <= fEndZ[0] - fInsideHalfCarTolerance)
      {
        areaOutside = true;
      }
    }
    else if (p.z() > fEndZ[1] - fInsideHalfCarTolerance)
    {
      areaBoundary = true;
      if (p.z() >= fEndZ[1] + fInsideHalfCarTolerance)
      {
        areaOutside = true;
      }
    }

    if (areaOutside)
    {
      outerState = kOutside;
    }
    else if (areaBoundary || outerDistance <= halftol)
    {
      outerState = kSurface;
    }
  }

  if (outerState == kOutside)
  {
    return kOutside;
  }

  const auto* innerHype = static_cast<G4TwistTubsHypeSide*>(fInnerHype);
  const G4double innerRho = innerHype->GetRhoAtPZ(p);
  const G4double innerDistance = pointRho - innerRho;  // +ve: inside
  if (innerDistance < -halftol)
  {
    return kOutside;
  }
  if (outerState == kSurface)
  {
    return kSurface;
  }
  if (innerDistance <= halftol)
  {
    return kSurface;
  }

  return kInside;
}

//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4TwistedTubs::SurfaceNormal(const G4ThreeVector& p) const
{
  //
  // return the normal unit vector to the Hyperbolical Surface at a point
  // p on (or nearly on) the surface
  //
  // Which of the three or four surfaces are we closest to?
  //

  G4double distance = kInfinity;

  G4VTwistSurface* surfaces[] = { fLatterTwisted, fFormerTwisted, fInnerHype,
                                  fOuterHype,     fLowerEndcap,   fUpperEndcap };

  G4ThreeVector xx;
  G4ThreeVector bestxx;
  G4int besti = -1;
  for (auto i = 0; i < 6; ++i)
  {
    G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
    if (tmpdistance < distance)
    {
      if (tmpdistance == 0.)
      {
        return surfaces[i]->GetNormal(xx, true);
      }
      distance = tmpdistance;
      bestxx = xx;
      besti = i;
    }
  }

  return surfaces[besti]->GetNormal(bestxx, true);
}

//=====================================================================
//* DistanceToIn (p, v) -----------------------------------------------

G4double G4TwistedTubs::DistanceToIn(const G4ThreeVector& p,
                                     const G4ThreeVector& v) const
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
    // If the particle is entering to the volume, return 0.
    //
    G4ThreeVector normal = SurfaceNormal(p);
    if (normal * v < 0)
    {
      return 0;
    }
  }

  // now, we can take smallest positive distance.

  // Initialize
  //
  G4double distance = kInfinity;

  // find intersections and choose nearest one.
  //
  G4VTwistSurface* surfaces[] = { fLowerEndcap,   fUpperEndcap, fLatterTwisted,
                                  fFormerTwisted, fInnerHype,   fOuterHype };

  G4ThreeVector xx;
  for (const auto& surface : surfaces)
  {
    G4double tmpdistance = surface->DistanceToIn(p, v, xx);
    if (tmpdistance < distance)
    {
      distance = tmpdistance;
    }
  }
  return distance;
}

//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4TwistedTubs::DistanceToIn(const G4ThreeVector& p) const
{
  // DistanceToIn(p):
  // Calculate distance to surface of shape from `outside',
  // allowing for tolerance

  //
  // Calculate DistanceToIn(p)
  //

  EInside currentside = Inside(p);

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
      // Initialize
      G4double distance = kInfinity;

      // find intersections and choose nearest one.
      G4VTwistSurface* surfaces[] = { fLowerEndcap,   fUpperEndcap, fLatterTwisted,
                                      fFormerTwisted, fInnerHype,   fOuterHype };

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
      G4Exception("G4TwistedTubs::DistanceToIn(p)", "GeomSolids0003", FatalException,
                  "Unknown point location!");
    }
  }  // switch end

  return kInfinity;
}

//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double G4TwistedTubs::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
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
    // If the particle is exiting from the volume, return 0.
    //
    G4ThreeVector normal = SurfaceNormal(p);
    if (normal * v > 0)
    {
      if (calcNorm)
      {
        *norm = normal;
        *validNorm = true;
      }
      return 0;
    }
  }

  // now, we can take smallest positive distance.

  // Initialize
  //
  G4double distance = kInfinity;

  // find intersections and choose nearest one.
  //
  G4VTwistSurface* surfaces[] = { fLatterTwisted, fFormerTwisted, fInnerHype,
                                  fOuterHype,     fLowerEndcap,   fUpperEndcap };

  G4int besti = -1;
  G4ThreeVector xx;
  G4ThreeVector bestxx;
  for (auto i = 0; i < 6; ++i)
  {
    G4double tmpdistance = surfaces[i]->DistanceToOut(p, v, xx);
    if (tmpdistance < distance)
    {
      distance = tmpdistance;
      bestxx = xx;
      besti = i;
    }
  }

  if (calcNorm)
  {
    if (besti != -1)
    {
      *norm = surfaces[besti]->GetNormal(bestxx, true);
      *validNorm = surfaces[besti]->IsValidNorm();
    }
  }

  return distance;
}

//=====================================================================
//* DistanceToOut (p) ----------------------------------------------

G4double G4TwistedTubs::DistanceToOut(const G4ThreeVector& p) const
{
  // DistanceToOut(p):
  // Calculate distance to surface of shape from `inside',
  // allowing for tolerance

  //
  // Calculate DistanceToOut(p)
  //

  EInside currentside = Inside(p);

  switch (currentside)
  {
    case (kOutside):
    {
    }
    case (kSurface):
    {
      return 0;
    }
    case (kInside):
    {
      // Initialize
      G4double distance = kInfinity;

      // find intersections and choose nearest one.
      G4VTwistSurface* surfaces[] = { fLatterTwisted, fFormerTwisted, fInnerHype,
                                      fOuterHype,     fLowerEndcap,   fUpperEndcap };

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
      G4Exception("G4TwistedTubs::DistanceToOut(p)", "GeomSolids0003", FatalException,
                  "Unknown point location!");
    }
  }  // switch end

  return 0.;
}

//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4TwistedTubs::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedTubs\n"
     << " Parameters: \n"
     << "    -ve end Z              : " << fEndZ[0] / mm << " mm \n"
     << "    +ve end Z              : " << fEndZ[1] / mm << " mm \n"
     << "    inner end radius(-ve z): " << fEndInnerRadius[0] / mm << " mm \n"
     << "    inner end radius(+ve z): " << fEndInnerRadius[1] / mm << " mm \n"
     << "    outer end radius(-ve z): " << fEndOuterRadius[0] / mm << " mm \n"
     << "    outer end radius(+ve z): " << fEndOuterRadius[1] / mm << " mm \n"
     << "    inner radius (z=0)     : " << fInnerRadius / mm << " mm \n"
     << "    outer radius (z=0)     : " << fOuterRadius / mm << " mm \n"
     << "    twisted angle          : " << fPhiTwist / degree << " degrees \n"
     << "    inner stereo angle     : " << fInnerStereo / degree << " degrees \n"
     << "    outer stereo angle     : " << fOuterStereo / degree << " degrees \n"
     << "    phi-width of a piece   : " << fDPhi / degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4TwistedTubs::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4TwistedTubs::GetExtent() const
{
  // Define the sides of the box into which the G4Tubs instance would fit.
  //
  G4ThreeVector pmin, pmax;
  BoundingLimits(pmin, pmax);
  return {pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmin.z(), pmax.z()};
}

//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4TwistedTubs::CreatePolyhedron() const
{
  // number of meshes
  //
  G4double absPhiTwist = std::abs(fPhiTwist);
  G4double dA = std::max(fDPhi, absPhiTwist);
  const G4int k = G4int(G4Polyhedron::GetNumberOfRotationSteps() * dA / twopi) + 2;
  const G4int n = G4int(G4Polyhedron::GetNumberOfRotationSteps() * absPhiTwist / twopi) + 2;

  const G4int nnodes = 4 * (k - 1) * (n - 2) + 2 * k * k;
  const G4int nfaces = 4 * (k - 1) * (n - 1) + 2 * (k - 1) * (k - 1);

  auto ph = new G4Polyhedron;
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];
  auto xyz = new G4double3[nnodes];  // number of nodes
  auto faces = new G4int4[nfaces];  // number of faces
  fLowerEndcap->GetFacets(k, k, xyz, faces, 0);
  fUpperEndcap->GetFacets(k, k, xyz, faces, 1);
  fInnerHype->GetFacets(k, n, xyz, faces, 2);
  fFormerTwisted->GetFacets(k, n, xyz, faces, 3);
  fOuterHype->GetFacets(k, n, xyz, faces, 4);
  fLatterTwisted->GetFacets(k, n, xyz, faces, 5);

  ph->createPolyhedron(nnodes, nfaces, xyz, faces);

  delete[] xyz;
  delete[] faces;

  return ph;
}

//=====================================================================
//* GetPolyhedron -----------------------------------------------------

G4Polyhedron* G4TwistedTubs::GetPolyhedron() const
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
//* CreateSurfaces ----------------------------------------------------

void G4TwistedTubs::CreateSurfaces()
{
  // create 6 surfaces of TwistedTub

  fLowerEndcap = new G4TwistTubsFlatSide("LowerEndcap", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                         fEndPhi, fEndZ, -1);

  fUpperEndcap = new G4TwistTubsFlatSide("UpperEndcap", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                         fEndPhi, fEndZ, 1);

  fLatterTwisted = new G4TwistTubsSide("LatterTwisted", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                       fEndPhi, fEndZ, fInnerRadius, fOuterRadius, fKappa, 1);
  fFormerTwisted = new G4TwistTubsSide("FormerTwisted", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                       fEndPhi, fEndZ, fInnerRadius, fOuterRadius, fKappa, -1);

  fInnerHype = new G4TwistTubsHypeSide("InnerHype", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                       fEndPhi, fEndZ, fInnerRadius, fOuterRadius, fKappa,
                                       fTanInnerStereo, fTanOuterStereo, -1);
  fOuterHype = new G4TwistTubsHypeSide("OuterHype", fEndInnerRadius, fEndOuterRadius, fDPhi,
                                       fEndPhi, fEndZ, fInnerRadius, fOuterRadius, fKappa,
                                       fTanInnerStereo, fTanOuterStereo, 1);

  G4int boundaryType = 0;
  fOuterHype->GetBoundaryParameters(G4VTwistSurface::sAxis0 & G4VTwistSurface::sAxisMin,
                                    fPhiBoundaryDirection[0], fPhiBoundaryOrigin[0], boundaryType);
  fOuterHype->GetBoundaryParameters(G4VTwistSurface::sAxis0 & G4VTwistSurface::sAxisMax,
                                    fPhiBoundaryDirection[1], fPhiBoundaryOrigin[1], boundaryType);

  const G4double angularTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  fInsideHalfCarTolerance = 0.5 * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4RotationMatrix toleranceRotation;
  fPhiToleranceRotation[1] = toleranceRotation.rotateZ(0.5 * angularTolerance);
  fPhiToleranceRotation[0] = toleranceRotation.rotateZ(-angularTolerance);

  // set neighbour surfaces
  //
  fLowerEndcap->SetNeighbours(fInnerHype, fLatterTwisted, fOuterHype, fFormerTwisted);
  fUpperEndcap->SetNeighbours(fInnerHype, fLatterTwisted, fOuterHype, fFormerTwisted);
  fLatterTwisted->SetNeighbours(fInnerHype, fLowerEndcap, fOuterHype, fUpperEndcap);
  fFormerTwisted->SetNeighbours(fInnerHype, fLowerEndcap, fOuterHype, fUpperEndcap);
  fInnerHype->SetNeighbours(fLatterTwisted, fLowerEndcap, fFormerTwisted, fUpperEndcap);
  fOuterHype->SetNeighbours(fLatterTwisted, fLowerEndcap, fFormerTwisted, fUpperEndcap);
}

//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTubs::GetEntityType() const
{
  return {"G4TwistedTubs"};
}

//=====================================================================
//* Clone -------------------------------------------------------------

G4VSolid* G4TwistedTubs::Clone() const
{
  return new G4TwistedTubs(*this);
}

//=====================================================================
//* GetLateralArea ----------------------------------------------------

G4double G4TwistedTubs::GetLateralArea(G4double a, G4double r, G4double z) const
{
  if (z == 0)
  {
    return 0.;
  }
  G4double h = std::abs(z);
  G4double area = h * a;
  if (std::abs(a - r) > kCarTolerance)
  {
    G4double aa = a * a;
    G4double hh = h * h;
    G4double rr = r * r;
    G4double cc = aa * hh / (rr - aa);
    G4double k = std::sqrt(aa + cc) / cc;
    G4double kh = k * h;
    area = 0.5 * a * (h * std::sqrt(1. + kh * kh) + std::asinh(kh) / k);
  }
  return GetDPhi() * area;
}

//=====================================================================
//* GetPhiCutArea -----------------------------------------------------

G4double G4TwistedTubs::GetPhiCutArea(G4double a, G4double r, G4double z) const
{
  if (GetDPhi() >= CLHEP::twopi || r <= 0 || z == 0)
  {
    return 0.;
  }
  G4double h = std::abs(z);
  G4double area = h * a;
  if (GetPhiTwist() > kCarTolerance)
  {
    G4double sinw = std::sin(0.5 * GetPhiTwist()) * h / GetZHalfLength();
    G4double p = sinw * r / h;
    G4double q = sinw * r / a;
    G4double pp = p * p;
    G4double qq = q * q;
    G4double pq = p * q;
    G4double sqroot = std::sqrt(pp + qq + 1);
    area =
      (pq * sqroot + 0.5 * p * (pp + 3.) * std::atanh(q / sqroot)
       + 0.5 * q * (qq + 3.) * std::atanh(p / sqroot) + std::atan(sqroot / (pq)) - CLHEP::halfpi)
      * h * a / (3. * pq);
  }
  return area;
}

//=====================================================================
//* GetCubicVolume ----------------------------------------------------

G4double G4TwistedTubs::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&twtubsMutex);
    G4double DPhi = GetDPhi();
    G4double Z0 = GetEndZ(0);
    G4double Z1 = GetEndZ(1);
    G4double Ain = GetInnerRadius();
    G4double Aout = GetOuterRadius();
    G4double R0in = GetEndInnerRadius(0);
    G4double R1in = GetEndInnerRadius(1);
    G4double R0out = GetEndOuterRadius(0);
    G4double R1out = GetEndOuterRadius(1);

    // V_hyperboloid = pi*h*(2*a*a + R*R)/3
    fCubicVolume = (2. * (Z1 - Z0) * (Aout + Ain) * (Aout - Ain)
                    + Z1 * (R1out + R1in) * (R1out - R1in) - Z0 * (R0out + R0in) * (R0out - R0in))
                   * DPhi / 6.;
    l.unlock();
  }
  return fCubicVolume;
}

//=====================================================================
//* GetSurfaceArea ----------------------------------------------------

G4double G4TwistedTubs::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&twtubsMutex);
    G4double dphi = GetDPhi();
    G4double Ainn = GetInnerRadius();
    G4double Aout = GetOuterRadius();
    G4double Rinn0 = GetEndInnerRadius(0);
    G4double Rout0 = GetEndOuterRadius(0);
    G4double Rinn1 = GetEndInnerRadius(1);
    G4double Rout1 = GetEndOuterRadius(1);
    G4double z0 = GetEndZ(0);
    G4double z1 = GetEndZ(1);

    G4double base0 = 0.5 * dphi * (Rout0 * Rout0 - Rinn0 * Rinn0);  // lower base
    G4double inner0 = GetLateralArea(Ainn, Rinn0, z0);  // lower inner surface
    G4double outer0 = GetLateralArea(Aout, Rout0, z0);  // lower outer surface
    G4double cut0 =  // lower phi cut
      GetPhiCutArea(Aout, Rout0, z0) - GetPhiCutArea(Ainn, Rinn0, z0);

    G4double base1 = base0;
    G4double inner1 = inner0;
    G4double outer1 = outer0;
    G4double cut1 = cut0;
    if (std::abs(z0) != std::abs(z1))
    {
      base1 = 0.5 * dphi * (Rout1 * Rout1 - Rinn1 * Rinn1);  // upper base
      inner1 = GetLateralArea(Ainn, Rinn1, z1);  // upper inner surface
      outer1 = GetLateralArea(Aout, Rout1, z1);  // upper outer surface
      cut1 =  // upper phi cut
        GetPhiCutArea(Aout, Rout1, z1) - GetPhiCutArea(Ainn, Rinn1, z1);
    }
    fSurfaceArea =
      base0 + base1
      + ((z0 * z1 < 0) ? (inner0 + inner1 + outer0 + outer1 + 2. * (cut0 + cut1))
                       : std::abs(inner0 - inner1 + outer0 - outer1 + 2. * (cut0 - cut1)));
    l.unlock();
  }
  return fSurfaceArea;
}

//=====================================================================
//* GetPointOnSurface -------------------------------------------------

G4ThreeVector G4TwistedTubs::GetPointOnSurface() const
{
  G4double z = (fEndZ[1] - fEndZ[0]) * G4QuickRand() + fEndZ[0];
  G4double phi, phimin, phimax;
  G4double x, xmin, xmax;
  G4double r, rmin, rmax;

  G4double a1 = fOuterHype->GetSurfaceArea();
  G4double a2 = fInnerHype->GetSurfaceArea();
  G4double a3 = fLatterTwisted->GetSurfaceArea();
  G4double a4 = fFormerTwisted->GetSurfaceArea();
  G4double a5 = fLowerEndcap->GetSurfaceArea();
  G4double a6 = fUpperEndcap->GetSurfaceArea();

  G4double chose = (a1 + a2 + a3 + a4 + a5 + a6) * G4QuickRand();

  if (chose < a1)
  {
    phimin = fOuterHype->GetBoundaryMin(z);
    phimax = fOuterHype->GetBoundaryMax(z);
    phi = (phimax - phimin) * G4QuickRand() + phimin;

    return fOuterHype->SurfacePoint(phi, z, true);
  }
  if ((chose >= a1) && (chose < a1 + a2))
  {
    phimin = fInnerHype->GetBoundaryMin(z);
    phimax = fInnerHype->GetBoundaryMax(z);
    phi = (phimax - phimin) * G4QuickRand() + phimin;

    return fInnerHype->SurfacePoint(phi, z, true);
  }
  if ((chose >= a1 + a2) && (chose < a1 + a2 + a3))
  {
    xmin = fLatterTwisted->GetBoundaryMin(z);
    xmax = fLatterTwisted->GetBoundaryMax(z);
    x = (xmax - xmin) * G4QuickRand() + xmin;

    return fLatterTwisted->SurfacePoint(x, z, true);
  }
  if ((chose >= a1 + a2 + a3) && (chose < a1 + a2 + a3 + a4))
  {
    xmin = fFormerTwisted->GetBoundaryMin(z);
    xmax = fFormerTwisted->GetBoundaryMax(z);
    x = (xmax - xmin) * G4QuickRand() + xmin;

    return fFormerTwisted->SurfacePoint(x, z, true);
  }
  if ((chose >= a1 + a2 + a3 + a4) && (chose < a1 + a2 + a3 + a4 + a5))
  {
    rmin = GetEndInnerRadius(0);
    rmax = GetEndOuterRadius(0);
    r = std::sqrt((sqr(rmax) - sqr(rmin)) * G4QuickRand() + sqr(rmin));

    phimin = fLowerEndcap->GetBoundaryMin(r);
    phimax = fLowerEndcap->GetBoundaryMax(r);
    phi = (phimax - phimin) * G4QuickRand() + phimin;

    return fLowerEndcap->SurfacePoint(phi, r, true);
  }

  rmin = GetEndInnerRadius(1);
  rmax = GetEndOuterRadius(1);
  r = rmin + (rmax - rmin) * std::sqrt(G4QuickRand());

  phimin = fUpperEndcap->GetBoundaryMin(r);
  phimax = fUpperEndcap->GetBoundaryMax(r);
  phi = (phimax - phimin) * G4QuickRand() + phimin;

  return fUpperEndcap->SurfacePoint(phi, r, true);
}
