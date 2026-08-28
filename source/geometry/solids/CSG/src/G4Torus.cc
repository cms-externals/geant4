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
// G4Torus implementation
//
// 30.10.96 V.Grichine: first implementation with G4Tubs elements in Fs
// 26.05.00 V.Grichine: added new fuctions developed by O.Cremonesi
// 31.08.00 E.Medernach: numerical computation of roots with bounding volume
// 11.01.01 E.Medernach: Use G4PolynomialSolver to find roots
// 03.05.05 V.Grichine: SurfaceNormal(p) according to J. Apostolakis proposal
// 25.08.05 O.Link: new methods for DistanceToIn/Out using JTPolynomialSolver
// 28.10.16 E.Tcherniaev: new CalculateExtent(); removed CreateRotatedVertices()
// 16.12.16 H.Burkhardt: use radius differences and hypot to improve precision
// --------------------------------------------------------------------

#include "G4Torus.hh"

#if !defined(G4GEOM_USE_UTORUS)

#  include "G4AffineTransform.hh"
#  include "G4AutoLock.hh"
#  include "G4BoundingEnvelope.hh"
#  include "G4GeomTools.hh"
#  include "G4GeometryTolerance.hh"
#  include "G4JTPolynomialSolver.hh"
#  include "G4Polyhedron.hh"
#  include "G4QuickRand.hh"
#  include "G4VGraphicsScene.hh"
#  include "G4VPVParameterisation.hh"
#  include "G4VoxelLimits.hh"

namespace
{
G4Mutex torusMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

// Private enums: Not for external use
namespace
{
// Used by distanceToOut
enum ESide
{
  kNull,
  kRMin,
  kRMax,
  kSPhi,
  kEPhi
};

// used by normal
enum ENorm
{
  kNRMin,
  kNRMax,
  kNSPhi,
  kNEPhi
};

G4int SolveCubic(G4double a, G4double b, G4double c, G4double* roots)
{
  constexpr G4double oneThird = 1. / 3.;
  const G4double sqrtThree = std::sqrt(3.);
  const G4double invSixSqrtThree = 1. / (6. * sqrtThree);
  const G4double p = b - a * a * oneThird;
  const G4double q = c - a * b * oneThird + 2. * a * a * a / 27.;
  G4double discriminant = 4. * p * p * p + 27. * q * q;

  if (discriminant >= 0.)
  {
    discriminant = std::sqrt(discriminant) * invSixSqrtThree;
    const G4double t = -0.5 * q + discriminant;
    const G4double u = 0.5 * q + discriminant;
    roots[0] = std::cbrt(t) - std::cbrt(u) - a * oneThird;
  }
  else
  {
    discriminant = std::sqrt(-discriminant);
    const G4double t = -0.5 * q;
    const G4double u = discriminant * invSixSqrtThree;
    const G4double magnitude = std::hypot(t, u);
    const G4double angle = std::acos(std::clamp(t / magnitude, -1., 1.));
    roots[0] = 2. * std::cbrt(magnitude) * std::cos(oneThird * angle) - a * oneThird;
  }

  const G4double t = roots[0] * roots[0] + a * roots[0] + b;
  const G4double u = a + roots[0];
  discriminant = u * u - 4. * t;
  if (discriminant < 0.) return 1;
  discriminant = std::sqrt(discriminant);
  roots[1] = 0.5 * (-u - discriminant);
  roots[2] = 0.5 * (-u + discriminant);
  return 3;
}

G4int SolveBiquadratic(G4double a, G4double e, G4double g, G4double* roots)
{
  G4double discriminant = e * e - 4. * g;
  if (discriminant < 0.) return 0;
  discriminant = std::sqrt(discriminant);
  G4int nRoots = 0;
  G4double h = 0.5 * (-e - discriminant);
  if (h >= 0.)
  {
    h = std::sqrt(h);
    roots[nRoots++] = -h - 0.25 * a;
    roots[nRoots++] = h - 0.25 * a;
  }
  h = 0.5 * (-e + discriminant);
  if (h >= 0.)
  {
    h = std::sqrt(h);
    roots[nRoots++] = -h - 0.25 * a;
    roots[nRoots++] = h - 0.25 * a;
  }
  std::sort(roots, roots + nRoots);
  return nRoots;
}

G4int SolveQuartic(G4double a, G4double b, G4double c, G4double d, G4double* roots)
{
  const G4double a2 = a * a;
  const G4double e = b - 3. * a2 / 8.;
  const G4double f = c + a * a2 / 8. - 0.5 * a * b;
  const G4double g = d - 3. * a2 * a2 / 256. + a2 * b / 16. - a * c / 4.;

  if (std::fabs(f) < 1.e-12)
  {
    const G4int nRoots = SolveBiquadratic(a, e, g, roots);
    if (nRoots > 0) return nRoots;
  }

  G4double cubicRoots[3];
  const G4int nCubic = SolveCubic(2. * e, e * e - 4. * g, -f * f, cubicRoots);
  G4double h = cubicRoots[0];
  if (nCubic == 1)
  {
    const G4double coefficientScale = std::max({1., std::fabs(e), std::fabs(g)});
    if (std::fabs(f) <= 1.e-3 * coefficientScale && h <= 1.e-8 * coefficientScale)
    {
      const G4int nRoots = SolveBiquadratic(a, e, g, roots);
      if (nRoots > 0) return nRoots;
    }
    if (h <= 0.) return 0;
  }
  else
  {
    for (G4int i = 0; i < 3; ++i)
    {
      if (cubicRoots[i] >= 0.)
      {
        h = cubicRoots[i];
        break;
      }
    }
    if (h <= 0.) return 0;
  }

  h = std::sqrt(h);
  const G4double j = 0.5 * (e + h * h - f / h);
  G4int nRoots = 0;
  G4double discriminant = h * h - 4. * j;
  if (discriminant >= 0.)
  {
    discriminant = std::sqrt(discriminant);
    roots[nRoots++] = 0.5 * (-h - discriminant) - 0.25 * a;
    roots[nRoots++] = 0.5 * (-h + discriminant) - 0.25 * a;
  }
  discriminant = h * h - 4. * g / j;
  if (discriminant >= 0.)
  {
    discriminant = std::sqrt(discriminant);
    roots[nRoots++] = 0.5 * (h - discriminant) - 0.25 * a;
    roots[nRoots++] = 0.5 * (h + discriminant) - 0.25 * a;
  }
  std::sort(roots, roots + nRoots);
  return nRoots;
}
}  // namespace

///////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4Torus::G4Torus(const G4String& pName, G4double pRmin, G4double pRmax, G4double pRtor,
                 G4double pSPhi, G4double pDPhi)
  : G4CSGSolid(pName)
{
  SetAllParameters(pRmin, pRmax, pRtor, pSPhi, pDPhi);
}

////////////////////////////////////////////////////////////////////////////
//
//

void G4Torus::SetAllParameters(G4double pRmin, G4double pRmax, G4double pRtor, G4double pSPhi,
                               G4double pDPhi)
{
  const G4double fEpsilon = 4.e-11;  // relative tolerance of radii

  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fRebuildPolyhedron = true;

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance = 0.5 * kCarTolerance;
  halfAngTolerance = 0.5 * kAngTolerance;

  if (pRtor >= pRmax + 1.e3 * kCarTolerance)  // Check swept radius, as in G4Cons
  {
    fRtor = pRtor;
  }
  else
  {
    std::ostringstream message;
    message << "Invalid swept radius for Solid: " << GetName() << G4endl
            << "        pRtor = " << pRtor << ", pRmax = " << pRmax;
    G4Exception("G4Torus::SetAllParameters()", "GeomSolids0002", FatalException, message);
  }

  // Check radii, as in G4Cons
  //
  if (pRmin < pRmax - 1.e2 * kCarTolerance && pRmin >= 0)
  {
    if (pRmin >= 1.e2 * kCarTolerance)
    {
      fRmin = pRmin;
    }
    else
    {
      fRmin = 0.0;
    }
    fRmax = pRmax;
  }
  else
  {
    std::ostringstream message;
    message << "Invalid values of radii for Solid: " << GetName() << G4endl
            << "        pRmin = " << pRmin << ", pRmax = " << pRmax;
    G4Exception("G4Torus::SetAllParameters()", "GeomSolids0002", FatalException, message);
  }

  // Relative tolerances
  //
  fRminTolerance = (fRmin) != 0.0 ? 0.5 * std::max(kRadTolerance, fEpsilon * (fRtor - fRmin)) : 0;
  fRmaxTolerance = 0.5 * std::max(kRadTolerance, fEpsilon * (fRtor + fRmax));

  // Check angles
  //
  if (pDPhi >= twopi)
  {
    fDPhi = twopi;
  }
  else
  {
    if (pDPhi > 0)
    {
      fDPhi = pDPhi;
    }
    else
    {
      std::ostringstream message;
      message << "Invalid Z delta-Phi for Solid: " << GetName() << G4endl
              << "        pDPhi = " << pDPhi;
      G4Exception("G4Torus::SetAllParameters()", "GeomSolids0002", FatalException, message);
    }
  }

  // Ensure psphi in 0-2PI or -2PI-0 range if shape crosses 0
  //
  fSPhi = pSPhi;

  if (fSPhi < 0)
  {
    fSPhi = twopi - std::fmod(std::fabs(fSPhi), twopi);
  }
  else
  {
    fSPhi = std::fmod(fSPhi, twopi);
  }

  if (fSPhi + fDPhi > twopi)
  {
    fSPhi -= twopi;
  }

  const G4double ePhi = fSPhi + fDPhi;
  const G4double cPhi = fSPhi + 0.5 * fDPhi;
  sinSPhi = std::sin(fSPhi);
  cosSPhi = std::cos(fSPhi);
  sinEPhi = std::sin(ePhi);
  cosEPhi = std::cos(ePhi);
  sinCPhi = std::sin(cPhi);
  cosCPhi = std::cos(cPhi);
  cosHDPhi = std::cos(0.5 * fDPhi);
  cosHDPhiIT = std::cos(0.5 * fDPhi - halfAngTolerance);
  cosHDPhiOT = std::cos(0.5 * fDPhi + halfAngTolerance);
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Torus::G4Torus(__void__& a) : G4CSGSolid(a) {}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Torus& G4Torus::operator=(const G4Torus& rhs)
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
  fRmin = rhs.fRmin;
  fRmax = rhs.fRmax;
  fRtor = rhs.fRtor;
  fSPhi = rhs.fSPhi;
  fDPhi = rhs.fDPhi;
  fRminTolerance = rhs.fRminTolerance;
  fRmaxTolerance = rhs.fRmaxTolerance;
  kRadTolerance = rhs.kRadTolerance;
  kAngTolerance = rhs.kAngTolerance;
  halfCarTolerance = rhs.halfCarTolerance;
  halfAngTolerance = rhs.halfAngTolerance;
  sinSPhi = rhs.sinSPhi;
  cosSPhi = rhs.cosSPhi;
  sinEPhi = rhs.sinEPhi;
  cosEPhi = rhs.cosEPhi;
  sinCPhi = rhs.sinCPhi;
  cosCPhi = rhs.cosCPhi;
  cosHDPhi = rhs.cosHDPhi;
  cosHDPhiIT = rhs.cosHDPhiIT;
  cosHDPhiOT = rhs.cosHDPhiOT;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Torus::ComputeDimensions(G4VPVParameterisation* p, const G4int n,
                                const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this, n, pRep);
}

////////////////////////////////////////////////////////////////////////////////
//
// Calculate the real roots to torus surface.
// Returns negative solutions as well.

G4int G4Torus::TorusRootsJT(const G4ThreeVector& p, const G4ThreeVector& v, G4double r,
                            std::array<G4double, 4>& roots) const
{
  G4int nRoots = 0;
  G4double c[5], srd[4], si[4];

  G4double Rtor2 = fRtor * fRtor, r2 = r * r;

  G4double pDotV = p.x() * v.x() + p.y() * v.y() + p.z() * v.z();
  G4double pRad2 = p.x() * p.x() + p.y() * p.y() + p.z() * p.z();

  G4double d = pRad2 - Rtor2;
  c[0] = 1.0;
  c[1] = 4 * pDotV;
  c[2] = 2 * ((d + 2 * pDotV * pDotV - r2) + 2 * Rtor2 * v.z() * v.z());
  c[3] = 4 * (pDotV * (d - r2) + 2 * Rtor2 * p.z() * v.z());
  c[4] = (d - r2) * (d - r2) + 4 * Rtor2 * (p.z() * p.z() - r2);

  G4JTPolynomialSolver torusEq;
  const G4int num = torusEq.FindRoots(c, 4, srd, si);

  for (G4int i = 0; i < num; ++i)
  {
    if (si[i] == 0.)
    {
      roots[nRoots++] = srd[i];
    }  // store real roots
  }

  std::sort(roots.begin(), roots.begin() + nRoots);  // sorting with <
  return nRoots;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate the real roots using a scale-normalised analytical quartic.
// Newton polishing and a backward-error check protect the fast path; the
// Jenkins-Traub solver is retained for ill-conditioned cases.

G4int G4Torus::TorusRoots(const G4ThreeVector& p, const G4ThreeVector& v, G4double r,
                          std::array<G4double, 4>& roots) const
{
  const G4double scale =
    std::max({fRtor + r, std::fabs(p.x()), std::fabs(p.y()), std::fabs(p.z()), 1. * mm});
  const G4ThreeVector scaledP = p / scale;
  const G4double Rtor = fRtor / scale;
  const G4double rs = r / scale;
  const G4double Rtor2 = Rtor * Rtor;
  const G4double r2 = rs * rs;
  const G4double pDotV = scaledP * v;
  const G4double pRad2 = scaledP.mag2();
  const G4double d = pRad2 - Rtor2;
  G4double c[5] = {1., 4. * pDotV, 2. * (d + 2. * pDotV * pDotV - r2 + 2. * Rtor2 * v.z() * v.z()),
                   4. * (pDotV * (d - r2) + 2. * Rtor2 * scaledP.z() * v.z()),
                   (d - r2) * (d - r2) + 4. * Rtor2 * (scaledP.z() * scaledP.z() - r2)};
  std::array<G4double, 4> solution;
  const G4int nSolutions = SolveQuartic(c[1], c[2], c[3], c[4], solution.data());

  G4int nRoots = 0;
  G4bool valid = true;
  for (G4int i = 0; i < nSolutions; ++i)
  {
    G4double root = solution[i];
    if (!std::isfinite(root))
    {
      valid = false;
      break;
    }

    for (G4int iteration = 0; iteration < 2; ++iteration)
    {
      const G4double value = (((root + c[1]) * root + c[2]) * root + c[3]) * root + c[4];
      const G4double derivative = ((4. * root + 3. * c[1]) * root + 2. * c[2]) * root + c[3];
      if (derivative == 0.) break;
      root -= value / derivative;
    }

    const G4double absRoot = std::fabs(root);
    const G4double value = (((root + c[1]) * root + c[2]) * root + c[3]) * root + c[4];
    const G4double magnitude =
      (((absRoot + std::fabs(c[1])) * absRoot + std::fabs(c[2])) * absRoot + std::fabs(c[3]))
        * absRoot
      + std::fabs(c[4]);
    if (!std::isfinite(root) || std::fabs(value) > 256. * DBL_EPSILON * magnitude)
    {
      valid = false;
      break;
    }
    roots[nRoots++] = root * scale;
  }

  if (!valid || (nRoots & 1) != 0) return TorusRootsJT(p, v, r, roots);
  std::sort(roots.begin(), roots.begin() + nRoots);
  return nRoots;
}

//////////////////////////////////////////////////////////////////////////////
//
// Interface for DistanceToIn and DistanceToOut.
// Calls TorusRoots and returns the smallest possible distance to
// the surface.
// Attention: Difference in DistanceToIn/Out for points p on the surface.

G4double G4Torus::SolveNumeric(const G4ThreeVector& p, const G4ThreeVector& v, G4double r,
                               G4bool IsDistanceToIn) const
{
  G4double tmin = kInfinity;
  G4double t, scal;

  // calculate the distances to the intersections with the Torus
  // from a given point p and direction v.
  //
  std::array<G4double, 4> roots;
  const G4int nRoots = TorusRoots(p, v, r, roots);

  G4ThreeVector ptmp;

  // determine the smallest non-negative solution
  //
  for (G4int k = 0; k < nRoots; ++k)
  {
    t = roots[k];

    if (t < -halfCarTolerance)
    {
      continue;
    }  // skip negative roots

    ptmp = p + t * v;  // calculate the position of the proposed intersection

    // We have to verify if this root is inside the region between
    // fSPhi and fSPhi + fDPhi
    //
    const G4double rhoi = std::hypot(ptmp.x(), ptmp.y());
    if ((fDPhi == twopi) || (rhoi == 0.)
        || ((ptmp.x() * cosCPhi + ptmp.y() * sinCPhi) / rhoi >= cosHDPhiOT))
    {
      // check if P is on the surface, and called from DistanceToIn
      // DistanceToIn has to return 0.0 if particle is going inside the solid

      if (IsDistanceToIn)
      {
        if (std::fabs(t) < halfCarTolerance)
        {
          // compute scalar product at position p : v.n
          // ( n taken from SurfaceNormal, not normalized )

          scal = v
                 * G4ThreeVector(p.x() * (1 - fRtor / std::hypot(p.x(), p.y())),
                                 p.y() * (1 - fRtor / std::hypot(p.x(), p.y())), p.z());

          // change sign in case of inner radius
          //
          if (r == GetRmin())
          {
            scal = -scal;
          }
          if (scal < 0)
          {
            return 0.0;
          }
        }
      }

      // check if P is on the surface, and called from DistanceToOut
      // DistanceToIn has to return 0.0 if particle is leaving the solid

      if (!IsDistanceToIn)
      {
        if (std::fabs(t) < halfCarTolerance)
        {
          // compute scalar product at position p : v.n
          //
          scal = v
                 * G4ThreeVector(p.x() * (1 - fRtor / std::hypot(p.x(), p.y())),
                                 p.y() * (1 - fRtor / std::hypot(p.x(), p.y())), p.z());

          // change sign in case of inner radius
          //
          if (r == GetRmin())
          {
            scal = -scal;
          }
          if (scal > 0)
          {
            return 0.0;
          }
        }
      }

      // check if distance is larger than 1/2 kCarTolerance
      //
      if (t > halfCarTolerance)
      {
        tmin = t;
        return tmin;
      }
    }
  }

  return tmin;
}

/////////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Torus::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double rmax = GetRmax();
  G4double rtor = GetRtor();
  G4double rint = rtor - rmax;
  G4double rext = rtor + rmax;
  G4double dz = rmax;

  // Find bounding box
  //
  if (GetDPhi() >= twopi)
  {
    pMin.set(-rext, -rext, -dz);
    pMax.set(rext, rext, dz);
  }
  else
  {
    G4TwoVector vmin, vmax;
    G4GeomTools::DiskExtent(rint, rext, GetSinStartPhi(), GetCosStartPhi(), GetSinEndPhi(),
                            GetCosEndPhi(), vmin, vmax);
    pMin.set(vmin.x(), vmin.y(), -dz);
    pMax.set(vmax.x(), vmax.y(), dz);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " << GetName() << " !"
            << "\npMin = " << pMin << "\npMax = " << pMax;
    G4Exception("G4Torus::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Torus::CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                const G4AffineTransform& pTransform, G4double& pMin,
                                G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Get bounding box
  BoundingLimits(bmin, bmax);

  // Check bounding box
  G4BoundingEnvelope bbox(bmin, bmax);
#  ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
#  endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis, pVoxelLimit, pTransform, pMin, pMax))
  {
    return exist = pMin < pMax;
  }

  // Get parameters of the solid
  G4double rmin = GetRmin();
  G4double rmax = GetRmax();
  G4double rtor = GetRtor();
  G4double dphi = GetDPhi();
  G4double sinStart = GetSinStartPhi();
  G4double cosStart = GetCosStartPhi();
  G4double sinEnd = GetSinEndPhi();
  G4double cosEnd = GetCosEndPhi();
  G4double rint = rtor - rmax;
  G4double rext = rtor + rmax;

  // Find bounding envelope and calculate extent
  //
  static const G4int NPHI = 24;  // number of steps for whole torus
  static const G4int NDISK = 16;  // number of steps for disk
  static const G4double sinHalfDisk = std::sin(pi / NDISK);
  static const G4double cosHalfDisk = std::cos(pi / NDISK);
  static const G4double sinStepDisk = 2. * sinHalfDisk * cosHalfDisk;
  static const G4double cosStepDisk = 1. - 2. * sinHalfDisk * sinHalfDisk;

  G4double astep = (360 / NPHI) * deg;  // max angle for one slice in phi
  G4int kphi = (dphi <= astep) ? 1 : (G4int)((dphi - deg) / astep) + 1;
  G4double ang = dphi / kphi;

  G4double sinHalf = std::sin(0.5 * ang);
  G4double cosHalf = std::cos(0.5 * ang);
  G4double sinStep = 2. * sinHalf * cosHalf;
  G4double cosStep = 1. - 2. * sinHalf * sinHalf;

  // define vectors for bounding envelope
  G4ThreeVectorList pols[NDISK + 1];
  for (auto& pol : pols)
  {
    pol.resize(4);
  }

  std::vector<const G4ThreeVectorList*> polygons;
  polygons.resize(NDISK + 1);
  for (G4int k = 0; k < NDISK + 1; ++k)
  {
    polygons[k] = &pols[k];
  }

  // set internal and external reference circles
  G4TwoVector rzmin[NDISK];
  G4TwoVector rzmax[NDISK];

  if ((rtor - rmin * sinHalfDisk) / cosHalf > (rtor + rmin * sinHalfDisk))
  {
    rmin = 0;
  }
  rmax /= cosHalfDisk;
  G4double sinCurDisk = sinHalfDisk;
  G4double cosCurDisk = cosHalfDisk;
  for (G4int k = 0; k < NDISK; ++k)
  {
    G4double rmincur = rtor + rmin * cosCurDisk;
    if (cosCurDisk < 0 && rmin > 0)
    {
      rmincur /= cosHalf;
    }
    rzmin[k].set(rmincur, rmin * sinCurDisk);

    G4double rmaxcur = rtor + rmax * cosCurDisk;
    if (cosCurDisk > 0)
    {
      rmaxcur /= cosHalf;
    }
    rzmax[k].set(rmaxcur, rmax * sinCurDisk);

    G4double sinTmpDisk = sinCurDisk;
    sinCurDisk = sinCurDisk * cosStepDisk + cosCurDisk * sinStepDisk;
    cosCurDisk = cosCurDisk * cosStepDisk - sinTmpDisk * sinStepDisk;
  }

  // Loop along slices in Phi. The extent is calculated as cumulative
  // extent of the slices
  pMin = kInfinity;
  pMax = -kInfinity;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);
  G4double sinCur1 = 0, cosCur1 = 0, sinCur2 = 0, cosCur2 = 0;
  for (G4int i = 0; i < kphi + 1; ++i)
  {
    if (i == 0)
    {
      sinCur1 = sinStart;
      cosCur1 = cosStart;
      sinCur2 = sinCur1 * cosHalf + cosCur1 * sinHalf;
      cosCur2 = cosCur1 * cosHalf - sinCur1 * sinHalf;
    }
    else
    {
      sinCur1 = sinCur2;
      cosCur1 = cosCur2;
      sinCur2 = (i == kphi) ? sinEnd : sinCur1 * cosStep + cosCur1 * sinStep;
      cosCur2 = (i == kphi) ? cosEnd : cosCur1 * cosStep - sinCur1 * sinStep;
    }
    for (G4int k = 0; k < NDISK; ++k)
    {
      G4double r1 = rzmin[k].x(), r2 = rzmax[k].x();
      G4double z1 = rzmin[k].y(), z2 = rzmax[k].y();
      pols[k][0].set(r1 * cosCur1, r1 * sinCur1, z1);
      pols[k][1].set(r2 * cosCur1, r2 * sinCur1, z2);
      pols[k][2].set(r2 * cosCur2, r2 * sinCur2, z2);
      pols[k][3].set(r1 * cosCur2, r1 * sinCur2, z1);
    }
    pols[NDISK] = pols[0];

    // get bounding box of current slice
    G4TwoVector vmin, vmax;
    G4GeomTools::DiskExtent(rint, rext, sinCur1, cosCur1, sinCur2, cosCur2, vmin, vmax);
    bmin.setX(vmin.x());
    bmin.setY(vmin.y());
    bmax.setX(vmax.x());
    bmax.setY(vmax.y());

    // set bounding envelope for current slice and adjust extent
    G4double emin, emax;
    G4BoundingEnvelope benv(bmin, bmax, polygons);
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
      break;
    }  // max possible extent
  }
  return (pMin < pMax);
}

//////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Torus::Inside(const G4ThreeVector& p) const
{
  G4double r, pt2, tolRMin, tolRMax;

  EInside in = kOutside;

  // General precals
  //
  r = std::hypot(p.x(), p.y());
  pt2 = p.z() * p.z() + (r - fRtor) * (r - fRtor);

  (fRmin != 0.0) ? (tolRMin = fRmin + fRminTolerance) : (tolRMin = 0);

  tolRMax = fRmax - fRmaxTolerance;

  if (pt2 >= tolRMin * tolRMin && pt2 <= tolRMax * tolRMax)
  {
    if (fDPhi == twopi || pt2 == 0)  // on torus swept axis
    {
      in = kInside;
    }
    else
    {
      // Try inner tolerant phi boundaries (=>inside)
      // if not inside, try outer tolerant phi boundaries

      const G4double cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / r;
      if (cosPsi >= cosHDPhiIT)
      {
        in = kInside;
      }
      else if (cosPsi >= cosHDPhiOT)
      {
        in = kSurface;
      }
    }
  }
  else  // Try generous boundaries
  {
    tolRMin = fRmin - fRminTolerance;
    tolRMax = fRmax + fRmaxTolerance;

    if (tolRMin < 0)
    {
      tolRMin = 0;
    }

    if ((pt2 >= tolRMin * tolRMin) && (pt2 <= tolRMax * tolRMax))
    {
      if ((fDPhi == twopi) || (pt2 == 0))  // Continuous in phi or on z-axis
      {
        in = kSurface;
      }
      else  // Try outer tolerant phi boundaries only
      {
        const G4double cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / r;
        if (cosPsi >= cosHDPhiOT)
        {
          in = kSurface;
        }
      }
    }
  }
  return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Torus::SurfaceNormal(const G4ThreeVector& p) const
{
  G4int noSurfaces = 0;
  G4double rho, pt;
  G4double distRMin = kInfinity;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;

  // To cope with precision loss
  //
  const G4double delta = std::max(10.0 * kCarTolerance, 1.0e-8 * (fRtor + fRmax));
  const G4double dAngle = 10.0 * kAngTolerance;

  G4ThreeVector nR, nPs, nPe;
  G4ThreeVector norm, sumnorm(0., 0., 0.);

  rho = std::hypot(p.x(), p.y());
  pt = std::hypot(p.z(), rho - fRtor);

  G4double distRMax = std::fabs(pt - fRmax);
  if (fRmin != 0.0)
  {
    distRMin = std::fabs(pt - fRmin);
  }

  if (rho > delta && pt != 0.0)
  {
    G4double redFactor = (rho - fRtor) / rho;
    nR = G4ThreeVector(p.x() * redFactor,  // p.x()*(1.-fRtor/rho),
                       p.y() * redFactor,  // p.y()*(1.-fRtor/rho),
                       p.z());
    nR *= 1.0 / pt;
  }

  if (fDPhi < twopi)
  {
    if (rho != 0.0)
    {
      if (p.x() * cosSPhi + p.y() * sinSPhi >= 0.)
      {
        distSPhi = std::fabs(p.x() * sinSPhi - p.y() * cosSPhi) / rho;
      }
      if (p.x() * cosEPhi + p.y() * sinEPhi >= 0.)
      {
        distEPhi = std::fabs(p.x() * sinEPhi - p.y() * cosEPhi) / rho;
      }
    }
    nPs = G4ThreeVector(sinSPhi, -cosSPhi, 0);
    nPe = G4ThreeVector(-sinEPhi, cosEPhi, 0);
  }
  if (distRMax <= delta)
  {
    ++noSurfaces;
    sumnorm += nR;
  }
  else if ((fRmin != 0.0) && (distRMin <= delta))  // Must not be on both Outer and Inner
  {
    ++noSurfaces;
    sumnorm -= nR;
  }

  //  To be on one of the 'phi' surfaces,
  //  it must be within the 'tube' - with tolerance

  if ((fDPhi < twopi) && (fRmin - delta <= pt) && (pt <= (fRmax + delta)))
  {
    if (distSPhi <= dAngle)
    {
      ++noSurfaces;
      sumnorm += nPs;
    }
    if (distEPhi <= dAngle)
    {
      ++noSurfaces;
      sumnorm += nPe;
    }
  }
  if (noSurfaces == 0)
  {
#  ifdef G4CSGDEBUG
    G4ExceptionDescription ed;
    ed.precision(16);

    EInside inIt = Inside(p);

    if (inIt != kSurface)
    {
      ed << " ERROR>  Surface Normal was called for Torus,"
         << " with point not on surface." << G4endl;
    }
    else
    {
      ed << " ERROR>  Surface Normal has not found a surface, "
         << " despite the point being on the surface. " << G4endl;
    }

    if (inIt != kInside)
    {
      ed << " Safety (Dist To In)  = " << DistanceToIn(p) << G4endl;
    }
    if (inIt != kOutside)
    {
      ed << " Safety (Dist to Out) = " << DistanceToOut(p) << G4endl;
    }
    ed << " Coordinates of point : " << p << G4endl;
    ed << " Parameters  of solid : " << G4endl << *this << G4endl;

    if (inIt == kSurface)
    {
      G4Exception("G4Torus::SurfaceNormal(p)", "GeomSolids1002", JustWarning, ed,
                  "Failing to find normal, even though point is on surface!");
    }
    else
    {
      static const char* NameInside[3] = {"Inside", "Surface", "Outside"};
      ed << "  The point is " << NameInside[inIt] << " the solid. " << G4endl;
      G4Exception("G4Torus::SurfaceNormal(p)", "GeomSolids1002", JustWarning, ed,
                  "Point p is not on surface !?");
    }
#  endif
    norm = ApproxSurfaceNormal(p);
  }
  else if (noSurfaces == 1)
  {
    norm = sumnorm;
  }
  else
  {
    norm = sumnorm.unit();
  }

  return norm;
}

//////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Torus::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  ENorm side;
  G4ThreeVector norm;
  G4double rho, pt;
  G4double distRMin, distRMax, distSPhi, distEPhi, distMin;

  rho = std::hypot(p.x(), p.y());
  pt = std::hypot(p.z(), rho - fRtor);

#  ifdef G4CSGDEBUG
  G4cout << " G4Torus::ApproximateSurfaceNormal called for point " << p << G4endl;
#  endif

  distRMax = std::fabs(pt - fRmax);

  if (fRmin != 0.0)  // First minimum radius
  {
    distRMin = std::fabs(pt - fRmin);

    if (distRMin < distRMax)
    {
      distMin = distRMin;
      side = kNRMin;
    }
    else
    {
      distMin = distRMax;
      side = kNRMax;
    }
  }
  else
  {
    distMin = distRMax;
    side = kNRMax;
  }
  if ((fDPhi < twopi) && (rho != 0.0))
  {
    distSPhi = kInfinity;
    distEPhi = kInfinity;
    if (p.x() * cosSPhi + p.y() * sinSPhi >= 0.)
    {
      distSPhi = std::fabs(p.x() * sinSPhi - p.y() * cosSPhi);
    }
    if (p.x() * cosEPhi + p.y() * sinEPhi >= 0.)
    {
      distEPhi = std::fabs(p.x() * sinEPhi - p.y() * cosEPhi);
    }

    if (distSPhi < distEPhi)  // Find new minimum
    {
      if (distSPhi < distMin)
      {
        side = kNSPhi;
      }
    }
    else
    {
      if (distEPhi < distMin)
      {
        side = kNEPhi;
      }
    }
  }
  if (rho == 0.0)
  {
    norm = G4ThreeVector(fRtor, 0, p.z()).unit();
  }  // x,y = 0
  if (pt == 0.0)
  {
    norm = p.unit();
  }  // rho = Rtor, z = 0
  switch (side)
  {
    case kNRMin:  // Inner radius
      if (rho == 0.0 || pt == 0.0) break;
      norm = G4ThreeVector(-p.x() * (1 - fRtor / rho) / pt, -p.y() * (1 - fRtor / rho) / pt,
                           -p.z() / pt);
      break;
    case kNRMax:  // Outer radius
      if (rho == 0.0 || pt == 0.0) break;
      norm =
        G4ThreeVector(p.x() * (1 - fRtor / rho) / pt, p.y() * (1 - fRtor / rho) / pt, p.z() / pt);
      break;
    case kNSPhi:
      norm = G4ThreeVector(sinSPhi, -cosSPhi, 0);
      break;
    case kNEPhi:
      norm = G4ThreeVector(-sinEPhi, cosEPhi, 0);
      break;
    default:  // Should never reach this case ...
      DumpInfo();
      G4Exception("G4Torus::ApproxSurfaceNormal()", "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break;
  }
  return norm;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points

G4double G4Torus::DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const
{
  // Get bounding box of full torus
  //
  G4double boxDx = fRtor + fRmax;
  G4double boxDy = boxDx;
  G4double boxDz = fRmax;
  G4double boxMax = boxDx;
  G4double boxMin = boxDz;

  // Check if point is traveling away
  //
  G4double distX = std::abs(p.x()) - boxDx;
  G4double distY = std::abs(p.y()) - boxDy;
  G4double distZ = std::abs(p.z()) - boxDz;
  if (distX >= -halfCarTolerance && p.x() * v.x() >= 0)
  {
    return kInfinity;
  }
  if (distY >= -halfCarTolerance && p.y() * v.y() >= 0)
  {
    return kInfinity;
  }
  if (distZ >= -halfCarTolerance && p.z() * v.z() >= 0)
  {
    return kInfinity;
  }

  // Calculate safety distance to bounding box
  // If point is too far, move it closer and calculate distance
  //
  G4double Dmax = 32 * boxMax;
  G4double safe = std::max(std::max(distX, distY), distZ);
  if (safe > Dmax)
  {
    G4double dist = safe - 1.e-8 * safe - boxMin;  // stay outside after the move
    dist += DistanceToIn(p + dist * v, v);
    return (dist >= kInfinity) ? kInfinity : dist;
  }

  // Find intersection with torus
  //
  G4double snxt = kInfinity, sphi = kInfinity;  // snxt = default return value

  // Precalculated trig for phi intersections - used by r,z intersections to
  //                                            check validity

  const G4bool seg = fDPhi < twopi;

  G4double tolORMin2;  // `generous' radii squared
  G4double tolORMax2;

  G4double Dist, xi, yi, zi, rhoi, it2;  // Intersection point variables

  G4double Comp;
  if (fRmin > fRminTolerance)  // Calculate tolerant rmin and rmax
  {
    tolORMin2 = (fRmin - fRminTolerance) * (fRmin - fRminTolerance);
  }
  else
  {
    tolORMin2 = 0;
  }
  tolORMax2 = (fRmax + fRmaxTolerance) * (fRmax + fRmaxTolerance);

  // Intersection with Rmax (possible return) and Rmin (must also check phi)

  snxt = SolveNumeric(p, v, fRmax, true);

  if (fRmin != 0.0)  // Possible Rmin intersection
  {
    const G4double distanceToInner = SolveNumeric(p, v, fRmin, true);
    if (distanceToInner < snxt)
    {
      snxt = distanceToInner;
    }
  }

  //
  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?

  if (seg)
  {
    Comp = v.x() * sinSPhi - v.y() * cosSPhi;  // Component in outwards
                                               // normal direction
    if (Comp < 0)
    {
      Dist = (p.y() * cosSPhi - p.x() * sinSPhi);

      if (Dist < halfCarTolerance)
      {
        sphi = Dist / Comp;
        if (sphi < snxt)
        {
          if (sphi < 0)
          {
            sphi = 0;
          }

          xi = p.x() + sphi * v.x();
          yi = p.y() + sphi * v.y();
          zi = p.z() + sphi * v.z();
          rhoi = std::hypot(xi, yi);
          it2 = zi * zi + (rhoi - fRtor) * (rhoi - fRtor);

          if (it2 >= tolORMin2 && it2 <= tolORMax2)
          {
            // r intersection is good - check intersecting
            // with correct half-plane
            //
            if ((yi * cosCPhi - xi * sinCPhi) <= 0)
            {
              snxt = sphi;
            }
          }
        }
      }
    }
    Comp = -(v.x() * sinEPhi - v.y() * cosEPhi);

    if (Comp < 0)  // Component in outwards normal dirn
    {
      Dist = -(p.y() * cosEPhi - p.x() * sinEPhi);

      if (Dist < halfCarTolerance)
      {
        sphi = Dist / Comp;

        if (sphi < snxt)
        {
          if (sphi < 0)
          {
            sphi = 0;
          }

          xi = p.x() + sphi * v.x();
          yi = p.y() + sphi * v.y();
          zi = p.z() + sphi * v.z();
          rhoi = std::hypot(xi, yi);
          it2 = zi * zi + (rhoi - fRtor) * (rhoi - fRtor);

          if (it2 >= tolORMin2 && it2 <= tolORMax2)
          {
            // z and r intersections good - check intersecting
            // with correct half-plane
            //
            if ((yi * cosCPhi - xi * sinCPhi) >= 0)
            {
              snxt = sphi;
            }
          }
        }
      }
    }
  }
  if (snxt < halfCarTolerance)
  {
    snxt = 0.0;
  }

  return snxt;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Torus::DistanceToIn(const G4ThreeVector& p) const
{
  G4double safePhi;

  const G4double rho = std::hypot(p.x(), p.y());
  const G4double pt = std::hypot(p.z(), rho - fRtor);
  G4double safe = std::max(fRmin - pt, pt - fRmax);

  if (fDPhi < twopi && (rho != 0.0))
  {
    const G4double cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / rho;

    if (cosPsi < cosHDPhi)  // Point lies outside phi range
    {
      if ((p.y() * cosCPhi - p.x() * sinCPhi) <= 0)
      {
        safePhi = std::fabs(p.x() * sinSPhi - p.y() * cosSPhi);
      }
      else
      {
        safePhi = std::fabs(p.x() * sinEPhi - p.y() * cosEPhi);
      }
      safe = std::max(safe, safePhi);
    }
  }
  return std::max(0., safe);
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection
//

G4double G4Torus::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                                const G4bool calcNorm, G4bool* validNorm, G4ThreeVector* n) const
{
  ESide side = kNull, sidephi = kNull;
  G4double snxt = kInfinity, sphi;

  // Vars for phi intersection
  //
  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, zi;

  // Radial intersections and general precalculations

  G4double rho = std::hypot(p.x(), p.y());
  G4double pt = hypot(p.z(), rho - fRtor);

  G4double pDotV = p.x() * v.x() + p.y() * v.y() + p.z() * v.z();

  G4double tolRMax = fRmax - fRmaxTolerance;

  G4double vDotNmax = pDotV - fRtor * (v.x() * p.x() + v.y() * p.y()) / rho;
  G4double pDotxyNmax = (1 - fRtor / rho);

  if ((pt * pt > tolRMax * tolRMax) && (vDotNmax >= 0))
  {
    // On tolerant boundary & heading outwards (or perpendicular to) outer
    // radial surface -> leaving immediately with *n for really convex part
    // only

    if (calcNorm && (pDotxyNmax >= -2. * fRmaxTolerance))
    {
      *n =
        G4ThreeVector(p.x() * (1 - fRtor / rho) / pt, p.y() * (1 - fRtor / rho) / pt, p.z() / pt);
      *validNorm = true;
    }

    return snxt = 0;  // Leaving by Rmax immediately
  }

  snxt = SolveNumeric(p, v, fRmax, false);
  side = kRMax;

  // rmin

  if (fRmin != 0.0)
  {
    G4double tolRMin = fRmin + fRminTolerance;

    if ((pt * pt < tolRMin * tolRMin) && (vDotNmax < 0))
    {
      if (calcNorm)
      {
        *validNorm = false;
      }  // Concave surface of the torus
      return snxt = 0;  // Leaving by Rmin immediately
    }

    const G4double distanceToInner = SolveNumeric(p, v, fRmin, false);
    if (distanceToInner < snxt)
    {
      snxt = distanceToInner;
      side = kRMin;
    }
  }

  if (fDPhi < twopi)  // Phi Intersections
  {
    const G4double vRho = std::hypot(v.x(), v.y());
    const G4bool directionInPhi =
      (vRho == 0.) || ((v.x() * cosCPhi + v.y() * sinCPhi) / vRho >= cosHDPhiOT);

    if ((p.x() != 0.0) || (p.y() != 0.0))  // Check if on z axis (rho not needed later)
    {
      pDistS = p.x() * sinSPhi - p.y() * cosSPhi;  // pDist -ve when inside
      pDistE = -p.x() * sinEPhi + p.y() * cosEPhi;

      // Comp -ve when in direction of outwards normal
      //
      compS = -sinSPhi * v.x() + cosSPhi * v.y();
      compE = sinEPhi * v.x() - cosEPhi * v.y();
      sidephi = kNull;

      if (((fDPhi <= pi) && ((pDistS <= halfCarTolerance) && (pDistE <= halfCarTolerance)))
          || ((fDPhi > pi) && ((pDistS <= halfCarTolerance) || (pDistE <= halfCarTolerance))))
      {
        // Inside both phi *full* planes

        if (compS < 0)
        {
          sphi = pDistS / compS;

          if (sphi >= -halfCarTolerance)
          {
            xi = p.x() + sphi * v.x();
            yi = p.y() + sphi * v.y();

            // Check intersecting with correct half-plane
            // (if not -> no intersect)
            //
            if ((std::fabs(xi) <= kCarTolerance) && (std::fabs(yi) <= kCarTolerance))
            {
              sidephi = kSPhi;
              if (directionInPhi)
              {
                sphi = kInfinity;
              }
            }
            else if (yi * cosCPhi - xi * sinCPhi >= 0)
            {
              sphi = kInfinity;
            }
            else
            {
              sidephi = kSPhi;
            }
          }
          else
          {
            sphi = kInfinity;
          }
        }
        else
        {
          sphi = kInfinity;
        }

        if (compE < 0)
        {
          sphi2 = pDistE / compE;

          // Only check further if < starting phi intersection
          //
          if ((sphi2 > -kCarTolerance) && (sphi2 < sphi))
          {
            xi = p.x() + sphi2 * v.x();
            yi = p.y() + sphi2 * v.y();

            if ((std::fabs(xi) <= kCarTolerance) && (std::fabs(yi) <= kCarTolerance))
            {
              // Leaving via ending phi
              //
              if (!directionInPhi)
              {
                sidephi = kEPhi;
                sphi = sphi2;
              }
            }
            else  // Check intersecting with correct half-plane
            {
              if ((yi * cosCPhi - xi * sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi;
                sphi = sphi2;
              }
            }
          }
        }
      }
      else
      {
        sphi = kInfinity;
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if (directionInPhi)
      {
        sphi = kInfinity;
      }
      else
      {
        sidephi = kSPhi;  // arbitrary
        sphi = 0;
      }
    }

    // Order intersections

    if (sphi < snxt)
    {
      snxt = sphi;
      side = sidephi;
    }
  }

  G4double rhoi, it, iDotxyNmax;
  // Note: by numerical computation we know where the ray hits the torus
  // So I propose to return the side where the ray hits

  if (calcNorm)
  {
    switch (side)
    {
      case kRMax:  // n is unit vector
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        zi = p.z() + snxt * v.z();
        rhoi = std::hypot(xi, yi);
        it = hypot(zi, rhoi - fRtor);

        iDotxyNmax = (1 - fRtor / rhoi);
        if (iDotxyNmax >= -2. * fRmaxTolerance)  // really convex part of Rmax
        {
          *n = G4ThreeVector(xi * (1 - fRtor / rhoi) / it, yi * (1 - fRtor / rhoi) / it, zi / it);
          *validNorm = true;
        }
        else
        {
          *validNorm = false;  // concave-convex part of Rmax
        }
        break;

      case kRMin:
        *validNorm = false;  // Rmin is concave or concave-convex
        break;

      case kSPhi:
        if (fDPhi <= pi)
        {
          *n = G4ThreeVector(sinSPhi, -cosSPhi, 0);
          *validNorm = true;
        }
        else
        {
          *validNorm = false;
        }
        break;

      case kEPhi:
        if (fDPhi <= pi)
        {
          *n = G4ThreeVector(-sinEPhi, cosEPhi, 0);
          *validNorm = true;
        }
        else
        {
          *validNorm = false;
        }
        break;

      default:

        // It seems we go here from time to time ...

        G4cout << G4endl;
        DumpInfo();
        std::ostringstream message;
        G4long oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid." << G4endl
                << "Position:" << G4endl << G4endl << "p.x() = " << p.x() / mm << " mm" << G4endl
                << "p.y() = " << p.y() / mm << " mm" << G4endl << "p.z() = " << p.z() / mm << " mm"
                << G4endl << G4endl << "Direction:" << G4endl << G4endl << "v.x() = " << v.x()
                << G4endl << "v.y() = " << v.y() << G4endl << "v.z() = " << v.z() << G4endl
                << G4endl << "Proposed distance :" << G4endl << G4endl << "snxt = " << snxt / mm
                << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Torus::DistanceToOut(p,v,..)", "GeomSolids1002", JustWarning, message);
        break;
    }
  }
  if (snxt < halfCarTolerance)
  {
    snxt = 0;
  }

  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Torus::DistanceToOut(const G4ThreeVector& p) const
{
  G4double safePhi;

  const G4double rho = std::hypot(p.x(), p.y());
  const G4double pt = std::hypot(p.z(), rho - fRtor);

#  ifdef G4CSGDEBUG
  if (Inside(p) == kOutside)
  {
    G4long oldprc = G4cout.precision(16);
    G4cout << G4endl;
    DumpInfo();
    G4cout << "Position:" << G4endl << G4endl;
    G4cout << "p.x() = " << p.x() / mm << " mm" << G4endl;
    G4cout << "p.y() = " << p.y() / mm << " mm" << G4endl;
    G4cout << "p.z() = " << p.z() / mm << " mm" << G4endl << G4endl;
    G4cout.precision(oldprc);
    G4Exception("G4Torus::DistanceToOut(p)", "GeomSolids1002", JustWarning,
                "Point p is outside !?");
  }
#  endif

  G4double safe = (fRmin != 0.0) ? std::min(pt - fRmin, fRmax - pt) : fRmax - pt;

  // Check if phi divided, Calc distances closest phi plane
  //
  if (fDPhi < twopi)  // Above/below central phi of Torus?
  {
    if ((p.y() * cosCPhi - p.x() * sinCPhi) <= 0)
    {
      safePhi = -(p.x() * sinSPhi - p.y() * cosSPhi);
    }
    else
    {
      safePhi = p.x() * sinEPhi - p.y() * cosEPhi;
    }
    safe = std::min(safe, safePhi);
  }
  return std::max(0., safe);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

G4GeometryType G4Torus::GetEntityType() const
{
  return {"G4Torus"};
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Torus::Clone() const
{
  return new G4Torus(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Torus::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Torus\n"
     << " Parameters: \n"
     << "    inner radius: " << fRmin / mm << " mm \n"
     << "    outer radius: " << fRmax / mm << " mm \n"
     << "    swept radius: " << fRtor / mm << " mm \n"
     << "    starting phi: " << fSPhi / degree << " degrees \n"
     << "    delta phi   : " << fDPhi / degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Torus::GetPointOnSurface() const
{
  G4double rrmin = fRmin * fRmin;
  G4double rrmax = fRmax * fRmax;
  G4double smax = twopi * fRtor * fDPhi * fRmax;
  G4double smin = twopi * fRtor * fDPhi * fRmin;
  G4double sphi = (fDPhi >= twopi) ? 0. : pi * (rrmax - rrmin);

  G4double u = G4QuickRand();
  G4double v = twopi * G4QuickRand();
  G4double select = (smax + smin + 2. * sphi) * G4QuickRand();
  G4double phi, r, ds;
  if (select < 2. * sphi)
  {
    // phi cut
    phi = fSPhi + fDPhi * (G4double)(select < sphi);
    r = std::sqrt(rrmin + (rrmax - rrmin) * u);
    ds = fRtor + r * std::cos(v);
  }
  else
  {
    // toroidal surface (rejection sampling)
    phi = fSPhi + fDPhi * u;
    r = (select < 2. * sphi + smax) ? fRmax : fRmin;
    ds = fRtor + r * std::cos(v);
    for (auto i = 0; i < 10; ++i)
    {
      if ((fRtor + r) * G4QuickRand() < ds)
      {
        break;
      }
      v = twopi * G4QuickRand();
      ds = fRtor + r * std::cos(v);
    }
  }
  return {ds * std::cos(phi), ds * std::sin(phi), r * std::sin(v)};
}

////////////////////////////////////////////////////////////////////////////
//
// GetCubicVolume

G4double G4Torus::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4AutoLock l(&torusMutex);
    fCubicVolume = fDPhi * CLHEP::pi * fRtor * (fRmax * fRmax - fRmin * fRmin);
    l.unlock();
  }
  return fCubicVolume;
}

////////////////////////////////////////////////////////////////////////////
//
// GetSurfaceArea

G4double G4Torus::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4AutoLock l(&torusMutex);
    fSurfaceArea = fDPhi * CLHEP::twopi * fRtor * (fRmax + fRmin);
    if (fDPhi < CLHEP::twopi)
    {
      fSurfaceArea = fSurfaceArea + CLHEP::twopi * (fRmax * fRmax - fRmin * fRmin);
    }
    l.unlock();
  }
  return fSurfaceArea;
}

///////////////////////////////////////////////////////////////////////
//
// Visualisation Functions

void G4Torus::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4Polyhedron* G4Torus::CreatePolyhedron() const
{
  return new G4PolyhedronTorus(fRmin, fRmax, fRtor, fSPhi, fDPhi);
}

#endif  // !defined(G4GEOM_USE_UTORUS)
