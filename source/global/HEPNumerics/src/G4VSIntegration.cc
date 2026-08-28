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
// G4VSIntegration
//
// Created 03.03.2025 V.Ivanchenko
//
// --------------------------------------------------------------------

#include "G4VSIntegration.hh"
#include "Randomize.hh"

namespace
{
  const G4double numLimit = 1.e-8;
}

void G4VSIntegration::InitialiseIntegrator(G4double acc, G4double f1, G4double f2, G4double de,
                                           G4double dmin, G4double dmax)
{
  if (acc > numLimit && acc < 0.2)
  {
    fAcc = acc;
  }
  if (f1 > 0.01 && f1 < 1.0)
  {
    fFactor1 = f1;
  }
  if (f2 > 1.0 && f2 < 5.0)
  {
    fFactor2 = f2;
  }
  if (de > 0.0)
  {
    fDelta = de;
  }
  fMinDelta = (dmin <= fDelta && dmin > 0.0) ? dmin : fDelta;
  fMaxDelta = (dmax > fDelta) ? dmax : fDelta;
  if (fVerbose > 2)
  {
    G4cout << "### G4VSIntegration::InitialiseIntegrator: "
           << "fAcc=" << fAcc << " fFact1=" << fFactor1 << " fFact2=" << fFactor2
           << " dE=" << fDelta << " dEmin=" << fMinDelta << " dEmax=" << fMaxDelta << G4endl;
  }
}

G4double G4VSIntegration::ComputeIntegral(const G4double emin, const G4double emax)
{
  G4double res = 0.0;
  if (emin >= emax)
  {
    return res;
  }
  fEmin = emin;
  fEmax = emax;

  // preparing smart binning
  G4int nbin = G4lrint((emax - emin) / fDelta) + 1;
  nbin = std::max(nbin, 6);
  G4double edelta = (emax - emin) / static_cast<G4double>(nbin);

  // prepare integration
  G4double x(emin), y(0.0);
  fPmax = ProbabilityDensityFunction(x);
  G4double problast = fPmax;

  // for some distributions it may happens that there is a local maximum
  // closed to maximal energy
  fE3 = fEmax - 0.02 * (fEmax - fEmin);
  fP3 = ProbabilityDensityFunction(fE3);
#ifdef G4VERBOSE
  if (fVerbose > 1)
  {
    G4cout << "### G4VSIntegration::ComputeIntegral: "
           << "Pmax=" << fPmax << " Emin=" << emin << " Emax=" << emax << " dE=" << edelta
           << " nbin=" << nbin << G4endl;
  }
#endif

  fE1 = fE2 = emax;
  fP1 = fP2 = 0.0;
  G4bool endpoint = false;
  x += edelta;
  // integration is performed over energy interval
  // it may be stopped if at a step the addition 
  G4int nn = G4lrint((emax - emin) / fMinDelta) + 1;
  for (G4int i=0; i < nn; ++i)
  {
    // the last point 
    if (x >= emax)
    {
      edelta += emax - x;
      x = emax;
      endpoint = true;
    }
    y = ProbabilityDensityFunction(x);
#ifdef G4VERBOSE
    if (fVerbose > 2)
    {
      G4cout << "    " << i << ".  E=" << x << "  prob=" << y << " Edel=" << edelta << G4endl;
    }
#endif
    // search for function maximum
    if (0.0 == fP1) {
      if (y >= fPmax) {
        fPmax = y;
        // end of the 1st energy interval
      } else if (y > 0.0 && y < fFactor1 * fPmax && !endpoint) {
        fE1 = x;
        fP1 = y;
      }
    // search for function maximum in the 2nd energy interval
    } else if (0.0 == fP2) {
      if (y >= fP1) {
        fP1 = y;
        // define end of the 2nd energy interval
      } else if (y > 0.0 && y < fFactor1 * fP1 && !endpoint) {
        fE2 = x;
        fP2 = y;
      }
    // check maximum in the 3d area
    } else if (y > fP2) {
      fP2 = y;
    }

    G4double del = (y + problast) * edelta * 0.5;
    res += del;

    // end of the loop condition when energy is at the maximum value
    // or upper estimation of remaining integral is very small
    if ((std::max(y, fP3) * (fEmax - x) < fAcc * res) || endpoint)
    {
      break;
    }
    problast = y;

    // smart next step definition
    if (del < res)
    {
      // step close to the maximum
      if (del > 0.8 * res)
      {
        edelta = std::max(0.7 * edelta, fMinDelta);
      }
      // step at the tail
      else if (del < 0.1 * res)
      {
        edelta = std::min(1.5 * edelta, fMaxDelta);
      }
    }
    x += edelta;
  }
#ifdef G4VERBOSE
  if (fVerbose > 1)
  {
    G4cout << "### G4VSIntegration:" << ModelName() << " " 
           << "I=" << res << " E1=" << fE1 << " E2=" << fE2 << " Pmax=" << fPmax
           << " P1=" << fP1 << " P2=" << fP2 << G4endl;
  }
#endif
  return res;
}

G4double G4VSIntegration::SampleValue()
{
  // should never happen
  if (fEmin >= fEmax)
  {
    return fEmin;
  }

  // added margine for non precise definition of a maximum
  G4double Q0 = fPmax;
  G4double Q1 = fP1;
  G4double Q2 = fP2;

  // check maximum in each energy interval
  if (fE3 >= fE2)
  {
    G4double E3 = fEmax - 0.5 * (fEmax - fE2);
    G4double Q4 = ProbabilityDensityFunction(E3);
    Q2 = std::max(std::max(Q2, fP3), Q4);
  }
  else if (fE3 >= fE1)
  {
    Q1 = std::max(Q1, fP3);
  }
  else
  {
    Q0 = std::max(Q0, fP3);
  }
  Q0 *= fFactor2;
  Q1 *= fFactor2;
  Q2 *= fFactor2;

  // 3 energy intervals
  G4double del1 = fE1 - fEmin;
  G4double del2 = fE2 - fE1;
  G4double del3 = fEmax - fE2;

  // upper limit for area in all energy intervals
  G4double p1 = del1 * Q0;
  G4double p2 = del2 * Q1;
  G4double p3 = del3 * Q2;

  // if the 3d energy interval is very small it should not be considered
  if (p3 > 0.0 && p3 < numLimit * p2)
  {
    p3 = 0.0;
    if (fE3 >= fE1)
    {
      Q1 = std::max(Q1, Q2);
    }
    del2 = fEmax - fE1;
    p2 = del2 * Q1;
  }

  // sampling in 3 energy intervals, some probabilities may be zero
  G4double sum = p1 + p2 + p3;
  G4double p12 = p1 + p2;

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  const G4int nmax = 100000;
  G4double e, gmax, gg;

  // main sampling loop
  G4double rndmarray[3];
  for (G4int n = 0; n < nmax; ++n)
  {
    rndm->flatArray(3, rndmarray);
    G4double p = sum * rndmarray[0];
    if (p <= p1)
    {
      gmax = Q0;
      e = del1 * rndmarray[1] + fEmin;
    }
    else if (p <= p12)
    {
      gmax = Q1;
      e = del2 * rndmarray[1] + fE1;
    }
    else
    {
      gmax = Q2;
      e = del3 * rndmarray[1] + fE2;
    }
    gg = ProbabilityDensityFunction(e);
    if ((gg > gmax || n >= nmax) && fVerbose > 0)
    {
      ++fnWarn;
      if (fnWarn <= fWarnLimit)
      {
        G4cout << "### G4VSIntegration::SampleValue() for " << ModelName() 
               << " n=" << n << " prob=" << gg << " gmax=" << gmax << G4endl;
        G4cout << "    E=" << e << " Emin=" << fEmin << " Emax=" << fEmax << " E1=" << fE1
               << " E2=" << fE2 << " F0=" << Q0 << " F1=" << Q1 << " F2=" << Q2
               << G4endl;
      }
      if (fnWarn == fWarnLimit)
      {
        G4cout << "### G4VSIntegration warnings are closed" << G4endl;
      }
    }
    if (gmax * rndmarray[2] <= gg)
    {
#ifdef G4VERBOSE
      if (fVerbose > 1)
      {
        G4cout << "### G4VSIntegration::SampleValue for " << ModelName() << " E=" << e
               << " Ntry=" << n << " Emin=" << fEmin << " Emax=" << fEmax << G4endl;
      }
#endif
      return e;
    }
  }
  // if sampling not converged, then sample uniformly in the 1st energy interval
  e = fEmin + rndm->flat() * del1;
#ifdef G4VERBOSE
  if (fVerbose > 1)
  {
    G4cout << "### G4VSIntegration::SampleValue for " << ModelName() << " E=" << e
           << " Ntry=" << nmax << " Emin=" << fEmin << " Emax=" << fEmax << G4endl;
  }
#endif
  return e;
}

const G4String& G4VSIntegration::ModelName() const
{
  return dummy;
}
