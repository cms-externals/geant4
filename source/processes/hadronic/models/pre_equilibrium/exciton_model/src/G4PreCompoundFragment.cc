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
//
// J. M. Quesada (August 2008).  
// Based  on previous work by V. Lara
//
// Modified:
// 06.09.2008 JMQ Also external choice has been added for:
//               - superimposed Coulomb barrier (if useSICB=true) 
// 20.08.2010 V.Ivanchenko cleanup
//

#include "G4PreCompoundFragment.hh"
#include "G4KalbachCrossSection.hh"
#include "G4ChatterjeeCrossSection.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4InterfaceToXS.hh"
#include "G4IsotopeList.hh"
#include "Randomize.hh"

namespace {
  const G4double eden = 1.0/CLHEP::MeV;
  const G4double fact = 0.25;   // for selection of region of the peak
  const G4double fsmall = 0.02; // tolerance for edges of the spectrum
  const G4double fmaj = 1.05;   // tolerance for majoranta
}

G4PreCompoundFragment::G4PreCompoundFragment(const G4ParticleDefinition* p,
					     G4VCoulombBarrier* aCoulBarrier)
  : G4VPreCompoundFragment(p, aCoulBarrier)
{}

G4double G4PreCompoundFragment::CalcEmissionProbability(const G4Fragment& fr)
{
  theEmissionProbability = (Initialize(fr)) ?
    IntegrateEmissionProbability(theMinKinEnergy, theMaxKinEnergy, fr) : 0.0;
  /*  
  G4cout << "## G4PreCompoundFragment::CalcEmisProb "
         << "Zf= " << fr.GetZ_asInt()
	 << " Af= " << fr.GetA_asInt()
	 << " Elow= " << theMinKinEnergy
	 << " Eup= " << theMaxKinEnergy
	 << " prob= " << theEmissionProbability
	 << " index=" << index << " Z=" << theZ << " A=" << theA
	 << G4endl;
  */
  return theEmissionProbability;
}

G4double 
G4PreCompoundFragment::IntegrateEmissionProbability(G4double low, G4double up,
                                                    const G4Fragment& fr)
{
  G4double res = 0.0;
  if (low >= up) { return res; }
  emin = low;
  emax = up; 
  G4double edelta = (up - low);
  G4int nbin = G4lrint(edelta*eden) + 1;
  nbin = std::max(nbin, 4);
  edelta /= static_cast<G4double>(nbin);

  G4double x(emin), y(0.0);
  const G4double edelmicro = edelta*fsmall;
  probmax = ProbabilityDistributionFunction(x + edelmicro, fr); 
  G4double problast = probmax;
  if (verbose > 1) {
    G4cout << "### G4PreCompoundFragment::IntegrateEmissionProbability: "
	   << "probmax=" << probmax << " Emin=" << emin
	   << " Emax=" << emax << " nbin=" << nbin 
	   << G4endl;
  }
  // fE1 - energy at half maximum after the peak
  fE1 = fE2 = emax;
  fP1 = fP2 = fP3 = 0.0;

  // edelmicro is used to avoid numerical problems
  fE3 = emax - edelmicro; 
  G4bool endpoint = false;
  for (G4int i=0; i<=nbin; ++i) {
    x += edelta;
    if(x >= fE3) { 
      edelta += emax - x;
      x = fE3;
      endpoint = true;
    }
    y = ProbabilityDistributionFunction(x, fr);
    if (verbose > 2) { 
      G4cout << "    " << i << ".  E= " << x << "  prob= " << y
	     << " Edel= " << edelta << G4endl;
    } 
    if (y >= probmax) {
      probmax = y;
    } else if (!endpoint) {
      if (0.0 == fP1) {
	if (y < fact*probmax) {
	  fE1 = x;
	  fP1 = y;
	}
      } else if (0.0 == fP2 && y < fact*fP1) {
	fE2 = x;
	fP2 = y;
      }
    }
    
    G4double del = (y + problast)*edelta*0.5;
    res += del;
    // end of the loop
    if (del < accuracy*res || endpoint) { break; }
    problast = y;
  }

  /*
  if(pVerbose > 1) { 
    G4cout << " Probability= " << pProbability << " probmax= " 
           << probmax << " emin=" << emin << " emax=" << emax 
	   << " E1=" << fE1 << " E2=" << fE2 << G4endl; 
  }
  */
  return res;
}

G4double G4PreCompoundFragment::CrossSection(G4double ekin)
{
  /*
  G4cout << "G4PreCompoundFragment::CrossSection OPTxs=" << OPTxs << " E=" << ekin
	 << " resZ=" << theResZ << " resA=" << theResA << " index=" << index
	 << " fXSection:" << fXSection << G4endl;
  */
  // compute power once
  if (OPTxs > 1 && 0 < index && theResA != lastA) {
    lastA = theResA;
    muu = G4KalbachCrossSection::ComputePowerParameter(lastA, index);
  }
  if (OPTxs == 0) { 
    recentXS = GetOpt0(ekin);
  } else if (OPTxs == 1) {
    G4int Z = std::min(theResZ, ZMAXNUCLEARDATA);
    //G4double e = std::max(ekin, lowEnergyLimitMeV[Z]);
    recentXS = fXSection->GetElementCrossSection(ekin, Z)/CLHEP::millibarn;

  } else if (OPTxs == 2) { 
    recentXS = G4ChatterjeeCrossSection::ComputeCrossSection(ekin, 
                                                             theCoulombBarrier, 
							     theResA13, muu, 
							     index, theZ, theResA); 

  } else { 
    recentXS = G4KalbachCrossSection::ComputeCrossSection(ekin, theCoulombBarrier, 
						          theResA13, muu, index,
						          theZ, theA, theResA);
  }
  return recentXS;
}  

G4double G4PreCompoundFragment::GetOpt0(G4double ekin) const
// OPT=0 : Dostrovski's cross section
{
  G4double r0 = theParameters->GetR0()*theResA13;
  // cross section is now given in mb (r0 is in mm) for the sake of consistency
  // with the rest of the options
  return 1.e+25*CLHEP::pi*r0*r0*theResA13*GetAlpha()*(1.0 + GetBeta()/ekin);
}

G4double G4PreCompoundFragment::SampleKineticEnergy(const G4Fragment& fr) 
{
  probmax *= fmaj;

  // two regions with flat and one with exponential majorant 
  G4double b = 0.0;
  G4double p3 = 0.0;

  if (fP1 > 0.0) {
    // exclude 2d and 3d areas from sampling
    if (fP1 > probmax*fact) {
      fP2 = fP1 = 0.0;
      fE1 = fE2 = emax;
      // 2d area is considered
    } else if (fP2 > 0.0 && fE2 < emax) {
      fP3 = 2*ProbabilityDistributionFunction(fE3, fr);
      // exclude 3d area from sampling
      if (fP3 > fP2*fact || fP3 == 0.0) {
	fP3 = 0.0;
	fE2 = emax;
      } else {
	b = G4Log(fP2/fP3)/(emax - fE2);
	p3 = (fP2 - fP3)/b;
      }
    }
  }
  G4double p1 = (fE1 - emin)*probmax;
  G4double p2 = (fE2 - fE1)*fP1;
  G4double sum = p1 + p2 + p3;
  G4double del1 = (fE1 - emin)*sum/p1;
  G4double del2 = (p2 > 0.0) ? (fE2 - fE1)/p2 : 0.0;

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  const G4int nmax = 1000;
  G4double ekin, gmax, gg;
  G4int n = 0;
  do {
    ++n;
    G4double q = rndm->flat();
    G4double p = sum*q;
    if (p <= p1) {
      gmax = probmax;
      ekin = del1*q + emin;
    } else if (p <= p1 + p2) {
      gmax = fP1;
      ekin = del2*(p - p1) + fE1;
    } else {
      G4double x = 1.0 - rndm->flat()*(1.0 - fP3/fP2);
      ekin = fE2 - G4Log(x)/b;
      gmax = fP2*x;
    }
    gg = ProbabilityDistributionFunction(ekin, fr);
    if ((gg > gmax || n > nmax) && verbose > 1) {
      G4cout << "### G4PreCompoundFragment::SampleKineticEnergy for Z= " << theZ 
             << " A= " << theA << " p1=" << p1 << " p2="
	     << p2 << " p3=" << p3
             << "\n    Warning n= " << n
	     << " prob/gmax=" << gg/gmax 
	     << " prob=" << gg << " gmax=" << gmax << " probmax=" << probmax 
	     << "\n    Ekin= " << ekin << " Emin= " << emin
	     << " Emax= " << emax << G4endl;
    }
  } while(gmax*rndm->flat() > gg && n < nmax);
  return ekin;
}

