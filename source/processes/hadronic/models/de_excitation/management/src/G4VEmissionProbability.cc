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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

namespace {
  const G4double edeltamin = 0.1*CLHEP::MeV;
  const G4double edeltamax = 2*CLHEP::MeV;
  const G4double fact = 0.25;   // for selection of region of the peak
  const G4double fsmall = 0.02; // tolerance for edges of the spectrum
  const G4double fmaj = 1.05;   // tolerance for majoranta
}

G4VEmissionProbability::G4VEmissionProbability(G4int Z, G4int A)
  : pVerbose(1), theZ(Z), theA(A), elimit(CLHEP::MeV)
{
  pNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  pG4pow = G4Pow::GetInstance();
  if(A > 0) { pEvapMass = G4NucleiProperties::GetNuclearMass(theA, theZ); }
  G4DeexPrecoParameters* param = pNuclearLevelData->GetParameters();
  OPTxs = param->GetDeexModelType();
}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = pNuclearLevelData->GetParameters();
  pVerbose = param->GetVerbose();
  fFD = param->GetDiscreteExcitationFlag();
  pTolerance = param->GetMinExcitation();
  pWidth = param->GetNuclearLevelWidth();
}

void G4VEmissionProbability::ResetIntegrator(G4double de, G4double eps)
{
  if(de > 0.0)  { elimit = de; }
  if(eps > 0.0) { accuracy = eps; }
}

G4double G4VEmissionProbability::EmissionProbability(const G4Fragment&, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::ComputeProbability(G4double, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::IntegrateProbability(G4double elow, 
                                                      G4double ehigh, 
                                                      G4double cb)
{
  pProbability = 0.0;
  if(elow >= ehigh) { return pProbability; }

  emin = elow;
  emax = ehigh;
  eCoulomb = cb;

  G4int nbin = G4lrint((emax - emin)/elimit) + 1;
  nbin = std::max(nbin, 4);

  // providing smart binning
  G4double edelta = (emax - emin)/static_cast<G4double>(nbin);

  G4double x(emin), y(0.0);
  const G4double edelmicro = edelta*fsmall;
  probmax = ComputeProbability(x + edelmicro, eCoulomb);
  G4double problast = probmax;
  if(pVerbose > 1) {
    G4cout << "### G4VEmissionProbability::IntegrateProbability: "
	   << "probmax=" << probmax << " Emin=" << emin
	   << " Emax=" << emax << " QB=" << cb << " nbin=" << nbin 
	   << G4endl;
  }
  // fE1 - energy at half maximum after the peak
  fE1 = fE2 = emax;
  fP1 = fP2 = fP3 = 0.0;

  // edelmicro is used to avoid numerical problems
  fE3 = emax - edelmicro; 
  G4bool endpoint = false;
  for(G4int i=0; i<=nbin; ++i) {
    x += edelta;
    if(x >= fE3) { 
      edelta += emax - x;
      x = fE3;
      endpoint = true;
    }
    y = ComputeProbability(x, eCoulomb);
    if(pVerbose > 2) { 
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
    pProbability += del;
    // end of the loop
    if(del < accuracy*pProbability || endpoint) { break; }
    problast = y;

    // smart step definition
    if(del != pProbability && del > 0.8*pProbability && 
       0.7*edelta > edeltamin) { 
      edelta *= 0.7;
    } else if(del < 0.1*pProbability && 1.5*edelta < edeltamax) { 
      edelta *= 1.5;
    }
  }

  if(pVerbose > 1) { 
    G4cout << " Probability= " << pProbability << " probmax= " 
           << probmax << " emin=" << emin << " emax=" << emax 
	   << " E1=" << fE1 << " E2=" << fE2 << G4endl; 
  }
  return pProbability;
}

G4double G4VEmissionProbability::SampleEnergy()
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
      fP3 = 2*ComputeProbability(fE3, eCoulomb);
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

  if(pVerbose > 1) {
    G4cout << "### G4VEmissionProbability::SampleEnergy: " 
	   << " Emin= " << emin << " Emax= " << emax 
           << "/n    E1=" << fE1 << " p1=" << p1 
	   << " probmax=" << probmax << " P2=" << fP2 << G4endl;
  }

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
    gg = ComputeProbability(ekin, eCoulomb);
    if(pVerbose > 2) {
      G4cout << "    " << n
	     << ". prob= " << gg << " probmax= " << probmax
	     << " Ekin= " << ekin << G4endl;
    }
    if((gg > gmax || n > nmax) && pVerbose > 1) {
      G4cout << "### G4VEmissionProbability::SampleEnergy for Z= " << theZ 
             << " A= " << theA << " Eex(MeV)=" << fExc << " p1=" << p1
             << "\n    Warning n= " << n
	     << " prob/gmax=" << gg/gmax 
	     << " prob=" << gg << " gmax=" << gmax << " probmax=" << probmax 
	     << "\n    Ekin= " << ekin << " Emin= " << emin
	     << " Emax= " << emax << G4endl;
    }
  } while(gmax*rndm->flat() > gg && n < nmax);
  G4double enew = FindRecoilExcitation(ekin);
  if(pVerbose > 1) {
    G4cout << "### SampleEnergy: Efinal= " 
	   << enew << " E=" << ekin << "  Eexc=" << fExcRes << G4endl;
  }
  return enew;
}

G4double G4VEmissionProbability::FindRecoilExcitation(const G4double e)
{
  G4double mass = pEvapMass + fExc;
    
  G4double m02 = pMass*pMass;
  G4double m12 = mass*mass;
  G4double m22 = pResMass*pResMass;
  G4double mres = std::sqrt(m02 + m12 - 2.*pMass*(mass + e));

  fExcRes = mres - pResMass;

  if(pVerbose > 1) {
    G4cout << "### FindRecoilExcitation for resZ= " 
           << resZ << " resA= " << resA 
           << " evaporated Z= " << theZ << " A= " << theA
	   << " Ekin= " << e << " Eexc= " << fExcRes << G4endl;
  }

  // residual nucleus is in the ground state
  if(fExcRes < pTolerance) {
    fExcRes = 0.0;
    return std::max(0.5*(m02 + m12 - m22)/pMass - mass, 0.0);
  }
  if(!fFD) { return e; }
 
  // select final state excitation
  auto lManager = pNuclearLevelData->GetLevelManager(resZ, resA);
  if(nullptr == lManager) { return e; }

  // levels are not known
  if(fExcRes > lManager->MaxLevelEnergy() + pTolerance) { return e; }

  // find level
  std::size_t idx = lManager->NearestLevelIndex(fExcRes);
  auto level = lManager->GetLevel(idx); 

  // unstable level
  if (level->GetTimeGamma() == 0.0) { return e; }

  // is possible to use level energy?
  G4double elevel = lManager->LevelEnergy(idx);
  if (std::abs(elevel - fExcRes) > pWidth || pMass < mass + pResMass + elevel) { 
    return e;
  }

  // long-lived level
  G4double massR = pResMass + elevel;
  G4double mr2 = massR*massR;
  fExcRes = elevel;
  return std::max(0.5*(m02 + m12 - mr2)/pMass - mass, 0.0);
}
