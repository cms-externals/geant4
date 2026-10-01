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
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4UnstableFragmentBreakUp
//
//      Author:        V.Ivanchenko
//
//      Creation date: 11 May 2010
//
// -------------------------------------------------------------------
//

#include "G4UnstableFragmentBreakUp.hh"

#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4LevelManager.hh"
#include "G4LorentzVector.hh"
#include "G4NuclearLevelData.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"

const G4int G4UnstableFragmentBreakUp::Zfr[] = {0, 1, 1, 1, 2, 2};
const G4int G4UnstableFragmentBreakUp::Afr[] = {1, 1, 2, 3, 3, 4};
G4double G4UnstableFragmentBreakUp::masses[] = {0., 0., 0., 0., 0., 0.};

namespace
{
constexpr G4double dmlimit = 5 * CLHEP::keV;
}

G4UnstableFragmentBreakUp::G4UnstableFragmentBreakUp()
{
  fLevelData = G4NuclearLevelData::GetInstance();
  if (0.0 == masses[0])
  {
    for (G4int i = 0; i < 6; ++i)
    {
      masses[i] = G4NucleiProperties::GetNuclearMass(Afr[i], Zfr[i]);
    }
  }
  for (G4int i = 0; i < 6; ++i)
  {
    prob[i] = mrec[i] = 0.0;
  }
  fSecID = G4PhysicsModelCatalog::GetModelID("model_G4UnstableFragmentBreakUp");
}

G4bool G4UnstableFragmentBreakUp::BreakUpChain(G4FragmentVector* results, G4Fragment* nucleus)
{
  // G4cout << "G4UnstableFragmentBreakUp::EmittedFragment" << G4endl;
  G4Fragment* frag = nullptr;

  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();

  G4LorentzVector lv = nucleus->GetMomentum();
  G4double time = nucleus->GetCreationTime();

  // look for the decay channel with normal masses
  // without Coulomb barrier and paring corrections
  if (fVerbose > 1)
  {
    G4cout << "#Unstable decay " << " Z= " << Z << " A= " << A
           << " Eex(MeV)= " << nucleus->GetExcitationEnergy() << G4endl;
  }
  G4double mass = lv.mag();
  G4double sum = 0.0;
  G4double dmass = 0.0;

  // only small energy non-concervation is allowed
  G4double maxprob = -0.5 * CLHEP::MeV;
  G4int idx = -1;
  for (G4int i = 0; i < 6; ++i)
  {
    G4int Zres = Z - Zfr[i];
    G4int Ares = A - Afr[i];
    G4double m1 = 0.0;
    G4double p = 0.0;
    if (Zres >= 0 && Ares >= Zres && Ares >= Afr[i])
    {
      m1 = G4NucleiProperties::GetNuclearMass(Ares, Zres);
      // decay probability is proportional to free energy
      p = mass - masses[i] - m1;
      if (p > maxprob)
      {
        maxprob = p;
        idx = i;
      }
      p = std::max(p, 0.0);
    }
    sum += p;
    prob[i] = sum;
    mrec[i] = m1;
  }
  // there is no decay channel
  if (-1 == idx)
  {
    return false;
  }

  // decay channel may be considered
  dmass = maxprob;
  if (sum > 0.0)
  {
    sum *= G4UniformRand();
    for (G4int i = 0; i < 6; ++i)
    {
      if (sum <= prob[i] && mrec[i] > 0.0)
      {
        idx = i;
        dmass = mass - masses[i] - mrec[i];
        break;
      }
    }
  }
  G4int Zres = Z - Zfr[idx];
  G4int Ares = A - Afr[idx];
  if (Zres < 0 || Ares < Zres)
  {
    return false;
  }

  if (dmass > 0.0)
  {
    // check if final particles are "stable" light mesons
    for (G4int i = 0; i < 6; ++i)
    {
      if (Zres == Zfr[i] && Ares == Afr[i])
      {
        dmass = 0.0;
        break;
      }
    }
    // randomize excitation energy
    if (dmass > 0.0)
    {
      G4double x;
      do
      {
        x = G4RandGauss::shoot(0.0, 0.25);
      } while (std::abs(x) > 1.0);
      dmass *= 0.5 * (1.0 + x);
    }
  }
  else
  {
    dmass = 0.0;
  }
  // 1 - light ion
  G4double mass1 = masses[idx];
  G4double mass2 = mrec[idx] + dmass;
  /*
  G4cout << "idx=" << idx << " Z=" << Z << " A=" << A << " Zres=" << Zres
         << " Ares=" << Ares << " dmass=" << dmass
   << " delm=" << mass - mass1 - mass2 << G4endl;
  G4cout << "    mass=" << mass << " mass1=" << mass1 << " mass2=" << mass2
         << G4endl;
  */
  // correct primary 4-vector if the reaction is impossible kinematically
  G4double m0 = mass1 + mass2;
  if (mass < m0)
  {
    G4double e0 = lv.e() - mass + m0 + 2 * dmlimit;
    mass = m0 + dmlimit;
    G4ThreeVector v0 = lv.vect().unit();
    G4double p0 = std::sqrt((e0 - mass) * (e0 + mass));
    lv.set(p0 * v0, e0);
  }

  // compute 4-momentum of light fragment
  G4double e1 = 0.5 * ((mass - mass2) * (mass + mass2) + mass1 * mass1) / mass;
  e1 = std::max(e1, mass1);
  G4double mom1 = std::sqrt((e1 - mass1) * (e1 + mass1));
  G4ThreeVector bst = lv.boostVector();
  G4ThreeVector v1 = G4RandomDirection();
  G4LorentzVector lv1 = G4LorentzVector(v1 * mom1, e1);
  lv1.boost(bst);
  frag = new G4Fragment(Afr[idx], Zfr[idx], lv1);
  frag->SetCreationTime(time);
  frag->SetCreatorModelID(fSecID);
  results->push_back(frag);

  // residual
  lv -= lv1;
  nucleus->SetZAandMomentum(lv, Zres, Ares);
  nucleus->SetCreatorModelID(fSecID);

  return true;
}

G4double G4UnstableFragmentBreakUp::GetEmissionProbability(G4Fragment*)
{
  return 0.0;
}
