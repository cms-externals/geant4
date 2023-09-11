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
//
// GEANT4 Class file
//
//
// File name:    G4ChargeExchangeXS
//
//

#include "G4ChargeExchangeXS.hh"
#include "G4DynamicParticle.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4HadronicParameters.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"  
#include "G4Pow.hh"

#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

G4ChargeExchangeXS::G4ChargeExchangeXS() 
{
  if(verboseLevel > 0){
    G4cout  << "G4ChargeExchangeXS::G4ChargeExchangeXS" << G4endl;
  }
  SetForAllAtomsAndEnergies(true);

  g4calc = G4Pow::GetInstance();
}

// Print the information of this .cc file
void G4ChargeExchangeXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ChargeExchangeXS calculates charge exchange cross section for "
          << "pi+, pi-, K+, K-, KL\n";
}

G4bool G4ChargeExchangeXS::IsIsoApplicable(const G4DynamicParticle*,
                                           G4int, G4int,
                                           const G4Element*, const G4Material*)
{
  return true;
}

G4double
G4ChargeExchangeXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
				       G4int Z, G4int A,
				       const G4Isotope*, const G4Element*,
				       const G4Material*)  
{
  G4double result = 0.0;
  auto part = aParticle->GetDefinition();
  G4int pdg = part->GetPDGEncoding();   

 // Get or calculate the nucleus mass, particle mass,particle kinetic energy 
 // and particle total energy 
  G4double tM = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double pE = aParticle->GetTotalEnergy();
  G4double pM = part->GetPDGMass(); 

  // Calculate the momentum of the bombarding particles and convert it to GeV/c^2 unit
  G4double p_momentum = std::sqrt(pE*pE - pM*pM)/(CLHEP::GeV);

  // Calculate s(lorentz invariant)
  G4double lorentz_s = tM*tM + 2*tM*pE +  pM*pM;

  // For unit conversion 
  const G4double inv1e7 = 1e-7;
  const G4double inv1e30 = 1e-30;

  // Calculate the charge exchange process cross section of different production of particles
  G4double Xsc_pi0 = (122)*g4calc->powA((lorentz_s*inv1e7), -1.23)*inv1e30*CLHEP::cm2;
  G4double Xsc_eta = (31)*g4calc->powA((lorentz_s*inv1e7), -1.53)*inv1e30*CLHEP::cm2;
  G4double Xsc_Kbar = (56.3)*g4calc->powA((p_momentum/10), -1.60)*inv1e30*CLHEP::cm2;

  // pi- + p -> sum of (pi0 + eta) + n 
  if (pdg == -211) {
    // The approximation of Glauber-Gribov formula -> extend it from interation with proton to different nucleon 
    // The cross section prop to Z^(2/3) 
    // The g4calc->powA(A,-beta_prime_pi*G4Log(A)) term is added for absorption of secondat particles 
    G4double Xsc_pi0_Z = Xsc_pi0*g4calc->Z23(Z)*g4calc->powA(A,-beta_prime_pi*G4Log(A));
    G4double Xsc_eta_Z = Xsc_eta*g4calc->Z23(Z)*g4calc->powA(A,-beta_prime_eta*G4Log(A));
    result = (Xsc_pi0_Z + Xsc_eta_Z)/(CLHEP::millibarn);
  }

  // pi+ + n -> sum of (pi0 + eta) + p
  else if (pdg == 211) {
    G4double Xsc_pi0_Z = Xsc_pi0*g4calc->Z23(A-Z)*g4calc->powA(A,-beta_prime_pi*G4Log(A));
    G4double Xsc_eta_Z = Xsc_eta*g4calc->Z23(A-Z)*g4calc->powA(A,-beta_prime_eta*G4Log(A));
    result = (Xsc_pi0_Z + Xsc_eta_Z)/(CLHEP::millibarn);
    
  }

  // K- + p -> Kbar + n
  else if (pdg == -321){
    result = g4calc->Z23(Z)*(Xsc_Kbar)/(CLHEP::millibarn);
  }

  // K+ + n -> Kbar + p
  else if (pdg == 321) {
    result = g4calc->Z23(A-Z)*(Xsc_Kbar)/(CLHEP::millibarn);
  }

  // KL 
  else if (pdg == 130) {
    // Cross section of K-long = 0.5*(Cross section of K+ + Cross section of K-)
    result = 0.5*(g4calc->Z23(Z) + g4calc->Z23(A-Z))*(Xsc_Kbar)/(CLHEP::millibarn);
  }
  
  return result;
}

const G4ParticleDefinition*
G4ChargeExchangeXS::SampleSecondaryType(const G4ParticleDefinition* part, 
				       G4int Z, G4int A, G4double ekin)
{

  const G4ParticleDefinition* pd = nullptr;
  G4int pdg = part->GetPDGEncoding();  

  // Get or calculate the nucleus mass, particle mass and the particle total energy 
  G4double tM = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double pM = part->GetPDGMass(); 
  G4double pE = pM + ekin; 

  // Calculate s(lorentz invariant)
  G4double lorentz_s = tM*tM + 2*tM*pE +  pM*pM;

  // For unit conversion
  const G4double inv1e7 = 1e-7;
  const G4double inv1e30 = 1e-30;
  
  // Calculate the charge exchange process cross section of different production of particles
  G4double Xsc_pi0 = (122)*g4calc->powA((lorentz_s*inv1e7), -1.23)*inv1e30*CLHEP::cm2;
  G4double Xsc_eta = (31)*g4calc->powA((lorentz_s*inv1e7), -1.53)*inv1e30*CLHEP::cm2;
 
  // pi- + p /  pi+ + n  
  // For pi+ case, it is the isotopic invariant. 
  if (std::abs(pdg) == 211){

    // The approximation of Glauber-Gribov formula -> extend it from interation with proton to different nucleon
    G4double Xsc_pi0_Z = Xsc_pi0*g4calc->Z23(Z)*g4calc->powA(A,-beta_prime_pi*G4Log(A));
    G4double Xsc_eta_Z = Xsc_eta*g4calc->Z23(Z)*g4calc->powA(A,-beta_prime_eta*G4Log(A));
    G4double Xsc_pi0_pb = Xsc_pi0_Z/(Xsc_pi0_Z+Xsc_eta_Z);

    if (G4UniformRand()>Xsc_pi0_pb){
      pd = G4Eta::Eta(); 
    }
    else{
      pd = G4PionZero::PionZero();
    }
  }

  // K- + p /  K+ + n 
  // Equal opportunity of producing k-short and k-long
  else if (std::abs(pdg) == 321){
    if (G4UniformRand()>0.5){
      pd = G4KaonZeroLong::KaonZeroLong(); 
    }
    else{
      pd = G4KaonZeroShort::KaonZeroShort();
    }
  }

  // KL + atom 
  else if (std::abs(pdg) == 130){
    G4double prob = G4double(Z)/G4double(A);
    if (G4UniformRand()>prob){
      pd = G4KaonMinus::KaonMinus();
    }
    else{
      pd = G4KaonPlus::KaonPlus();
    }
  }

  return pd;
} 
