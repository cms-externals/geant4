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
// by V. Lara

#include "G4PreCompoundIon.hh"
#include "G4PreCompoundParameters.hh"
//#include "G4EvaporationLevelDensityParameter.hh"


G4double G4PreCompoundIon::
ProbabilityDistributionFunction(const G4double eKin, 
				const G4Fragment& aFragment)
{
  if ( !IsItPossible(aFragment) ) return 0.0;
  
  const G4double r0 = G4PreCompoundParameters::GetAddress()->Getr0();

  G4double U = aFragment.GetExcitationEnergy();
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  G4double N = P + H;
  
  // g = (6.0/pi2)*a*A
  //  G4EvaporationLevelDensityParameter theLDP;
  G4double g0 = (6.0/pi2)*aFragment.GetA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(static_cast<G4int>(aFragment.GetA()),static_cast<G4int>(aFragment.GetZ()),U);
  G4double g1 = (6.0/pi2)*GetRestA() * 
    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(GetRestA(),GetRestZ(),U);
  //no longer used:  G4double gj = (6.0/pi2)*GetA() *
  //  -----"-----    G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  //    theLDP.LevelDensityParameter(GetA(),GetZ(),U);

  G4double A0 = ((P*P+H*H+P-H)/4.0 - H/2.0)/g0;
  G4double A1 = std::max(0.0,(A0*g0 + GetA()*(GetA()-2.0*P-1.0)/4.0)/g1);
  G4double Aj = GetA()*(GetA()+1.0)/4.0;

  G4double E0 = std::max(0.0,U - A0);
  if (E0 == 0.0) return 0.0;
  //  G4double E1 = std::max(0.0,GetMaximalKineticEnergy() + GetCoulombBarrier() - eKin - A1);
  G4double E1 = (std::max(0.0,GetMaximalKineticEnergy() - eKin - A1))/g1; //JMQ fix
  //  G4double Ej = std::max(0.0,U + GetBindingEnergy() - Aj);
  G4double Ej = std::max(0.0,eKin + GetBindingEnergy() - Aj); //JMQ fix

  G4double pA = (3.0/4.0) * std::sqrt(std::max(0.0, 2.0/(GetReducedMass()*(eKin+GetBindingEnergy()))))*
      GetAlpha()*(eKin+GetBeta())/(r0*std::pow(GetRestA(),1.0/3.0)) *
      CoalescenceFactor(aFragment.GetA()) * FactorialFactor(N,P);

  G4double pB = std::pow((g1*E1)/(g0*E0),N-GetA()-1.0)*(g1/g0);

  //  G4double pC = std::pow((gj*Ej)/(g0*E0),GetA()-1.0)*(gj/g0)/E0;
  G4double pC = std::pow((g1*Ej)/(g0*E0),GetA()-1.0)*(g1/g0)/E0; //JMQ fix

  G4double Probability = pA * pB * pC * GetRj(aFragment.GetNumberOfParticles(), aFragment.GetNumberOfCharged());

  return Probability;
}