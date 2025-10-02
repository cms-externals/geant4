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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modification: 13.08.2025 V.Ivanchenko rewrite

#include <numeric>
#include "G4StatMFMicroCanonical.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4Pow.hh"

namespace
{
  constexpr G4int fMaxMultiplicity = 4;
}
  
// constructor
G4StatMFMicroCanonical::G4StatMFMicroCanonical() 
{
  fPartitionManagerVector.reserve(fMaxMultiplicity);
}

// destructor
G4StatMFMicroCanonical::~G4StatMFMicroCanonical() 
{
  if (!fPartitionManagerVector.empty()) {
    for (auto const & p : fPartitionManagerVector) { delete p; }
  }
}

void G4StatMFMicroCanonical::Initialise(const G4Fragment& theFragment) 
{
  fPartitionManagerVector.clear();
  // Excitation Energy 
  G4double U = theFragment.GetExcitationEnergy();

  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  G4double x = 1.0 - (2*Z)/G4double(A);
  G4Pow* g4calc = G4Pow::GetInstance();
    
  // Configuration temperature
  G4double TConfiguration = std::sqrt(8.0*U/G4double(A));
  
  // Free internal energy at Temperature T = 0
  pFreeInternalE0 = A*( 
			// Volume term (for T = 0)
			-G4StatMFParameters::GetE0() +  
			// Symmetry term
			G4StatMFParameters::GetGamma0()*x*x 
			) + 
    // Surface term (for T = 0)
    G4StatMFParameters::GetBeta0()*g4calc->Z23(A) + 
    // Coulomb term 
    CLHEP::elm_coupling*0.6*(Z*Z)/(G4StatMFParameters::Getr0()*g4calc->Z13(A));

  //G4cout << "Tconf=" <<  TConfiguration << " freeE=" << pFreeInternalE0 << G4endl;
    
  // Mean breakup multiplicity
  pMeanMultiplicity = 0.0;
  
  // Mean channel temperature
  pMeanTemperature = 0.0;
  
  // Mean channel entropy
  pMeanEntropy = 0.0;
  
  // Calculate entropy of compound nucleus
  G4double SCompoundNucleus = CalcEntropyOfCompoundNucleus(theFragment,TConfiguration);
  
  // Statistical weight of compound nucleus
  fWCompoundNucleus = 1.0; 
  
  // Statistical weight
  G4double W = fWCompoundNucleus;
    
  // Maximal fragment multiplicity allowed in direct simulation  
  for (G4int im = 2; im <= fMaxMultiplicity; ++im) {
    auto ptr = new G4StatMFMicroManager(theFragment,im,pFreeInternalE0,SCompoundNucleus);
    fPartitionManagerVector.push_back(ptr);
    W += ptr->GetProbability();
  }

  // Normalization of statistical weights
  for (auto & ptr : fPartitionManagerVector) {
    ptr->Normalize(W);
    pMeanMultiplicity += ptr->GetMeanMultiplicity();
    pMeanTemperature += ptr->GetMeanTemperature();
    pMeanEntropy += ptr->GetMeanEntropy();
  }

  fWCompoundNucleus /= W;
  
  pMeanMultiplicity += fWCompoundNucleus;
  pMeanTemperature += TConfiguration * fWCompoundNucleus;
  pMeanEntropy += SCompoundNucleus * fWCompoundNucleus;
}

G4double G4StatMFMicroCanonical::CalcFreeInternalEnergy(const G4Fragment& theFragment,
							G4double T)
{
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  G4double A13 = G4Pow::GetInstance()->Z13(A);
  
  G4double InvLevelDensityPar = G4StatMFParameters::GetEpsilon0()
    *(1.0 + 3.0/G4double(A-1));
  
  G4double VolumeTerm = (-G4StatMFParameters::GetE0()+T*T/InvLevelDensityPar)*A;
  
  G4double SymmetryTerm = G4StatMFParameters::GetGamma0()
    *(A - 2*Z)*(A - 2*Z)/G4double(A);
  
  G4double SurfaceTerm = (G4StatMFParameters::Beta(T)
			  - T*G4StatMFParameters::DBetaDT(T))*A13*A13;
  
  G4double CoulombTerm = elm_coupling*0.6*Z*Z/(G4StatMFParameters::Getr0()*A13);

  G4double sum = VolumeTerm + SymmetryTerm + SurfaceTerm + CoulombTerm;
  //G4cout << "G4StatMFMicroCanonical::CalcFreeInternalEnergy " << sum
  //	 << " " << VolumeTerm << " " << SymmetryTerm << " " << SurfaceTerm
  //	 <<  " " << CoulombTerm << G4endl;
  return sum;
}

G4double 
G4StatMFMicroCanonical::CalcEntropyOfCompoundNucleus(const G4Fragment & theFragment,
						     G4double & TConf)
  // Calculates Temperature and Entropy of compound nucleus
{
  G4int A = theFragment.GetA_asInt();
  G4double U = theFragment.GetExcitationEnergy();
  G4double A13 = G4Pow::GetInstance()->Z13(A);
  
  G4double Ta = std::max(std::sqrt(U/(0.125*A)),0.0012*MeV); 
  G4double Tb = Ta;
  
  G4double ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Ta);
  G4double Da = (U+pFreeInternalE0-ECompoundNucleus)/U;
  G4double Db = 0.0;
    
  G4double InvLevelDensity = CalcInvLevelDensity(A);
  
  // bracketing the solution
  if (Da == 0.0) {
    TConf = Ta;
    return 2*Ta*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Ta)*A13*A13;
  } else if (Da < 0.0) {
    do {
      Tb -= 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+pFreeInternalE0-ECompoundNucleus)/U;
    } while (Db < 0.0);
  } else {
    do {
      Tb += 0.5*Tb;
      ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tb);
      Db = (U+pFreeInternalE0-ECompoundNucleus)/U;
    } while (Db > 0.0);
  }
  
  G4double eps = 1.0e-14 * std::abs(Tb-Ta);
  
  for (G4int i = 0; i < 1000; i++) {
    G4double Tc = (Ta+Tb)*0.5;
    if (std::abs(Ta-Tb) <= eps) {
      TConf = Tc;
      return 2*Tc*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Tc)*A13*A13;
    }
    ECompoundNucleus = CalcFreeInternalEnergy(theFragment,Tc);
    G4double Dc = (U+pFreeInternalE0-ECompoundNucleus)/U;
    
    if (Dc == 0.0) {
      TConf = Tc;
      return 2*Tc*A/InvLevelDensity - G4StatMFParameters::DBetaDT(Tc)*A13*A13;
    }
    
    if (Da*Dc < 0.0) {
      Tb = Tc;
      Db = Dc;
    } else {
      Ta = Tc;
      Da = Dc;
    } 
  }
  /*
  G4cout << "G4StatMFMicroCanoncal::CalcEntropyOfCompoundNucleus: "
	 << " faild calculate the temperature" << G4endl;
  */
  
  return 0.0;
}

G4StatMFChannel*  G4StatMFMicroCanonical::ChooseAandZ(const G4Fragment& theFragment)
{
  // Choice of fragment atomic numbers and charges 
  // We choose a multiplicity (1,2,3,...) and then a channel
  G4double RandNumber = G4UniformRand();
  G4int A = theFragment.GetA_asInt();
  G4int Z = theFragment.GetZ_asInt();
  
  if (RandNumber < fWCompoundNucleus) { 
	
    G4StatMFChannel * aChannel = new G4StatMFChannel;
    aChannel->CreateFragment(A, Z);
    return aChannel;
	
  } else {
    G4double AccumWeight = fWCompoundNucleus;
    for (auto & ptr : fPartitionManagerVector) {
      AccumWeight += ptr->GetProbability();
      if (RandNumber <= AccumWeight) {
	return ptr->ChooseChannel(A, Z, pMeanTemperature);
      }
    }
  }

  return nullptr;
}

G4double G4StatMFMicroCanonical::CalcInvLevelDensity(G4int anA)
{
  G4double res = 0.0;
  if (anA > 1) {
    res = G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(anA - 1.0));
  }
  return res;
}
