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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name: G4ChargeExchangeXS
//
 
// Class Description:
// This is a class for charge exchange hadronic cross section. 
// Class Description - End

#ifndef G4ChargeExchangeXS_h
#define G4ChargeExchangeXS_h 

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include <vector>

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Isotope;
class G4Element;
class G4Material;
class G4Pow;

class G4ChargeExchangeXS final : public G4VCrossSectionDataSet
{
public: 

  G4ChargeExchangeXS();

  ~G4ChargeExchangeXS() override = default;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element*, const G4Material*) override;

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso,
                              const G4Element* elm,
                              const G4Material* mat) override;

  void CrossSectionDescription(std::ostream&) const override;

  const G4ParticleDefinition* SampleSecondaryType(const G4ParticleDefinition* part, 
				       G4int Z, G4int A, G4double ekin);

  G4Pow* g4calc;

  G4ChargeExchangeXS & operator=(const G4ChargeExchangeXS &right) = delete;
  G4ChargeExchangeXS(const G4ChargeExchangeXS&) = delete;
  
  // beta_prime value for calculation of cross section of change exchange interactions 
  // with different nucleon 
  const G4double beta_prime_pi = 0.0410;
  const G4double beta_prime_eta = 0.0402;

};

#endif
