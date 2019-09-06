//---------------------------------------------------------------------------
//
// ClassName:  G4HadronicAbsorptionFritiofWithBinaryCascade
//
// Author:     Alberto Ribon
//
// Date:       July 2019
//
// Modified:
//
// Class Description:
//
// Intermediate base class for hadronic absorption at rest using
// Fritiof (FTF) coupled with Binary Cascade (BIC) for the nuclear
// de-excitation. 
//
// Note: it is applicable only for anti_proton and anti_neutron,
//       but not for light anti-ion because the propagate method
//       for nucleus-nucleus interaction:
//         G4VIntraNuclearTransportModel::Propagate()
//       is not implemented.
//
//---------------------------------------------------------------------------

#ifndef G4HadronicAbsorptionFritiofWithBinaryCascade_h
#define G4HadronicAbsorptionFritiofWithBinaryCascade_h 1

#include "globals.hh"
#include "G4HadronStoppingProcess.hh"

class G4ParticleDefinition;
class G4LundStringFragmentation;
class G4ExcitedStringDecay;


class G4HadronicAbsorptionFritiofWithBinaryCascade : public G4HadronStoppingProcess { 

public:
  G4HadronicAbsorptionFritiofWithBinaryCascade( G4ParticleDefinition* pdef = 0 ); 
  virtual ~G4HadronicAbsorptionFritiofWithBinaryCascade();
  
  G4bool IsApplicable( const G4ParticleDefinition& );

  void ProcessDescription( std::ostream& outFile ) const;

private:
  // hide assignment operator as private 
  G4HadronicAbsorptionFritiofWithBinaryCascade& operator=( const G4HadronicAbsorptionFritiofWithBinaryCascade& );
  G4HadronicAbsorptionFritiofWithBinaryCascade( const G4HadronicAbsorptionFritiofWithBinaryCascade& );
  
  G4ParticleDefinition* pdefApplicable;

  G4LundStringFragmentation* theLund;
  G4ExcitedStringDecay* theStringDecay;

};

#endif

