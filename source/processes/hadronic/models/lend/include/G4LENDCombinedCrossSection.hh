#ifndef G4LENDCombinedCrossSection_h
#define G4LENDCombinedCrossSection_h 1

// Class Description
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End
// 170912 First implementation done by T. Koi (SLAC/EPP)

#include "G4LENDCrossSection.hh"

class G4LENDElasticCrossSection;
class G4LENDInelasticCrossSection;
class G4LENDCaptureCrossSection;
class G4LENDFissionCrossSection;
class G4HadProjectile;

class G4LENDCombinedCrossSection : public G4LENDCrossSection
{

   public:
      G4LENDCombinedCrossSection( G4ParticleDefinition* pd );
      ~G4LENDCombinedCrossSection(){;};

      void BuildPhysicsTable( const G4ParticleDefinition& );

      G4double GetIsoCrossSection( const G4DynamicParticle*, G4int /*Z*/, G4int /*A*/ ,
                                   const G4Isotope* , const G4Element* , const G4Material* );

      G4int SelectChannel( const G4DynamicParticle*, G4int /*Z*/, G4int /*A*/ ,
                           const G4Isotope* , const G4Element* , const G4Material* );

   private:
      G4LENDElasticCrossSection* elasticXS;
      G4LENDInelasticCrossSection* inelasticXS;
      G4LENDCaptureCrossSection* captureXS;
      G4LENDFissionCrossSection* fissionXS;

};
#endif
