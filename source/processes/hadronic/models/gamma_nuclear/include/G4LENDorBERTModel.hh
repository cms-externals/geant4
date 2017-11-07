#ifndef G4LENDorBERTModel_h
#define G4LENDorBERTModel_h 1

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

#include "G4LENDModel.hh"

class G4LENDCombinedModel;
class G4CascadeInterface;

class G4LENDorBERTModel : public G4LENDModel
{
   public: 
      G4LENDorBERTModel( G4ParticleDefinition* pd );
      ~G4LENDorBERTModel(){;};

      void BuildPhysicsTable(const G4ParticleDefinition&);

      G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );
     
   private: 
      G4LENDCombinedModel* lend;
      G4CascadeInterface* bert;
};

#endif
