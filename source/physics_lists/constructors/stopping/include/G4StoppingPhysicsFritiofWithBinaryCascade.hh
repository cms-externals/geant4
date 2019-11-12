//---------------------------------------------------------------------------
//
// ClassName:  G4StoppingPhysicsFritiofWithBinaryCascade
//
// Author:     Alberto Ribon
//
// Date:       July 2019
//
// Modified:  
//
// Class Description:
//
// This class provides the nuclear capture at rest of negatively charged
// particles, using: Bertini for pi-, K-, Sigma-, Xi-, and Omega-; 
//                   Fritiof/Binary for anti-proton and anti-neutron;
//                   Fritiof/Precompound for anti-lambda, anti-sigma0, 
//                   anti-sigma+, anti-xi0 and anti-nuclei;
//                   another model for mu-.
//
//----------------------------------------------------------------------------

#ifndef G4StoppingPhysicsFritiofWithBinaryCascade_h
#define G4StoppingPhysicsFritiofWithBinaryCascade_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"


class G4StoppingPhysicsFritiofWithBinaryCascade : public G4VPhysicsConstructor {
  public: 
    G4StoppingPhysicsFritiofWithBinaryCascade( G4int ver = 1 );
    G4StoppingPhysicsFritiofWithBinaryCascade( const G4String& name, G4int ver = 1,
                                               G4bool UseMuonMinusCapture = true );
    virtual ~G4StoppingPhysicsFritiofWithBinaryCascade();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  private:
    G4int  verbose;
    static G4ThreadLocal G4bool wasActivated;
    G4bool useMuonMinusCapture;
};


#endif
