#ifndef G4AntiOmegabMinusInelasticProcess_h
#define G4AntiOmegabMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_omega_b- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiOmegabMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiOmegabMinusInelasticProcess( const G4String& processName = "anti_omega_b-Inelastic" );
    ~G4AntiOmegabMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

