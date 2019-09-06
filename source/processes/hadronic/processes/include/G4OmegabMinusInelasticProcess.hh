#ifndef G4OmegabMinusInelasticProcess_h
#define G4OmegabMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for omega_b- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4OmegabMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4OmegabMinusInelasticProcess( const G4String& processName = "omega_b-Inelastic" );
    ~G4OmegabMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

