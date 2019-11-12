#ifndef G4OmegacZeroInelasticProcess_h
#define G4OmegacZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for omega_c0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4OmegacZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4OmegacZeroInelasticProcess( const G4String& processName = "omega_c0Inelastic" );
    ~G4OmegacZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

