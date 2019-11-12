#ifndef G4AntiOmegacZeroInelasticProcess_h
#define G4AntiOmegacZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_omega_c0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiOmegacZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiOmegacZeroInelasticProcess( const G4String& processName = "anti_omega_c0Inelastic" );
    ~G4AntiOmegacZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

