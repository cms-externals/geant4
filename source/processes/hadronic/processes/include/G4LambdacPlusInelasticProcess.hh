#ifndef G4LambdacPlusInelasticProcess_h
#define G4LambdacPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for lambda_c+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4LambdacPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4LambdacPlusInelasticProcess( const G4String& processName = "lamba_c+Inelastic" );
    ~G4LambdacPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

