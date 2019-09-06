#ifndef G4AntiLambdacPlusInelasticProcess_h
#define G4AntiLambdacPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_lambda_c+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiLambdacPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiLambdacPlusInelasticProcess( const G4String& processName = "anti_lamba_c+Inelastic" );
    ~G4AntiLambdacPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

