#ifndef G4AntiLambdabInelasticProcess_h
#define G4AntiLambdabInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_lambda_b Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiLambdabInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiLambdabInelasticProcess( const G4String& processName = "anti_lamba_bInelastic" );
    ~G4AntiLambdabInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

