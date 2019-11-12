#ifndef G4LambdabInelasticProcess_h
#define G4LambdabInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for lambda_b Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4LambdabInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4LambdabInelasticProcess( const G4String& processName = "lamba_bInelastic" );
    ~G4LambdabInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

