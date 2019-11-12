#ifndef G4XicZeroInelasticProcess_h
#define G4XicZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for xi_c0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4XicZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4XicZeroInelasticProcess( const G4String& processName = "xi_c0Inelastic" );
    ~G4XicZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

