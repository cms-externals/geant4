#ifndef G4XicPlusInelasticProcess_h
#define G4XicPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for xi_c+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4XicPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4XicPlusInelasticProcess( const G4String& processName = "xi_c+Inelastic" );
    ~G4XicPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

