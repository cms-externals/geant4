#ifndef G4XibZeroInelasticProcess_h
#define G4XibZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for xi_b0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4XibZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4XibZeroInelasticProcess( const G4String& processName = "xi_b0Inelastic" );
    ~G4XibZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

