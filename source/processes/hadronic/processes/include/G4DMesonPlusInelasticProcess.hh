#ifndef G4DMesonPlusInelasticProcess_h
#define G4DMesonPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for D+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4DMesonPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4DMesonPlusInelasticProcess( const G4String& processName = "D+Inelastic" );
    ~G4DMesonPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

