#ifndef G4BcMesonPlusInelasticProcess_h
#define G4BcMesonPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for Bc+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BcMesonPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BcMesonPlusInelasticProcess( const G4String& processName = "Bc+Inelastic" );
    ~G4BcMesonPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

