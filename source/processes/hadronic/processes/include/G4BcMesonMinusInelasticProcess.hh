#ifndef G4BcMesonMinusInelasticProcess_h
#define G4BcMesonMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for Bc- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BcMesonMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BcMesonMinusInelasticProcess( const G4String& processName = "Bc-Inelastic" );
    ~G4BcMesonMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

