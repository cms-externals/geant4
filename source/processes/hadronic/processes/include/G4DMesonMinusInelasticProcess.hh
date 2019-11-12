#ifndef G4DMesonMinusInelasticProcess_h
#define G4DMesonMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for D- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4DMesonMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4DMesonMinusInelasticProcess( const G4String& processName = "D-Inelastic" );
    ~G4DMesonMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

