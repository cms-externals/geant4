#ifndef G4DMesonZeroInelasticProcess_h
#define G4DMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for D0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4DMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4DMesonZeroInelasticProcess( const G4String& processName = "D0Inelastic" );
    ~G4DMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

