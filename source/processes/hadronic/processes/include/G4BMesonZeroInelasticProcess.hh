#ifndef G4BMesonZeroInelasticProcess_h
#define G4BMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for B0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BMesonZeroInelasticProcess( const G4String& processName = "B0Inelastic" );
    ~G4BMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

