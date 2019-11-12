#ifndef G4BMesonMinusInelasticProcess_h
#define G4BMesonMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for B- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BMesonMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BMesonMinusInelasticProcess( const G4String& processName = "B-Inelastic" );
    ~G4BMesonMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

