#ifndef G4BMesonPlusInelasticProcess_h
#define G4BMesonPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for B+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BMesonPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BMesonPlusInelasticProcess( const G4String& processName = "B+Inelastic" );
    ~G4BMesonPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

