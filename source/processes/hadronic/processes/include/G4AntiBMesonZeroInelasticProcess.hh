#ifndef G4AntiBMesonZeroInelasticProcess_h
#define G4AntiBMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_B0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiBMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiBMesonZeroInelasticProcess( const G4String& processName = "anti_B0Inelastic" );
    ~G4AntiBMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

