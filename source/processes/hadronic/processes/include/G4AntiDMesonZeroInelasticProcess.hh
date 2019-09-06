#ifndef G4AntiDMesonZeroInelasticProcess_h
#define G4AntiDMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_D0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiDMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiDMesonZeroInelasticProcess( const G4String& processName = "anti_D0Inelastic" );
    ~G4AntiDMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

