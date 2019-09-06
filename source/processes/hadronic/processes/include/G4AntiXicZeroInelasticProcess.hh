#ifndef G4AntiXicZeroInelasticProcess_h
#define G4AntiXicZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_xi_c0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiXicZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiXicZeroInelasticProcess( const G4String& processName = "anti_xi_c0Inelastic" );
    ~G4AntiXicZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

