#ifndef G4AntiXicPlusInelasticProcess_h
#define G4AntiXicPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_xi_c+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiXicPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiXicPlusInelasticProcess( const G4String& processName = "anti_xi_c+Inelastic" );
    ~G4AntiXicPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

