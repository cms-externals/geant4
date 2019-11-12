#ifndef G4AntiXibZeroInelasticProcess_h
#define G4AntiXibZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_xi_b0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiXibZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiXibZeroInelasticProcess( const G4String& processName = "anti_xi_b0Inelastic" );
    ~G4AntiXibZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

