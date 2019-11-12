#ifndef G4AntiXibMinusInelasticProcess_h
#define G4AntiXibMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_xi_b- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiXibMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiXibMinusInelasticProcess( const G4String& processName = "anti_xi_b-Inelastic" );
    ~G4AntiXibMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

