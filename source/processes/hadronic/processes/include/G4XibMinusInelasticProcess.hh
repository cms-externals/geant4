#ifndef G4XibMinusInelasticProcess_h
#define G4XibMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for xi_b- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4XibMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4XibMinusInelasticProcess( const G4String& processName = "xi_b-Inelastic" );
    ~G4XibMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

