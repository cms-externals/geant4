#ifndef G4DsMesonMinusInelasticProcess_h
#define G4DsMesonMinusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for Ds- Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4DsMesonMinusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4DsMesonMinusInelasticProcess( const G4String& processName = "Ds-Inelastic" );
    ~G4DsMesonMinusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

