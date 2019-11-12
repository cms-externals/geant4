#ifndef G4BsMesonZeroInelasticProcess_h
#define G4BsMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for Bs0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4BsMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4BsMesonZeroInelasticProcess( const G4String& processName = "Bs0Inelastic" );
    ~G4BsMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

