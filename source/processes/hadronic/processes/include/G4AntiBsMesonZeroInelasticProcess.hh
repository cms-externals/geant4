#ifndef G4AntiBsMesonZeroInelasticProcess_h
#define G4AntiBsMesonZeroInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for anti_Bs0 Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4AntiBsMesonZeroInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4AntiBsMesonZeroInelasticProcess( const G4String& processName = "anti_Bs0Inelastic" );
    ~G4AntiBsMesonZeroInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

