#ifndef G4DsMesonPlusInelasticProcess_h
#define G4DsMesonPlusInelasticProcess_h 1

// Author: Alberto Ribon (CERN) 2019
 
// Class Description
// Process for Ds+ Inelastic scattering 
// Class Description - End

#include "G4HadronInelasticProcess.hh"
 

class G4DsMesonPlusInelasticProcess : public G4HadronInelasticProcess {
  public:  
    G4DsMesonPlusInelasticProcess( const G4String& processName = "Ds+Inelastic" );
    ~G4DsMesonPlusInelasticProcess() {}
    virtual void ProcessDescription( std::ostream& outFile ) const;
};

#endif

