#include "G4BcMesonPlusInelasticProcess.hh"
#include "G4BcMesonPlus.hh"
#include <iostream>


G4BcMesonPlusInelasticProcess::G4BcMesonPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BcMesonPlus::Definition() ) {}


void G4BcMesonPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BcMesonPlusInelasticProcess handles the inelastic scattering of\n" 
          << "Bc+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

