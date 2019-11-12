#include "G4DMesonPlusInelasticProcess.hh"
#include "G4DMesonPlus.hh"
#include <iostream>


G4DMesonPlusInelasticProcess::G4DMesonPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4DMesonPlus::Definition() ) {}


void G4DMesonPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4DMesonPlusInelasticProcess handles the inelastic scattering of\n" 
          << "D+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

