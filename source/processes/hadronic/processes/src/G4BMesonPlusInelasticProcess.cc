#include "G4BMesonPlusInelasticProcess.hh"
#include "G4BMesonPlus.hh"
#include <iostream>


G4BMesonPlusInelasticProcess::G4BMesonPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BMesonPlus::Definition() ) {}


void G4BMesonPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BMesonPlusInelasticProcess handles the inelastic scattering of\n" 
          << "B+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

