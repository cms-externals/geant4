#include "G4BMesonMinusInelasticProcess.hh"
#include "G4BMesonMinus.hh"
#include <iostream>


G4BMesonMinusInelasticProcess::G4BMesonMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BMesonMinus::Definition() ) {}


void G4BMesonMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BMesonMinusInelasticProcess handles the inelastic scattering of\n" 
          << "B- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

