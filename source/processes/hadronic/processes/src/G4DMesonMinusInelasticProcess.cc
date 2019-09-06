#include "G4DMesonMinusInelasticProcess.hh"
#include "G4DMesonMinus.hh"
#include <iostream>


G4DMesonMinusInelasticProcess::G4DMesonMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4DMesonMinus::Definition() ) {}


void G4DMesonMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4DMesonMinusInelasticProcess handles the inelastic scattering of\n" 
          << "D- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

