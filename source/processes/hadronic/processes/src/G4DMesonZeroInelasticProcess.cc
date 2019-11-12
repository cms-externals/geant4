#include "G4DMesonZeroInelasticProcess.hh"
#include "G4DMesonZero.hh"
#include <iostream>


G4DMesonZeroInelasticProcess::G4DMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4DMesonZero::Definition() ) {}


void G4DMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4DMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "D0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

