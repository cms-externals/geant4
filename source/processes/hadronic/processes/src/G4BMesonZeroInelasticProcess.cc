#include "G4BMesonZeroInelasticProcess.hh"
#include "G4BMesonZero.hh"
#include <iostream>


G4BMesonZeroInelasticProcess::G4BMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BMesonZero::Definition() ) {}


void G4BMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "B0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

