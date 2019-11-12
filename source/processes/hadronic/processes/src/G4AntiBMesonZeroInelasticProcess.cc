#include "G4AntiBMesonZeroInelasticProcess.hh"
#include "G4AntiBMesonZero.hh"
#include <iostream>


G4AntiBMesonZeroInelasticProcess::G4AntiBMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiBMesonZero::Definition() ) {}


void G4AntiBMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiBMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_B0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

