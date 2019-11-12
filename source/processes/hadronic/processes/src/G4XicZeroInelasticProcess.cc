#include "G4XicZeroInelasticProcess.hh"
#include "G4XicZero.hh"
#include <iostream>


G4XicZeroInelasticProcess::G4XicZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4XicZero::Definition() ) {}


void G4XicZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4XicZeroInelasticProcess handles the inelastic scattering of\n" 
          << "xi_c0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

