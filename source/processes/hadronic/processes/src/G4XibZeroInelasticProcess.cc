#include "G4XibZeroInelasticProcess.hh"
#include "G4XibZero.hh"
#include <iostream>


G4XibZeroInelasticProcess::G4XibZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4XibZero::Definition() ) {}


void G4XibZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4XibZeroInelasticProcess handles the inelastic scattering of\n" 
          << "xi_b0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

