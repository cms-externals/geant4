#include "G4OmegacZeroInelasticProcess.hh"
#include "G4OmegacZero.hh"
#include <iostream>


G4OmegacZeroInelasticProcess::G4OmegacZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4OmegacZero::Definition() ) {}


void G4OmegacZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4OmegacZeroInelasticProcess handles the inelastic scattering of\n" 
          << "omega_c0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

