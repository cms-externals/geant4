#include "G4AntiOmegacZeroInelasticProcess.hh"
#include "G4AntiOmegacZero.hh"
#include <iostream>


G4AntiOmegacZeroInelasticProcess::G4AntiOmegacZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiOmegacZero::Definition() ) {}


void G4AntiOmegacZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiOmegacZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_omega_c0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

