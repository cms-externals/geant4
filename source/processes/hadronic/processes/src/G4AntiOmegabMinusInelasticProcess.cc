#include "G4AntiOmegabMinusInelasticProcess.hh"
#include "G4AntiOmegabMinus.hh"
#include <iostream>


G4AntiOmegabMinusInelasticProcess::G4AntiOmegabMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiOmegabMinus::Definition() ) {}


void G4AntiOmegabMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiOmegabMinusInelasticProcess handles the inelastic scattering of\n" 
          << "anti_omega_b- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

