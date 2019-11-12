#include "G4OmegabMinusInelasticProcess.hh"
#include "G4OmegabMinus.hh"
#include <iostream>


G4OmegabMinusInelasticProcess::G4OmegabMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4OmegabMinus::Definition() ) {}


void G4OmegabMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4OmegabMinusInelasticProcess handles the inelastic scattering of\n" 
          << "omega_b- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

