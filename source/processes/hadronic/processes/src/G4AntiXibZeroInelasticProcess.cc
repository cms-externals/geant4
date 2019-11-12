#include "G4AntiXibZeroInelasticProcess.hh"
#include "G4AntiXibZero.hh"
#include <iostream>


G4AntiXibZeroInelasticProcess::G4AntiXibZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiXibZero::Definition() ) {}


void G4AntiXibZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiXibZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_xi_b0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

