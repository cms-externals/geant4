#include "G4AntiXicZeroInelasticProcess.hh"
#include "G4AntiXicZero.hh"
#include <iostream>


G4AntiXicZeroInelasticProcess::G4AntiXicZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiXicZero::Definition() ) {}


void G4AntiXicZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiXicZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_xi_c0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

