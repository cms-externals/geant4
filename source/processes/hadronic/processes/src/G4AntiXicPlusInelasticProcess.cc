#include "G4AntiXicPlusInelasticProcess.hh"
#include "G4AntiXicPlus.hh"
#include <iostream>


G4AntiXicPlusInelasticProcess::G4AntiXicPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiXicPlus::Definition() ) {}


void G4AntiXicPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiXicPlusInelasticProcess handles the inelastic scattering of\n" 
          << "anti_xi_c+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

