#include "G4XicPlusInelasticProcess.hh"
#include "G4XicPlus.hh"
#include <iostream>


G4XicPlusInelasticProcess::G4XicPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4XicPlus::Definition() ) {}


void G4XicPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4XicPlusInelasticProcess handles the inelastic scattering of\n" 
          << "xi_c+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

