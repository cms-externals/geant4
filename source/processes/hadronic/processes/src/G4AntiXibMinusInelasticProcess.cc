#include "G4AntiXibMinusInelasticProcess.hh"
#include "G4AntiXibMinus.hh"
#include <iostream>


G4AntiXibMinusInelasticProcess::G4AntiXibMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiXibMinus::Definition() ) {}


void G4AntiXibMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiXibMinusInelasticProcess handles the inelastic scattering of\n" 
          << "anti_xi_b- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

