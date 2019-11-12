#include "G4BcMesonMinusInelasticProcess.hh"
#include "G4BcMesonMinus.hh"
#include <iostream>


G4BcMesonMinusInelasticProcess::G4BcMesonMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BcMesonMinus::Definition() ) {}


void G4BcMesonMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BcMesonMinusInelasticProcess handles the inelastic scattering of\n" 
          << "Bc- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

