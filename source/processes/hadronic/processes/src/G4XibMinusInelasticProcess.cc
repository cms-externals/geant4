#include "G4XibMinusInelasticProcess.hh"
#include "G4XibMinus.hh"
#include <iostream>


G4XibMinusInelasticProcess::G4XibMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4XibMinus::Definition() ) {}


void G4XibMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4XibMinusInelasticProcess handles the inelastic scattering of\n" 
          << "xi_b- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

