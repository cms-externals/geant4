#include "G4DsMesonMinusInelasticProcess.hh"
#include "G4DsMesonMinus.hh"
#include <iostream>


G4DsMesonMinusInelasticProcess::G4DsMesonMinusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4DsMesonMinus::Definition() ) {}


void G4DsMesonMinusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4DsMesonMinusInelasticProcess handles the inelastic scattering of\n" 
          << "Ds- from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

