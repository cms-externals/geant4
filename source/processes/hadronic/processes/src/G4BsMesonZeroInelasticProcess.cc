#include "G4BsMesonZeroInelasticProcess.hh"
#include "G4BsMesonZero.hh"
#include <iostream>


G4BsMesonZeroInelasticProcess::G4BsMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4BsMesonZero::Definition() ) {}


void G4BsMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4BsMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "Bs0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

