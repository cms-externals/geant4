#include "G4AntiBsMesonZeroInelasticProcess.hh"
#include "G4AntiBsMesonZero.hh"
#include <iostream>


G4AntiBsMesonZeroInelasticProcess::G4AntiBsMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiBsMesonZero::Definition() ) {}


void G4AntiBsMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiBsMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_Bs0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

