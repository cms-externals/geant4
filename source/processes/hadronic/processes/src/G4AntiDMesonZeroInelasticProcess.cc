#include "G4AntiDMesonZeroInelasticProcess.hh"
#include "G4AntiDMesonZero.hh"
#include <iostream>


G4AntiDMesonZeroInelasticProcess::G4AntiDMesonZeroInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiDMesonZero::Definition() ) {}


void G4AntiDMesonZeroInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiDMesonZeroInelasticProcess handles the inelastic scattering of\n" 
          << "anti_D0 from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

