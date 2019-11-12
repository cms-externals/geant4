#include "G4AntiLambdacPlusInelasticProcess.hh"
#include "G4AntiLambdacPlus.hh"
#include <iostream>


G4AntiLambdacPlusInelasticProcess::G4AntiLambdacPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiLambdacPlus::Definition() ) {}


void G4AntiLambdacPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiLambdacPlusInelasticProcess handles the inelastic scattering of\n" 
          << "anti_lambda_c+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

