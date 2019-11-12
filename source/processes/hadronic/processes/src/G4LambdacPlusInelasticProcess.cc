#include "G4LambdacPlusInelasticProcess.hh"
#include "G4LambdacPlus.hh"
#include <iostream>


G4LambdacPlusInelasticProcess::G4LambdacPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4LambdacPlus::Definition() ) {}


void G4LambdacPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4LambdacPlusInelasticProcess handles the inelastic scattering of\n" 
          << "lambda_c+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

