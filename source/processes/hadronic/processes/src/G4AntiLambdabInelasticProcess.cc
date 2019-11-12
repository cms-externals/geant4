#include "G4AntiLambdabInelasticProcess.hh"
#include "G4AntiLambdab.hh"
#include <iostream>


G4AntiLambdabInelasticProcess::G4AntiLambdabInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4AntiLambdab::Definition() ) {}


void G4AntiLambdabInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4AntiLambdabInelasticProcess handles the inelastic scattering of\n" 
          << "anti_lambda_b from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

