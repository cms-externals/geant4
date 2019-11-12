#include "G4LambdabInelasticProcess.hh"
#include "G4Lambdab.hh"
#include <iostream>


G4LambdabInelasticProcess::G4LambdabInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4Lambdab::Definition() ) {}


void G4LambdabInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4LambdabInelasticProcess handles the inelastic scattering of\n" 
          << "lambda_b from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

