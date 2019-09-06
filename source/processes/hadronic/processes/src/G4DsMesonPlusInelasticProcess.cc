#include "G4DsMesonPlusInelasticProcess.hh"
#include "G4DsMesonPlus.hh"
#include <iostream>


G4DsMesonPlusInelasticProcess::G4DsMesonPlusInelasticProcess( const G4String& name )
 :G4HadronInelasticProcess( name, G4DsMesonPlus::Definition() ) {}


void G4DsMesonPlusInelasticProcess::ProcessDescription( std::ostream& outFile ) const {
  outFile << "G4DsMesonPlusInelasticProcess handles the inelastic scattering of\n" 
          << "Ds+ from nuclei by invoking one or more hadronic models and one\n"
          << "or more hadronic cross section sets.\n";
}

