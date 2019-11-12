/// \file BiasingOperator.cc
/// \brief Implementation of the BiasingOperator class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BiasingOperator.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "BiasingOperation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BiasingOperator::BiasingOperator() : G4VBiasingOperator( "BiasingOperator" ) {
  fBiasingOperation = new BiasingOperation( "BiasingOperation" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BiasingOperator::AddParticle( G4String particleName ) {
  const G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle( particleName );
  if ( particle == 0 ) {
    G4ExceptionDescription ed;
    ed << "Particle `" << particleName << "' not found !" << G4endl;
    G4Exception( "BiasingOperator::AddParticle(...)", "BiasError", JustWarning, ed );
    return;
  }
  fParticlesToBias.push_back( particle );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation* BiasingOperator::
ProposeFinalStateBiasingOperation( const G4Track* , 
                                   const G4BiasingProcessInterface* callingProcess ) {
  // Apply the biasing operation only for inelastic processes of:
  // proton, neutron, pion+ and pion-
  if ( callingProcess  &&  callingProcess->GetWrappedProcess()  &&  
       ( callingProcess->GetWrappedProcess()->GetProcessName() == "protonInelastic"  ||
         callingProcess->GetWrappedProcess()->GetProcessName() == "neutronInelastic" || 
         callingProcess->GetWrappedProcess()->GetProcessName() == "pi+Inelastic"     || 
         callingProcess->GetWrappedProcess()->GetProcessName() == "pi-Inelastic" ) ) {
    return fBiasingOperation;
  } else {  
    return 0; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
