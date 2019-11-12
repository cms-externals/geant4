/// \file BiasingOperation.hh
/// \brief Definition of the BiasingOperation class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BiasingOperation_hh
#define BiasingOperation_hh 1

#include "G4VBiasingOperation.hh"

class G4ProtonInelasticProcess;
class G4NeutronInelasticProcess;
class G4PionPlusInelasticProcess;
class G4PionMinusInelasticProcess;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BiasingOperation : public G4VBiasingOperation {
  // The biasing operation implemented in this class is indeed a "trick" to 
  // use FTFP+INCLXX instead of FTFP+BERT for determining the final-state of
  // proton, neutron, pion+, pion- inelastic interactions happening in one
  // particular logical volume, Tracking_region, where the biasing is applied.
  public:
    BiasingOperation( G4String name );
    virtual ~BiasingOperation();
    virtual G4VParticleChange* ApplyFinalStateBiasing( const G4BiasingProcessInterface*, 
                                                       const G4Track*, const G4Step*, G4bool& );
    // Unused :
    virtual const G4VBiasingInteractionLaw* 
      ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, 
                                             G4ForceCondition& ) { return 0; }
    virtual G4double 
      DistanceToApplyOperation( const G4Track*, G4double, G4ForceCondition* ) { return DBL_MAX; }
    virtual G4VParticleChange*
      GenerateBiasingFinalState( const G4Track*, const G4Step* ) { return 0; }

  private:
    G4ProtonInelasticProcess*    fProtonInelasticProcess;
    G4NeutronInelasticProcess*   fNeutronInelasticProcess;
    G4PionPlusInelasticProcess*  fPionPlusInelasticProcess;
    G4PionMinusInelasticProcess* fPionMinusInelasticProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

