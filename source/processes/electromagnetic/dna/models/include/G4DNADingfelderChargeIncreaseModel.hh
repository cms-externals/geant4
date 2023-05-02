//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#ifndef G4DNADingfelderChargeIncreaseModel_h
#define G4DNADingfelderChargeIncreaseModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4NistManager.hh"

class G4DNADingfelderChargeIncreaseModel : public G4VEmModel
{
public:

  explicit G4DNADingfelderChargeIncreaseModel(const G4ParticleDefinition* p = nullptr,
		          const G4String& nam = "DNADingfelderChargeIncreaseModel");
  ~G4DNADingfelderChargeIncreaseModel() override = default;
  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  G4double CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax) override;
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;
  inline void SelectStationary(G4bool input);
protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;
private:
  // Water density table
  const std::vector<G4double>* fpMolWaterDensity = nullptr;

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  G4bool isInitialised = false, statCode = false;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  G4int verboseLevel = 0;
  
  // Partial cross section

  G4double PartialCrossSection(const G4double& energy, const G4int& level, const G4ParticleDefinition* particle);

  G4double Sum(const G4double& energy, const G4ParticleDefinition* particle);

  G4int RandomSelect(const G4double& energy, const G4ParticleDefinition* particle);
  
  G4int numberOfPartialCrossSections[2] = {0}; // 2 is the particle type index
  G4double f0[2][2] = {{0, 0},{0, 0}};
  G4double a0[2][2] = {{0, 0},{0, 0}};
  G4double a1[2][2] = {{0, 0},{0, 0}};
  G4double b0[2][2] = {{0, 0},{0, 0}};
  G4double b1[2][2] = {{0, 0},{0, 0}};
  G4double c0[2][2] = {{0, 0},{0, 0}};
  G4double d0[2][2] = {{0, 0},{0, 0}};
  G4double x0[2][2] = {{0, 0},{0, 0}};
  G4double x1[2][2] = {{0, 0},{0, 0}};

  // Final state
  G4int NumberOfFinalStates(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
  G4ParticleDefinition* OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
  G4double IncomingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition, G4int finalStateIndex);
  G4DNADingfelderChargeIncreaseModel & operator=(const  G4DNADingfelderChargeIncreaseModel &right)= delete;
  G4DNADingfelderChargeIncreaseModel(const  G4DNADingfelderChargeIncreaseModel&) = delete;
  // Reusable particle definitions
  G4ParticleDefinition* hydrogenDef = nullptr;
  G4ParticleDefinition* alphaPlusPlusDef = nullptr;
  G4ParticleDefinition* alphaPlusDef = nullptr;
  G4ParticleDefinition* heliumDef = nullptr;
  
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4DNADingfelderChargeIncreaseModel::SelectStationary (G4bool input)
{ 
  statCode = input;
}		 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

