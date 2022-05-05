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
// G4TransportationWithMsc
//
// Class Description:
//
// It is a generic process of transportation with multiple scattering included
// in the step limitation and propagation.
//
// Original author: Jonas Hahnfeld, 2022

#ifndef G4TrasportationWithMsc_h
#define G4TrasportationWithMsc_h 1

#include "G4Transportation.hh"

#include <vector>

class G4EmModelManager;
class G4LossTableManager;
class G4ParticleDefinition;
class G4Region;
class G4VMscModel;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4TransportationWithMsc : public G4Transportation
{
 public:
  enum class ScatteringType
  {
    MultipleScattering,
  };

  explicit G4TransportationWithMsc(ScatteringType type, G4int verbosity = 0);

  ~G4TransportationWithMsc() override;

  void AddMscModel(G4VMscModel* mscModel, G4int order = 0,
                   const G4Region* region = nullptr);

 public:
  void PreparePhysicsTable(const G4ParticleDefinition& part) override;
  void BuildPhysicsTable(const G4ParticleDefinition& part) override;

  void StartTracking(G4Track* track) override;

  G4double AlongStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4double currentMinimumStep, G4double& proposedSafety,
    G4GPILSelection* selection) override;

 private:
  const ScatteringType fType;

  G4LossTableManager* fEmManager;
  G4EmModelManager* fModelManager;
  const G4ParticleDefinition* fFirstParticle = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
