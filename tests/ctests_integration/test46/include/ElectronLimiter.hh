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
#ifndef ElectronLimiter_h
#define ElectronLimiter_h 1

// V.Ivanchenko 2013/10/19
// step limiter and killer for e+,e- and other charged particles


#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4ParticleChangeForGamma.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Step;
class G4Track;
class G4Region;
class G4ParticleDefinition;

class ElectronLimiter : public G4VEmProcess {
public:
  explicit ElectronLimiter(const G4ParticleDefinition *);

  ~ElectronLimiter() override;

  G4bool IsApplicable(const G4ParticleDefinition &) override;

  void InitialiseProcess(const G4ParticleDefinition *) override;

  void StartTracking(G4Track *) override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                                G4double previousStepSize,
                                                G4ForceCondition *condition) override;

  G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &) override;

  inline void SetTrackingCutEcal(G4double cut, G4double fac, G4double rms)
  { limitEcal = cut; factEcal = fac; rmsEcal = rms; };

  inline void SetTrackingCutHcal(G4double cut, G4double fac, G4double rms)
  { limitHcal = cut; factHcal = fac; rmsHcal = rms; };

  inline const G4ParticleDefinition* GetParticle() const
  { return particle; }

private:

  const G4ParticleDefinition* particle;

  const G4Region* regionEcal;
  const G4Region* regionHcal;

  G4double limitEcal;
  G4double factEcal;
  G4double rmsEcal;

  G4double limitHcal;
  G4double factHcal;
  G4double rmsHcal;

  G4bool insideEcal;
};

#endif
