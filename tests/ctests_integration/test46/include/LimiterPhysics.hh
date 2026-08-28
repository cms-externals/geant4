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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef LimiterPhysics_h
#  define LimiterPhysics_h 1

#  include "G4VPhysicsConstructor.hh"
#  include "globals.hh"

#  include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LimiterPhysics : public G4VPhysicsConstructor
{
  public:

    explicit LimiterPhysics(const G4String& name = "limiter");
    ~LimiterPhysics() override;

    void ConstructParticle() override;

    void ConstructProcess() override;

    void SetParameters(const G4String& part, const G4String& reg, G4double lim, G4double fac,
                       G4double rms);

  private:

    std::vector<G4String> fPart;
    std::vector<G4String> fRegion;
    std::vector<G4double> fLimit;
    std::vector<G4double> fFactor;
    std::vector<G4double> fRMS;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
