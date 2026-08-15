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

#ifndef TEST_HH
#define TEST_HH

#include "globals.hh"

#include <ostream>

class Test
{
  public:

    Test(const G4String& physQuantity, const G4String& category, G4double refValue,
         G4double refRelError);
    ~Test();

    void SetSimulationValue(G4double value);

    void SetMetaData(G4double energy, const G4String& particle, const G4String& material)
    {
      fPrimaryEnergy = energy;
      fParticleName = particle;
      fMaterialName = material;
    }

    G4String GetPhysQuantity() { return fQuantity; }

    std::ostream& Print(std::ostream& os) const;

    G4bool Passed() { return fPassed; }

  private:

    const G4String fQuantity;
    const G4String fUnitCategory;

    const G4double fReferenceValue;
    const G4double fReferenceRelError;

    G4double fComputedValue;
    G4double fDiffPercent;
    G4int fNmbSigmas;
    G4bool fPassed;

    G4double fPrimaryEnergy;
    G4String fParticleName;
    G4String fMaterialName;
};

inline std::ostream& operator<<(std::ostream& os, const Test& test)
{
  return test.Print(os);
}

#endif
