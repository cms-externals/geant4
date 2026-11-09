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

#include "G4DNAIonisation.hh"

#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4LEPTSIonisationModel.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4SystemOfUnits.hh"

// SEB
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAIonisation::G4DNAIonisation(const G4String& processName, G4ProcessType type)
  : G4VEmProcess(processName, type)
{
  SetProcessSubType(fLowEnergyIonisation);
  SetBuildTableFlag(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNAIonisation::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAIonisation::InitialiseProcess(const G4ParticleDefinition* p)
{
  if (!isInitialised)
  {
    isInitialised = true;

    G4String name = p->GetParticleName();
    G4EmParameters* param = G4EmParameters::Instance();
    const G4double emaxElectronDNA = param->MaxDNAElectronEnergy();
    const G4double emaxProtonDNA = param->MaxDNAProtonEnergy();
    const G4double emaxIonDNA = param->MaxDNAIonEnergyPerNucleon();

    if (name == "e-")
    {
      if (EmModel() == nullptr)
      {
        auto born = new G4DNABornIonisationModel();
        SetEmModel(born);
        born->SetHighEnergyLimit(emaxElectronDNA);
      }
      AddEmModel(1, EmModel());
    }
    else if (name == "e+")
    {
      if (EmModel() == nullptr)
      {
        auto lepts = new G4LEPTSIonisationModel();
        SetEmModel(lepts);
        lepts->SetLowEnergyLimit(CLHEP::eV);
        lepts->SetHighEnergyLimit(emaxElectronDNA);
      }
      AddEmModel(1, EmModel());
    }
    else if (name == "proton")
    {
      if (EmModel(0) == nullptr)  // MK : Is this a reliable test ? VI: it is useful
      {
        G4double elim = 500 * CLHEP::keV;
        auto rudd = new G4DNARuddIonisationExtendedModel();
        rudd->SetHighEnergyLimit(elim);
        SetEmModel(rudd);

        auto born = new G4DNABornIonisationModel();
        born->SetLowEnergyLimit(elim);
        born->SetHighEnergyLimit(emaxProtonDNA);
        SetEmModel(born);
      }

      AddEmModel(1, EmModel());
      if (EmModel(1) != nullptr) AddEmModel(2, EmModel(1));
    }
    else
    {
      if (EmModel() == nullptr)
      {
        auto rudd = new G4DNARuddIonisationExtendedModel();
        SetEmModel(rudd);
        rudd->SetHighEnergyLimit(emaxIonDNA);
      }
      AddEmModel(1, EmModel());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNAIonisation::PrintInfo()
{
  if (EmModel(1) != nullptr)
  {
    G4cout << " Total cross sections computed from " << EmModel(0)->GetName() << " and "
           << EmModel(1)->GetName() << " models" << G4endl;
  }
  else
  {
    G4cout << " Total cross sections computed from " << EmModel()->GetName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
