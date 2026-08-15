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
// PhysListEmStandardISS.cc
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandardISS.hh"

#include "G4BuilderType.hh"
#include "G4ComptonScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4EmParameters.hh"
#include "G4GammaConversion.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LossTableManager.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuIonisation.hh"
#include "G4MuPairProduction.hh"
#include "G4MuPairProductionModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RayleighScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eSingleCoulombScatteringModel.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hCoulombScatteringModel.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4hPairProduction.hh"
#include "G4hPairProductionModel.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardISS::PhysListEmStandardISS(const G4String& name, G4double Th, const G4String& Mod)
  : G4VPhysicsConstructor(name), th(Th), model(Mod)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(1);
  param->SetLowestElectronEnergy(10 * eV);
  param->SetMscThetaLimit(0.0);
  param->SetFluo(true);
  param->SetAuger(true);
  param->SetPixe(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardISS::~PhysListEmStandardISS() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardISS::ConstructProcess()
{
  G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();
  G4hBremsstrahlung* pib = new G4hBremsstrahlung();
  G4hPairProduction* pip = new G4hPairProduction();
  G4hBremsstrahlung* kb = new G4hBremsstrahlung();
  G4hPairProduction* kp = new G4hPairProduction();

  // muon & hadron scattering
  G4CoulombScattering* muss = new G4CoulombScattering(false);
  muss->SetEmModel(new G4hCoulombScatteringModel());
  G4CoulombScattering* piss = new G4CoulombScattering(false);
  piss->SetEmModel(new G4hCoulombScatteringModel());
  G4CoulombScattering* kss = new G4CoulombScattering(false);
  kss->SetEmModel(new G4hCoulombScatteringModel());

  // Add standard EM Processes
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for (const auto& particleName : partList.PartNames())
  {
    G4ParticleDefinition* particle = table->FindParticle(particleName);
    if (!particle)
    {
      continue;
    }
    if (particleName == "gamma")
    {
      G4ComptonScattering* cs = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());

      G4PhotoElectricEffect* pee = new G4PhotoElectricEffect();
      pee->SetEmModel(new G4LivermorePhotoElectricModel());

      ph->RegisterProcess(cs, particle);
      ph->RegisterProcess(pee, particle);
      ph->RegisterProcess(new G4GammaConversion(), particle);
      ph->RegisterProcess(new G4RayleighScattering(), particle);
    }
    else if (particleName == "e-")
    {
      G4CoulombScattering* ss = new G4CoulombScattering(false);
      G4eSingleCoulombScatteringModel* mod = new G4eSingleCoulombScatteringModel();

      mod->SetLowEnergyLimit(1. * keV);
      mod->SetRecoilThreshold(th);
      mod->SetXSectionModel(model);
      ss->SetEmModel(mod);

      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(ss, particle);
    }
    else if (particleName == "e+")
    {
      G4CoulombScattering* ss = new G4CoulombScattering(false);

      ph->RegisterProcess(new G4eIonisation(), particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
      ph->RegisterProcess(ss, particle);
    }
    else if (particleName == "mu+" || particleName == "mu-")
    {
      ph->RegisterProcess(new G4MuIonisation(), particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);
      ph->RegisterProcess(muss, particle);
    }
    else if (particleName == "alpha" || particleName == "He3" || particleName == "deuteron"
             || particleName == "triton" || particleName == "proton"
             || particleName == "GenericIon")
    {
      ph->RegisterProcess(new G4ionIonisation(), particle);
      G4CoulombScattering* cs = new G4CoulombScattering(false);
      G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
      mod->SetLowEnergyLimit(100. * keV);
      mod->SetRecoilThreshold(th);
      cs->SetEmModel(mod);
      ph->RegisterProcess(cs, particle);
    }
    else if (particleName == "GenericIon")
    {
      ph->RegisterProcess(new G4ionIonisation(), particle);
      ph->RegisterProcess(new G4CoulombScattering(false), particle);
    }
    else if (particleName == "pi+" || particleName == "pi-")
    {
      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pib, particle);
      ph->RegisterProcess(pip, particle);
      ph->RegisterProcess(piss, particle);
    }
    else if (particleName == "kaon+" || particleName == "kaon-")
    {
      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(kb, particle);
      ph->RegisterProcess(kp, particle);
      ph->RegisterProcess(kss, particle);
    }
    else if (particleName == "anti_proton")
    {
      G4CoulombScattering* pss = new G4CoulombScattering(false);
      pss->SetEmModel(new G4hCoulombScatteringModel());

      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pss, particle);
    }
    else if (particleName == "B+" || particleName == "B-" || particleName == "D+"
             || particleName == "D-" || particleName == "Ds+" || particleName == "Ds-"
             || particleName == "anti_He3" || particleName == "anti_alpha"
             || particleName == "anti_deuteron" || particleName == "anti_lambda_c+"
             || particleName == "anti_omega-" || particleName == "anti_sigma_c+"
             || particleName == "anti_sigma_c++" || particleName == "anti_sigma+"
             || particleName == "anti_sigma-" || particleName == "anti_triton"
             || particleName == "anti_xi_c+" || particleName == "anti_xi-"
             || particleName == "lambda_c+" || particleName == "omega-"
             || particleName == "sigma_c+" || particleName == "sigma_c++"
             || particleName == "sigma+" || particleName == "sigma-" || particleName == "tau+"
             || particleName == "tau-" || particleName == "xi_c+" || particleName == "xi-")
    {
      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(new G4CoulombScattering(false), particle);
    }
  }

  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
