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
//
#include "Tst25PhysicsList.hh"

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4HadronicParameters.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MesonConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <iomanip>

Tst25PhysicsList::Tst25PhysicsList() : G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst25PhysicsList::~Tst25PhysicsList() {}

void Tst25PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();
}

void Tst25PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst25PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuPairProduction.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"

void Tst25PhysicsList::ConstructEM()
{
  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();
  while ((*myParticleIterator)())
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma")
    {
      // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
    }
    else if (particleName == "e-")
    {
      // electron
      //  Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);
    }
    else if (particleName == "e+")
    {
      // positron
      //  Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);
    }
    else if (particleName == "mu+" || particleName == "mu-")
    {
      // muon
      //  Construct processes for muon+/-
      pmanager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung(), -1, -1, 3);
      pmanager->AddProcess(new G4MuPairProduction(), -1, -1, 4);
    }
    else if (particleName == "GenericIon")
    {
      pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation(), -1, 2, 2);
    }
    else
    {
      if ((particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino")
          && (!particle->IsShortLived()))
      {
        // all others charged particles except geantino
        pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4hIonisation(), -1, 2, 2);
      }
    }
  }
}

// Hadron Processes
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"

// Elastic scattering model
#include "G4HadronElastic.hh"

// Neutron capture
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronRadCapture.hh"

// Light ion model
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4VCrossSectionDataSet.hh"

// Geneator models
#include "G4CascadeInterface.hh"
#include "G4CompetitiveFission.hh"
#include "G4ExcitationHandler.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4StringModel.hh"
#include "G4TheoFSGenerator.hh"

// ConstructHad()

void Tst25PhysicsList::ConstructHad()
{
  // this will be the model class for high energies
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;

  // Evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  theHandler->SetMinEForMultiFrag(3 * MeV);

  // Pre equilibrium stage
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);

  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);

  // High energy parts
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4FTFModel;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19 * GeV);
  theTheoModel->SetMaxEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());

  // FTFP (for low and medium energy anti-baryons)
  G4TheoFSGenerator* ftfp_anti = new G4TheoFSGenerator("FTFP");
  G4FTFModel* ftfModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
  ftfModel->SetFragmentationModel(stringDecay);
  G4GeneratorPrecompoundInterface* preInt = new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  preInt->SetDeExcitation(preEquilib);

  ftfp_anti->SetTransport(preInt);
  ftfp_anti->SetHighEnergyGenerator(ftfModel);
  ftfp_anti->SetMinEnergy(0.0);
  ftfp_anti->SetMaxEnergy(20.0 * GeV);

  G4CascadeInterface* theCasc = new G4CascadeInterface;
  theCasc->SetMinEnergy(65 * MeV);
  thePreEquilib->SetMaxEnergy(70 * MeV);

  G4CascadeInterface* bertini = new G4CascadeInterface;
  bertini->SetMaxEnergy(20 * GeV);

  G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction;
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(20 * GeV);
  G4VCrossSectionDataSet* ionXS = new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc);

  G4VLongitudinalStringDecay* theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  // done with the generator model (most of the above is also available as default)
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();
  while ((*myParticleIterator)())
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4PionPlus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "pi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4PionMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonPlus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0S")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroShort::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0L")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroLong::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Proton::Definition());
      theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCasc);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiProton::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "neutron")
    {
      // elastic scattering
      pmanager->AddDiscreteProcess(theElasticProcess);

      // inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Neutron::Definition());
      theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCasc);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

      // capture
      G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
      G4NeutronRadCapture* captureModel = new G4NeutronRadCapture;
      theCaptureProcess->RegisterMe(captureModel);
      G4NeutronCaptureXS* theCaptureXS = new G4NeutronCaptureXS;
      theCaptureProcess->AddDataSet(theCaptureXS);
      pmanager->AddDiscreteProcess(theCaptureProcess);
    }
    else if (particleName == "anti_neutron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiNeutron::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Lambda::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiLambda::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4SigmaPlus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4SigmaMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiSigmaPlus::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiSigmaMinus::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4XiZero::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4XiMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4OmegaMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiXiZero::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiXiMinus::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiOmegaMinus::Definition());
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "deuteron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Deuteron::Definition());
      theInelasticProcess->AddDataSet(ionXS);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "triton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Triton::Definition());
      theInelasticProcess->AddDataSet(ionXS);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "alpha")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Alpha::Definition());
      theInelasticProcess->AddDataSet(ionXS);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void Tst25PhysicsList::ConstructLeptHad() {}

#include "G4Decay.hh"
void Tst25PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();
  while ((*myParticleIterator)())
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle))
    {
      pmanager->AddProcess(theDecayProcess);
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void Tst25PhysicsList::SetCuts()
{
  if (verboseLevel > 0)
  {
    G4cout << "Tst25PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue / mm << " (mm)" << G4endl;
  }
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}
