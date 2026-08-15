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

#include "Tst33PhysicsList.hh"

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ComptonScattering.hh"
#include "G4Decay.hh"
#include "G4GammaConversion.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicParameters.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MesonConstructor.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuIonisation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuPairProduction.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "globals.hh"

#include <iomanip>

// Elastic models
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElastic.hh"

// Inelastic models
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4NeutronRadCapture.hh"
#include "G4PreCompoundModel.hh"
#include "G4TheoFSGenerator.hh"

// Cross sections
#include "G4AntiNuclElastic.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4VCrossSectionDataSet.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4ParallelWorldProcess.hh"

Tst33PhysicsList::Tst33PhysicsList() : G4VUserPhysicsList()
{
  paraWorldName.clear();
  SetVerboseLevel(1);
}

Tst33PhysicsList::~Tst33PhysicsList()
{
  paraWorldName.clear();
}

void Tst33PhysicsList::ConstructParticle()
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

void Tst33PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructProcess()
{
  AddTransportation();
  AddScoringProcess();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

void Tst33PhysicsList::ConstructEM()
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
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, -1);
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 1);
      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 2);
    }
    else if (particleName == "e+")
    {
      // positron
      //  Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, -1);

      pmanager->AddProcess(new G4eIonisation(), -1, 2, 1);
      pmanager->AddProcess(new G4eBremsstrahlung(), -1, -1, 2);
      pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 3);
    }
    else if (particleName == "mu+" || particleName == "mu-")
    {
      // muon
      //  Construct processes for muon+
      pmanager->AddProcess(new G4MuMultipleScattering(), -1, 1, -1);
      pmanager->AddProcess(new G4MuIonisation(), -1, 2, 1);
      pmanager->AddProcess(new G4MuBremsstrahlung(), -1, -1, 2);
      pmanager->AddProcess(new G4MuPairProduction(), -1, -1, 3);
    }
    else if (particleName == "GenericIon")
    {
      pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, -1);
      pmanager->AddProcess(new G4hIonisation(), -1, 2, 1);
    }
    else
    {
      if ((particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino")
          && (!particle->IsShortLived()))
      {
        // all others charged particles except geantino
        pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, -1);
        pmanager->AddProcess(new G4hIonisation(), -1, 2, 1);
      }
    }
  }
}

// Hadron Processes

void Tst33PhysicsList::ConstructHad()
{
  // Elastic models
  const G4double elastic_elimitPi = 1.0 * GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy(elastic_elimitPi);
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();
  elastic_he->SetMinEnergy(elastic_elimitPi);

  // Inelastic scattering
  const G4double theFTFMin0 = 0.0 * GeV;
  const G4double theFTFMin1 = 4.0 * GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  const G4double theBERTMin = 0.0 * GeV;
  const G4double theBERTMax = 5.0 * GeV;

  G4FTFModel* theStringModel = new G4FTFModel;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface(thePreEquilib);

  G4TheoFSGenerator* theFTFModel0 = new G4TheoFSGenerator("FTFP");
  theFTFModel0->SetHighEnergyGenerator(theStringModel);
  theFTFModel0->SetTransport(theCascade);
  theFTFModel0->SetMinEnergy(theFTFMin0);
  theFTFModel0->SetMaxEnergy(theFTFMax);

  G4TheoFSGenerator* theFTFModel1 = new G4TheoFSGenerator("FTFP");
  theFTFModel1->SetHighEnergyGenerator(theStringModel);
  theFTFModel1->SetTransport(theCascade);
  theFTFModel1->SetMinEnergy(theFTFMin1);
  theFTFModel1->SetMaxEnergy(theFTFMax);

  G4CascadeInterface* theBERTModel = new G4CascadeInterface;
  theBERTModel->SetMinEnergy(theBERTMin);
  theBERTModel->SetMaxEnergy(theBERTMax);

  G4HadronicInteraction* theIonBC = new G4BinaryLightIonReaction();
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(5 * GeV);

  G4VCrossSectionDataSet* thePiPlusData = new G4BGGPionInelasticXS(G4PionPlus::Definition());
  G4VCrossSectionDataSet* thePiMinusData = new G4BGGPionInelasticXS(G4PionMinus::Definition());

  G4VCrossSectionDataSet* theAntiNucleonData =
    new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS);

  G4ComponentGGHadronNucleusXsc* ggXsec = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet* theGGEl = new G4CrossSectionElastic(ggXsec);

  G4ComponentGGNuclNuclXsc* ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet* theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);
  G4VCrossSectionDataSet* theGGNuclNuclEl = new G4CrossSectionElastic(ggNuclNuclXsec);

  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();
  while ((*myParticleIterator)())
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGPionElasticXS(particle));
      theElasticProcess->RegisterMe(elastic_lhep1);
      theElasticProcess->RegisterMe(elastic_he);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4PionPlus::Definition());
      theInelasticProcess->AddDataSet(thePiPlusData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "pi-")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGPionElasticXS(particle));
      theElasticProcess->RegisterMe(elastic_lhep1);
      theElasticProcess->RegisterMe(elastic_he);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4PionMinus::Definition());
      theInelasticProcess->AddDataSet(thePiMinusData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // Absorption
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4PionMinus::Definition()),
                               ordDefault);
    }

    else if (particleName == "kaon+")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonPlus::Definition());
      theInelasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsKaonPlusInelasticXS::Default_Name()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "kaon0S")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroShort::Definition());
      theInelasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsKaonZeroInelasticXS::Default_Name()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "kaon0L")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroLong::Definition());
      theInelasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsKaonZeroInelasticXS::Default_Name()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "kaon-")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonMinus::Definition());
      theInelasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsKaonMinusInelasticXS::Default_Name()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // Absorption
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4KaonMinus::Definition()),
                               ordDefault);
    }

    else if (particleName == "proton")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsProtonElasticXS::Default_Name()));
      theElasticProcess->RegisterMe(elastic_chip);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Proton::Definition());
      theInelasticProcess->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "anti_proton")
    {
      // Elastic scattering
      const G4double elastic_elimitAntiNuc = 100.0 * CLHEP::MeV;
      G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
      elastic_anuc->SetMinEnergy(elastic_elimitAntiNuc);
      G4CrossSectionElastic* elastic_anucxs =
        new G4CrossSectionElastic(elastic_anuc->GetComponentCrossSection());
      G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
      elastic_lhep2->SetMaxEnergy(elastic_elimitAntiNuc);
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(elastic_anucxs);
      theElasticProcess->RegisterMe(elastic_lhep2);
      theElasticProcess->RegisterMe(elastic_anuc);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiProton::Definition());
      theInelasticProcess->AddDataSet(theAntiNucleonData);
      theInelasticProcess->RegisterMe(theFTFModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // Absorption
      pmanager->AddRestProcess(new G4HadronicAbsorptionFritiof(G4AntiProton::Definition()),
                               ordDefault);
    }

    else if (particleName == "neutron")
    {
      // elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(
        G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
          G4ChipsNeutronElasticXS::Default_Name()));
      theElasticProcess->RegisterMe(elastic_chip);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Neutron::Definition());
      theInelasticProcess->AddDataSet(new G4BGGNucleonInelasticXS(G4Neutron::Neutron()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // capture
      G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
      theCaptureProcess->AddDataSet(new G4NeutronCaptureXS());
      G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);
    }

    else if (particleName == "anti_neutron")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering (include annihilation on-fly)
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4AntiNeutron::Definition());
      theInelasticProcess->AddDataSet(theAntiNucleonData);
      theInelasticProcess->RegisterMe(theFTFModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "deuteron")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNuclNuclEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Deuteron::Definition());
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "triton")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNuclNuclEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Triton::Definition());
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }

    else if (particleName == "alpha")
    {
      // Elastic scattering
      G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(theGGNuclNuclEl);
      theElasticProcess->RegisterMe(elastic_lhep0);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Alpha::Definition());
      theInelasticProcess->AddDataSet(theGGNuclNuclData);
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void Tst33PhysicsList::ConstructLeptHad()
{
  ;
}

void Tst33PhysicsList::ConstructGeneral()
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

void Tst33PhysicsList::SetCuts()
{
  if (verboseLevel > 0)
  {
    G4cout << "Tst33PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue / mm << " (mm)" << G4endl;
  }
  //   "G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}

void Tst33PhysicsList::AddScoringProcess()
{
  G4int npw = paraWorldName.size();
  for (G4int i = 0; i < npw; i++)
  {
    G4ParallelWorldProcess* theParallelWorldProcess = new G4ParallelWorldProcess("ParaWorldProc");
    theParallelWorldProcess->SetParallelWorld(paraWorldName[i]);

    auto myParticleIterator = GetParticleIterator();
    myParticleIterator->reset();
    while ((*myParticleIterator)())
    {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(theParallelWorldProcess);
      if (theParallelWorldProcess->IsAtRestRequired(particle))
      {
        pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9900);
      }
      pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
      pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9900);
    }
  }
}
