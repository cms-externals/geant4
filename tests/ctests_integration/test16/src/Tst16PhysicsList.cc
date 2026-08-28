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
//

#include "Tst16PhysicsList.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MesonConstructor.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <iomanip>

Tst16PhysicsList::Tst16PhysicsList() : G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst16PhysicsList::~Tst16PhysicsList() {}

void Tst16PhysicsList::ConstructParticle()
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

void Tst16PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst16PhysicsList::ConstructProcess()
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

void Tst16PhysicsList::ConstructEM()
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
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung = new G4eBremsstrahlung();
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep, 2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);
    }
    else if (particleName == "e+")
    {
      // positron
      //  Construct processes for positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation = new G4eplusAnnihilation();
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep, 2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
    }
    else if (particleName == "mu+" || particleName == "mu-")
    {
      // muon
      //  Construct processes for muon+
      G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
      G4VProcess* aBremsstrahlung = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction = new G4MuPairProduction();
      G4VProcess* anIonisation = new G4MuIonisation();
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep, 2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(aBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(aPairProduction, idxPostStep, 4);
    }
    else if (particleName == "GenericIon")
    {
      G4VProcess* aionIonization = new G4hIonisation;
      G4VProcess* aMultipleScattering = new G4hMultipleScattering();
      pmanager->AddProcess(aionIonization);
      pmanager->AddProcess(aMultipleScattering);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(aionIonization, idxAlongStep, 2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(aionIonization, idxPostStep, 2);
    }
    else if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0)
             && (particle->GetParticleName() != "chargedgeantino"))
    {
      // all others charged particles except geantino
      G4VProcess* aMultipleScattering = new G4hMultipleScattering();
      G4VProcess* anIonisation = new G4hIonisation();
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep, 2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
    }
  }
}

// Hadron Processes
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

// Models
#include "G4BinaryLightIonReaction.hh"
#include "G4CascadeInterface.hh"
#include "G4HadronElastic.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"

//
// ConstructHad()
//
// Makes discrete physics processes for the hadrons
//

void Tst16PhysicsList::ConstructHad()
{
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  // Inelastic hadronic model for hadrons
  G4CascadeInterface* bertini = new G4CascadeInterface;

  // Inelastic hadronic model for ions
  G4BinaryLightIonReaction* binaryCascade = new G4BinaryLightIonReaction;
  binaryCascade->SetMinEnergy(0.0);
  binaryCascade->SetMaxEnergy(10 * GeV);

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
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "pi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4PionMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4PionMinus::Definition()),
                               ordDefault);
    }
    else if (particleName == "kaon+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonPlus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0S")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroShort::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0L")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonZeroLong::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4KaonMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4KaonMinus::Definition()),
                               ordDefault);
    }
    else if (particleName == "proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Proton::Definition());
      theInelasticProcess->RegisterMe(bertini);
      // now the cross-sections.
      G4VCrossSectionDataSet* theProtonData1 = new G4BGGNucleonInelasticXS(G4Proton::Proton());
      theInelasticProcess->AddDataSet(theProtonData1);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "neutron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Neutron::Definition());
      theInelasticProcess->RegisterMe(bertini);
      // now the cross-sections.
      G4VCrossSectionDataSet* theNeutronData1 = new G4NeutronInelasticXS;
      theInelasticProcess->AddDataSet(theNeutronData1);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Lambda::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4SigmaPlus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4SigmaMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4XiZero::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4XiMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4OmegaMinus::Definition());
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "deuteron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Deuteron::Definition());
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "triton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Triton::Definition());
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "alpha")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
        new G4HadronInelasticProcess("inelastic", G4Alpha::Definition());
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void Tst16PhysicsList::ConstructLeptHad()
{
  ;
}

#include "G4Decay.hh"
void Tst16PhysicsList::ConstructGeneral()
{
  // Add Decay Process
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
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void Tst16PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}
