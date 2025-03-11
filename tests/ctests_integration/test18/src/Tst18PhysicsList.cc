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

#include "Tst18PhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"

#include "G4IonTable.hh"
#include "G4Ions.hh"

#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmParameters.hh"

// Hadron Processes
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

// Models
#include "G4HadronElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"

// Stopping processes
#include "G4HadronStoppingProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"
#include "G4HadronicAbsorptionFritiof.hh"


Tst18PhysicsList::Tst18PhysicsList():  G4VUserPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*m;

  SetVerboseLevel(1);

  fEmPhys = new G4EmStandardPhysics(1);
  fDecayPhys = new G4DecayPhysics(1);
  fRadDecayPhys = new G4RadioactiveDecayPhysics(1);
  G4EmParameters::Instance()->SetVerbose(1);
  G4EmParameters::Instance()->SetAuger(true);
}

Tst18PhysicsList::~Tst18PhysicsList()
{
}

void Tst18PhysicsList::ConstructParticle()
{
  fEmPhys->ConstructParticle();
  fDecayPhys->ConstructParticle();
  fRadDecayPhys->ConstructParticle();
}

void Tst18PhysicsList::ConstructProcess()
{
  AddTransportation();

  fEmPhys->ConstructProcess();
  fDecayPhys->ConstructProcess();
  fRadDecayPhys->ConstructProcess();

  // Elastic model and process
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  // Inelastic hadronic model for hadrons
  G4CascadeInterface* bertini = new G4CascadeInterface;

  // Inelastic hadronic model for ions
  G4BinaryLightIonReaction* binaryCascade = new G4BinaryLightIonReaction;
  binaryCascade->SetMinEnergy(0.0);
  binaryCascade->SetMaxEnergy(110*MeV);

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while ((*myParticleIterator)()) {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4PionPlus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4PionMinus::Definition()), ordDefault);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4HadronicAbsorptionBertini(G4KaonMinus::Definition()), ordDefault);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Lambda::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4SigmaPlus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4SigmaMinus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4XiZero::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4XiMinus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
	new G4HadronInelasticProcess( "inelastic", G4OmegaMinus::Definition() );
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

