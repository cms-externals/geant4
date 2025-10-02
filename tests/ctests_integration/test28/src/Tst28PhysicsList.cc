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

#include <iomanip>                

#include "Tst28PhysicsList.hh"
#include "globals.hh"
#include "G4ios.hh"
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
#include "G4HadronicParameters.hh"


Tst28PhysicsList::Tst28PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst28PhysicsList::~Tst28PhysicsList()
{
}

void Tst28PhysicsList::ConstructParticle()
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

void Tst28PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst28PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst28PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst28PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst28PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst28PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst28PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void Tst28PhysicsList::ConstructEM()
{
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,-1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,1);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,2);
  
    } else if (particleName == "e+") {
      // Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,-1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,1);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,2);      
      pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,3);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      // Construct processes for muon+
      pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,-1);
      pmanager->AddProcess(new G4MuIonisation(),-1,2,1);
      pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,2);
      pmanager->AddProcess(new G4MuPairProduction(),-1,-1,3);       
     
    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,-1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,1);
 
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")&&
          (!particle->IsShortLived()) ) {
       // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,-1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,1);  
      }
    }
  }
}

// Hadron Processes
#include "G4HadronElasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4HadronInelasticProcess.hh"

// Low energy Models
#include "G4HadronElastic.hh"

// Neutron capture model and cross sections
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"

// Geneator models
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4CompetitiveFission.hh"
#include "G4GeneratorPrecompoundInterface.hh"

#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// -- bc
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"

#include "G4CascadeInterface.hh"

#include "G4EMDissociation.hh"
#include "G4EMDissociationCrossSection.hh"
#include "G4WilsonAbrasionModel.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"

void Tst28PhysicsList::ConstructHad()
{
  G4ComponentGGHadronNucleusXsc* gg = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet* xs_el = new G4CrossSectionElastic(gg);
  G4VCrossSectionDataSet* xs_in = new G4CrossSectionInelastic(gg);
  G4ComponentGGNuclNuclXsc* ggnn = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet* xs_el_nn = new G4CrossSectionElastic(ggnn);
  G4VCrossSectionDataSet* xs_in_nn = new G4CrossSectionInelastic(ggnn);

  // Build QGSP
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
       
  // Evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  //theHandler->SetMinEForMultiFrag(3*MeV);
	
  // Pre equilibrium stage 
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);

  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);  
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19*GeV);
  theTheoModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  
  // FTFP (for low and medium energy anti-baryons)
  G4TheoFSGenerator* ftfp_anti = new G4TheoFSGenerator("FTFP");
  G4FTFModel* stringModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecay =
     new G4ExcitedStringDecay(new G4LundStringFragmentation);
  stringModel->SetFragmentationModel(stringDecay);
  G4GeneratorPrecompoundInterface* cascade = new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  cascade->SetDeExcitation(preEquilib);

  // Medium energy models  
  G4BinaryCascade* theBC = new G4BinaryCascade;
  G4BinaryLightIonReaction* theIonBC= new G4BinaryLightIonReaction;
  // theIonBC->SetMinEnergy(1*MeV);
  theIonBC->SetMinEnergy(0.0);
  theIonBC->SetMaxEnergy(20*GeV);

  G4CascadeInterface* bertini = new G4CascadeInterface;
  bertini->SetMinEnergy(0.0);
  bertini->SetMaxEnergy(20*GeV);

  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
  theElasticProcess->AddDataSet(xs_el);
  theElasticProcess->RegisterMe(new G4HadronElastic());
  G4HadronElasticProcess* theElasticProcess1 = new G4HadronElasticProcess();
  theElasticProcess1->AddDataSet(xs_el);
  theElasticProcess1->RegisterMe(new G4HadronElastic());
  G4HadronElasticProcess* theElasticProcess2 = new G4HadronElasticProcess();
  theElasticProcess2->AddDataSet(xs_el_nn);
  theElasticProcess2->RegisterMe(new G4HadronElastic());

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
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(theBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {   
      // elastic scattering
      pmanager->AddDiscreteProcess(theElasticProcess1);

      // inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(theBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

      // capture
      G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
      G4NeutronRadCapture* captureModel = new G4NeutronRadCapture;
      theCaptureProcess->RegisterMe(captureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiNeutron::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Lambda::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiLambda::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4SigmaPlus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4SigmaMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiSigmaPlus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiSigmaMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4XiZero::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4XiMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
	new G4HadronInelasticProcess( "inelastic", G4OmegaMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
 
    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiXiZero::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4AntiXiMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess =
	new G4HadronInelasticProcess( "inelastic", G4AntiOmegaMinus::Definition() );
      theInelasticProcess->AddDataSet(xs_in);
      theInelasticProcess->RegisterMe(ftfp_anti);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess2);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
      theInelasticProcess->AddDataSet(xs_in_nn);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess2);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
      theInelasticProcess->AddDataSet(xs_in_nn);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess2);
      G4HadronInelasticProcess* theInelasticProcess = 
	new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );
      theInelasticProcess->AddDataSet(xs_in_nn);
      theInelasticProcess->RegisterMe(theIonBC);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "GenericIon") {
      G4HadronInelasticProcess* theInelasticProcess =
	new G4HadronInelasticProcess( "inelastic", G4GenericIon::Definition() );
      theInelasticProcess->AddDataSet(xs_in_nn);
      G4BinaryLightIonReaction* theGenIonBC= new G4BinaryLightIonReaction;
      theGenIonBC->SetMinEnergy(0*MeV);
      theGenIonBC->SetMaxEnergy(10*GeV);
      G4WilsonAbrasionModel* theGenIonWil= new G4WilsonAbrasionModel;
      theGenIonWil->SetMinEnergy(0*MeV);
      theGenIonWil->SetMaxEnergy(10*GeV);
      // theInelasticProcess->RegisterMe(theGenIonBC);
      theInelasticProcess->RegisterMe(theGenIonWil);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
	 
      // em dissociation code
      G4HadronInelasticProcess* theEMDiss =
	new G4HadronInelasticProcess( "EMdissociation", G4GenericIon::Definition() );
      G4EMDissociation* theEMDModel = new G4EMDissociation;
      theEMDModel->SetMinEnergy(0);
      theEMDModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
      G4EMDissociationCrossSection* theEMDXsec = new G4EMDissociationCrossSection;
      theEMDiss->RegisterMe(theEMDModel);
      theEMDiss->AddDataSet(theEMDXsec);
      pmanager->AddDiscreteProcess(theEMDiss);
    }
  }
}

void Tst28PhysicsList::ConstructLeptHad()
{}

#include "G4Decay.hh"
void Tst28PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

