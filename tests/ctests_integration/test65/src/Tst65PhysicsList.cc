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
// 080901 Add dump neutron Cross Section
//        Add Thermal Scattering by T. Koi
// 091118 Change multiple scattering processes to particle dedicated by T. Koi
// 110906 Migrate to new interface "hadr-man-V09-04-10"
//        From Process::GetMicroscopicCrossSection
//        To Process::GetElementCrossSection
//
#include <iomanip>                

#include "Tst65PhysicsList.hh"
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


Tst65PhysicsList::Tst65PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst65PhysicsList::~Tst65PhysicsList()
{
}

void Tst65PhysicsList::ConstructParticle()
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

void Tst65PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst65PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst65PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst65PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst65PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst65PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst65PhysicsList::ConstructProcess()
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

//#include "G4MultipleScattering.hh"
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

void Tst65PhysicsList::ConstructEM()
{
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
     pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
     
     pmanager->AddProcess(new G4eIonisation(),-1,2,2);
     pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2);
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")&&
          (!particle->IsShortLived()) ) { 
     // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
     }
    }
  }
}

// Hadron Processes
#include "G4HadronElasticProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4HadronInelasticProcess.hh"

// Low-energy Models

#include "G4HadronElastic.hh"
#include "G4LFission.hh"
#include "G4NeutronRadCapture.hh"

//#include "G4LEPionPlusInelastic.hh"
//#include "G4LEPionMinusInelastic.hh"
//#include "G4LEKaonPlusInelastic.hh"
//#include "G4LEKaonZeroSInelastic.hh"
//#include "G4LEKaonZeroLInelastic.hh"
//#include "G4LEKaonMinusInelastic.hh"
//#include "G4LEProtonInelastic.hh"
//#include "G4LEAntiProtonInelastic.hh"
//#include "G4LENeutronInelastic.hh"
//#include "G4LEAntiNeutronInelastic.hh"
//#include "G4LELambdaInelastic.hh"
//#include "G4LEAntiLambdaInelastic.hh"
//#include "G4LESigmaPlusInelastic.hh"
//#include "G4LESigmaMinusInelastic.hh"
//#include "G4LEAntiSigmaPlusInelastic.hh"
//#include "G4LEAntiSigmaMinusInelastic.hh"
//#include "G4LEXiZeroInelastic.hh"
//#include "G4LEXiMinusInelastic.hh"
//#include "G4LEAntiXiZeroInelastic.hh"
//#include "G4LEAntiXiMinusInelastic.hh"
//#include "G4LEDeuteronInelastic.hh"
//#include "G4LETritonInelastic.hh"
//#include "G4LEAlphaInelastic.hh"
//#include "G4LEOmegaMinusInelastic.hh"
//#include "G4LEAntiOmegaMinusInelastic.hh"

// -- low energy neutron models
//#include "G4ParticleHPCapture.hh"
//#include "G4ParticleHPCaptureData.hh"
//#include "G4ParticleHPInelastic.hh"
//#include "G4ParticleHPInelasticData.hh"
//#include "G4ParticleHPFission.hh"
//#include "G4ParticleHPFissionData.hh"
//#include "G4ParticleHPElastic.hh"
//#include "G4ParticleHPElasticData.hh"
#include "G4LENDElasticCrossSection.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4LENDFissionCrossSection.hh"
#include "G4LENDElastic.hh"
#include "G4LENDCapture.hh"
#include "G4LENDInelastic.hh"
#include "G4LENDFission.hh"

#include "G4CrossSectionDataStore.hh"

//#include "G4ParticleHPThermalScattering.hh"
//#include "G4ParticleHPThermalScatteringData.hh"
//
// ConstructHad()
//
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering, Inelastic scattering,
// Fission (for neutron only), and Capture (neutron).
//
// F.W.Jones  06-JUL-1998
//

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4LENDCombinedCrossSection.hh"
#include "G4LENDCombinedModel.hh"
#include "G4LENDorBERTModel.hh"
#include "G4HadronicProcess.hh"
#include "G4ZeroXS.hh"
#include "G4PhotoNuclearCrossSection.hh"

void Tst65PhysicsList::ConstructHad()
{

//Prepairing models

//Most hadrons
//Bertini at low energies, then FTFP
//theFTFPwtBERT + theBertini
//
   G4TheoFSGenerator * theFTFPwtBERT;
   G4TheoFSGenerator * theFTFP;
   G4PreCompoundModel * thePreEquilib;
   G4ExcitationHandler * theHandler;
   G4GeneratorPrecompoundInterface * theCascade;
   G4FTFModel * theStringModel;
   G4ExcitedStringDecay * theStringDecay;
   G4LundStringFragmentation * theLund;
   G4CascadeInterface * theBertini;

   theFTFPwtBERT = new G4TheoFSGenerator("FTFP");
   
   theFTFPwtBERT->SetMinEnergy( 2.*GeV );
   theFTFPwtBERT->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
 
   theStringModel = new G4FTFModel;
   theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);
 
   theCascade = new G4GeneratorPrecompoundInterface;
   thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
   theCascade->SetDeExcitation(thePreEquilib);  
 
   theFTFPwtBERT->SetTransport(theCascade);
   theFTFPwtBERT->SetHighEnergyGenerator(theStringModel);
   
   theBertini = new G4CascadeInterface;
   theBertini->SetMinEnergy( 0.*GeV );
   theBertini->SetMaxEnergy( 6.*GeV );
 
//AntiHyperons:
//Use FTFP for full energy range
//theFTFP
// 
   theFTFP = new G4TheoFSGenerator("FTFP");
   theFTFP->SetMinEnergy( 0.*GeV );
   theFTFP->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
   theFTFP->SetTransport(theCascade);
   theFTFP->SetHighEnergyGenerator(theStringModel);

//Ions
//Binary Cascade + FTFP
//theIonBC + theFTFPforIon
//
   G4ExcitationHandler* handler = new G4ExcitationHandler();
   G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);

   G4TheoFSGenerator * theFTFPforIon;

// Binary Cascade
   G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction(thePreCompound);
   theIonBC->SetMinEnergy(0.0);
   theIonBC->SetMaxEnergy(4*GeV);

// FTFP
   theFTFPforIon = theFTFPwtBERT;

//Prepairing models END

   G4HadronElasticProcess* theElasticProcess = 
                                    new G4HadronElasticProcess;
   G4HadronElastic* theElasticModel = new G4HadronElastic;
   theElasticProcess->RegisterMe(theElasticModel);

   auto myParticleIterator=GetParticleIterator();
   myParticleIterator->reset();
   while ((*myParticleIterator)()) {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "gamma") {

      G4HadronInelasticProcess* thePhotoNuclearProcess =
	new G4HadronInelasticProcess( "photonNuclear", G4Gamma::Definition() );
      thePhotoNuclearProcess->AddDataSet( new G4PhotoNuclearCrossSection );
      G4LENDCombinedCrossSection* endlInelasticXS = new G4LENDCombinedCrossSection( particle );
      //endlInelasticXS->ChangeDefaultEvaluation( "ENDF.B-VII.0" );
      endlInelasticXS->AllowNaturalAbundanceTarget();
      //endlInelasticXS->AllowAnyCandidateTarget();
      endlInelasticXS->DumpLENDTargetInfo( true );
      //G4LENDCombinedModel* endlInelasticFS = new G4LENDCombinedModel( particle );
      G4LENDorBERTModel* endlInelasticFS = new G4LENDorBERTModel( particle );
      //endlInelasticFS->ChangeDefaultEvaluation( "ENDF.B-VII.0" );
      endlInelasticFS->AllowNaturalAbundanceTarget();
      //endlInelasticFS->AllowAnyCandidateTarget();
      endlInelasticFS->DumpLENDTargetInfo( true );
      //
      G4CascadeInterface* theGammaReaction = new G4CascadeInterface;
      theGammaReaction->SetMinEnergy(20*MeV);

      
      thePhotoNuclearProcess->AddDataSet( endlInelasticXS );
      thePhotoNuclearProcess->RegisterMe( endlInelasticFS );
      thePhotoNuclearProcess->RegisterMe( theGammaReaction );
      pmanager->AddDiscreteProcess( thePhotoNuclearProcess );

/*
      G4HadronicProcess* thePhotoFissionProcess = new G4HadronicProcess( "photoFission", fFission );
      G4LENDFissionCrossSection* endlFissionXS = new G4LENDFissionCrossSection( particle );
      G4LENDFission* endlFissionFS = new G4LENDFission( particle );
      thePhotoFissionProcess->AddDataSet( new G4ZeroXS );
      thePhotoFissionProcess->AddDataSet( endlFissionXS );
      thePhotoFissionProcess->RegisterMe( endlFissionFS );
      pmanager->AddDiscreteProcess( thePhotoFissionProcess );

      G4HadronicProcess* thePhotoCaptureProcess = new G4HadronicProcess( "photoCapture", fCapture );
      G4LENDCaptureCrossSection* endlCaptureXS = new G4LENDCaptureCrossSection( particle );
      G4LENDCapture* endlCaptureFS = new G4LENDCapture( particle );
      thePhotoCaptureProcess->AddDataSet( new G4ZeroXS );
      thePhotoCaptureProcess->AddDataSet( endlCaptureXS );
      thePhotoCaptureProcess->RegisterMe( endlCaptureFS );
      pmanager->AddDiscreteProcess( thePhotoCaptureProcess );
*/

      }
      else if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4PionPlus::Definition() );
         //G4LEPionPlusInelastic* theInelasticModel = 
         //                       new G4LEPionPlusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
         //G4LEPionMinusInelastic* theInelasticModel = 
         //                       new G4LEPionMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
         //G4LEKaonPlusInelastic* theInelasticModel = new G4LEKaonPlusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
         //G4LEKaonZeroSInelastic* theInelasticModel = 
         //                    new G4LEKaonZeroSInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
         //G4LEKaonZeroLInelastic* theInelasticModel = 
         //                    new G4LEKaonZeroLInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
         //G4LEKaonMinusInelastic* theInelasticModel = 
         //                        new G4LEKaonMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
         //G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
         //G4LEAntiProtonInelastic* theInelasticModel = 
         //                      new G4LEAntiProtonInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         
         // elastic scattering
         G4HadronElasticProcess* theElasticProcess1 = 
                                    new G4HadronElasticProcess;
         G4HadronElastic* theElasticModel1 = new G4HadronElastic;
         G4LENDElastic * theElasticNeutron = new G4LENDElastic(G4Neutron::Neutron());
         theElasticNeutron->AllowNaturalAbundanceTarget();
         theElasticNeutron->DumpLENDTargetInfo( true );
         theElasticProcess1->RegisterMe(theElasticModel1);
         theElasticModel1->SetMinEnergy(19*MeV);
         //theElasticNeutron->SetMinEnergy(4*eV);
         theElasticProcess1->RegisterMe(theElasticNeutron);
         G4LENDElasticCrossSection* theNeutronData = new G4LENDElasticCrossSection(G4Neutron::Neutron());
         theNeutronData->AllowNaturalAbundanceTarget();
         theNeutronData->DumpLENDTargetInfo( true );
         theElasticProcess1->AddDataSet(theNeutronData);

//080901 TK add Thermal Scattering 
         //G4ParticleHPThermalScattering * theThermal = new G4ParticleHPThermalScattering;
         //theElasticProcess1->RegisterMe(theThermal);
         //G4ParticleHPThermalScatteringData * theThermalData = new G4ParticleHPThermalScatteringData;
         //theElasticProcess1->AddDataSet(theThermalData);
         pmanager->AddDiscreteProcess(theElasticProcess1);
         
          // inelastic scattering
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
         //G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
         //theInelasticModel->SetMinEnergy(19*MeV);
         //theInelasticProcess->RegisterMe(theInelasticModel);
         
         G4CascadeInterface* theBertiniForN = new G4CascadeInterface;
         theBertiniForN->SetMinEnergy( 19.*MeV );
         theBertiniForN->SetMaxEnergy( 6.*GeV );
         
         theInelasticProcess->RegisterMe(theBertiniForN);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         G4LENDInelastic * theLENeutronInelasticModel = new G4LENDInelastic(G4Neutron::Neutron());
         theLENeutronInelasticModel->AllowNaturalAbundanceTarget();
         theLENeutronInelasticModel->DumpLENDTargetInfo( true );
         theInelasticProcess->RegisterMe(theLENeutronInelasticModel);
         G4LENDInelasticCrossSection * theNeutronData1 = new G4LENDInelasticCrossSection(G4Neutron::Neutron());
         theNeutronData1->AllowNaturalAbundanceTarget();
         theNeutronData1->DumpLENDTargetInfo( true );
         theInelasticProcess->AddDataSet(theNeutronData1);
         pmanager->AddDiscreteProcess(theInelasticProcess);
         
          // fission
         G4NeutronFissionProcess* theFissionProcess =
                                    new G4NeutronFissionProcess;
         G4LFission* theFissionModel = new G4LFission;
         theFissionModel->SetMinEnergy(19*MeV);
         theFissionProcess->RegisterMe(theFissionModel);
         G4LENDFission * theLENeutronFissionModel = new G4LENDFission(G4Neutron::Neutron());
         theLENeutronFissionModel->DumpLENDTargetInfo( true );
         theFissionProcess->RegisterMe(theLENeutronFissionModel);
         G4LENDFissionCrossSection * theNeutronData2 = new G4LENDFissionCrossSection(G4Neutron::Neutron());
         theNeutronData2->DumpLENDTargetInfo( true );
         theFissionProcess->AddDataSet(theNeutronData2);
         pmanager->AddDiscreteProcess(theFissionProcess);
         
         // capture
         G4NeutronCaptureProcess* theCaptureProcess =
                                    new G4NeutronCaptureProcess;
         G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
         theCaptureModel->SetMinEnergy(19*MeV);
         theCaptureProcess->RegisterMe(theCaptureModel);
         G4LENDCapture * theLENeutronCaptureModel = new G4LENDCapture(G4Neutron::Neutron());
         theLENeutronCaptureModel->AllowNaturalAbundanceTarget();
         theLENeutronCaptureModel->DumpLENDTargetInfo( true );
         theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
         G4LENDCaptureCrossSection * theNeutronData3 = new G4LENDCaptureCrossSection(G4Neutron::Neutron());
         theNeutronData3->AllowNaturalAbundanceTarget();
         theNeutronData3->DumpLENDTargetInfo( true );
         theCaptureProcess->AddDataSet(theNeutronData3);
         pmanager->AddDiscreteProcess(theCaptureProcess);


//080901 TK add 

//         G4cout << G4endl;
//         G4cout << "Cross Section Dump Through HPData->DumpPhysicsTable() " << G4endl;
//         G4cout << "NeutronHP Elastic Cross Sections " << G4endl;
//         theNeutronData->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "NeutronHP Inelastic Cross Sections " << G4endl;
//         theNeutronData1->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "NeutronHP Capture Cross Sections " << G4endl;
//         theNeutronData2->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "NeutronHP Fission Cross Sections " << G4endl;
//         theNeutronData3->DumpPhysicsTable( *G4Neutron::Neutron() );

/*
         G4cout << G4endl;
         G4cout << "Cross Section Dump Through Process::GetElementCrossSection " << G4endl;
         G4int nmat = G4Material::GetNumberOfMaterials();
         std::map< G4int , G4double > tmp_map;
         for ( G4int imat = 0 ; imat < nmat ; imat ++ )
         {
            const G4Material* mat = (*(G4Material::GetMaterialTable()))[ imat ];
            G4int ne = mat->GetNumberOfElements();
            for ( G4int iele = 0 ; iele < ne ; iele ++ )
            {

               // 110906 TK migaration to "hadr-man-V09-04-10"
               //const G4Element* ele = (*(G4Element::GetElementTable()))[ iele ];
               const G4Element* ele = mat->GetElement( iele );
               std::pair < G4int , G4double > apair( ele->GetIndex() , mat->GetTemperature());
               if ( !(tmp_map.insert ( apair )).second ) continue; 

               //G4cout << ele->GetName() << ": " << mat->GetTemperature()/kelvin << " kelvin " << mat->GetName() << G4endl;
               G4cout << ele->GetName() << ": " << mat->GetTemperature()/kelvin << " kelvin " << G4endl;
               G4cout << "Energy[eV]" << '\t' << "Elastic[mb]" << '\t' << "Inelastic[mb]" << '\t' << "Capture[mb]" << '\t' << "Fission[mb]" << G4endl;

               for ( G4int ie = 1 ; ie < 150 ; ie++ )
               {
                  G4double e = 1.0e-5*eV*std::pow( 10.0 , 1.0*ie/10 );  
                  G4DynamicParticle* dp = new G4DynamicParticle( G4Neutron::Neutron(), G4ThreeVector(1,0,0), e);
                  G4cout.precision(7);
                  G4cout << std::scientific << e/eV 
                                    << '\t' << ((G4HadronicProcess*)theElasticProcess1) ->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theInelasticProcess)->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theCaptureProcess)  ->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theFissionProcess)  ->GetElementCrossSection( dp , ele , mat )/millibarn
                                            << G4endl;
               } 
               // 1 - 20 MeV
//               for ( G4int ie = 1 ; ie < 20 ; ie++ )
//               {
//                  G4DynamicParticle* dp = new G4DynamicParticle( G4Neutron::Neutron(), G4ThreeVector(1,0,0), ie*MeV);
//                  G4cout << ie*MeV/eV   
//                         << " " << ((G4HadronicProcess*)theElasticProcess1)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theInelasticProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theCaptureProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theFissionProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << G4endl;
//               }

            }
         }
         G4cout << G4endl;
*/

      }  
      else if (particleName == "anti_neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiNeutron::Definition() );
         //G4LEAntiNeutronInelastic* theInelasticModel = 
         //                      new G4LEAntiNeutronInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Lambda::Definition() );
         //G4LELambdaInelastic* theInelasticModel = new G4LELambdaInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiLambda::Definition() );
         //G4LEAntiLambdaInelastic* theInelasticModel = 
         //                       new G4LEAntiLambdaInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaPlus::Definition() );
         //G4LESigmaPlusInelastic* theInelasticModel = 
         //                        new G4LESigmaPlusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaMinus::Definition() );
         //G4LESigmaMinusInelastic* theInelasticModel = 
         //                        new G4LESigmaMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaPlus::Definition() );
         //G4LEAntiSigmaPlusInelastic* theInelasticModel = 
         //                        new G4LEAntiSigmaPlusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaMinus::Definition() );
         //G4LEAntiSigmaMinusInelastic* theInelasticModel = 
         //                        new G4LEAntiSigmaMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiZero::Definition() );
         //G4LEXiZeroInelastic* theInelasticModel = 
         //                        new G4LEXiZeroInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiMinus::Definition() );
         //G4LEXiMinusInelastic* theInelasticModel = 
         //                        new G4LEXiMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiZero::Definition() );
         //G4LEAntiXiZeroInelastic* theInelasticModel = 
         //                        new G4LEAntiXiZeroInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiMinus::Definition() );
         //G4LEAntiXiMinusInelastic* theInelasticModel = 
         //                        new G4LEAntiXiMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
         //G4LEDeuteronInelastic* theInelasticModel = 
         //                        new G4LEDeuteronInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
         //G4LETritonInelastic* theInelasticModel = 
         //                        new G4LETritonInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );
         //G4LEAlphaInelastic* theInelasticModel = 
         //                        new G4LEAlphaInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4OmegaMinus::Definition() );
         //G4LEOmegaMinusInelastic* theInelasticModel = 
         //                        new G4LEOmegaMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiOmegaMinus::Definition() );
         //G4LEAntiOmegaMinusInelastic* theInelasticModel = 
         //                        new G4LEAntiOmegaMinusInelastic;
         //theInelasticProcess->RegisterMe(theInelasticModel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
   }
}

void Tst65PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void Tst65PhysicsList::ConstructGeneral()
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

void Tst65PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Tst65PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
