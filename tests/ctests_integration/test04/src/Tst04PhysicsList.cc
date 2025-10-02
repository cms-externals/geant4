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
// 080901 Add dump neutron Cross Section
//        Add Thermal Scattering by T. Koi
// 091118 Change multiple scattering processes to particle dedicated by T. Koi
// 110906 Migrate to new interface "hadr-man-V09-04-10"
//        From Process::GetMicroscopicCrossSection
//        To Process::GetElementCrossSection
//
#include <iomanip>                

#include "Tst04PhysicsList.hh"
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


Tst04PhysicsList::Tst04PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst04PhysicsList::~Tst04PhysicsList()
{
}

void Tst04PhysicsList::ConstructParticle()
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

void Tst04PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst04PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst04PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst04PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst04PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst04PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst04PhysicsList::ConstructProcess()
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

void Tst04PhysicsList::ConstructEM()
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

// Low energy neutron models
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPFission.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4CrossSectionDataStore.hh"

#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"

#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

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

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4NuclNuclDiffuseElastic.hh"
#include "G4CrossSectionElastic.hh"
#include "G4HadronicParameters.hh"

void Tst04PhysicsList::ConstructHad()
{

  //Prepairing models

  //Most hadrons
  //Bertini at low energies, then FTFP
  //theFTFPwtBERT + theBertini
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
   theIonBC->SetMinEnergy(19.9*MeV);
   theIonBC->SetMaxEnergy(4*GeV);

   // FTFP
   theFTFPforIon = theFTFPwtBERT;

   //Hadron Cross Section 
   G4VComponentCrossSection* theGGXS = new G4ComponentGGHadronNucleusXsc();
   G4VCrossSectionDataSet* theHNInel = new G4CrossSectionInelastic( theGGXS );
   G4VCrossSectionDataSet* theHNEl   = new G4CrossSectionElastic(theGGXS);

   G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
   theElasticProcess->AddDataSet(theHNEl);
   G4HadronElastic* theElasticModel = new G4HadronElastic;
   theElasticProcess->RegisterMe(theElasticModel);

   G4HadronElasticProcess* theElasticProcess1 = 
                                    new G4HadronElasticProcess;
   theElasticProcess1->AddDataSet(theHNEl);

   //Ion Inelasitc Cross Section 
   G4VComponentCrossSection* theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc();
   G4VCrossSectionDataSet* theNuclNuclData = new G4CrossSectionInelastic( theGGNuclNuclXS );
   G4VCrossSectionDataSet* ionElasticXSDataSet = new G4CrossSectionElastic(theGGNuclNuclXS);

   //Ion Elasitc Cross Section and Model 
   G4HadronElasticProcess* ionElasticProcess = new G4HadronElasticProcess("ionElastic");
   //Cross Section
   ionElasticXSDataSet->SetMinKinEnergy(0.0);
   ionElasticProcess->AddDataSet(ionElasticXSDataSet);
   //Model
   G4NuclNuclDiffuseElastic* ionElastic = new G4NuclNuclDiffuseElastic;
   ionElastic->SetMinEnergy(0.0);
   ionElasticProcess->RegisterMe(ionElastic);
   

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
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         //theInelasticProcess->RegisterMe(theBertini);
         G4CascadeInterface* theBertiniForP = new G4CascadeInterface;
         theBertiniForP->SetMinEnergy( 19.*MeV );
         theBertiniForP->SetMaxEnergy( 6.*GeV );
         theInelasticProcess->RegisterMe(theBertiniForP);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         G4ParticleHPInelasticData* hpInelasticXS = new G4ParticleHPInelasticData( G4Proton::Proton() );
         theInelasticProcess->AddDataSet( hpInelasticXS );
         G4ParticleHPInelastic* hpInelasticFS = new G4ParticleHPInelastic( G4Proton::Proton() , "ParticleHPInelastic" );
         theInelasticProcess->RegisterMe( hpInelasticFS );
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         
          // elastic scattering
         G4HadronElastic* theElasticModel1 = new G4HadronElastic;
         G4ParticleHPElastic * theElasticNeutron = new G4ParticleHPElastic;
         theElasticProcess1->RegisterMe(theElasticModel1);
         theElasticModel1->SetMinEnergy(19*MeV);
         theElasticNeutron->SetMinEnergy(4*eV);
         theElasticProcess1->RegisterMe(theElasticNeutron);
         G4ParticleHPElasticData * theNeutronData = new G4ParticleHPElasticData;
         theElasticProcess1->AddDataSet(theNeutronData);

//080901 TK add Thermal Scattering 
         G4ParticleHPThermalScattering * theThermal = new G4ParticleHPThermalScattering;
         theElasticProcess1->RegisterMe(theThermal);
         G4ParticleHPThermalScatteringData * theThermalData = new G4ParticleHPThermalScatteringData;
         theElasticProcess1->AddDataSet(theThermalData);
         pmanager->AddDiscreteProcess(theElasticProcess1);
         
          // inelastic scattering
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         
         G4CascadeInterface* theBertiniForN = new G4CascadeInterface;
         theBertiniForN->SetMinEnergy( 19.*MeV );
         theBertiniForN->SetMaxEnergy( 6.*GeV );
         
         theInelasticProcess->RegisterMe(theBertiniForN);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);

         G4ParticleHPInelastic * theLENeutronInelasticModel = new G4ParticleHPInelastic(G4Neutron::Neutron(),"ParticleHPInelastic");
         G4ParticleHPInelasticData * theNeutronData1 = new G4ParticleHPInelasticData(G4Neutron::Neutron());
         theInelasticProcess->RegisterMe(theLENeutronInelasticModel);
         theInelasticProcess->AddDataSet(theNeutronData1);

         pmanager->AddDiscreteProcess(theInelasticProcess);
         
          // fission
         G4NeutronFissionProcess* theFissionProcess =
                                    new G4NeutronFissionProcess;
         G4LFission* theFissionModel = new G4LFission;
         theFissionModel->SetMinEnergy(19*MeV);
         theFissionProcess->RegisterMe(theFissionModel);
         G4ParticleHPFission * theLENeutronFissionModel = new G4ParticleHPFission;
         theFissionProcess->RegisterMe(theLENeutronFissionModel);
         G4ParticleHPFissionData * theNeutronData2 = new G4ParticleHPFissionData;
         theFissionProcess->AddDataSet(theNeutronData2);
         pmanager->AddDiscreteProcess(theFissionProcess);
         
         // capture
         G4NeutronCaptureProcess* theCaptureProcess =
                                    new G4NeutronCaptureProcess;
         G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
         theCaptureModel->SetMinEnergy(19*MeV);
         theCaptureProcess->RegisterMe(theCaptureModel);
         G4ParticleHPCapture * theLENeutronCaptureModel = new G4ParticleHPCapture;
         theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
         G4ParticleHPCaptureData * theNeutronData3 = new G4ParticleHPCaptureData;
         theCaptureProcess->AddDataSet(theNeutronData3);
         pmanager->AddDiscreteProcess(theCaptureProcess);


//080901 TK add 

//         G4cout << G4endl;
//         G4cout << "Cross Section Dump Through HPData->DumpPhysicsTable() " << G4endl;
//         G4cout << "ParticleHP Elastic Cross Sections " << G4endl;
//         theNeutronData->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "ParticleHP Inelastic Cross Sections " << G4endl;
//         theNeutronData1->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "ParticleHP Capture Cross Sections " << G4endl;
//         theNeutronData2->DumpPhysicsTable( *G4Neutron::Neutron() );
//         G4cout << "ParticleHP Fission Cross Sections " << G4endl;
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
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Lambda::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiLambda::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaPlus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaPlus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiZero::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiZero::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(ionElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         theInelasticProcess->AddDataSet(theNuclNuclData);
         G4ParticleHPInelasticData* hpInelasticXS = new G4ParticleHPInelasticData( G4Deuteron::Deuteron() );
         theInelasticProcess->AddDataSet( hpInelasticXS );
         G4ParticleHPInelastic* hpInelasticFS = new G4ParticleHPInelastic( G4Deuteron::Deuteron() , "ParticleHPInelastic" );
         theInelasticProcess->RegisterMe( hpInelasticFS );
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(ionElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         theInelasticProcess->AddDataSet(theNuclNuclData);
         G4ParticleHPInelasticData* hpInelasticXS = new G4ParticleHPInelasticData( G4Triton::Triton() );
         theInelasticProcess->AddDataSet( hpInelasticXS );
         G4ParticleHPInelastic* hpInelasticFS = new G4ParticleHPInelastic( G4Triton::Triton() , "ParticleHPInelastic" );
         theInelasticProcess->RegisterMe( hpInelasticFS );
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "He3") {
         pmanager->AddDiscreteProcess(ionElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4He3::Definition() );
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         theInelasticProcess->AddDataSet(theNuclNuclData);
         G4ParticleHPInelasticData* hpInelasticXS = new G4ParticleHPInelasticData( G4He3::He3() );
         theInelasticProcess->AddDataSet( hpInelasticXS );
         G4ParticleHPInelastic* hpInelasticFS = new G4ParticleHPInelastic( G4He3::He3() , "ParticleHPInelastic" );
         theInelasticProcess->RegisterMe( hpInelasticFS );
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(ionElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         theInelasticProcess->AddDataSet(theNuclNuclData);
         G4ParticleHPInelasticData* hpInelasticXS = new G4ParticleHPInelasticData( G4Alpha::Alpha() );
         theInelasticProcess->AddDataSet( hpInelasticXS );
         G4ParticleHPInelastic* hpInelasticFS = new G4ParticleHPInelastic( G4Alpha::Alpha() , "ParticleHPInelastic" );
         theInelasticProcess->RegisterMe( hpInelasticFS );
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4OmegaMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiOmegaMinus::Definition() );
         theInelasticProcess->AddDataSet(theHNInel);
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
   }
}

void Tst04PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void Tst04PhysicsList::ConstructGeneral()
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

void Tst04PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Tst04PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
