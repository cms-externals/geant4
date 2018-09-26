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
// $Id: G4HadronPhysicsShielding.cc 111436 2018-08-10 07:40:29Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2010 Tatsumi Koi, Gunter Folger
//   created from G4HadronPhysicsFTFP_BERT
//
// Modified:
//
// 2014.08.05 K.L.Genser added provisions for modifing the Bertini to
//            FTF transition energy region
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsShielding.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4ParticleHPBGGNucleonInelasticXS.hh"
#include "G4ParticleHPJENDLHEInelasticData.hh"
#include "G4ParticleHPInelasticData.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsShielding);

G4ThreadLocal G4HadronPhysicsShielding::ThreadPrivate* 
G4HadronPhysicsShielding::tpdata = 0;

G4HadronPhysicsShielding::G4HadronPhysicsShielding( G4int )
    :  G4VPhysicsConstructor("hInelastic Shielding")
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(9.5*GeV)
    , maxBertiniEnergy_(9.9*GeV)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name) 
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(9.5*GeV)
    , maxBertiniEnergy_(9.9*GeV)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}

G4HadronPhysicsShielding::G4HadronPhysicsShielding(const G4String& name,
                                G4int /*verbose*/, G4double minFTFPEnergy, G4double maxBertiniEnergy)
    :  G4VPhysicsConstructor(name)
    , useLEND_(false)
    , evaluation_()
    , minFTFPEnergy_(minFTFPEnergy)
    , maxBertiniEnergy_(maxBertiniEnergy)
    , minNonHPNeutronEnergy_(19.9*MeV)
{}

#include "G4NeutronLENDBuilder.hh"
void G4HadronPhysicsShielding::CreateModels()
{
  G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)

  tpdata->theNeutrons=new G4NeutronBuilder( true ); // Fission on
  tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(quasiElasticFTF);
  tpdata->theFTFPNeutron->SetMinEnergy(minFTFPEnergy_);
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron);
  tpdata->theNeutrons->RegisterMe(tpdata->theBertiniNeutron=new G4BertiniNeutronBuilder);
  tpdata->theBertiniNeutron->SetMinEnergy(minNonHPNeutronEnergy_);
  tpdata->theBertiniNeutron->SetMaxEnergy(maxBertiniEnergy_);
  //tpdata->theNeutrons->RegisterMe(tpdata->theHPNeutron=new G4NeutronPHPBuilder);

  if ( useLEND_ != true )
     tpdata->theNeutrons->RegisterMe(tpdata->theLENeutron=new G4NeutronPHPBuilder);
  else
  {
     tpdata->theNeutrons->RegisterMe(tpdata->theLENeutron=new G4NeutronLENDBuilder(evaluation_));
  }

  tpdata->thePro=new G4ProtonBuilder;
  tpdata->theFTFPPro=new G4FTFPProtonBuilder(quasiElasticFTF);
  tpdata->theFTFPPro->SetMinEnergy(minFTFPEnergy_);
  tpdata->thePro->RegisterMe(tpdata->theFTFPPro);
  tpdata->thePro->RegisterMe(tpdata->theBertiniPro=new G4BertiniProtonBuilder);
  tpdata->theBertiniPro->SetMaxEnergy(maxBertiniEnergy_);

  tpdata->thePiK=new G4PiKBuilder;
  tpdata->theFTFPPiK=new G4FTFPPiKBuilder(quasiElasticFTF);
  tpdata->theFTFPPiK->SetMinEnergy(minFTFPEnergy_);
  tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK);
  tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK=new G4BertiniPiKBuilder);
  tpdata->theBertiniPiK->SetMaxEnergy(maxBertiniEnergy_);

  tpdata->theHyperon=new G4HyperonFTFPBuilder;
    
  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

G4HadronPhysicsShielding::~G4HadronPhysicsShielding()
{
  if (!tpdata) return;

  delete tpdata->theNeutrons;
  delete tpdata->theBertiniNeutron;
  delete tpdata->theFTFPNeutron;
  //delete tpdata->theHPNeutron;
  delete tpdata->theLENeutron;
    
  delete tpdata->thePiK;
  delete tpdata->theBertiniPiK;
  delete tpdata->theFTFPPiK;
    
  delete tpdata->thePro;
  delete tpdata->theBertiniPro;
  delete tpdata->theFTFPPro;    
    
  delete tpdata->theHyperon;
  delete tpdata->theAntiBaryon;
  delete tpdata->theFTFPAntiBaryon;

  delete tpdata; tpdata=0;
}

void G4HadronPhysicsShielding::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
}

#include "G4ProcessManager.hh"
void G4HadronPhysicsShielding::ConstructProcess()
{
  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();

  tpdata->thePro->Build();
  tpdata->theNeutrons->Build();
    
  G4ParticleHPJENDLHEInelasticData* theNeutronHPJENDLHEInelastic=new G4ParticleHPJENDLHEInelasticData;
  //Register the G4ParticleHPJENDLHEInelasticData as the 2nd priority.
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->GetCrossSectionDataStore()->AddDataSet(theNeutronHPJENDLHEInelastic,1);
    
  tpdata->thePiK->Build();

  // --- Kaons ---
  G4ComponentGGHadronNucleusXsc* xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);

  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (capture) {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
    theNeutronRadCapture->SetMinEnergy( minNonHPNeutronEnergy_ ); 
    capture->RegisterMe( theNeutronRadCapture );
  }
  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess(neutron);
  if (fission) {
    G4LFission* theNeutronLEPFission = new G4LFission();
    theNeutronLEPFission->SetMinEnergy( minNonHPNeutronEnergy_ );
    theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
    fission->RegisterMe( theNeutronLEPFission );
  }
}
