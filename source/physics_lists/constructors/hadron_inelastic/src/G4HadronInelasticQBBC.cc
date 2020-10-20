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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronInelasticQBBC
//
// Author: 2 October 2009 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronInelasticQBBC.hh"

#include "G4SystemOfUnits.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronicInteraction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

#include "G4ParticleInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionInelastic.hh"

#include "G4QGSBuilder.hh"
#include "G4FTFBuilder.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4NeutronRadCapture.hh"

#include "G4PreCompoundModel.hh"
#include "G4HadronicInteractionRegistry.hh"

#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronInelasticQBBC);

G4HadronInelasticQBBC::G4HadronInelasticQBBC(G4int ver) 
  : G4VHadronPhysics("hInelasticQBBC"),verbose(ver)
{}

G4HadronInelasticQBBC::G4HadronInelasticQBBC(const G4String&, G4int ver, 
    G4bool, G4bool,G4bool, G4bool, G4bool) : G4HadronInelasticQBBC(ver)
{}

G4HadronInelasticQBBC::~G4HadronInelasticQBBC()
{}

void G4HadronInelasticQBBC::ConstructProcess()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // configure models
  const G4double eminFtf = param->GetMinEnergyTransitionFTF_Cascade();
  const G4double emaxBert = param->GetMaxEnergyTransitionFTF_Cascade();
  const G4double emaxBertPions = 12*CLHEP::GeV;
  const G4double emax = param->GetMaxEnergy();
  if(verbose > 1) {
    G4cout << "### HadronInelasticQBBC Construct Process Emin(FTFP)= " 
           << eminFtf/CLHEP::GeV << " GeV; Emax(BERT)= " << emaxBert/CLHEP::GeV
           << " GeV; Emax(BERTpions)= " << emaxBertPions/CLHEP::GeV 
           << " GeV; Emax= " << emax/CLHEP::GeV << " GeV" << G4endl;
  }

  // PreCompound and Evaporation models are instantiated here
  G4PreCompoundModel* thePreCompound = nullptr;
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  thePreCompound = static_cast<G4PreCompoundModel*>(p);
  if(!thePreCompound) { thePreCompound = new G4PreCompoundModel(); }
 
  //G4HadronicInteraction* theQGSP = 
  //  BuildModel(new G4QGSBuilder("QGSP",thePreCompound,true,false),12.5*GeV,emax);
  G4HadronicInteraction* theFTFP = 
    BuildModel(new G4FTFBuilder("FTFP",thePreCompound),eminFtf,emax);

  G4CascadeInterface* casc = new G4CascadeInterface();
  casc->usePreCompoundDeexcitation();
  G4HadronicInteraction* theBERT = NewModel(casc,1.0*GeV,emaxBert);

  casc = new G4CascadeInterface();
  casc->usePreCompoundDeexcitation();
  G4HadronicInteraction* theBERT1 = NewModel(casc,0.0*GeV,emaxBertPions);

  //G4cout << "G4HadronInelasticQBBC::ConstructProcess new Binary"<< G4endl;
  G4BinaryCascade* bic = new G4BinaryCascade(thePreCompound);
  G4HadronicInteraction* theBIC = NewModel(bic,0.0,1.5*GeV);

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronicProcess* hp = 
    new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  hp->AddDataSet(new G4ParticleInelasticXS(particle));
  hp->RegisterMe(theFTFP);
  hp->RegisterMe(theBERT);
  hp->RegisterMe(theBIC);
  ph->RegisterProcess(hp, particle);
  if( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );

  // n
  particle = G4Neutron::Neutron();
  hp = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  hp->AddDataSet(new G4NeutronInelasticXS());
  hp->RegisterMe(theFTFP);
  hp->RegisterMe(theBERT);
  hp->RegisterMe(theBIC);
  ph->RegisterProcess(hp, particle);
  if( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
       
  hp = new G4HadronCaptureProcess("nCapture");
  hp->RegisterMe(new G4NeutronRadCapture());
  hp->AddDataSet(new G4NeutronCaptureXS());
  ph->RegisterProcess(hp, particle);

  // pi+
  particle = G4PionPlus::PionPlus();
  hp = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  hp->AddDataSet(new G4BGGPionInelasticXS(particle));
  hp->RegisterMe(theFTFP);
  hp->RegisterMe(theBERT1);
  ph->RegisterProcess(hp, particle);
  if( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // pi-
  particle = G4PionMinus::PionMinus();
  hp = new G4HadronInelasticProcess( particle->GetParticleName()+"Inelastic", particle );
  hp->AddDataSet(new G4BGGPionInelasticXS(particle));
  hp->RegisterMe(theFTFP);
  hp->RegisterMe(theBERT1);
  ph->RegisterProcess(hp, particle);
  if( useFactorXS ) hp->MultiplyCrossSectionBy( param->XSFactorPionInelastic() );

  // kaons
  G4HadronicBuilder::BuildKaonsFTFP_BERT();

  // high energy particles
  if( emax > param->EnergyThresholdForHeavyHadrons() ) {

    // pbar, nbar, anti light ions
    G4HadronicBuilder::BuildAntiLightIonsFTFP();

    // hyperons
    G4HadronicBuilder::BuildHyperonsFTFP_BERT();

    // b-, c- baryons and mesons
    if( param->EnableBCParticles() ) {
      G4HadronicBuilder::BuildBCHadronsFTFP_BERT();
    }
  }
}
