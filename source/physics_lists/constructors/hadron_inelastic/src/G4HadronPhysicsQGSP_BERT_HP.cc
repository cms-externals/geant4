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
// ClassName:   G4HadronPhysicsQGSP_BERT_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 15.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 20.06.2006 G.Folger: Bertini applies to Kaons, i.e. use SetMinEnergy
//            instead of SetMinPionEnergy
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP
//
//----------------------------------------------------------------------------
//
#include "G4HadronPhysicsQGSP_BERT_HP.hh"

#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4HadronicParameters.hh"
#include "G4LFission.hh"
#include "G4NeutronBuilder.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronHPCaptureXS.hh"
#include "G4NeutronPHPBuilder.hh"
#include "G4NeutronRadCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysListUtil.hh"
#include "G4ProcessManager.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <iomanip>

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BERT_HP);

G4HadronPhysicsQGSP_BERT_HP::G4HadronPhysicsQGSP_BERT_HP(G4int verb)
  : G4HadronPhysicsQGSP_BERT_HP("hInelastic QGSP_BERT_HP")
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BERT_HP::G4HadronPhysicsQGSP_BERT_HP(const G4String& name, G4bool)
  : G4HadronPhysicsQGSP_BERT(name)
{
  minBERT_neutron = 19.9 * CLHEP::MeV;
  G4HadronicParameters::Instance();
}

void G4HadronPhysicsQGSP_BERT_HP::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  G4NeutronBuilder neu(true);  // Fission on

  G4QGSPNeutronBuilder qgs(QuasiElasticQGS);
  qgs.SetMinEnergy(minQGSP_neutron);
  neu.RegisterMe(&qgs);

  G4FTFPNeutronBuilder ftf(QuasiElasticFTF);
  ftf.SetMinEnergy(minFTFP_neutron);
  ftf.SetMaxEnergy(maxFTFP_neutron);
  neu.RegisterMe(&ftf);

  G4BertiniNeutronBuilder bert;
  bert.SetMinEnergy(minBERT_neutron);
  bert.SetMaxEnergy(maxBERT_neutron);
  neu.RegisterMe(&bert);

  G4NeutronPHPBuilder hp;
  neu.RegisterMe(&hp);

  // build all models
  neu.Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(neutron);
  if(nullptr != inel && useFactorXS)
  {
    inel->MultiplyCrossSectionBy(param->XSFactorNucleonInelastic());
  }

  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
  if (nullptr != capture)
  {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture;
    theNeutronRadCapture->SetMinEnergy(minBERT_neutron);
    capture->RegisterMe(theNeutronRadCapture);
  }

  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess( neutron );
  if (nullptr != fission)
  {
    G4LFission* theNeutronLEPFission = new G4LFission;
    theNeutronLEPFission->SetMinEnergy(minBERT_neutron);
    theNeutronLEPFission->SetMaxEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());
    fission->RegisterMe(theNeutronLEPFission);
  }
}
