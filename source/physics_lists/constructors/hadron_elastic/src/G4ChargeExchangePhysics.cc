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
// ClassName:   G4ChargeExchangePhysics
//
// Author: 19 November 2008 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4ChargeExchangePhysics.hh"

#include "G4ChargeExchangeXS.hh"
#include "G4ChargeExchange.hh"

#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronInelasticProcess.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4ChargeExchangePhysics);

G4ChargeExchangePhysics::G4ChargeExchangePhysics(G4int ver)
  : G4VPhysicsConstructor("chargeExchange")
{
  // because it is an addition, the type of this constructor is 0
  G4HadronicParameters::Instance()->SetVerboseLevel(ver);
  if(ver > 1) G4cout << "### ChargeExchangePhysics" << G4endl;
}

void G4ChargeExchangePhysics::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

void G4ChargeExchangePhysics::ConstructProcess()
{
  auto xs = new G4ChargeExchangeXS();
  auto model = new G4ChargeExchange(xs);

  if (G4HadronicParameters::Instance()->GetVerboseLevel() > 1) {
    G4cout << "### ChargeExchangePhysics Construct Processes with the model <" 
	   << model->GetModelName() << "> and x-section <" 
	   << xs->GetName() << ">" << G4endl;
  }

  // pi-
  G4ParticleDefinition* part = G4PionMinus::PionMinus(); 
  auto proc =
    new G4HadronInelasticProcess(part->GetParticleName()+"ChargeEx", part);
  proc->AddDataSet( xs );
  proc->RegisterMe( model );
  G4ProcessManager* pman = part->GetProcessManager();
  pman->AddDiscreteProcess(proc);

  // pi+
  part = G4PionPlus::PionPlus(); 
  proc = new G4HadronInelasticProcess(part->GetParticleName()+"ChargeEx", part);
  proc->AddDataSet( xs );
  proc->RegisterMe( model );
  pman = part->GetProcessManager();
  pman->AddDiscreteProcess(proc);

  // kaon-
  part = G4KaonMinus::KaonMinus(); 
  proc = new G4HadronInelasticProcess(part->GetParticleName()+"ChargeEx", part);
  proc->AddDataSet( xs );
  proc->RegisterMe( model );
  pman = part->GetProcessManager();
  pman->AddDiscreteProcess(proc);


  // kaon+
  part = G4KaonPlus::KaonPlus(); 
  proc = new G4HadronInelasticProcess(part->GetParticleName()+"ChargeEx", part);
  proc->AddDataSet( xs );
  proc->RegisterMe( model );
  pman = part->GetProcessManager();
  pman->AddDiscreteProcess(proc);

  // KL
  part = G4KaonZeroLong::KaonZeroLong();
  proc = new G4HadronInelasticProcess(part->GetParticleName()+"ChargeEx", part);
  proc->AddDataSet( xs );
  proc->RegisterMe( model );
  pman = part->GetProcessManager();
  pman->AddDiscreteProcess(proc);
}
