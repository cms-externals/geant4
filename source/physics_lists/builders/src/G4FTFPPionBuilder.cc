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
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#include "G4FTFPPionBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4HadronicParameters.hh"


G4FTFPPionBuilder::
G4FTFPPionBuilder(G4bool quasiElastic) 
{
  theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  theModel = new G4TheoFSGenerator("FTFP");

  G4FTFModel* theStringModel = new G4FTFModel();
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());

  G4GeneratorPrecompoundInterface* theCascade = 
    new G4GeneratorPrecompoundInterface();

  theModel->SetHighEnergyGenerator(theStringModel);
  if (quasiElastic) {
     theModel->SetQuasiElasticChannel(new G4QuasiElasticChannel());
  } 

  theModel->SetTransport(theCascade);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
}

G4FTFPPionBuilder::~G4FTFPPionBuilder() 
{
}

void G4FTFPPionBuilder::
Build(G4HadronInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  if ( aP->GetParticleDefinition() == G4PionPlus::Definition() ) { 
    aP->AddDataSet( new G4BGGPionInelasticXS( G4PionPlus::Definition() ) );
  } else if ( aP->GetParticleDefinition() == G4PionMinus::Definition() ) { 
    aP->AddDataSet( new G4BGGPionInelasticXS( G4PionMinus::Definition() ) );
  }
  aP->RegisterMe(theModel);
}
