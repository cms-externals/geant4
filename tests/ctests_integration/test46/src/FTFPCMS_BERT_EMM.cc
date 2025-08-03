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
#include "FTFPCMS_BERT_EMM.hh"
#include "CMSEmStandardPhysicsHcal.hh"
#include "CMSHadronPhysicsFTFP_BERT.hh"

#include "G4SystemOfUnits.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"

FTFPCMS_BERT_EMM::FTFPCMS_BERT_EMM(G4int ver)
{

  defaultCutValue = 0.7*mm;  
  SetVerboseLevel(ver);
  
  G4cout << "You are using " << "FTFP_BERT_EMM "
	 << G4endl;

  // EM Physics
  RegisterPhysics( new CMSEmStandardPhysicsHcal(ver));

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver));

  //  RegisterPhysics( new G4HadronPhysicsFTFP_BERT(ver)); 
  RegisterPhysics( new CMSHadronPhysicsFTFP_BERT(ver)); 
  
  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));
}

FTFPCMS_BERT_EMM::~FTFPCMS_BERT_EMM()
{}

void FTFPCMS_BERT_EMM::SetCuts()
{
  SetCutsWithDefault();   
}
