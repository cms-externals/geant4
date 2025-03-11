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
#include "TstPrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"

#include "TstReader.hh"
#include "TstTarget.hh"

#include "G4RunManager.hh"
// --> obsolete --> #include "G4HadronCrossSections.hh"
// #include "G4CrossSectionDataStore.hh"
#include "G4VCrossSectionDataSet.hh"
// --> obsolete --> #include "G4HadronInelasticDataSet.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
//
#include "G4HadronicProcessStore.hh"

#include "G4RunManagerKernel.hh"

#ifdef G4_USE_FLUKA
// interface to FLUKA.CERN, if specified
#include "FLUKAInelasticScatteringXS.hh"
#endif

TstPrimaryGeneratorAction::~TstPrimaryGeneratorAction()
{

   if ( fPartGun ) delete fPartGun;

}

void TstPrimaryGeneratorAction::InitBeam( TstReader* pset )
{

   if ( fIsInit ) return;
   
   fIsInit=true;
   
   fConfigPtr=pset;

   G4int NPart = 1; 
   fPartGun = new G4ParticleGun( NPart );
      
   G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* partDef = partTable->FindParticle(pset->GetBeamParticle()); 
   // G4ParticleDefinition* partDef = (G4ParticleTable::GetParticleTable())->FindParticle(pset->GetBeamParticle());
   G4double partMass = partDef->GetPDGMass();
   G4double partMom = pset->GetBeamMomentum();
   G4double partEnergy = std::sqrt( partMom*partMom + partMass*partMass ); // this is total energy
   
   fPartGun->SetParticleDefinition( partDef );
   fPartGun->SetParticleEnergy( partEnergy-partMass ); // this has to be kinetic energy
   fPartGun->SetParticlePosition( pset->GetPosition() ); // it's already in mm (from TstReader)

   // in principle, this should be calculated on the event-by-event basis
   // in case beam direction varies
   
   const TstTarget* target = dynamic_cast<const TstTarget*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   
   const G4Element* elm = target->GetCurrentMaterial()->GetElement(0);
   G4int A = (G4int)(elm->GetN()+0.5);
   G4int Z = (G4int)(elm->GetZ()+0.5);
   G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);
   
   G4DynamicParticle dParticle( partDef, pset->GetDirection(), partEnergy-partMass); // 3rd arg has to be EKin  
   
   fLabV.setX(0.);
   fLabV.setY(0.);
   fLabV.setZ( std::sqrt( partEnergy*(partEnergy+2.0*partMass) )/GeV );
   fLabV.setT( (partEnergy+partMass+amass)/GeV );
   fLabP.setX(0.);
   fLabP.setY(0.);
   fLabP.setZ( std::sqrt(partEnergy*(partEnergy+2.0*partMass))/GeV );
   fLabP.setT( (partEnergy+partMass+G4Proton::Proton()->GetPDGMass())/GeV );

   /* this basically defaults to the Gheisha XSec while there may be more modern versions
   fXSecOnTarget = (G4HadronCrossSections::Instance())->GetInelasticCrossSection( &dParticle, Z, A );
   */

/*
   // Check if PhysList is there
   // ... BUT it also needs to be initialized (BuildPhysicsTable, etc.),
   // otherwise the store returns ZERO;
   // however, so far can't see any clear way how to check for that
   // even if there's a flag fIsPhysicsTableBuilt in G4VUserPhysicsList... 
   //
   if ( G4RunManagerKernel::GetRunManagerKernel()->GetPhysicsList() != nullptr )
   {
      std::cout << " PhysicsList is there !" << std::endl;
      // If so, extract XS from the store
      //
      G4HadronicProcessStore* const store = G4HadronicProcessStore::Instance();      
      fXSecOnTarget = store->GetInelasticCrossSectionPerAtom(partDef, 
                                                             dParticle.GetKineticEnergy(), 
						             elm, 
							     target->GetCurrentMaterial());
      
      return;
   }
*/      

   // This comes to play if PhysList is NOT there and/or 
   // if the physics table isn't properly built (yet) 
   //
   G4VCrossSectionDataSet* cs = 0;

   if ( pset->GetPhysics().find("FLUKA") == std::string::npos )
   {
      if ( partDef == G4Proton::Definition() || partDef == G4Neutron::Definition() ) 
      {
         cs = new G4BGGNucleonInelasticXS( partDef );
      } 
      else if ( partDef == G4PionPlus::Definition() || partDef == G4PionMinus::Definition() ) 
      {
         cs = new G4BGGPionInelasticXS( partDef );
      } 
      else if ( partDef->GetBaryonNumber() > 1 ) 
      {  // Ions
         cs = new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc );
      } 
      else if ( partDef == G4AntiProton::Definition() ||
	       partDef == G4AntiNeutron::Definition() ||
	       partDef == G4AntiDeuteron::Definition() ||
	       partDef == G4AntiTriton::Definition() ||
	       partDef == G4AntiHe3::Definition() ||
	       partDef == G4AntiAlpha::Definition() ) 
      {
         cs = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS ); 
      } 
      else 
      {
         cs = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
      }
   } // end check if no-FLUKA
#ifdef G4_USE_FLUKA
   else 
   {
      cs = new FLUKAInelasticScatteringXS();
   }
#endif      

   if ( cs ) 
   {
      cs->BuildPhysicsTable(*partDef);
      fXSecOnTarget = cs->GetCrossSection( &dParticle, elm );
   }
      
/* move up
   fLabV.setX(0.);
   fLabV.setY(0.);
   fLabV.setZ( std::sqrt( partEnergy*(partEnergy+2.0*partMass) )/GeV );
   fLabV.setT( (partEnergy+partMass+amass)/GeV );
   fLabP.setX(0.);
   fLabP.setY(0.);
   fLabP.setZ( std::sqrt(partEnergy*(partEnergy+2.0*partMass))/GeV );
   fLabP.setT( (partEnergy+partMass+G4Proton::Proton()->GetPDGMass())/GeV );
*/   
   return;

}

void TstPrimaryGeneratorAction::GeneratePrimaries( G4Event* evt ) 
{

   // assert( fConfigPtr );

   // G4ThreeVector dir( fConfigPtr->GetDirection() );
   fPartGun->SetParticleMomentumDirection( fConfigPtr->GetDirection()  );
  
   fPartGun->GeneratePrimaryVertex( evt );
     
   return;

}


