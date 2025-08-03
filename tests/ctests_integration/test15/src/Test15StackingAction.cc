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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// StackingAction program
// --------------------------------------------------------------

#include "Test15StackingAction.hh"

#include "Test15DetectorConstruction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4HadronicProcess.hh"
#include "G4AnalysisManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Test15EventAction.hh"

Test15StackingAction::Test15StackingAction(Test15EventAction* eventAction)
  : evtAction(eventAction) {

  fNumber_newtracks = 0;
  fNeutron = 0;
  fProton = 0;
  fDeuteron = 0;
  fOther = 0;
  // messenger defaults
  killGammasFlag  = 0;

  // global geometry navigator
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();

}


Test15StackingAction::~Test15StackingAction() {
  
}


void Test15StackingAction::NewStage() {;}

    
void Test15StackingAction::PrepareNewEvent() {;}




//x2018 void Test15StackingAction::printCrossSection(G4ParticleDefinition * particleType,
//        const G4Element* element)
void Test15StackingAction::printCrossSection(G4ParticleDefinition * ,
       const G4Element* )
{

  /*
     G4ProcessManager* processManager=particleType->GetProcessManager();

     G4ProcessVector* processes= processManager->GetProcessList();
     G4cout << "Number of processes " << processes->entries() << G4endl;
     for ( G4int aProc=0; aProc < processes->size(); aProc++ )
     {
       G4cout << " process " << (*processes)[aProc]->GetProcessName()
               << ", type " << (*processes)[aProc]->GetProcessType()
               << G4endl;
       if (   (*processes)[aProc]->GetProcessType() ==4
           && (*processes)[aProc]->GetProcessName() != "PionMinusAbsorptionAtRest"
           && (*processes)[aProc]->GetProcessName() != "KaonMinusAbsorption"
           && (*processes)[aProc]->GetProcessName() != "AntiProtonAnnihilationAtRest"
          )
       {
          G4HadronicProcess * hproc=static_cast<G4HadronicProcess *>((*processes)[aProc]);
           G4double energy=0.;
          G4double mass=particleType->GetPDGMass();
          for ( G4double mom=1.*MeV, energy=std::sqrt(sqr(mass)+sqr(mom)); energy < 1.*TeV; mom *=1.01,
energy=std::sqrt(sqr(mass)+sqr(mom)))
          {
            const G4DynamicParticle * dynParticle=new G4DynamicParticle(particleType, G4ThreeVector(1,0,0), energy-mass);
              G4double cross_sec = hproc->GetMicroscopicCrossSection(dynParticle, element, 293.);

              G4double xs_lhep = (G4HadronCrossSections::Instance())->
                GetInelasticCrossSection(dynParticle, element);

              G4cout << "mat, part KE, p, X " << element->GetName()
           << "    " << particleType->GetParticleName()
           << "    " << (*processes)[aProc]->GetProcessName()
                   << "    " << (energy-mass)/GeV
           << "    " << mom/GeV
           << "    " << cross_sec *1000./ barn
           << "    " << xs_lhep*1000./barn
           << G4endl;
             delete dynParticle;
          }
       }
     }
  */
 return;
}


G4ClassificationOfNewTrack
Test15StackingAction::ClassifyNewTrack(const G4Track * aTrack) {

  // G4ClassificationOfNewTrack classification_w = fWaiting;
  // return classification_w;
  fNumber_newtracks++;
  if(fNumber_newtracks%20000 == 0) G4cout << " current number of tracks: " << fNumber_newtracks << G4endl
					 << " Neutron: " << fNeutron << " Proton: " << fProton 
					 << " Deuteron: " << fDeuteron << " Other: " << fOther << G4endl;
  //  G4cout << "ClassifyNewTrack----------------" << G4endl;

  G4ParticleDefinition* particleType = aTrack->GetDefinition();

  G4String particleName;
  G4String leadName;
  G4int particle_type = -9;

  if ( particleType == G4Neutron::NeutronDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 1;
    fNeutron++;
    evtAction->AddToNeutronStack();
    particleName = "neutron";
  } else if ( particleType == G4Proton::ProtonDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 2;
    fProton++;
    particleName = "proton";
  } else if ( particleType == G4Deuteron::DeuteronDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 3;
    particleName = "deuteron";
  } else if ( particleType == G4Triton::TritonDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 4;
    particleName = "triton";
  } else if ( particleType == G4Gamma::GammaDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 5;
    particleName = "gamma";
  } else if ( particleType == G4Electron::ElectronDefinition() && aTrack->GetTrackID() != 1 ) {
    particleName = "e-";
    particle_type = 6;
  } else if ( particleType == G4Positron::PositronDefinition() && aTrack->GetTrackID() != 1 ) {
    particle_type = 7;
    particleName = "e+";
  } else if ( particleType->GetAtomicMass() != 208 && particleType->GetAtomicNumber() == 82 && aTrack->GetTrackID() != 1 ) {
    leadName = "lead";
    particle_type = 8;
    //     G4cout << G4endl;
    //     G4cout << " lead produced!!!!!!!!! " << G4endl;
    //     G4cout << G4endl;
    if ( particleType->GetAtomicMass() == 204 && particleType->GetAtomicNumber() == 82 && aTrack->GetTrackID() != 1 ) {
      particleName = "pb204";
    } else if ( particleType->GetAtomicMass() == 206 && particleType->GetAtomicNumber() == 82 && aTrack->GetTrackID() != 1 ) {
      particleName = "pb206";
    } else if ( particleType->GetAtomicMass() == 207 && particleType->GetAtomicNumber() == 82 && aTrack->GetTrackID() != 1 ) {
      particleName = "pb207";
    } else if ( particleType->GetAtomicMass() == 208 && particleType->GetAtomicNumber() == 82 && aTrack->GetTrackID() != 1 ) {
      particleName = "pb208";
    }
  } else if ( particleType->GetAtomicMass() == 2 && particleType->GetAtomicNumber() == 1 && aTrack->GetTrackID() != 1 ) {
    particle_type = 9;
    fDeuteron++;
    particleName = "deuteron";
  } else if ( particleType->GetAtomicMass() == 3 && particleType->GetAtomicNumber() == 1 && aTrack->GetTrackID() != 1 ) {
    particle_type = 10;
    particleName = "triton";
  } else {
    if(particleType->GetAtomicMass()>41) particle_type = 11;
    else particle_type = 12;
    fOther++;
    particleName = "other";
    static G4int counter = 0;
    counter++;
    if(counter < 50) {
      G4cout << " particle type: " << particleType->GetParticleType() << G4endl;
      G4cout << " A: " << particleType->GetAtomicMass() << G4endl;
      G4cout << " Z: " << particleType->GetAtomicNumber() << G4endl;
      if(particleType->GetAtomicMass() < 200 && particleType->GetAtomicMass() > 10) {
	G4cout << " found a problem? " << G4endl;
	//   getchar();
      }
    }
  }
    //    G4double Tproton=aTrack->GetKineticEnergy();
  // G4double KineticEnergy=aTrack->GetKineticEnergy();
  // G4double Momentum=aTrack->GetMomentum().mag();

  // G4int atomicNumber = particleType->GetAtomicNumber();
  // G4int atomicMass = particleType->GetAtomicMass();
  // G4double mass=particleType->GetPDGMass();

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if(particle_type != 6) analysisManager->FillH1(6,particle_type);
  // analysisManager->FillNtupleDColumn(2,0, energy);
  // analysisManager->AddNtupleRow(2);  
  

  //  G4cout << " Got here1 ! " << G4endl;

  // still needed? 27/9/15:

  // if(aTrack->GetTrackID() !=1) {

  //   //  const G4VProcess* proc=aTrack->GetCreatorProcess();
  //   //  G4String ProcName=proc->GetProcessName();
  
  //   Test15AnalysisManager* analyse =  Test15AnalysisManager::getInstance();
  //   analyse->analyseStack(KineticEnergy, particleName, Momentum, leadName, atomicNumber, atomicMass, mass);

  // }

  //  G4cout << " Got here2 ! " << G4endl;

  /*
  if ( particleType == G4Gamma::GammaDefinition() && aTrack->GetTrackID() != 1 ) {
    const G4VProcess* proc=aTrack->GetCreatorProcess();
    G4double Tgamma=aTrack->GetKineticEnergy();
    G4String ProcName=proc->GetProcessName();

    G4cout << " Gamma process: " << ProcName << " and energy in MeV: " << Tgamma/MeV << G4endl;

  }
  */


  //xxx
  /*
  G4ParticleDefinition * aproton = G4Proton::Proton();
  G4ParticleDefinition * aneutron = G4Neutron::Neutron();
  G4ParticleDefinition * apiplus = G4PionPlus::PionPlus();
  G4ParticleDefinition * apiminus = G4PionMinus::PionMinus();

  G4Material* material=0;
  const G4Element* element=0;
  G4VPhysicalVolume* pvolume=0;
  G4LogicalVolume* lvolume=0;
  G4String volumeName="";
  const G4VTouchable* touchable=aTrack->GetTouchable();
  if ( touchable ) pvolume=touchable->GetVolume();
  if ( pvolume ) {
     volumeName=pvolume->GetName();
     lvolume=pvolume->GetLogicalVolume();
  }
  if ( lvolume ) material=lvolume->GetMaterial();

  static std::vector<G4String> elementsDone;
  if ( material )
  {
      for (G4int ielem=0; ielem< material->GetNumberOfElements(); ielem++)
      {
        element=material->GetElement(ielem);

        if ( element )
        {
           G4String name=element->GetName();
           G4cout << " __Element__ = " << name << G4endl;;
           G4bool isDone=false;
           G4cout << " elements done size = " << elementsDone.size() << G4endl;
           for (int i=0; i<elementsDone.size(); i++)
           {
              isDone= isDone || elementsDone[i] == name;
              G4cout << " ... " << i << " " <<  elementsDone[i] << " - " << name << G4endl;
           }
           if ( ! isDone )
           {
       elementsDone.push_back(name);
       printCrossSection( aproton, element);
       printCrossSection( aneutron, element);
       printCrossSection( apiplus, element);
       printCrossSection( apiminus, element);
       printCrossSection( G4KaonPlus::KaonPlus(), element);
       printCrossSection( G4KaonMinus::KaonMinus(), element);
       printCrossSection( G4AntiProton::AntiProton(), element);
           }

         } else { G4cout << "element not found" << G4endl;}
       }
  }
*/
  //xxx

// to go in stacking action to see which process created which process created which particle - eg. 600 keV protons in lead from 2.5GeV protons...
//

  /*
  G4ParticleDefinition* particleType = aTrack->GetDefinition();

  if(aTrack->GetTrackID() == 1) energy = 0.0;
  
  //  if ( particleType == G4Proton::ProtonDefinition() && aTrack->GetTrackID() != 1 ) {
  if ( particleType == G4Neutron::NeutronDefinition() && aTrack->GetTrackID() != 1 ) {
    const G4VProcess* proc=aTrack->GetCreatorProcess();
    //    G4double Tproton=aTrack->GetKineticEnergy();
    G4double Tneutron=aTrack->GetKineticEnergy();
    energy += Tneutron;

    G4String ProcName=proc->GetProcessName();

    //    if ( Tproton < 0.575 and Tproton > 0.570*keV ) {
    //    G4cout << " Energy: " << G4BestUnit(Tneutron,"Energy") << G4endl;
    if ( Tneutron < 0.10 and Tneutron > 0.05 ) {
      //  && ProcName!="Decay"
      //   && ProcName!="eBrem"
      //  ) {
      
      // G4cout << " gamma CreatorProc " << proc << " " <<proc->GetProcessName()<< G4endl;
      G4String volumeName("_Unknown_Vol_");
      G4String materialName("_Unknown_Mat_");
      G4Material* material=0;
      G4VPhysicalVolume* pvolume=0;
      G4LogicalVolume* lvolume=0;
      const G4VTouchable* touchable=aTrack->GetTouchable();
      if ( touchable ) pvolume=touchable->GetVolume();
      if ( pvolume ) {
	volumeName=pvolume->GetName();
	lvolume=pvolume->GetLogicalVolume();
      }
      if ( lvolume ) material=lvolume->GetMaterial();
      if ( material ) materialName= material->GetName();
      G4cout << ProcName << " 1 neutron E(GeV) "<< aTrack->GetTotalEnergy()/GeV <<
	" "
	     << materialName << " "
	     << volumeName << G4endl;
      G4cout << ProcName << " 2 neutron kinetic in keV: " << Tneutron/keV << G4endl;
      G4cout << ProcName << " 3 parent ID: " << aTrack->GetParentID() << " track id: " << aTrack->GetTrackID() << " sum energy: " << energy << G4endl;
      G4cout << ProcName << " 4 neutron position(m): " << aTrack->GetPosition()/m
	     << G4endl;
      G4cout << ProcName << " 5 neutron direction " <<
	aTrack->GetMomentumDirection() << G4endl;
    }
    
  }

  */
  G4ClassificationOfNewTrack classification = fWaiting;

  return classification;

}
