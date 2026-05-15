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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//
// History
//   January 23rd, 2023 : first implementation
//	 March 22nd, 2026   : Modified by Ilker Parmaksiz , Opticks photon collection
//
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
#include "SteppingAction.hh"
// Cats Headers
#include "ConfigurationManager.hh"
#include "Event.hh"
#include "PhotonHit.hh"
#include "PhotonSD.hh"
//G4 Headers
#include "G4AutoLock.hh"
#include "G4Cerenkov.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Scintillation.hh"
#include "G4Step.hh"

#include <G4String.hh>
#include <G4Types.hh>
#include <G4VHit.hh>
#include <G4VHitsCollection.hh>
#ifdef WITH_G4OPTICKS
#include "SEvt.hh"
#include "U4.hh"
#include "G4CXOpticks.hh"
#endif

SteppingAction::SteppingAction() {}
SteppingAction::~SteppingAction() {}
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  auto aTrack = aStep->GetTrack();
  (void) aTrack;	
#ifdef WITH_G4OPTICKS
  auto configMng = ConfigurationManager::getInstance();

  if (configMng->isEnable_opticks())
  {
    G4int fNumPhotons = 0;  // number of scintillation photons this step
    G4SteppingManager* fpSteppingManager =
      G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
    G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();
    if (stepStatus != fAtRestDoItProc)
    {
      G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
      size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
      for (size_t i3 = 0; i3 < MAXofPostStepLoops; i3++)
      {
        if ((*procPost)[i3]->GetProcessName() == "Cerenkov")
        {
          const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
          G4double charge = aParticle->GetDefinition()->GetPDGCharge();
          const G4Material* aMaterial = aTrack->GetMaterial();
          G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
          G4MaterialPropertyVector* Rindex = MPT->GetProperty(kRINDEX);
          G4Cerenkov* proc = (G4Cerenkov*)(*procPost)[i3];
          fNumPhotons = proc->GetNumPhotons();
          Photoncounter += fNumPhotons;
          if (fNumPhotons > 0)
          {
            G4double Pmin = Rindex->Energy(0);
            G4double Pmax = Rindex->GetMaxEnergy();
            G4double nMax = Rindex->GetMaxValue();
            G4double beta1 = aStep->GetPreStepPoint()->GetBeta();
            G4double beta2 = aStep->GetPostStepPoint()->GetBeta();
            G4double beta = (beta1 + beta2) * 0.5;
            G4double BetaInverse = 1. / beta;
            G4double maxCos = BetaInverse / nMax;
            G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
            G4double MeanNumberOfPhotons1 =
              proc->GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
            G4double MeanNumberOfPhotons2 =
              proc->GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);
            U4::CollectGenstep_G4Cerenkov_modified(aTrack, aStep, fNumPhotons, BetaInverse, Pmin,
                                                   Pmax, maxCos, maxSin2, MeanNumberOfPhotons1,
                                                   MeanNumberOfPhotons2);

            GenStepcounter++;
          }
        }
        if ((*procPost)[i3]->GetProcessName() == "Scintillation")
        {
          G4Scintillation* proc1 = (G4Scintillation*)(*procPost)[i3];
          fNumPhotons = proc1->GetNumPhotons();
          Photoncounter += fNumPhotons;
          G4double timeconst = 0.0;
          if (fNumPhotons > 0)
          {
            const G4Material* aMaterial = aTrack->GetMaterial();
            G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
            timeconst = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
            U4::CollectGenstep_DsG4Scintillation_r4695(aTrack, aStep, fNumPhotons, 1, timeconst);
            GenStepcounter++;
          }
        }
      }


      if (Photoncounter > SEventConfig::MaxPhoton())
      {

        G4int num_PhotonCollected = SEvt::GetNumPhotonCollected(0);
        G4int num_hits = 0;
        if (ConfigurationManager::getInstance()->isEnable_verbose())
        {
          G4int inum_genstep = SEvt::GetNumGenstepFromGenstep(0);
          G4int inum_photon = SEvt::GetNumPhotonFromGenstep(0);
          G4int num_PhotonGenstepMax = SEvt::GetNumPhotonGenstepMax(0);
          G4cout << "-------------------------------------------------------------------"<< G4endl;
          G4cout << "SteppingAction: GenStepcounter:           " << GenStepcounter << G4endl;
          G4cout << "SteppingAction: PhotonCounter:            " << Photoncounter << G4endl;
          G4cout << "SteppingAction: GetNumPhotonFromGenstep:  " << inum_photon << G4endl;
          G4cout << "SteppingAction: GetNumGenstepFromGenstep: " << inum_genstep << G4endl;
          G4cout << "SteppingAction: GetNumPhotonCollected:    " << num_PhotonCollected<< G4endl;
          G4cout << "SteppingAction: GetNumPhotonGenstepMax:   " << num_PhotonGenstepMax<< G4endl;
        }

        G4RunManager* rm = G4RunManager::GetRunManager();
        const G4Event* event = rm->GetCurrentEvent();
        G4int eventid = event->GetEventID();
        if (num_PhotonCollected)
        {

          G4CXOpticks::Get()->simulate(eventid, false);
          cudaDeviceSynchronize();
          num_hits = SEvt::GetNumHit(0);

          if (num_hits > 0)
          {
            G4HCtable* hctable = G4SDManager::GetSDMpointer()->GetHCtable();
            for (G4int i = 0; i < hctable->entries(); ++i)
            {
              G4String sdn = hctable->GetSDname(i);

              size_t found = sdn.find("PhotonDetector");
              if (found != G4String::npos)
              {
                PhotonSD* aSD = (PhotonSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdn);
                aSD->AddOpticksHits();
              }
            }
          }
        }
        if (ConfigurationManager::getInstance()->isEnable_verbose())
        {
          G4cout << "SteppingAction: NumHits:   "                << num_hits<< G4endl;
          G4cout << "-------------------------------------------------------------------"<<G4endl;
        }
        SteppingAction::ResetPhotoncounter();
        SteppingAction::ResetGenStepcounter();
        if (num_PhotonCollected>0) G4CXOpticks::Get()->reset(eventid);
      }
    }
  }
#endif

}
