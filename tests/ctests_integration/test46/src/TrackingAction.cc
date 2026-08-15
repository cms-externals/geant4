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
/// \file electromagnetic/TestEm9/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

//---------------------------------------------------------------------------
//
// ClassName:   TrackingAction
//
// Description: Implementation file for MC truth.
//
// Author:      V.Ivanchenko 17/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TrackingAction.hh"

#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackingAction::TrackingAction() : G4UserTrackingAction(), fHisto(HistoManager::GetPointer()) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackingAction::~TrackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  if (0 == aTrack->GetParentID())
  {
    fHisto->SetBeamEnergy(aTrack->GetKineticEnergy());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  if (0 == aTrack->GetParentID())
  {
    const G4Step* step = aTrack->GetStep();
    const G4StepPoint* poststep = step->GetPostStepPoint();
    G4int n_neut = 0;
    G4int n_char = 0;
    G4double z = poststep->GetPosition().z();
    // G4cout << "### end of track " <<  poststep->GetProcessDefinedStep()->GetProcessSubType()
    //	   << G4endl;
    if (poststep->GetProcessDefinedStep()->GetProcessSubType() == 121)
    {
      const std::vector<const G4Track*>* tracks = step->GetSecondaryInCurrentStep();
      G4int nn = tracks->size();
      const G4double elim = 100 * MeV;
      for (G4int i = 0; i < nn; ++i)
      {
        /*
              G4cout << i << ". E(GeV)= " << (*tracks)[i]->GetKineticEnergy()/GeV
               << "  id= " << (*tracks)[i]->GetDefinition()->GetPDGEncoding()
               << " " << (*tracks)[i]->GetDefinition()->GetParticleName() << G4endl;
        */
        if ((*tracks)[i]->GetKineticEnergy() > elim)
        {
          G4int id = (*tracks)[i]->GetDefinition()->GetPDGEncoding();
          if (id == 22)
          {
            ++n_neut;
          }
          else if (id == 111)
          {
            n_neut += 2;
          }
          else if (0.0 != (*tracks)[i]->GetDefinition()->GetPDGCharge())
          {
            ++n_char;
          }
        }
      }
    }
    // G4cout << "  Nchar= " << n_char << " n_neut= " << n_neut << G4endl;
    fHisto->SetFirstInteraction(z, n_char, n_neut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
