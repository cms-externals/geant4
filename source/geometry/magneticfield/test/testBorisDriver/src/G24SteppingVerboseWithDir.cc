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
// Author : J. Apostolakis,  (EP/SFT CERN)  2019-21
//
//  Refines the existing G4SteppingVerbose, adding:
//   - direction of particle
//   - safety (at step's start point)

#include "G24SteppingVerboseWithDir.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

static const unsigned int wp = 10;  // width of position number
static const int precE = 8;  // precision digits of Energy

using std::setprecision;
using std::setw;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G24SteppingVerboseWithDir::G24SteppingVerboseWithDir() : G4SteppingVerbose() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G24SteppingVerboseWithDir::~G24SteppingVerboseWithDir() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G24SteppingVerboseWithDir::Banner()
{
  G4cout << G4endl;
  G4cout << setw(5) << "Step#" << " " << setw(wp) << "X (mm)" << " " << setw(wp) << "Y (mm)" << " "
         << setw(wp) << "Z (mm)" << " " << setw(11) << "Ek (MeV)" << " " << setw(11)
         << "E0-Ek (MeV)" << " " << setw(8) << "dE (MeV)" << " " << setw(7) << "Step(mm)" << " "
         << setw(12) << "Track(mm)" << " " << setw(7) << "Safety(mm)" << " | " << setw(11)
         << "Dir.x " << " " << setw(11) << "Dir.y " << " " << setw(11) << " Dir.z " << " | "
         << setw(11) << "1- Dir.y " << " | " << setw(15) << "Material" << setw(14) << "Process"
         << setw(15) << "NextVolu" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G24SteppingVerboseWithDir::StepInfo()
{
  CopyState();

  G4int prec = G4cout.precision(12);

  if (verboseLevel >= 1)
  {
    if (verboseLevel >= 4) VerboseTrack();
    if (verboseLevel >= 3)
    {
      Banner();
    }

    G4cout << setw(5) << fTrack->GetCurrentStepNumber() << " ";
    G4cout << std::fixed << setprecision(6) << setw(wp) << fTrack->GetPosition().x() / CLHEP::mm
           << " " << setw(wp) << fTrack->GetPosition().y() / CLHEP::mm
           << " "
           //   << std::defaultfloat // Small cheat --- it's zero here!
           << setw(wp) << fTrack->GetPosition().z() / CLHEP::mm << " ";
    G4cout << std::defaultfloat << setprecision(precE) << setw(11)
           << fTrack->GetKineticEnergy() / CLHEP::MeV << " " << setprecision(precE / 2) << setw(11)
           << (fTrack->GetKineticEnergy() - fInitialTrackKinE) / CLHEP::MeV << " "
           << setprecision(3) << setw(8) << fStep->GetTotalEnergyDeposit() / CLHEP::MeV << " ";
    G4cout << setprecision(5) << setw(7) << fStep->GetStepLength() / CLHEP::mm << " ";
    G4cout.precision(6);
    G4cout << setw(12) << fTrack->GetTrackLength() / CLHEP::mm << " ";

    G4cout.precision(5);
    G4StepPoint* preStepPt = fTrack->GetStep()->GetPreStepPoint();
    if (preStepPt != nullptr)
    {
      G4double safety = preStepPt->GetSafety();
      G4cout << setw(8) << ((safety > 1.0e-9) ? safety : 0.0) << " ";
    }
    else
    {
      G4cout << setw(8) << "-0.0_(N/A)" << " ";
    }

    G4cout << "| ";
    G4cout << setw(11) << fTrack->GetMomentumDirection().x() << " " << setw(11)
           << fTrack->GetMomentumDirection().y() << " " << setw(11)
           << fTrack->GetMomentumDirection().z() << " | " << setw(11)
           << 1.0 - fTrack->GetMomentumDirection().y() << " | ";

    G4cout << setw(15) << fTrack->GetMaterial()->GetName() << " ";
    G4VProcess const* pds = fStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4cout << setw(14) << ((pds != nullptr) ? pds->GetProcessName() : G4String("User Limit"))
           << " ";
    // if( fStepStatus != fWorldBoundary){
    if (fTrack->GetNextVolume() != nullptr)
    {
      G4cout << setw(10) << fTrack->GetNextVolume()->GetName();
      G4cout << setw(4) << fTrack->GetNextVolume()->GetCopyNo();
    }
    else
    {
      G4cout << setw(10) << "OutOfWorld";
    }

    G4cout << " ";

    G4cout << G4endl;

    if (verboseLevel == 2)
    {
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt + fN2ndariesAlongStepDoIt + fN2ndariesPostStepDoIt;
      if (tN2ndariesTot > 0)
      {
        G4cout << "    :----- List of 2ndaries - "
               << "#SpawnInStep=" << setw(3) << tN2ndariesTot << "(Rest=" << setw(2)
               << fN2ndariesAtRestDoIt << ",Along=" << setw(2) << fN2ndariesAlongStepDoIt
               << ",Post=" << setw(2) << fN2ndariesPostStepDoIt << "), "
               << "#SpawnTotal=" << setw(3) << (*fSecondary).size() << " ---------------" << G4endl;

        for (size_t lp1 = (*fSecondary).size() - tN2ndariesTot; lp1 < (*fSecondary).size(); lp1++)
        {
          G4cout << "    : " << setw(8) << (*fSecondary)[lp1]->GetPosition().x() / CLHEP::mm << " "
                 << setw(8) << (*fSecondary)[lp1]->GetPosition().y() / CLHEP::mm << " " << setw(8)
                 << (*fSecondary)[lp1]->GetPosition().z() / CLHEP::mm << " " << setw(8)
                 << (*fSecondary)[lp1]->GetKineticEnergy() / CLHEP::MeV << " " << setw(10)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName() << " ";
          G4cout << G4endl;
        }

        G4cout << "    :-----------------------------"
               << "----------------------------------"
               << "-- EndOf2ndaries Info ---------------" << G4endl;
      }
    }
  }
  G4cout.precision(prec);
  // G4cout<< "exit G24SteppingVerboseWithDir::StepInfo   " <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G24SteppingVerboseWithDir::TrackingStarted()
{
  CopyState();
  fInitialTrackKinE = fTrack->GetKineticEnergy();
  G4int prec = G4cout.precision(6);
  if (verboseLevel > 0)
  {
    Banner();

    G4cout << setw(5) << fTrack->GetCurrentStepNumber() << " " << setw(wp)
           << fTrack->GetPosition().x() / CLHEP::mm << " " << setw(wp)
           << fTrack->GetPosition().y() / CLHEP::mm << " " << setw(wp)
           << fTrack->GetPosition().z() / CLHEP::mm << " "  // << " E_0 = "
           << setprecision(precE) << setw(11) << fTrack->GetKineticEnergy() / CLHEP::MeV << " "
           << setw(11) << (fTrack->GetKineticEnergy() - fInitialTrackKinE) / CLHEP::MeV << " ";
    // << " ( = " << setw(8) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy") << " ) "
    G4cout.precision(3);
    G4cout << setw(8) << fStep->GetTotalEnergyDeposit() / CLHEP::MeV << " " << setw(7)
           << fStep->GetStepLength() / CLHEP::mm << " ";
    G4cout.precision(6);
    G4cout << setw(12) << fTrack->GetTrackLength() / CLHEP::mm << " ";

    G4Step const* step = fTrack->GetStep();
    G4StepPoint const* preStepPt = nullptr;
    if (step != nullptr) preStepPt = fTrack->GetStep()->GetPreStepPoint();
    if (preStepPt != nullptr)
    {
      G4double safety = preStepPt->GetSafety();
      G4cout << setw(8) << ((safety > 1.0e-9) ? safety : 0.0) << " ";
    }
    else
    {
      G4cout << setw(8) << "N/A" << " ";
    }
    G4cout << "| ";
    G4cout << setw(11) << fTrack->GetMomentumDirection().x() << " " << setw(11)
           << fTrack->GetMomentumDirection().y() << " " << setw(11)
           << fTrack->GetMomentumDirection().z() << " | " << setw(11)
           << 1.0 - fTrack->GetMomentumDirection().y() << " | ";

    // Material and volume name
    if (fTrack->GetVolume())
    {
      G4Material* material = fTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
      if (material != nullptr)
      {
        G4cout << setw(15) << material->GetName() << " ";
      }
      else
      {
        G4cout << setw(15) << "No-Material" << " ";
      }
      G4cout << setw(14) << "initStep" << " ";
      G4cout << setw(10) << fTrack->GetVolume()->GetName();
      G4cout << setw(4) << fTrack->GetVolume()->GetCopyNo() << " (initial Volume)";
    }
    else
    {
      G4cout << setw(10) << "No-Material" << " ";
      G4cout << setw(14) << "initStep" << " ";
      G4cout << setw(10) << "OutOfWorld" << " ";
    }
    // fTrack->GetMaterial() fails at track start !!!!
    G4cout << G4endl;
  }
  G4cout.precision(prec);
}
