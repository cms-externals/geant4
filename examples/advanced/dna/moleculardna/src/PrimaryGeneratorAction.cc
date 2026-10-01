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
/// \file PrimaryGeneratorAction.cc
/// \brief Primary Generator for Molecular DNA simulation

#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "PSFReader.hh"    //for GetNextPrimary

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  fpGPS = new G4GeneralParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fpGPS;
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  if (PSFReader::Instance().GetPhaseSpacePattern().empty()) {
    // If no PSF is given, use the macro-defined source
    fpGPS->GeneratePrimaryVertex(anEvent);
  }
  // For PSF
  else
  {
    // If the PSF has not been loaded yet
    // it should have been either in the main or in the BeginOfRunAction
    if (!PSFReader::Instance().IsLoaded() || !PSFReader::Instance().GetNumberOfEntries()) {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "PSF is not yet loaded !";
      G4Exception("PrimaryGeneratorAction" "GeneratePrimaries", "PSF is not yet loaded",
                  FatalException, exceptionDescription);
    } else // PSF has been loaded
    {
      anEvent->AddPrimaryVertex(PSFReader::Instance().GetNextPrimary());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
