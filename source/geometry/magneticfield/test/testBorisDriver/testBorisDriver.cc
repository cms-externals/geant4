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
#include "G4MuonMinus.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManagerFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "FTFP_BERT.hh"
#include "G24ActionInitialization.hh"
#include "G24DetectorConstruction.hh"
#include "G24PhysicsList.hh"
#include "G24SteppingVerboseWithDir.hh"
#include <stdlib.h>

#include <iostream>

// Arguments (argv) of main
// 1. minEpsilon    =  requested accuracy of trajectory
// 2. beta          =  beta (velocity / c_light) of test particle
// 3. step_fraction =  fraction used a factor to reduce the 2*pi (estimated) circumferance
//                  ( we could in future change it to an integer 'steps' and calculate step_fraction
//                  = 1 / steps )

int main(int argc, char* argv[])
{
  G4double minEpsilon = 1.0e-4;
  G4double beta = 0.1;  // default
  G4double step_fraction = 1.0;
  G4double mass_muon_c2 = 105.6583755 * CLHEP::MeV;
  G4double mass_muon_kg =
    ((mass_muon_c2) / c_squared) / kg;  // G4MuonMinus::MuonMinusDefinition()->GetPDGMass();

  G4double bval_si = 1.0;  // CLHEP::electron_mass_c2/mass_muon_c2; //tesla

  G4VSteppingVerbose::SetInstance(new G24SteppingVerboseWithDir);
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  // char *ptr;
  if (argc > 1)
  {
    minEpsilon = strtod(argv[1], nullptr);
    G4cout << " - argv[1] : minEpsilon    = " << minEpsilon << G4endl;
  }
  else
  {
    G4cout << " - Default : minEpsilon    = " << minEpsilon << " (change it as argv[1] )" << G4endl;
  }
  if (argc > 2)
  {
    beta = strtod(argv[2], nullptr);
    G4cout << " - argv[2] : beta          = " << beta << G4endl;
  }
  else
  {
    G4cout << " - Default : beta          = " << beta << "       (change it as argv[2] )" << G4endl;
  }
  if (argc > 3)
  {
    step_fraction = strtod(argv[3], nullptr);
    G4cout << " - argv[3] : step_fraction = " << step_fraction << G4endl;
  }
  else
  {
    G4cout << " - Default : step_fraction = " << step_fraction << " change it as argv[3] "
           << G4endl;
  }

  G24DetectorConstruction* detectorConstructor = new G24DetectorConstruction;
  detectorConstructor->minEpsilon = minEpsilon;
  detectorConstructor->step_fraction = step_fraction;

  // constexpr G4double mass_e_si = (electron_mass_c2/c_squared)/CLHEP::kg;;  //9.108e-31;
  G4cout << "Mass of Muon " << mass_muon_kg << " kg " << G4endl;

  G4double velocity_si = beta * CLHEP::c_light / CLHEP::m * CLHEP::second;
  constexpr G4double charge_si = CLHEP::e_SI;
  G4double gamma = 1 / std::sqrt(1 - sqr(beta));

  G4double radius_si = (gamma * mass_muon_kg * velocity_si) / (charge_si * bval_si);
  G4cout << "Beta: " << beta << G4endl;
  G4cout << "Radius: " << radius_si << " m " << G4endl;
  G4cout << "step_fraction " << step_fraction << G4endl;
  const G4double circum_elec = 2 * CLHEP::pi * radius_si * m;

  runManager->SetUserInitialization(detectorConstructor);
  detectorConstructor->circumference = circum_elec;
  detectorConstructor->bz_si = bval_si;
  // detectorConstructor->energy = 0.5*mass_e*velocity*velocity;

  // runManager->SetUserInitialization(new G24PhysicsList);
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  auto actionini = new G24ActionInitialization;
  actionini->fbeta = beta;
  actionini->mass_c2 = mass_muon_c2;
  runManager->SetUserInitialization(actionini);

  // auto visManager = new G4VisExecutive;
  // visManager->Initialize();

  runManager->Initialize();
  G4UImanager* UI = G4UImanager::GetUIpointer();

  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // start a run
  int numberOfEvent = 1;
  runManager->BeamOn(numberOfEvent);

  // job termination
  // delete visManager;
  delete runManager;

  return 0;
}
