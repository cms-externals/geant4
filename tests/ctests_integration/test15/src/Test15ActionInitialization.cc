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
/// \file biasing/Test15/src/Test15ActionInitialization.cc
/// \brief Implementation of the Test15ActionInitialization class
//
//
//

#include "Test15ActionInitialization.hh"
#include "Test15PrimaryGeneratorAction.hh"
#include "Test15RunAction.hh"
#include "Test15SteppingAction.hh"
#include "Test15StackingAction.hh"
#include "Test15DetectorConstruction.hh"
#include "Test15EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


Test15ActionInitialization::Test15ActionInitialization()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Test15ActionInitialization::~Test15ActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Test15ActionInitialization::BuildForMaster() const
{
  SetUserAction(new Test15RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Test15ActionInitialization::Build() const
{

  // set user action classes
  Test15PrimaryGeneratorAction* Test15Generator = new Test15PrimaryGeneratorAction;
  SetUserAction(Test15Generator);
  SetUserAction(new Test15RunAction);
  // Test15EventAction* eventAction = new Test15EventAction(Test15Generator);
  Test15EventAction* eventAction = new Test15EventAction();
  SetUserAction(eventAction);
  SetUserAction(new Test15SteppingAction(eventAction));
  SetUserAction(new Test15StackingAction(eventAction));

  // // RunAction is inherited by EventAction for output filenames - will all
  // // change when implement proper analysis manager?
  // Test15RunAction* Test15Run = new Test15RunAction;
  // runManager->SetUserAction(Test15Run);
  // Test15EventAction* eventAction = new Test15EventAction(Test15Run,Test15Generator);

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
