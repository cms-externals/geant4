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
/// \file Tst202ActionInitialization.cc
/// \brief Implementation of the Tst202ActionInitialization class

#include "Tst202ActionInitialization.hh"

#include "Tst202PrimaryGeneratorAction.hh"
// #include "Tst202RunAction.hh"
// #include "Tst202EventAction.hh"
// #include "Tst202SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst202ActionInitialization::Tst202ActionInitialization() : G4VUserActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst202ActionInitialization::~Tst202ActionInitialization() {}

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// void Tst202ActionInitialization::BuildForMaster() const
//{
//  Tst202RunAction* runAction = new Tst202RunAction;
//  SetUserAction(runAction);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst202ActionInitialization::Build() const
{
  SetUserAction(new Tst202PrimaryGeneratorAction);

  //  Tst202RunAction* runAction = new Tst202RunAction;
  //  SetUserAction(runAction);
  //
  //  Tst202EventAction* eventAction = new Tst202EventAction(runAction);
  //  SetUserAction(eventAction);
  //
  //  SetUserAction(new Tst202SteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
