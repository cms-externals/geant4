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
// SteppingAction header
// --------------------------------------------------------------

#ifndef Test15SteppingAction_h
#define Test15SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <map>

class Test15EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test15SteppingAction : public G4UserSteppingAction
{
  public:
    Test15SteppingAction(Test15EventAction*);
    virtual ~Test15SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  public:
    G4double GetPrimaryEnergy() {return startEnergy;};
    G4double GetPrimaryTime() {return startTime;};

  private:

  G4int number_shells;

  G4double outer_radius[26];
  G4double inner_radius[26];
  G4double shell_outer_radius;
  G4double shell_inner_radius;


    G4double                    startEnergy;
    G4double                    startTime;

    Test15EventAction*    evtAction;  //pointer to event action

    G4bool flag;

    std::map<G4int,G4double,std::less<G4int> > parent_energy;
    std::map<G4int,G4String,std::less<G4int> > parent_particle;
    std::map<G4int,G4int,std::less<G4int> > parent_particleID;

    G4int number_generations;

};

#endif
