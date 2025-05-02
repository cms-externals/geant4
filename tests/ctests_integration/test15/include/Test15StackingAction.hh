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
// StackingAction header
// --------------------------------------------------------------

#ifndef Test15StackingAction_H
#define Test15StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ParticleDefinition.hh"


class G4Navigator;
class G4Track;
class G4Element;
class G4ParticleDefinition;

class Test15EventAction;

class Test15StackingAction : public G4UserStackingAction {

  public:
    Test15StackingAction(Test15EventAction*);
    virtual ~Test15StackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();

    void printCrossSection(G4ParticleDefinition * particleType, const G4Element* element);

  private:

  G4int fNumber_newtracks;
  G4int fNeutron;
  G4int fProton;
  G4int fDeuteron;
  G4int fOther;
    G4bool killGammasFlag;

    // G4double energy;

    G4Navigator* gNavigator; 

  Test15EventAction* evtAction;

  public:
    inline void SetKillGammasFlag(G4bool val)     {killGammasFlag  = val;};

};

#endif

