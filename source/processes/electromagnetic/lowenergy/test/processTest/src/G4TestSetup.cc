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
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4TestSetup.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include "G4DynamicParticle.hh"

#include "G4UnitsTable.hh"
#include "Randomize.hh"


G4TestSetup::G4TestSetup(const G4Material* aMaterial,
			 const G4ParticleDefinition* def,
			 G4double minEnergy, G4double  maxEnergy)
  : eMin(minEnergy), eMax(maxEnergy)
{
  // const_cast needed because geometry and tracking objects require
  // non-const materials and particle definitions in their constructors (why?)
  part = const_cast<G4ParticleDefinition*>(def);
  material = const_cast<G4Material*>(aMaterial);
  track = 0;
  step = new G4Step;
  physicalFrame = 0;
}

G4TestSetup:: ~G4TestSetup()
{
  delete physicalFrame;
  physicalFrame = 0;
  delete step;
  step = 0;
  delete track;
  track = 0;
}

  void G4TestSetup::makeGeometry()
{
  G4double dimX = 10.*cm;
  G4double dimY = 10.*cm;
  G4double dimZ = 10.*cm;
  
  G4Box box("Box",dimX, dimY, dimZ);
  
  G4LogicalVolume logicalVol(&box,material,"LogicalBox", 0, 0, 0);
  logicalVol.SetMaterial(material); 
  
  physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
				    "PhysicalBox",&logicalVol,0,false,0);

}

const G4Track* G4TestSetup::makeTrack()
{
  G4double energy = eMin + (eMax - eMin) * G4UniformRand();

  if (track == 0)
    {   
      G4double initX = 0.; 
      G4double initY = 0.; 
      G4double initZ = 1.;
      G4ParticleMomentum direction(initX,initY,initZ);
      G4DynamicParticle* dynamicPart = new G4DynamicParticle(part,direction,energy);
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      
      track = new G4Track(dynamicPart,time,position);
      
      // do I really need this?     
      G4GRSVolume* touche = new G4GRSVolume(physicalFrame,0,position);   
      track->SetTouchable(touche);
    }
  else
    {
      track->SetKineticEnergy(energy);
    }

  G4double e = track->GetKineticEnergy();
  G4cout << "Track energy = " << e/MeV << " MeV" << G4endl;

  return track;
}

const G4Step* G4TestSetup::makeStep()
{
  step->SetTrack(track);

  G4ThreeVector aPosition(0.,0.,0.);
  G4ThreeVector newPosition(0.,0.,1.*mm);
  
  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);
  G4StepPoint* newPoint = new G4StepPoint();
  newPoint->SetPosition(newPosition);
  newPoint->SetMaterial(material);
  newPoint->SetSafety(safety);
  step->SetPostStepPoint(newPoint);
  step->SetStepLength(1*mm);

  track->SetStep(step); 
  
return step;
}

