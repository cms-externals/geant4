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

#include "ExN01PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e-"));
  ////First Test =No Intersection
  G4double momen=(0.3*2.*0.70001)*GeV;
  G4double kinEnergy=std::sqrt(momen*momen+0.000510999*0.000510999*GeV);
  particleGun->SetParticleEnergy(kinEnergy);
  particleGun->SetParticlePosition(G4ThreeVector(-70.0001*cm, 0.0*cm,-0.0*cm));
  /////Second Test=Intersection But Almost Not !!!Very slow step but not in Propagator
  //G4double momen=(0.3*2.*0.50001)*GeV;
  // G4double kinEnergy=std::sqrt(momen*momen+0.000510999*0.000510999*GeV);
  //  particleGun->SetParticleEnergy(kinEnergy);
  //  particleGun->SetParticlePosition(G4ThreeVector(-70.0001*cm, 0.0*cm,-0.0*cm));
  /////Third Test=No Intersection But One "Difficult" Point
  //Test simple 
  //   G4double momen=(0.3*2.*0.20001)*GeV;
  //   G4double kinEnergy=std::sqrt(momen*momen+0.000510999*0.000510999*GeV);
  //  particleGun->SetParticleEnergy(kinEnergy);
  //  particleGun->SetParticlePosition(G4ThreeVector(-70.0001*cm, 0.0*cm,-0.0*cm));
  
    G4cout<<"Energy="<<kinEnergy/GeV<<"  Radius="<<momen/0.6/mm<<G4endl;
   
}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4int i = anEvent->GetEventID() % 3;
  G4ThreeVector v(0.0,1.0,0.01);
  // switch(i)
  //{
  //  case 0:
  //    break;
  //  case 1:
  //    v.setY(0.1);
  //    break;
  //  case 2:
  //    v.setZ(0.1);
  //    break;
  //}
  particleGun->SetParticleMomentumDirection(v);
  particleGun->GeneratePrimaryVertex(anEvent);
}


