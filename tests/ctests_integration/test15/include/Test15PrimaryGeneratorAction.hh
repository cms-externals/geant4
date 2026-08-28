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
// PrimaryGeneratorAction header
// --------------------------------------------------------------

#ifndef Test15PrimaryGenerator_h
#  define Test15PrimaryGenerator_h 1

#  include "G4VUserPrimaryGeneratorAction.hh"
#  include "globals.hh"

class G4GeneralParticleSource;

class G4ParticleGun;
class G4Event;

class Test15PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    Test15PrimaryGeneratorAction();
    virtual ~Test15PrimaryGeneratorAction();

  public:

    virtual void GeneratePrimaries(G4Event* anEvent);

  private:

    G4GeneralParticleSource* particleGun;

  private:

    //      const long* seeds;
    long seeds[2];
    G4double energy_pri;

  public:

    const long* GetEventSeeds() { return seeds; };
    G4double GetEnergyPrimary() { return energy_pri; };
};

#endif
