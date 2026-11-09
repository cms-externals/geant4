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
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//
// This physics list is taken from the underground_physics example with small
// modifications. Its physics is registered as a constructor so additional
// constructors can be used through the modular physics list interface.
// The physics covered would be suitable for a low background
// experiment including the neutron_hp package
//
// Additional constructors must complement the built-in physics (for example,
// step limiters, biasing or parallel worlds). The built-in constructor has type
// zero: ReplacePhysics does not replace its EM, optical, hadronic or decay
// components. Adding overlapping physics, including G4OpticalPhysics, can
// register duplicate processes.
//
//
//
// PhysicsList header
// --------------------------------------------------------------

#ifndef LBE_HH
#define LBE_HH

#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4StoppingPhysics;  // This builder encapsulate stopping processes

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LBE : public G4VModularPhysicsList
{
  public:

    LBE(G4int ver = 1);
    virtual ~LBE();

    // delete copy constructor and assignment operator
    LBE(const LBE&) = delete;
    LBE& operator=(const LBE& right) = delete;

  public:

    virtual void SetCuts();

  protected:

    // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructOp();

    /*
    // these methods Construct all particles in each category
    virtual void ConstructAllBosons();
    virtual void ConstructAllLeptons();
    virtual void ConstructAllMesons();
    virtual void ConstructAllBaryons();
    virtual void ConstructAllIons();
    virtual void ConstructAllShortLiveds();
    */

    // Called after the modular list adds transportation, before LBE physics.
    // Subclasses customizing the former AddTransportation hook should use this
    // hook without adding transportation themselves.
    virtual void AddSpecialCuts();

  private:

    // Construct LBE physics under the modular list's construction lock.
    class Physics : public G4VPhysicsConstructor
    {
      public:

        explicit Physics(LBE* list);

        void ConstructParticle() override;
        void ConstructProcess() override;

      private:

        LBE* physicsList;
    };

    G4int VerboseLevel;
    G4int OpVerbLevel;

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    // uncomment if the code if corresponding code in LBE.cc is un-commented
    //  G4double cutForProton;
    //  G4double cutForAlpha;
    //  G4double cutForGenericIon;

    G4StoppingPhysics* stoppingPhysics;

    // these methods Construct particles
    void ConstructMyBosons();
    void ConstructMyLeptons();
    void ConstructMyMesons();
    void ConstructMyBaryons();
    void ConstructMyIons();
    void ConstructMyShortLiveds();
};

#endif
