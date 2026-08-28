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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Ilker Parmaksiz
//
// History
//   November 15th, 2025 : first implementation of GPU hits collection
//
// ********************************************************************
//
/// \file OpticksHitHandler.cc
/// \brief Implementation of the OpticksHitHandler class

#ifndef GDMLOPTICKS_OPTICKSHITHANDLER_HH
#define GDMLOPTICKS_OPTICKSHITHANDLER_HH
// Opticks Headers
#include "OpticksGenstep.h"
#include "OpticksPhoton.hh"
#include "SEvt.hh"
#include "QSim.hh"
#include "G4CXOpticks.hh"
// G4 Headers
#include "G4RunManager.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"

class OpticksHitHandler
{
  public:

    static OpticksHitHandler* getInstance()
    {
      G4AutoLock lock(&mtx);
      if (instance == nullptr)
      {
        instance = new OpticksHitHandler();
      }
      return instance;
    };

    struct OpticksHit
    {
        G4int hit_id;
        G4int sensor_id;
        G4double time;
        G4ThreeVector position;
        G4ThreeVector direction;
        G4ThreeVector polarization;
        G4double wavelength;
        G4int creationId;
        unsigned orient;
        unsigned boundary_flag;
        unsigned flag_mask;
        unsigned iindex;
    };

    void CollectHits();
    void AddHits();
    void SaveHits();
    G4int GetCreationProcessId(unsigned flagmask);
    std::vector<OpticksHit> GetHits();

  private:

    OpticksHitHandler() {};
    static OpticksHitHandler* instance;
    static G4Mutex mtx;
    std::vector<sphoton> sphotons;
    std::vector<OpticksHit> hits;
};

inline G4int OpticksHitHandler::GetCreationProcessId(unsigned flagmask)
{
  if (OpticksPhoton::HasCerenkovFlag(flagmask)) return 0;
  if (OpticksPhoton::HasScintillationFlag(flagmask)) return 1;
  return -1;
}

inline std::vector<OpticksHitHandler::OpticksHit> OpticksHitHandler::GetHits()
{
  std::vector<OpticksHit> hits_swapped;
  {
    G4AutoLock lock(&mtx);
    hits_swapped.swap(hits); // This creates a copy of the hits vector
  }
  return hits_swapped;
}
#endif  // GDMLOPTICKS_OPTICKSHITHANDLER_HH
