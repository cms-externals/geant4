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

// Opticks Hit Collection
// CaTS Headers
#include "Opticks/OpticksHitHandler.hh"

OpticksHitHandler* OpticksHitHandler::instance = nullptr;
G4Mutex OpticksHitHandler::mtx;

void OpticksHitHandler::CollectHits()
{
  // Collecting Opticks Photons
  G4AutoLock lock(&mtx);
  SEvt* sev = SEvt::Get_EGPU();
  sphoton::Get(sphotons, sev->getHit());
  auto run = G4RunManager::GetRunManager();
  G4int eventID = run->GetCurrentEvent()->GetEventID();
  G4ThreeVector position, direction, polarization;
  for (auto& hit : sphotons)
  {
    OpticksHit ohit = OpticksHit();
    position = G4ThreeVector(hit.pos.x, hit.pos.y, hit.pos.z);
    direction = G4ThreeVector(hit.mom.x, hit.mom.y, hit.mom.z);
    polarization = G4ThreeVector(hit.pol.x, hit.pol.y, hit.pol.z);
    ohit.hit_id = hit.index;
    ohit.sensor_id = hit.identity;
    ohit.position = position;
    ohit.direction = direction;
    ohit.polarization = polarization;
    ohit.time = hit.time;
    ohit.boundary_flag = hit.boundary();
    ohit.wavelength = hit.wavelength;
    ohit.orient = hit.orient();
    ohit.flag_mask = hit.flagmask;
    ohit.boundary_flag = hit.boundary();
    ohit.iindex = hit.iindex();
    ohit.creationId = GetCreationProcessId(hit.flagmask);
    hits.push_back(ohit);
  }

  // clear the hits
  sphotons.clear();
  sphotons.shrink_to_fit();
  G4CXOpticks::Get()->reset(eventID);
  QSim::Get()->reset(eventID);

}
