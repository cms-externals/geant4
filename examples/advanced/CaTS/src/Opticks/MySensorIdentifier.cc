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
//   November 15th, 2025 : first implementation of SensorIdentification Class for Opticks
//
// ********************************************************************
//
/// \file MySensorIdentifier.cc
/// \brief Implementation of the MySensorIdentifier class

//Cats Header
#include "Opticks/MySensorIdentifier.hh"
//G4 Headers
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

MySensorIdentifier::MySensorIdentifier(std::map<G4String, G4int>& ids) : fDetectIds(ids) {}
MySensorIdentifier::~MySensorIdentifier() {}
G4int MySensorIdentifier::getInstanceIdentity(const G4VPhysicalVolume* pv) const
{
  // For instanced geometry, just return the copy number
  // Not using
  return -1;
}

G4int MySensorIdentifier::getGlobalIdentity(const G4VPhysicalVolume* pv, const G4VPhysicalVolume* ppv)
{

  if (fDetectIds.size() != 0)
  {
    auto it = fDetectIds.find(pv->GetName());
    if (it != fDetectIds.end())
    {
      // Return the build detector id or just generate one
      // if the configured id is > 0 return it, otherwise generate a new positive id
      return (it->second > 0) ? it->second : ++ids;
    }
    return -1;
  }
  G4cout << " Could not find any detector IDs" << G4endl;

  return -1;
}

void MySensorIdentifier::setLevel(int _level) {}
