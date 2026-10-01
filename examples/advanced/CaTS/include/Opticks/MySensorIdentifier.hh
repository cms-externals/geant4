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

#ifndef GDMLOPTICKS_MYSENSORIDENTIFIER_HH
#define GDMLOPTICKS_MYSENSORIDENTIFIER_HH

//Opticks Headers
#include "U4SensorIdentifier.h"
//G4 Headers
#include "G4String.hh"
// C++ Headers
#include <map>
#include <string>
#include <vector>

class MySensorIdentifier : public U4SensorIdentifier
{
  public:

    MySensorIdentifier(std::map<G4String, G4int>& ids);
    ~MySensorIdentifier();
    virtual void setLevel(G4int _level);
    G4int getGlobalIdentity(const G4VPhysicalVolume*, const G4VPhysicalVolume*) override;
    G4int getInstanceIdentity(const G4VPhysicalVolume*) const override;

  private:

    std::map<G4String, G4int>& fDetectIds;
    G4int ids = 0;
};

#endif  // GDMLOPTICKS_MYSENSORIDENTIFIER_HH
