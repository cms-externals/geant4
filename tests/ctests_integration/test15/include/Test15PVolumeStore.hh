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
// GEANT4 tag
//
// ----------------------------------------------------------------------
// Class Test15PVolumeStore
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Test15PVolumeStore_hh
#  define Test15PVolumeStore_hh Test15PVolumeStore_hh

#  include "G4GeometryCell.hh"
#  include "G4GeometryCellComp.hh"
#  include "globals.hh"

#  include <set>

typedef std::set<G4GeometryCell, G4GeometryCellComp> Test15SetGeometryCell;

class Test15PVolumeStore
{
  public:

    Test15PVolumeStore();
    ~Test15PVolumeStore();

    void AddPVolume(const G4GeometryCell& cell);
    const G4VPhysicalVolume* GetPVolume(const G4String& name) const;
    G4int Size();
    G4String GetPNames() const;

  private:

    Test15SetGeometryCell fSetGeometryCell;
};

#endif
