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
// GEANT 4 class source file
//
// Test15PVolumeStore.cc
//
// ----------------------------------------------------------------------

#include "Test15PVolumeStore.hh"
#include <sstream>


#include "G4VPhysicalVolume.hh"

Test15PVolumeStore::Test15PVolumeStore(){}
Test15PVolumeStore::~Test15PVolumeStore(){}
  
void Test15PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  Test15SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  G4cout << " adding to pVolumeStore " << G4endl;
  if (it != fSetGeometryCell.end()) {
    G4cout << "Test15PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *Test15PVolumeStore::
GetPVolume(const G4String &name) const {
  const G4VPhysicalVolume *pvol = 0;
  for (Test15SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    G4cout << " Looking at volume: " << vol.GetName() << " compared to: " << name << G4endl;
    if (vol.GetName() == name) {
      pvol =  &vol;
    } 
  }
  if (!pvol) {
    G4cout << "Test15PVolumeStore::GetPVolume: no physical volume named: " 
	   << name << ", found" << G4endl;
  }
  return pvol;
}

G4String Test15PVolumeStore::GetPNames() const {
  G4String NameString;
  for (Test15SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    std::ostringstream os;
    os << vol.GetName() << "_" << it->GetReplicaNumber() 
       << "\n";
    G4String cellname = os.str();
    
    //    G4String cellname(vol.GetName());
    //    cellname += G4String("_");
    //    cellname += std::str(it->GetReplicaNumber());

    NameString += cellname;
  }
  return NameString;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int Test15PVolumeStore::Size() {
  return fSetGeometryCell.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
