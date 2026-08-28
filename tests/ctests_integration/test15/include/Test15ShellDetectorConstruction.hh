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

#ifndef Test15ShellDetectorConstruction_hh
#  define Test15ShellDetectorConstruction_hh Test15ShellDetectorConstruction_hh

#  include "G4GeometryCell.hh"
#  include "G4VUserParallelWorld.hh"
#  include "globals.hh"

#  include "Test15PVolumeStore.hh"

#  include <map>
#  include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VIStore;

class Test15ShellDetectorConstruction : public G4VUserParallelWorld
{
  public:

    Test15ShellDetectorConstruction(G4String worldName);
    ~Test15ShellDetectorConstruction();

    G4VPhysicalVolume& GetWorldVolumeAddress() const;
    G4VPhysicalVolume* GetWorldVolume();

    // G4VPhysicalVolume &GetShellVolumeAddress() const;

    const G4VPhysicalVolume& GetPhysicalVolumeByName(const G4String& name) const;
    G4String ListPhysNamesAsG4String();

    G4String GetCellName(G4int i);
    G4GeometryCell GetGeometryCell(G4int i);

    void SetSensitive();

    virtual void Construct();
    virtual void ConstructSD();

    G4VIStore* CreateImportanceStore();
    // create an importance store, caller is responsible for deleting it

    G4String int_to_string(int x);

  private:

    // void Construct();

    Test15PVolumeStore fPVolumeStore;

    G4VPhysicalVolume* ghostWorld;

    // G4LogicalVolume* shellLogical;
    // G4VPhysicalVolume* shellPhysical;
    G4LogicalVolume* radialLogical[26];
    G4VPhysicalVolume* radialPhysical[26];

    std::vector<G4LogicalVolume*> fLogicalVolumeVector;

    G4int number_shells;

    G4double outer_radius[26];
    G4double inner_radius[26];
    G4double shell_outer_radius;
    G4double shell_inner_radius;

    G4bool olap_test;
};

inline G4String Test15ShellDetectorConstruction::int_to_string(int i)
{
  std::ostringstream os;
  //   std::string o_string;
  if (!(os << i)) return "null";
  return os.str();
}

#endif
