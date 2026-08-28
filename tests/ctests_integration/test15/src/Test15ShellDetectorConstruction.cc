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

#include "Test15ShellDetectorConstruction.hh"

#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

#include <sstream>

// For Primitive Scorers
#include "G4MultiFunctionalDetector.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"
#include "G4SDManager.hh"
#include "G4SDParticleFilter.hh"

// for importance biasing
#include "G4IStore.hh"

// #include "Test15AnalysisManager.hh"

Test15ShellDetectorConstruction::Test15ShellDetectorConstruction(G4String worldName)
  : G4VUserParallelWorld(worldName), fLogicalVolumeVector(), olap_test(true)
{
  //  Construct();
  number_shells = 26;
  //  G4double shell_thickness = 10.0*mm;
  G4double shell_thickness = 2.0 * mm;

  shell_outer_radius = 457.0 * mm;
  shell_inner_radius = shell_outer_radius - shell_thickness;

  G4double radii_start[] = {200.0 * cm, 190.0 * cm, 185.0 * cm, 175.0 * cm, 165.0 * cm, 150.0 * cm,
                            140.0 * cm, 130.0 * cm, 120.0 * cm, 110.0 * cm, 100.0 * cm, 90.0 * cm,
                            80.0 * cm,  70.0 * cm,  60.0 * cm,  50.0 * cm,  45.7 * cm,  40.0 * cm,
                            30.0 * cm,  25.0 * cm,  20.0 * cm,  15.0 * cm,  10.0 * cm,  8.0 * cm,
                            5.0 * cm,   3.0 * cm};  // care not to double count 45.7cm shell

  for (G4int i = 0; i < number_shells; ++i)
  {
    outer_radius[i] = radii_start[i];
    inner_radius[i] = radii_start[i] - shell_thickness;
  }
}

Test15ShellDetectorConstruction::~Test15ShellDetectorConstruction()
{
  fLogicalVolumeVector.clear();
}

void Test15ShellDetectorConstruction::Construct()
{
  G4cout << " constructing parallel world " << G4endl;

  G4Material* dummyMat = 0;

  // GetWorld methods create a clone of the mass world to the parallel world (!)
  //  via the transportation manager
  ghostWorld = GetWorld();
  G4cout << " Test15ShellDetectorConstruction:: ghostWorldName = " << ghostWorld->GetName()
         << G4endl;
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  fLogicalVolumeVector.push_back(worldLogical);

  //  fPVolumeStore.AddPVolume(G4GeometryCell(*pWorldVolume, 0));
  fPVolumeStore.AddPVolume(G4GeometryCell(*ghostWorld, 0));

  // still needed? 27/09/15:
  // Test15AnalysisManager* analysis =  Test15AnalysisManager::getInstance();

  G4Colour red(1.0, 0.0, 0.0);
  G4VisAttributes* shell_vat = new G4VisAttributes(red);
  //  shell_log->SetVisAttributes(G4VisAttributes::Invisible);
  shell_vat->SetVisibility(true);
  shell_vat->SetForceSolid(true);

  // still needed? 27/09/15:
  // G4double outer_radius = analysis->GetShellOuterRadius();
  // G4double inner_radius = analysis->GetShellInnerRadius();

  // Don't double place shell at 45.7cm!!
  /*
  G4Sphere* shellSphere = new G4Sphere("shellSphere",shell_inner_radius,shell_outer_radius,0.*deg,
  360.*deg, 0.*deg, 180.*deg);
  // G4Box* shellSphere = new
  G4Box("shellSphere",0.5*shell_inner_radius,0.5*shell_outer_radius,0.5*shell_outer_radius);

  G4cout << " OUTER -shell is: " << shell_inner_radius << " outer: " << shell_outer_radius <<
  G4endl;

    // creating 18 slobs of 10 cm thicknes

    // logical parallel cells

  shellLogical =
    new G4LogicalVolume(shellSphere, dummyMat, "aShell_log");

  fLogicalVolumeVector.push_back(shellLogical);

  G4String name =  GetCellName(0);

  shellPhysical = new G4PVPlacement(0, G4ThreeVector(-0., -0., -0.), name, shellLogical, ghostWorld,
  false,0,olap_test);

  shellLogical->SetVisAttributes(shell_vat);

  G4GeometryCell cell(*shellPhysical, 0);

  G4cout << " adding pvolume within ShellDetectorConstruction " << G4endl;

  fPVolumeStore.AddPVolume(cell);

*/

  G4String name;  // =  GetCellName(0);

  // still needed? 27/09/15:
  // G4int number_shells = analysis->GetNumberShells();
  // linked above

  for (G4int i = 0; i < number_shells; ++i)
  {
    // still needed? 27/09/15:
    // G4double outer_radius = analysis->GetRadialOuterRadius(i);
    // G4double inner_radius = analysis->GetRadialInnerRadius(i);

    G4cout << " i-shell is: " << i << " inner_radius[i]: " << inner_radius[i]
           << " outer: " << outer_radius[i] << G4endl;
    G4Sphere* radialSphere = new G4Sphere("shellSphere", inner_radius[i], outer_radius[i], 0. * deg,
                                          360. * deg, 0. * deg, 180. * deg);
    // G4Box* radialSphere = new
    // G4Box("shellSphere",0.5*inner_radius[i],0.5*outer_radius[i],0.5*outer_radius);

    // creating 18 slobs of 10 cm thicknes

    // logical parallel cells

    G4String logical_name = "aRadial_log_" + int_to_string(i);

    std::ostringstream os;
    os << i;
    logical_name.append(os.str());

    G4cout << " logical name: " << logical_name << G4endl;
    radialLogical[i] = new G4LogicalVolume(radialSphere, dummyMat, logical_name);
    //    new G4LogicalVolume(radialSphere, dummyMat, "aRadial_log");

    fLogicalVolumeVector.push_back(radialLogical[i]);
    //    name = GetCellName(i+1);
    name = GetCellName(i);

    radialPhysical[i] = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), name, radialLogical[i],
                                          ghostWorld, false, i, olap_test);

    radialLogical[i]->SetVisAttributes(shell_vat);

    G4GeometryCell radcell(*radialPhysical[i], i);

    G4cout << " adding pvolume within ShellDetectorConstruction " << G4endl;

    fPVolumeStore.AddPVolume(radcell);
  }

  SetSensitive();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4VPhysicalVolume&
Test15ShellDetectorConstruction::GetPhysicalVolumeByName(const G4String& name) const
{
  return *fPVolumeStore.GetPVolume(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Test15ShellDetectorConstruction::ListPhysNamesAsG4String()
{
  G4String names(fPVolumeStore.GetPNames());
  return names;
}

G4VPhysicalVolume& Test15ShellDetectorConstruction::GetWorldVolumeAddress() const
{
  return *ghostWorld;
}

// G4VPhysicalVolume &Test15ShellDetectorConstruction::GetShellVolumeAddress() const{
//    return *shellPhysical;
// }

G4VPhysicalVolume* Test15ShellDetectorConstruction::GetWorldVolume()
{
  return ghostWorld;
}

G4String Test15ShellDetectorConstruction::GetCellName(G4int i)
{
  std::ostringstream os;
  os << "cell_";
  if (i < 10)
  {
    os << "0";
  }
  os << i;
  G4String name = os.str();
  G4cout << " cell name is : " << name << G4endl;
  return name;
}

G4GeometryCell Test15ShellDetectorConstruction::GetGeometryCell(G4int i)
{
  G4String name(GetCellName(i));
  const G4VPhysicalVolume* p = 0;
  p = fPVolumeStore.GetPVolume(name);
  if (p)
  {
    return G4GeometryCell(*p, 0);
  }
  else
  {
    G4cout << "Test15ShellDetectorConstruction::GetGeometryCell: couldn't get G4GeometryCell"
           << G4endl;
    return G4GeometryCell(*ghostWorld, -2);
  }
}

//--------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Test15ShellDetectorConstruction::SetSensitive()
{
  //  -------------------------------------------------
  //   The collection names of defined Primitives are
  //   0       ConcreteSD/Collisions
  //   1       ConcreteSD/CollWeight
  //   2       ConcreteSD/Population
  //   3       ConcreteSD/TrackEnter
  //   4       ConcreteSD/SL
  //   5       ConcreteSD/SLW
  //   6       ConcreteSD/SLWE
  //   7       ConcreteSD/SLW_V
  //   8       ConcreteSD/SLWE_V
  //  -------------------------------------------------

  // moved to ConstructSD() for MT compliance
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Test15ShellDetectorConstruction::ConstructSD()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //
  // Sensitive Detector Name
  G4String concreteSDname = "ConcreteSD";

  //------------------------
  // MultiFunctionalDetector
  //------------------------
  //
  // Define MultiFunctionalDetector with name.
  G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(concreteSDname);
  SDman->AddNewDetector(MFDet);  // Register SD to SDManager

  G4String fltName, particleName;
  G4SDParticleFilter* neutronFilter =
    new G4SDParticleFilter(fltName = "neutronFilter", particleName = "neutron");

  MFDet->SetFilter(neutronFilter);

  for (std::vector<G4LogicalVolume*>::iterator it = fLogicalVolumeVector.begin();
       it != fLogicalVolumeVector.end(); it++)
  {
    //      (*it)->SetSensitiveDetector(MFDet);
    SetSensitiveDetector((*it)->GetName(), MFDet);
  }

  G4String psName;
  G4PSNofCollision* scorer0 = new G4PSNofCollision(psName = "Collisions");
  MFDet->RegisterPrimitive(scorer0);

  G4PSNofCollision* scorer1 = new G4PSNofCollision(psName = "CollWeight");
  scorer1->Weighted(true);
  MFDet->RegisterPrimitive(scorer1);

  G4PSPopulation* scorer2 = new G4PSPopulation(psName = "Population");
  MFDet->RegisterPrimitive(scorer2);

  G4PSTrackCounter* scorer3 = new G4PSTrackCounter(psName = "TrackEnter", fCurrent_In);
  MFDet->RegisterPrimitive(scorer3);

  G4PSTrackLength* scorer4 = new G4PSTrackLength(psName = "SL");
  MFDet->RegisterPrimitive(scorer4);

  G4PSTrackLength* scorer5 = new G4PSTrackLength(psName = "SLW");
  scorer5->Weighted(true);
  MFDet->RegisterPrimitive(scorer5);

  G4PSTrackLength* scorer6 = new G4PSTrackLength(psName = "SLWE");
  scorer6->Weighted(true);
  scorer6->MultiplyKineticEnergy(true);
  MFDet->RegisterPrimitive(scorer6);

  G4PSTrackLength* scorer7 = new G4PSTrackLength(psName = "SLW_V");
  scorer7->Weighted(true);
  scorer7->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer7);

  G4PSTrackLength* scorer8 = new G4PSTrackLength(psName = "SLWE_V");
  scorer8->Weighted(true);
  scorer8->MultiplyKineticEnergy(true);
  scorer8->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer8);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VIStore* Test15ShellDetectorConstruction::CreateImportanceStore()
{
  G4cout << " Test15ShellDetectorConstruction:: Creating Importance Store " << G4endl;
  if (!fPVolumeStore.Size())
  {
    G4Exception("Test15ShellDetectorConstruction::CreateImportanceStore", "exampleTest15_0001",
                RunMustBeAborted, "no physical volumes created yet!");
  }

  // creating and filling the importance store

  //  G4IStore *istore = new G4IStore(*fWorldVolume);

  G4IStore* istore = G4IStore::GetInstance(GetName());

  G4GeometryCell gWorldVolumeCell(GetWorldVolumeAddress(), 0);

  G4double imp = 1;

  istore->AddImportanceGeometryCell(1, gWorldVolumeCell);

  // set importance values and create scorers

  // set importance values and create scorers
  //  G4int number_shells = analysis->GetNumberShells();
  G4int cell(26);
  for (cell = 0; cell < 26; cell++)
  {
    G4GeometryCell gCell = GetGeometryCell(cell);
    G4cout << " adding cell: " << cell << " replica: " << gCell.GetReplicaNumber()
           << " name: " << gCell.GetPhysicalVolume().GetName() << G4endl;
    imp = 1;
    //    G4double imp = std::pow(2.0,cell-1);
    // x    aIstore.AddImportanceGeometryCell(imp, gCell);
    istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), cell);
    // adding the standard G4CellScorer for 17 concrete cells
    //     if (cell<18) {
    //       b02store.AddG4CellScorer(gCell);
    //     }
  }

  return istore;

  //--------------------------------------------------------------------
}

//--------------------------------------------------------------------
