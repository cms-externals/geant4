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
/// \file RunAction.cc
/// \brief Implementation of the volumeTreeVisitor::RunAction class

#include "RunAction.hh"

#include "G4RunManager.hh"

#include "G4TouchableWalker.hh"
#include "G4TouchableCollector.hh"

#include "G4TransportationManager.hh"

#include "DetectorConstruction.hh"


namespace volumeTreeVisitor
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/ )
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  static G4RunManager* rm= G4RunManager::GetRunManager();
  const G4VUserDetectorConstruction* detectorCons= rm->GetUserDetectorConstruction();
  auto  myDetectorCons = static_cast<const volumeTreeVisitor::DetectorConstruction*>(detectorCons);
  const G4LogicalVolume* targetLV= myDetectorCons->GetTargetLogical();

  static G4int runNumber =0;
  runNumber++;
  G4cout << "Run number: " << runNumber << G4endl;

  G4bool verbose= true;

  if( runNumber <= 5) {
    G4cout << G4endl;
    G4cout << "RunAction -- BeginOfRunAction " << G4endl
           << "==============================" << G4endl << G4endl;

    G4int found = FindTouchablesAndReport(targetLV, verbose);
    G4cout << "Found " << found << " touchables in targetLV" << G4endl;

    for( G4int idx=0; idx<myDetectorCons->GetNChambers(); ++idx ) {
      const G4LogicalVolume* chamberLV= myDetectorCons->GetChamberLogical(idx);
      found = FindTouchablesAndReport(chamberLV, verbose);
      G4cout << "Found " << found << " touchables in chamberLV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int RunAction::FindTouchablesAndReport(const G4LogicalVolume* seekingLogVol, G4bool verbose)
{
  G4cout << "RunAction::FindTouchablesAndReport constructing G4TouchableCollector() with logicalVol = "
         << seekingLogVol << " name = " << seekingLogVol->GetName() << G4endl;

  G4TouchableCollector listBuilder(seekingLogVol);
  G4TouchableWalker    lbWalker(listBuilder);

  G4VPhysicalVolume* worldPV=
          G4TransportationManager::GetTransportationManager()
                                 ->GetNavigatorForTracking()
                                 ->GetWorldVolume();

  lbWalker.Walk(worldPV);

  auto listOfTouchables = listBuilder.GetListOfMatchingTouchables();

  const char* VolumeTypeNames[4] = { "kNormal",  "kReplica",  "kParameterised",  "kExternal" };

  G4cout << "Touchables collected - with matching LV are: "<< G4endl;
  for(auto item : listOfTouchables )
  {
    G4int depth= item->GetHistoryDepth();
    G4VPhysicalVolume* pv= item->GetVolume();
    G4VSolid* solid= item->GetSolid();
    G4int replicaNo= item->GetReplicaNumber();
    G4LogicalVolume *lv= pv->GetLogicalVolume();
    EVolume volumeType = pv->VolumeType();
    G4String volTypeName= G4String( VolumeTypeNames[volumeType]);

    G4String logVolName= ( lv != nullptr ? lv->GetName() : G4String("-Not-Available-") );

    G4cout << "Touchable top-volume: " // << Indent(depth)
           << "[PV] " << pv->GetName()
           << "  [ type= " << volTypeName << " ] "
           << "  repNum= " << replicaNo
           << "  pv-copy=" << pv->GetCopyNo()
           << " -- at depth= " << depth
           << "  logical volume= '" << logVolName << "' ";
   if( verbose )
    {
      G4String    solidName    = ( solid != nullptr ? solid->GetName() : G4String("-Not-Available-") );
      G4Material* material     = ( lv != nullptr ? lv->GetMaterial() : nullptr );
      G4String    materialName = ( material != nullptr ? lv->GetMaterial()->GetName() : G4String("-Not-Available-") );
      G4cout  << " solid = '" << solidName << "' "
              << " material = '" << materialName << "' ";
    }
    G4cout << G4endl;

    if(lv != seekingLogVol){
        G4ExceptionDescription desc;
        desc << "Mismatch in logical volume - does not agree with requested Target";
        desc << "  lv= " << lv << " name= " << lv->GetName()
           << " seekingLogVol= " << seekingLogVol << " name= " << seekingLogVol->GetName();
        G4Exception("RunAction::FindTouchablesAndReport", "User-code",
            JustWarning, desc);
    }
  }
  G4cout << "End of list."<<G4endl<<G4endl;

  G4int numberOfTouchables= listOfTouchables.size();
  return numberOfTouchables;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace volumeTreeVisitor
