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
// $Id: G4TouchableGlobalTransformScene.cc 66773 2013-01-12 14:48:08Z allison $
//
//
// John Allison  22nd August 2018
// An artificial scene to find the global transfromation of a physical volume.

#include "G4TouchableGlobalTransformScene.hh"

//#include "G4PhysicalVolumeModel.hh"

G4TouchableGlobalTransformScene::G4TouchableGlobalTransformScene
(G4PhysicalVolumeModel* pPVModel,
 const G4ModelingParameters::PVNameCopyNoPath& requiredTouchable)
:fpPVModel(pPVModel)
,fRequiredTouchable(requiredTouchable)
,fpFoundGlobalTransform(nullptr)
{}

G4TouchableGlobalTransformScene::~G4TouchableGlobalTransformScene () {}

void G4TouchableGlobalTransformScene::ProcessVolume (const G4VSolid& /*solid*/) {

  const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>&
  fullPVPath = fpPVModel->GetFullPVPath();

  fpFoundGlobalTransform = nullptr;
  if (fRequiredTouchable.size() == fullPVPath.size()) {
    G4ModelingParameters::PVNameCopyNoPathConstIterator iNameCopyNo;
    std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>::const_iterator
    iPVNodeId;
    for (iNameCopyNo = fRequiredTouchable.begin(), iPVNodeId = fullPVPath.begin();
         iNameCopyNo != fRequiredTouchable.end();
         ++iNameCopyNo, ++iPVNodeId) {
      if (!(
            iNameCopyNo->GetName() ==
            iPVNodeId->GetPhysicalVolume()->GetName() &&
            iNameCopyNo->GetCopyNo() ==
            iPVNodeId->GetPhysicalVolume()->GetCopyNo()
            )) {
        break;
      }
    }
    if (iNameCopyNo == fRequiredTouchable.end()) {
      fpFoundGlobalTransform = fpPVModel->GetCurrentTransform();
      fpPVModel->Abort();  // No need to look further.
    }
  }
}
