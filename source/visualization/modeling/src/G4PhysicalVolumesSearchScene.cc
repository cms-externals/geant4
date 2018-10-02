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
// $Id: G4PhysicalVolumesSearchScene.cc 111780 2018-09-03 18:20:03Z allison $
//
// 
// John Allison  5th September 2018, based on G4PhysicalVolumeSearchScene
// An artificial scene to find physical volumes. Instead of returning the
// first occurence (G4PhysicalVolumeSearchScene) this class (note the extra
// 's' in the name of this class) returns a vector of all occurences.

#include "G4PhysicalVolumesSearchScene.hh"

void G4PhysicalVolumesSearchScene::ProcessVolume (const G4VSolid&)
{
  G4VPhysicalVolume* pCurrentPV = fpSearchVolumeModel->GetCurrentPV();
  const G4String& name = pCurrentPV->GetName();
  G4int copyNo = fpSearchVolumeModel->GetCurrentPVCopyNo();
  if (fRequiredPhysicalVolumeName == name) {
    if ((fRequiredCopyNo < 0 ||  // I.e., ignore negative request.
         fRequiredCopyNo == copyNo)) {
      auto path = fpSearchVolumeModel->GetFullPVPath();
      path.pop_back();  // Base node is one up from found node.
      fFindings.push_back
      (Findings
       (fpSearchVolumeModel->GetTopPhysicalVolume(),
        pCurrentPV,
        copyNo,
        fpSearchVolumeModel->GetCurrentDepth(),
        path,
        *fpCurrentObjectTransformation));
    }
  }
}
