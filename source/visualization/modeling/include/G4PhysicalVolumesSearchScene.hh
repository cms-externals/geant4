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
// $Id: G4PhysicalVolumesSearchScene.hh 111780 2018-09-03 18:20:03Z allison $
//
// 
// John Allison  5th September 2018, based on G4PhysicalVolumeSearchScene
// An artificial scene to find physical volumes. Instead of returning the
// first occurence (G4PhysicalVolumeSearchScene) this class (note the extra
// 's' in the name of this class) returns a vector of all occurences.

#ifndef G4PHYSICALVOLUMESSEARCHSCENE_HH
#define G4PHYSICALVOLUMESSEARCHSCENE_HH

#include "G4PseudoScene.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

class G4PhysicalVolumesSearchScene: public G4PseudoScene
{
public:

  G4PhysicalVolumesSearchScene
  (G4PhysicalVolumeModel* pSearchVolumeModel,       // usually a world
   const G4String&        requiredPhysicalVolumeName,
   G4int                  requiredCopyNo)
  : fpSearchVolumeModel           (pSearchVolumeModel)
  , fRequiredPhysicalVolumeName   (requiredPhysicalVolumeName)
  , fRequiredCopyNo               (requiredCopyNo)
  {}

  virtual ~G4PhysicalVolumesSearchScene () {}

  struct Findings
  {
    Findings
    (G4VPhysicalVolume* pSearchPV,
     G4VPhysicalVolume* pFoundPV,
     G4int foundPVCopyNo = 0,
     G4int foundDepth = 0,
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
     foundBasePVPath =
     std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>(),
     G4Transform3D foundObjectTransformation = G4Transform3D())
    : fpSearchPV(pSearchPV)
    , fpFoundPV(pFoundPV)
    , fFoundPVCopyNo(foundPVCopyNo)
    , fFoundDepth(foundDepth)
    , fFoundBasePVPath(foundBasePVPath)
    , fFoundObjectTransformation(foundObjectTransformation) {}
    G4VPhysicalVolume*   fpSearchPV;   // Searched physical volume.
    G4VPhysicalVolume*   fpFoundPV;    // Found physical volume.
    G4int                fFoundPVCopyNo;  // Found Copy number.
    G4int                fFoundDepth;  // Found depth.
    std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
    fFoundBasePVPath;    // Base path (e.g., empty for world volume)
    G4Transform3D        fFoundObjectTransformation;  // Found transformation.
  };

  const std::vector<Findings>& GetFindings() const
  {return fFindings;}

private:

  void ProcessVolume(const G4VSolid&);
  
  const G4PhysicalVolumeModel* fpSearchVolumeModel;
  G4String              fRequiredPhysicalVolumeName;
  G4int                 fRequiredCopyNo;
  std::vector<Findings> fFindings;
};

#endif
