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
// 
// John Allison  May 2021
//
// G4Mesh encapsulates and validates a nested parameterisation, which we
// call a "mesh". If a valid mesh cannot be created out of this
// G4VPhysicalVolume* (which will probably be most common), it will
// have a type "invalid". Then, usually, it may simply be destroyed.
// The overhead of an invalid attempt is expected to be small.

#include "G4Mesh.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4VNestedParameterisation.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Tet.hh"

std::map<G4int,G4String> G4Mesh::fEnumMap = {
  {invalid,"invalid"},
  {rectangle,"rectangle"},
  {cylinder,"cylinder"},
  {sphere,"sphere"},
  {tetrahedron,"tetrahedron"}
};

G4Mesh::G4Mesh (G4VPhysicalVolume* containerVolume,const G4Transform3D& transform)
: fpContainerVolume(containerVolume)
, fpParameterisedVolume(nullptr)
, fMeshType(invalid)
, fMeshDepth(0)
, fTransform(transform)
{
  if (fpContainerVolume == nullptr) return;
    
  G4VPhysicalVolume* pv0 = fpContainerVolume;
  G4VPhysicalVolume* pv1 = nullptr;
  G4VPhysicalVolume* pv2 = nullptr;
  G4VPhysicalVolume* pv3 = nullptr;
  G4LogicalVolume*   lv0 = pv0->GetLogicalVolume();
  G4LogicalVolume*   lv1 = nullptr;
  G4LogicalVolume*   lv2 = nullptr;
  G4LogicalVolume*   lv3 = nullptr;
  
  // Check if this is a container for a parameterisation
  G4bool isContainer = false;
  if (lv0->GetNoDaughters()) {
    fMeshDepth++;
    pv1 = lv0->GetDaughter(0);
    lv1 = pv1->GetLogicalVolume();
    auto pvParam1 = dynamic_cast<G4PVParameterised*>(pv1);
    if (pvParam1) {
      // A 1-deep parametrisation (so not nested)
      isContainer = true;
    } else {
      if (lv1->GetNoDaughters()) {
        fMeshDepth++;
        pv2 = lv1->GetDaughter(0);
        lv2 = pv2->GetLogicalVolume();
        auto pvParam2 = dynamic_cast<G4PVParameterised*>(pv2);
        if (pvParam2) {
          auto param2 = pvParam2->GetParameterisation();
          auto nestedParam2 = dynamic_cast<G4VNestedParameterisation*>(param2);
          if (nestedParam2) { // 2-deep nested mesh
            isContainer = true;
          }
        } else {
          if (lv2->GetNoDaughters()) {
            fMeshDepth++;
            pv3 = lv2->GetDaughter(0);
            lv3 = pv3->GetLogicalVolume();
            auto pvParam3 = dynamic_cast<G4PVParameterised*>(pv3);
            if (pvParam3) {
              auto param3 = pvParam3->GetParameterisation();
              auto nestedParam3 = dynamic_cast<G4VNestedParameterisation*>(param3);
              if (nestedParam3) {  // 3-deep nested mesh
                isContainer = true;
              }
            }
          }
        }
      }
    }
  }
  if (isContainer) {
    // Get type
    if (fMeshDepth == 1) {  // Could be a tetrahedral mesh
      G4VSolid* pEndSol = lv1->GetSolid ();
      if (dynamic_cast<G4Tet*>(pEndSol)) {
        fMeshType = tetrahedron;
        fpParameterisedVolume = pv1;
      }
    } else {  // Take type from container
      G4VSolid* pContainerSol = lv0->GetSolid ();
      if (dynamic_cast<G4Box*>(pContainerSol)) {
        fMeshType = rectangle;
        if (fMeshDepth == 3) {
          pv1->GetReplicationData
          (f3DRPs.fAxis1,f3DRPs.fNreplica1,f3DRPs.fWidth1,f3DRPs.fOffset1,f3DRPs.fConsuming1);
          pv2->GetReplicationData
          (f3DRPs.fAxis2,f3DRPs.fNreplica2,f3DRPs.fWidth2,f3DRPs.fOffset2,f3DRPs.fConsuming2);
          pv3->GetReplicationData
          (f3DRPs.fAxis3,f3DRPs.fNreplica3,f3DRPs.fWidth3,f3DRPs.fOffset3,f3DRPs.fConsuming3);
          auto pBox = static_cast<G4Box*>(lv3->GetSolid());
          f3DRPs.fHalfX = pBox->GetXHalfLength();
          f3DRPs.fHalfY = pBox->GetYHalfLength();
          f3DRPs.fHalfZ = pBox->GetZHalfLength();
        }
        fpParameterisedVolume = pv3;
      } else if (dynamic_cast<G4Tubs*>(pContainerSol)) {
        fMeshType = cylinder;
      } else if (dynamic_cast<G4Sphere*>(pContainerSol)) {
        fMeshType = sphere;
      }
    }
  }
}

G4Mesh::~G4Mesh () {}

std::ostream& operator << (std::ostream& os, const G4Mesh& mesh) {
  os << "G4Mesh: ";
  os << "\nContainer: " << mesh.GetContainerVolume()->GetName();
  const auto& map = mesh.GetEnumMap();
  const auto& typeEntry = map.find(mesh.GetMeshType());
  G4String type;
  if (typeEntry != map.end()) {
    type = typeEntry->second;
  } else {
    type = "unrecognised";
  }
  os << "\nType: " << type;
  os << "\nDepth: " << mesh.GetMeshDepth();
  os << "\nTranslation: " << mesh.GetTransform().getTranslation();
  os << "\nRotation: " << mesh.GetTransform().getRotation();
  if (mesh.GetMeshType() == G4Mesh::rectangle &&
      mesh.GetMeshDepth() == 3) {
    // Print ThreeDRectangleParameters
  }
  return os;
}
