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
// Tests that G4GeometryManager::CloseGeometry() itself - not just a direct
// call to a checker's CheckGeometry()/CheckAndMark() entry point - now runs
// the parameterisation-reassignment consistency pass by default, and that
// UseNestedParameterisationCheckAtClose(true) really dispatches through the
// navigation module's registered G4ParameterisationCheckVisitor rather than
// the flat G4ParameterisationChecker.
//
// One geometry, two independent parameterised branches:
//  - FlatSlab: an ordinary (non-nested) parameterisation whose
//    ComputeTransformation() reassigns its logical volume for one copy.
//    The flat checker exercises ComputeTransformation(), so this must be
//    caught by the DEFAULT pass, with no explicit enable call at all.
//  - NestedBox: nested inside a nested-hosting parameterisation, whose
//    ComputeMaterial() reassigns its logical volume for one copy. The flat
//    checker never probes ComputeMaterial() for a nested parameterisation
//    (it has no live parent touchable to offer), so this must NOT be caught
//    by the first (default, non-nested) CloseGeometry() pass, but MUST be
//    caught once UseNestedParameterisationCheckAtClose(true) is set and
//    CloseGeometry() is called again.
//
// Author: John Apostolakis (CERN), August 2026
// --------------------------------------------------------------------

#include <cassert>

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VNestedParameterisation.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

namespace
{
G4Material *darkMaterial, *brightMaterial, *defaultMaterial;
G4LogicalVolume* flatRebindLV = nullptr;
G4LogicalVolume* nestedRebindLV = nullptr;
G4bool doNestedRebind = false;  // toggled between the two CloseGeometry() passes

// Ordinary (non-nested) parameterisation: reassigns its logical volume for
// copy 0 unconditionally, from ComputeTransformation() - exactly the hook
// the flat G4ParameterisationChecker always exercises.
class FlatParam : public G4VPVParameterisation
{
  public:

    FlatParam(G4int numCopies, G4double xFullWidth, G4double yFullWidth, G4double zFullWidth)
      : fNumCopies(numCopies), fXfullWidth(xFullWidth), fYfullWidth(yFullWidth), fZfullWidth(zFullWidth)
    {}

    void ComputeTransformation(const G4int n, G4VPhysicalVolume* pRep) const override
    {
      pRep->SetTranslation(G4ThreeVector((n - ((fNumCopies - 1.0) / 2.)) * fXfullWidth, 0., 0.));
      if (n == 0)
      {
        pRep->SetLogicalVolume(flatRebindLV);
      }
    }

    void ComputeDimensions(G4Box& pBox, const G4int, const G4VPhysicalVolume*) const override
    {
      pBox.SetXHalfLength(fXfullWidth * 0.5);
      pBox.SetYHalfLength(fYfullWidth * 0.5);
      pBox.SetZHalfLength(fZfullWidth * 0.5);
    }

  private:

    G4int fNumCopies;
    G4double fXfullWidth, fYfullWidth, fZfullWidth;
};

// First level, hosting the nested parameterisation below it (plain, not
// itself reassigning anything).
class TopParam : public G4VPVParameterisation
{
  public:

    TopParam(G4int numRows, G4double xFullWidth, G4double yFullWidth, G4double zFullWidth)
      : fNumRows(numRows), fXfullWidth(xFullWidth), fYfullWidth(yFullWidth), fZfullWidth(zFullWidth)
    {}

    void ComputeTransformation(const G4int n, G4VPhysicalVolume* pRep) const override
    {
      pRep->SetTranslation(G4ThreeVector((n - ((fNumRows - 1.0) / 2.)) * fXfullWidth, 0., 0.));
    }

    void ComputeDimensions(G4Box& pBox, const G4int, const G4VPhysicalVolume*) const override
    {
      pBox.SetXHalfLength(fXfullWidth * 0.5);
      pBox.SetYHalfLength(fYfullWidth * 0.5);
      pBox.SetZHalfLength(fZfullWidth * 0.5);
    }

  private:

    G4int fNumRows;
    G4double fXfullWidth, fYfullWidth, fZfullWidth;
};

// Nested inside TopParam's copies. Reassigns its logical volume for copy 0
// from ComputeMaterial(), the hook the flat checker skips for nested
// parameterisations.
class NestedParam : public G4VNestedParameterisation
{
  public:

    NestedParam(G4int numCols, G4double yHalfWidth, G4double xHalfWidth, G4double zHalfWidth)
      : fNumCols(numCols),
        fYHalfWidth(yHalfWidth),
        fYFullWidth(yHalfWidth * 2.0),
        fXHalfWidth(xHalfWidth),
        fZHalfWidth(zHalfWidth)
    {}

    void ComputeTransformation(const G4int n, G4VPhysicalVolume* pRep) const override
    {
      pRep->SetTranslation(G4ThreeVector(0., (n - ((fNumCols - 1) / 2.)) * fYFullWidth, 0.));
      pRep->SetRotation(nullptr);
    }

    G4Material* ComputeMaterial(G4VPhysicalVolume* currentVol, const G4int repNo,
                                const G4VTouchable* parentTouch) override
    {
      if (parentTouch == nullptr)
      {
        G4Exception("NestedParam::ComputeMaterial()", "Null parent TouchHist", FatalException,
                   "Null pointer as parent touchable pointer.");
      }

      const G4int parentReplicaNo = parentTouch->GetReplicaNumber();
      const G4bool odd = ((repNo + parentReplicaNo) % 2) == 1;
      G4Material* material = odd ? darkMaterial : brightMaterial;

      if (doNestedRebind && (repNo == 0))
      {
        currentVol->SetLogicalVolume(nestedRebindLV);
        nestedRebindLV->SetMaterial(material);
      }
      else
      {
        currentVol->GetLogicalVolume()->SetMaterial(material);
      }
      return material;
    }

    G4int GetNumberOfMaterials() const override { return 2; }

    G4Material* GetMaterial(G4int idx) const override
    {
      return ((idx % 2) == 0) ? brightMaterial : darkMaterial;
    }

    void ComputeDimensions(G4Box& pBox, const G4int, const G4VPhysicalVolume*) const override
    {
      pBox.SetXHalfLength(fXHalfWidth);
      pBox.SetYHalfLength(fYHalfWidth);
      pBox.SetZHalfLength(fZHalfWidth);
    }

  private:

    G4int fNumCols;
    G4double fYHalfWidth, fYFullWidth;
    G4double fXHalfWidth, fZHalfWidth;
};

// Builds: World -> FlatMother(placement) -> FlatSlab (FlatParam)
//               -> NestedMother(placement) -> Slab (TopParam) -> Box (NestedParam)
G4VPhysicalVolume* BuildGeometry()
{
  G4double a, fractionmass, density;
  G4int z, ncomponents;

  auto* N = new G4Element("Nitrogen", "N", z = 7, a = 14.01 * g / mole);
  auto* O = new G4Element("Oxygen", "O", z = 8, a = 16.00 * g / mole);
  auto* Air = new G4Material("CGDPC_Air", density = 1.290 * mg / cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);
  defaultMaterial = Air;

  darkMaterial = new G4Material("CGDPC_Dark", z = 82, a = 207.19 * g / mole,
                                density = 11.35 * g / cm3);
  brightMaterial = new G4Material("CGDPC_Bright", z = 13, a = 26.98 * g / mole,
                                  density = 2.70 * g / cm3);
  auto* Cu = new G4Material("CGDPC_Cu", z = 29, a = 63.55 * g / mole, density = 8.96 * g / cm3);

  const G4double worldHalf = 500. * cm;
  auto* worldBox = new G4Box("CGDPC_World", worldHalf, worldHalf, worldHalf);
  auto* worldLog = new G4LogicalVolume(worldBox, defaultMaterial, "CGDPC_WorldLV");
  auto* worldPhys = new G4PVPlacement(nullptr, G4ThreeVector(), worldLog, "CGDPC_World",
                                      nullptr, false, 0);

  // --- Flat branch ---------------------------------------------------
  const G4int numFlat = 4;
  const G4double xFlatHalf = 10. * cm, yFlatHalf = 100. * cm, zFlatHalf = 100. * cm;

  auto* flatMotherBox = new G4Box("CGDPC_FlatMother", xFlatHalf * numFlat, yFlatHalf, zFlatHalf);
  auto* flatMotherLog =
    new G4LogicalVolume(flatMotherBox, defaultMaterial, "CGDPC_FlatMotherLV");
  new G4PVPlacement(nullptr, G4ThreeVector(-200. * cm, 0., 0.), flatMotherLog, "CGDPC_FlatMother",
                    worldLog, false, 0);

  auto* flatBox = new G4Box("CGDPC_FlatSlab", xFlatHalf, yFlatHalf, zFlatHalf);
  auto* flatLog = new G4LogicalVolume(flatBox, defaultMaterial, "CGDPC_FlatSlabLV");
  auto* flatParam = new FlatParam(numFlat, xFlatHalf * 2., yFlatHalf * 2., zFlatHalf * 2.);
  new G4PVParameterised("CGDPC_FlatSlab", flatLog, flatMotherLog, kXAxis, numFlat, flatParam);

  auto* flatRebindBox = new G4Box("CGDPC_FlatRebindBox", xFlatHalf, yFlatHalf, zFlatHalf);
  flatRebindLV = new G4LogicalVolume(flatRebindBox, Cu, "CGDPC_FlatRebindLV");

  // --- Nested branch ---------------------------------------------------
  const G4int numSlabs = 4;
  const G4double xSlabHalf = 10. * cm, yTopHalf = 100. * cm, zTopHalf = 100. * cm;

  auto* nestedMotherBox = new G4Box("CGDPC_NestedMother", xSlabHalf * numSlabs, yTopHalf, zTopHalf);
  auto* nestedMotherLog =
    new G4LogicalVolume(nestedMotherBox, defaultMaterial, "CGDPC_NestedMotherLV");
  new G4PVPlacement(nullptr, G4ThreeVector(200. * cm, 0., 0.), nestedMotherLog,
                    "CGDPC_NestedMother", worldLog, false, 0);

  auto* topBox = new G4Box("CGDPC_Top", xSlabHalf * numSlabs, yTopHalf, zTopHalf);
  auto* topLog = new G4LogicalVolume(topBox, defaultMaterial, "CGDPC_TopLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), topLog, "CGDPC_Top", nestedMotherLog, false, 0);

  auto* slabBox = new G4Box("CGDPC_Slab", xSlabHalf, yTopHalf, zTopHalf);
  auto* slabLog = new G4LogicalVolume(slabBox, defaultMaterial, "CGDPC_SlabLV");
  auto* topParam = new TopParam(numSlabs, xSlabHalf * 2., yTopHalf * 2., zTopHalf * 2.);
  new G4PVParameterised("CGDPC_Slab", slabLog, topLog, kXAxis, numSlabs, topParam);

  const G4int numBoxesY = 4;
  const G4double xBoxHalf = xSlabHalf, yBoxHalf = yTopHalf / numBoxesY, zBoxHalf = zTopHalf;

  auto* smallBox = new G4Box("CGDPC_SmallBox", xBoxHalf, yBoxHalf, zBoxHalf);
  auto* smallLog = new G4LogicalVolume(smallBox, defaultMaterial, "CGDPC_SmallLV");
  auto* nestedParam = new NestedParam(numBoxesY, yBoxHalf, xBoxHalf, zBoxHalf);
  new G4PVParameterised("CGDPC_Nested", smallLog, slabLog, kYAxis, numBoxesY, nestedParam);

  auto* nestedRebindBox = new G4Box("CGDPC_NestedRebindBox", xBoxHalf, yBoxHalf, zBoxHalf);
  nestedRebindLV = new G4LogicalVolume(nestedRebindBox, Cu, "CGDPC_NestedRebindLV");

  return worldPhys;
}
}  // namespace

int main()
{
  doNestedRebind = true;  // the nested branch reassigns from copy zero onward
  G4VPhysicalVolume* world = BuildGeometry();

  G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking()
    ->SetWorldVolume(world);

  assert(!flatRebindLV->IsMaterialMutable());
  assert(!nestedRebindLV->IsMaterialMutable());

  // --- Pass 1: default settings, no explicit enable call at all ---------
  // The basic (flat) pass must now run without CheckParameterisedReassignmentAtClose()
  // ever being called: it is on by default. It exercises ComputeTransformation(),
  // so it catches the FLAT branch's reassignment...
  G4GeometryManager::GetInstance()->CloseGeometry(false, false, world);

  assert(flatRebindLV->IsMaterialMutable());
  // ...but the nested option was never requested, so the flat checker's
  // documented gap applies: it never probes ComputeMaterial() for a nested
  // parameterisation, so it cannot see this reassignment.
  assert(!nestedRebindLV->IsMaterialMutable());

  // --- Pass 2: the higher-level 'ON' option selects the nested-aware pass
  G4GeometryManager::GetInstance()->OpenGeometry(world);
  G4GeometryManager::GetInstance()->UseNestedParameterisationCheckAtClose(true);
  G4GeometryManager::GetInstance()->CloseGeometry(false, false, world);

  assert(nestedRebindLV->IsMaterialMutable());

  G4GeometryManager::GetInstance()->UseNestedParameterisationCheckAtClose(false);
  G4GeometryManager::GetInstance()->OpenGeometry(world);

  G4cout << "testG4CloseGeometryDefaultParamCheck: PASSED" << G4endl;
  return 0;
}
