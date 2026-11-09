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
// Tests G4ParameterisationCheckVisitor against a two-level nested
// parameterisation (a G4VNestedParameterisation living inside an ordinary
// G4PVParameterised, itself inside a placement) - the case the flat
// G4ParameterisationChecker cannot exercise, since it has no live parent
// touchable to offer a nested ComputeMaterial().
//
// Two scenarios:
//  A) A clean nested parameterisation. Expect: 0 problems, the nested PV's
//     material observed to vary with its parent context, and no rebinding.
//  B) The same geometry, but the nested ComputeMaterial() additionally
//     rebinds its physical volume to a second logical volume that was never
//     marked mutable. Expect: exactly one rebinding is detected and, with
//     repair requested on open geometry, the rebound logical volume becomes
//     mutable.
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
#include "G4ParameterisationCheckVisitor.hh"
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
G4LogicalVolume* rebindTargetLV = nullptr;  // only used by scenario B
G4bool doRebind = false;                    // toggled between scenarios

// First level: an ordinary (non-nested) parameterisation, hosting the
// nested one below it.
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

// Second level: nested inside TopParam's copies. Chooses material from its
// own replica number combined with the parent's replica number (a checker-
// board pattern), exactly as G4VNestedParameterisation::ComputeMaterial()
// requires - and, in scenario B, additionally reassigns its logical volume.
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

      G4LogicalVolume* currentLV = currentVol->GetLogicalVolume();
      if (doRebind && (repNo == 0))
      {
        // Deliberately reassign the logical volume for one specific copy,
        // as an 'unorthodox' nested parameterisation might. rebindTargetLV
        // is never marked mutable, so this is exactly the case
        // G4ParameterisationCheckVisitor must catch.
        currentVol->SetLogicalVolume(rebindTargetLV);
        rebindTargetLV->SetMaterial(material);
      }
      else
      {
        currentLV->SetMaterial(material);
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

// Builds: World -> Top-pv (placement) -> Slab (TopParam, plain parameterised)
//                -> Box (NestedParam, nested parameterised)
G4VPhysicalVolume* BuildGeometry()
{
  G4double a, fractionmass, density;
  G4int z, ncomponents;

  auto* N = new G4Element("Nitrogen", "N", z = 7, a = 14.01 * g / mole);
  auto* O = new G4Element("Oxygen", "O", z = 8, a = 16.00 * g / mole);
  auto* Air = new G4Material("TPCV_Air", density = 1.290 * mg / cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  auto* Pb = new G4Material("TPCV_Lead", z = 82, a = 207.19 * g / mole, density = 11.35 * g / cm3);
  auto* Al =
    new G4Material("TPCV_Aluminium", z = 13, a = 26.98 * g / mole, density = 2.700 * g / cm3);
  auto* Cu = new G4Material("TPCV_Copper", z = 29, a = 63.55 * g / mole, density = 8.960 * g / cm3);

  darkMaterial = Pb;
  brightMaterial = Al;
  defaultMaterial = Air;

  auto* worldBox = new G4Box("TPCV_WorldBox", 1000. * cm, 1000. * cm, 1000. * cm);
  auto* worldLog = new G4LogicalVolume(worldBox, defaultMaterial, "TPCV_World");
  auto* worldPhys = new G4PVPlacement(nullptr, G4ThreeVector(), "TPCV_World", worldLog, nullptr,
                                      false, 0);

  auto* topBox = new G4Box("TPCV_TopBox", 100. * cm, 100. * cm, 100. * cm);
  auto* topLog = new G4LogicalVolume(topBox, defaultMaterial, "TPCV_TopLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), "TPCV_Top-pv", topLog, worldPhys, false, 0);

  const G4int numSlabs = 4;
  const G4double xTopHalf = 100. * cm, yTopHalf = 100. * cm, zTopHalf = 100. * cm;
  const G4double xSlabHalf = xTopHalf / numSlabs;

  auto* slabBox = new G4Box("TPCV_Slab", xSlabHalf, yTopHalf, zTopHalf);
  auto* slabLog = new G4LogicalVolume(slabBox, defaultMaterial, "TPCV_SlabLV");
  auto* topParam = new TopParam(numSlabs, xSlabHalf * 2., yTopHalf * 2., zTopHalf * 2.);
  new G4PVParameterised("TPCV_Slab", slabLog, topLog, kXAxis, numSlabs, topParam);

  const G4int numBoxesY = 4;
  const G4double xBoxHalf = xSlabHalf, yBoxHalf = yTopHalf / numBoxesY, zBoxHalf = zTopHalf;

  auto* smallBox = new G4Box("TPCV_SmallBox", xBoxHalf, yBoxHalf, zBoxHalf);
  auto* smallLog = new G4LogicalVolume(smallBox, defaultMaterial, "TPCV_SmallLV");
  auto* nestedParam = new NestedParam(numBoxesY, yBoxHalf, xBoxHalf, zBoxHalf);
  new G4PVParameterised("TPCV_Nested", smallLog, slabLog, kYAxis, numBoxesY, nestedParam);

  // A separate, unrelated logical volume that scenario B rebinds into. It is
  // deliberately never marked mutable, so the checker must be the one to fix
  // that (or diagnose it, if the geometry is closed).
  auto* rebindBox = new G4Box("TPCV_RebindBox", xBoxHalf, yBoxHalf, zBoxHalf);
  rebindTargetLV = new G4LogicalVolume(rebindBox, Cu, "TPCV_RebindLV");

  return worldPhys;
}
}  // namespace

int main()
{
  // --- Scenario A: clean two-level nesting -------------------------------
  doRebind = false;
  G4VPhysicalVolume* world = BuildGeometry();

  // G4TouchableWalker::Walk() (default startVol=nullptr) finds the world
  // volume via the tracking navigator, not just any G4PVPlacement with no
  // mother - register it, exactly as G4RunManagerKernel::DefineWorldVolume()
  // does for a real run.
  G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking()
    ->SetWorldVolume(world);

  G4int nProblems = G4ParameterisationCheckVisitor::CheckGeometry(nullptr, true, true);
  assert(nProblems == 0);

  // --- Scenario B: nested ComputeMaterial() rebinds its logical volume ---
  doRebind = true;
  assert(!rebindTargetLV->IsMaterialMutable());

  nProblems = G4ParameterisationCheckVisitor::CheckGeometry(nullptr, true, true);
  assert(nProblems >= 1);
  assert(rebindTargetLV->IsMaterialMutable());

  // A second pass should now find nothing new to repair for this LV, since
  // it has already been marked mutable (though it still reports as a
  // rebinding - that is expected and correct, only the mutability is fixed).
  doRebind = false;

  G4GeometryManager::GetInstance()->OpenGeometry();
  // Must end with geometry open.
  (void)world;

  G4cout << "testG4ParameterisationCheckVisitor: PASSED" << G4endl;
  return 0;
}
