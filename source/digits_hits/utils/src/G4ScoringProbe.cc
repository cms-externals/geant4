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
// G4ScoringProbe
// --------------------------------------------------------------------

#include "G4ScoringProbe.hh"

#include "G4AutoLock.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UImanager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoreColorMap.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

namespace
{
G4Mutex logvolmutex = G4MUTEX_INITIALIZER;
}

G4ScoringProbe::G4ScoringProbe(const G4String& lvName, G4double half_size, G4bool checkOverlap)
  : G4VScoringMesh(lvName),
    chkOverlap(checkOverlap),
    layeredMaterialName("none"),
    layeredMaterial(nullptr)
{
  fShape = MeshShape::probe;
  logVolName = lvName;
  probeSize = half_size;
  G4double hs[] = {half_size, half_size, half_size};
  SetSize(hs);
  G4int nBin[] = {1, 1, 1};
  SetNumberOfSegments(nBin);
  regName = lvName + "_region";
  if (G4Threading::IsMasterThread())
  {
    new G4Region(regName);
  }
  fDivisionAxisNames[0] = "x";
  fDivisionAxisNames[1] = "y";
  fDivisionAxisNames[2] = "z";
}

void G4ScoringProbe::List() const
{
  G4cout << "G4ScoringProbe : " << logVolName << G4endl;
  std::size_t np = posVec.size();
  for (std::size_t i = 0; i < np; ++i)
  {
    G4cout << " >> probe #" << i << " at " << posVec[i] << G4endl;
  }
  G4VScoringMesh::List();
}

void G4ScoringProbe::Draw(RunScore* map, G4VScoreColorMap* colorMap, G4int)
{
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager != nullptr)
  {
    // delete transient objects stored in the viewer
    if (refreshDraw)
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/clearTransients");
    }
    refreshDraw = true;

    // search max. & min. values for probes
    G4double min = DBL_MAX;
    G4double max = 0.;

    // get scored value
    std::vector<G4double> scored_vals;
    auto itr = map->GetMap()->begin();
    for (; itr != map->GetMap()->end(); itr++)
    {
      G4double val = (itr->second->sum_wx()) / fDrawUnitValue;
      scored_vals.push_back(val);
      if (min > val) min = val;
      if (max < val) max = val;
    }
    if (colorMap->IfFloatMinMax() && drawColorChart >= 0)
    {
      colorMap->SetMinMax(min, max);
    }

    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);

    pVisManager->BeginDraw();

    G4double c[4];
    std::size_t np = posVec.size();
    for (std::size_t i = 0; i < np; ++i)
    {
      colorMap->GetMapColor(scored_vals[i], c);
      att.SetColour(c[0], c[1], c[2]);

      G4Box aBox(logVolName + "_sc", probeSize, probeSize, probeSize);

      G4Polyhedron* poly = aBox.GetPolyhedron();

      G4Translate3D trans = G4Translate3D(posVec[i]);
      poly->Transform(trans);
      poly->SetVisAttributes(&att);
      pVisManager->Draw(*poly);
    }

    pVisManager->EndDraw();

    if (drawColorChart >= 0)
    {
      colorMap->SetOffset(drawColorChart);
      colorMap->SetPSUnit(fDrawUnit);
      colorMap->SetPSName(fWorldName, fDrawPSName);
      colorMap->DrawColorChart();
    }
    drawColorChart = 0;

    // refresh the viewer
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/refresh");
  }
}

void G4ScoringProbe::SetupGeometry(G4VPhysicalVolume* worldPhys)
{
  if (G4Threading::IsMasterThread())
  {
    auto worldLog = worldPhys->GetLogicalVolume();
    auto region = G4RegionStore::GetInstance()->GetRegion(regName);
    assert(region != nullptr);
    region->AddRootLogicalVolume(worldLog);
    region->SetWorld(worldPhys);

    auto boxSolid = new G4Box(logVolName + "_solid", probeSize, probeSize, probeSize);
    fMeshElementLogical = new G4LogicalVolume(boxSolid, layeredMaterial, logVolName + "_log");

    std::size_t np = posVec.size();
    for (std::size_t i = 0; i < np; ++i)
    {
      new G4PVPlacement(nullptr, posVec[i], fMeshElementLogical, logVolName + "_phy", worldLog,
                        false, (G4int)i, chkOverlap);
    }

    auto wisatt = new G4VisAttributes(G4Colour(.5, .5, .5));
    wisatt->SetVisibility(false);
    worldLog->SetVisAttributes(wisatt);
    auto visatt = new G4VisAttributes(G4Colour(.5, .5, .5));
    visatt->SetVisibility(true);
    fMeshElementLogical->SetVisAttributes(visatt);
  }
  else
  {
    G4AutoLock l(&logvolmutex);
    fMeshElementLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(logVolName, false);
    assert(fMeshElementLogical != nullptr);
    l.unlock();
  }

  fMeshElementLogical->SetSensitiveDetector(fMFD);
}

G4bool G4ScoringProbe::SetMaterial(G4String val)
{
  if (val == "none")
  {
    layeredMaterialName = val;
    layeredMassFlg = false;
    layeredMaterial = nullptr;
  }
  else
  {
    G4AutoLock l(&logvolmutex);
    auto mat = G4NistManager::Instance()->FindOrBuildMaterial(val);
    if (mat == nullptr)
    {
      return false;
    }
    layeredMaterialName = val;
    layeredMassFlg = true;
    layeredMaterial = mat;
    if (G4Threading::IsMasterThread())
    {
      auto region = G4RegionStore::GetInstance()->GetRegion(regName);
      assert(region != nullptr);
      region->UpdateMaterialList();
    }
    l.unlock();
  }
  return true;
}
