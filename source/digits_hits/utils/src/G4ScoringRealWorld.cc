#include "G4ScoringRealWorld.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Region.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"

#include "G4SystemOfUnits.hh"

G4ScoringRealWorld::G4ScoringRealWorld(G4String lvName)
  :G4VScoringMesh(lvName)
{
  fShape = MeshShape::realWorldLogVol;
  logVolName = lvName;
  G4double size[] = {0.,0.,0.};
  SetSize(size);
  G4int nBin[] = {1,1,1};
  SetNumberOfSegments(nBin);
}

G4ScoringRealWorld::~G4ScoringRealWorld()
{
}

void G4ScoringRealWorld::List() const
{
  G4cout << "G4ScoringRealWorld : " << logVolName << G4endl;
  G4VScoringMesh::List();
}

void G4ScoringRealWorld::SetupGeometry(G4VPhysicalVolume* ) 
{
  auto store = G4LogicalVolumeStore::GetInstance();
  auto itr = store->begin();
  for(;itr!=store->end();itr++)
  {
    if((*itr)->GetName()==logVolName)
    {
      fMeshElementLogical = (*itr);  // Logical volume to score
      G4int nb = 0;
      auto pvStore = G4PhysicalVolumeStore::GetInstance();
      auto pvItr = pvStore->begin();
      for(;pvItr!=pvStore->end();pvItr++)
      {
        if((*pvItr)->GetLogicalVolume()==(*itr))
        { nb += (*pvItr)->GetMultiplicity(); }
      }
      G4int nBin[] = {nb,1,1};
      SetNumberOfSegments(nBin);
      // check if this logical volume belongs to the real world
      auto region = (*itr)->GetRegion();
      if(region && !(region->IsInMassGeometry()))
      {
        G4ExceptionDescription ed;
        ed << "Logical Volume with name <" << logVolName << "> is not used in the mass world.";
        G4Exception("G4ScoringRealWorld","SWV0001",FatalException,ed);
      }  
      // set the sensitive detector
      fMeshElementLogical->SetSensitiveDetector(fMFD);
      return;
    }
  }
  G4ExceptionDescription ed;
  ed << "Logical Volume with name <" << logVolName << "> is not found";
  G4Exception("G4ScoringRealWorld","SWV0000",FatalException,ed);
}


