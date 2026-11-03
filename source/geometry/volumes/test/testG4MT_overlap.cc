#include "G4Box.hh"
#include "G4Geantino.hh"
#include "G4GeomTestVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4MTRunManager.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <cstdio>

class Detector final : public G4VUserDetectorConstruction
{
  public:

    G4VPhysicalVolume* Construct() override
    {
      auto* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* worldSolid = new G4Box("world", 100 * cm, 100 * cm, 100 * cm);
      auto* worldLogical = new G4LogicalVolume(worldSolid, material, "worldLogical");
      fWorld = new G4PVPlacement(nullptr, {}, worldLogical, "world", nullptr, false, 0);

      auto* childSolid = new G4Box("child", 10 * cm, 10 * cm, 10 * cm);
      auto* childLogical = new G4LogicalVolume(childSolid, material, "childLogical");
      new G4PVPlacement(nullptr, {-30 * cm, 0, 0}, childLogical, "child0", worldLogical, false, 0);
      new G4PVPlacement(nullptr, {0, 0, 0}, childLogical, "child1", worldLogical, false, 1);
      new G4PVPlacement(nullptr, {30 * cm, 0, 0}, childLogical, "child2", worldLogical, false, 2);
      return fWorld;
    }

    G4VPhysicalVolume* GetWorld() const { return fWorld; }

  private:

    G4VPhysicalVolume* fWorld = nullptr;
};

class Physics final : public G4VUserPhysicsList
{
  protected:

    void ConstructParticle() override { G4Geantino::GeantinoDefinition(); }
    void ConstructProcess() override { AddTransportation(); }
    void SetCuts() override { SetCutsWithDefault(); }
};

class Actions final : public G4VUserActionInitialization
{
  public:

    class Primary final : public G4VUserPrimaryGeneratorAction
    {
      public:

        void GeneratePrimaries(G4Event*) override {}
    };

    void Build() const override { SetUserAction(new Primary); }
};

int main()
{
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT, 2);
  if (dynamic_cast<G4MTRunManager*>(runManager) == nullptr)
  {
    return 10;
  }

  auto* detector = new Detector;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Physics);
  runManager->SetUserInitialization(new Actions);

  runManager->InitializeGeometry();
  std::fprintf(stderr, "MT check before worker initialization\n");
  G4GeomTestVolume before(detector->GetWorld(), 0., 2000, false);
  before.TestOverlapInTree();

  runManager->Initialize();
  std::fprintf(stderr, "MT check after worker initialization\n");
  G4GeomTestVolume after(detector->GetWorld(), 0., 2000, false);
  after.TestOverlapInTree();
  after.TestRecursiveOverlap();

  delete runManager;
  return 0;
}
