#include "G4Box.hh"
#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4GeomTestVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4MTRunManager.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPVParameterisation.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <cstdio>
#include <cstdlib>

class BoxParameterisation final : public G4VPVParameterisation
{
  public:

    explicit BoxParameterisation(G4bool overlapping) : fPitch(overlapping ? 12. * mm : 30. * mm) {}

    void ComputeTransformation(G4int copyNo, G4VPhysicalVolume* volume) const override
    {
      volume->SetTranslation({(copyNo - 2) * fPitch, 0., 0.});
      volume->SetRotation(nullptr);
    }

    void ComputeDimensions(G4Box& box, G4int copyNo, const G4VPhysicalVolume*) const override
    {
      box.SetXHalfLength((8. + copyNo) * mm);
      box.SetYHalfLength(8. * mm);
      box.SetZHalfLength(8. * mm);
    }

  private:

    G4double fPitch;
};

class Detector final : public G4VUserDetectorConstruction
{
  public:

    explicit Detector(G4bool overlapping) : fOverlapping(overlapping) {}

    G4VPhysicalVolume* Construct() override
    {
      auto* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      auto* worldLogical = new G4LogicalVolume(
        new G4Box("worldSolid", 100. * cm, 100. * cm, 100. * cm), material, "worldLogical");
      fWorld = new G4PVPlacement(nullptr, {}, worldLogical, "world", nullptr, false, 0);

      auto* containerLogical = new G4LogicalVolume(
        new G4Box("containerSolid", 20. * cm, 20. * cm, 20. * cm), material, "containerLogical");
      new G4PVPlacement(nullptr, {}, containerLogical, "container", worldLogical, false, 0);

      auto* cellLogical = new G4LogicalVolume(new G4Box("cellSolid", 8. * mm, 8. * mm, 8. * mm),
                                              material, "cellLogical");
      fParameterised = new G4PVParameterised("cell", cellLogical, containerLogical, kXAxis, 5,
                                             new BoxParameterisation(fOverlapping));

      auto* markerLogical = new G4LogicalVolume(new G4Box("markerSolid", 2. * cm, 2. * cm, 2. * cm),
                                                material, "markerLogical");
      for (G4int i = 0; i < 4; ++i)
      {
        new G4PVPlacement(nullptr, {-60. * cm + i * 30. * cm, 60. * cm, 0.}, markerLogical,
                          "marker", worldLogical, false, i);
      }
      return fWorld;
    }

    G4VPhysicalVolume* GetWorld() const { return fWorld; }
    G4PVParameterised* GetParameterised() const { return fParameterised; }

  private:

    G4bool fOverlapping;
    G4VPhysicalVolume* fWorld = nullptr;
    G4PVParameterised* fParameterised = nullptr;
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

int main(int argc, char** argv)
{
  const G4bool overlapping = argc > 1 && std::atoi(argv[1]) != 0;
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT, 2);
  if (dynamic_cast<G4MTRunManager*>(runManager) == nullptr)
  {
    return 10;
  }

  auto* detector = new Detector(overlapping);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Physics);
  runManager->SetUserInitialization(new Actions);
  runManager->InitializeGeometry();

  std::fprintf(stderr, "parameterised check before workers (%s)\n",
               overlapping ? "overlap" : "clear");
  G4GeomTestVolume before(detector->GetWorld(), 0., 1000, false);
  before.TestOverlapInTree();
  G4bool detectedBefore = detector->GetParameterised()->CheckOverlaps(1000, 0., false, 1);

  runManager->Initialize();
  std::fprintf(stderr, "parameterised check after workers (%s)\n",
               overlapping ? "overlap" : "clear");
  G4GeomTestVolume after(detector->GetWorld(), 0., 1000, false);
  after.TestOverlapInTree();
  after.TestRecursiveOverlap();
  G4bool detectedAfter = detector->GetParameterised()->CheckOverlaps(1000, 0., false, 1);

  delete runManager;
  if (detectedBefore != overlapping || detectedAfter != overlapping)
  {
    std::fprintf(stderr, "unexpected overlap result: before=%d after=%d\n", detectedBefore,
                 detectedAfter);
    return 20;
  }
  return 0;
}
