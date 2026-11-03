#include "G4Box.hh"
#include "G4GeomTestVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreadPool.hh"

#include "PTL/ThreadData.hh"

#include <cstdio>

int main()
{
  std::fprintf(stderr, "construct\n");
  auto* worldSolid = new G4Box("world", 100 * cm, 100 * cm, 100 * cm);
  auto* worldLogical = new G4LogicalVolume(worldSolid, nullptr, "world");
  auto* world = new G4PVPlacement(nullptr, {}, worldLogical, "world", nullptr, false, 0);

  auto* childSolid = new G4Box("child", 10 * cm, 10 * cm, 10 * cm);
  auto* childLogical = new G4LogicalVolume(childSolid, nullptr, "child");
  new G4PVPlacement(nullptr, {-30 * cm, 0, 0}, childLogical, "child0", worldLogical, false, 0);
  new G4PVPlacement(nullptr, {0, 0, 0}, childLogical, "child1", worldLogical, false, 1);
  new G4PVPlacement(nullptr, {30 * cm, 0, 0}, childLogical, "child2", worldLogical, false, 2);

  G4GeomTestVolume test(world, 0., 2000, false);
  std::fprintf(stderr, "direct\n");
  auto* placement = static_cast<G4PVPlacement*>(worldLogical->GetDaughter(0));
  auto overlap = placement->CheckOverlaps(2000, 0., false, 1);
  std::fprintf(stderr, "tree\n");
  test.TestOverlapInTree();
  std::fprintf(stderr, "recursive\n");
  G4ThreadPool::Config config;
  config.pool_size = 2;
  auto* previousThreadData = PTL::ThreadData::GetInstance();
  {
    G4ThreadPool pool(config);
    test.TestRecursiveOverlap();
    pool.destroy_threadpool();
  }
  auto* taskThreadData = PTL::ThreadData::GetInstance();
  if (taskThreadData != previousThreadData)
  {
    delete taskThreadData;
  }
  PTL::ThreadData::GetInstance() = previousThreadData;

  auto* overlapping =
    new G4PVPlacement(nullptr, {5 * cm, 0, 0}, childLogical, "overlapping", worldLogical, false, 3);
  if (!overlapping->CheckOverlaps(2000, 0., false, 1))
  {
    return 2;
  }

  std::fprintf(stderr, "done\n");
  return overlap ? 1 : 0;
}
