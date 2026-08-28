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
//-------------------------------------------------------------------//
// Unit test for G4Cache
//-------------------------------------------------------------------//
#include "G4Cache.hh"

#include <gtest/gtest.h>

#include <vector>

//---------------------------------------------------------------------------//
// "Fixtures" for tests
struct A
{
    A() = default;
    A(int _a) : a(_a) {}
    int a{0};
};

struct B
{
    int b{1};
};

G4ThreadFunReturnType myfunc(G4ThreadFunArgType val)
{
  auto& asc = *((G4MapCache<int, double>*)val);
  auto tid = std::hash<G4Pid_t>{}(G4Threading::G4GetPidId()) % 100000;

  asc[1] = tid;
  asc[2] = 2.2;
  asc[30] = 3.3;
  EXPECT_EQ(asc[1], tid) << "Wrong first element";
  EXPECT_EQ(asc[2], 2.2) << "Wrong second element on thread " << tid;
  EXPECT_EQ(asc[30], 3.3) << "Wrong third element on thread " << tid;

  G4Cache<B> bb;
  EXPECT_EQ(bb.Get().b, 1) << "Wrong Cache content for default initialization on thread " << tid;
  return nullptr;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(G4Cache, Threading)
{
#ifndef G4MULTITHREADED
  GTEST_SKIP() << "Test only valid in multithreaded mode";
#endif

  int nthreads = 2;  // num threads
  G4Thread tid[2];  // = new G4Thread[nthreads];

  G4MapCache<int, double> aSharedCache;
  aSharedCache.Insert(1, 0.1);
  aSharedCache.Insert(2, 0.2);
  for (int idx = 0; idx < nthreads; ++idx)
  {
    G4Thread* tr = &tid[idx];
    G4THREADCREATE(tr, myfunc, &(aSharedCache));
  }

  for (int idx = 0; idx < nthreads; ++idx)
  {
    G4THREADJOIN(tid[idx]);
  }

  // Passing shared cache over to threads should not
  // have modified value here
  EXPECT_EQ(aSharedCache.Size(), 2) << "Wrong cache size after passing to threads";
  EXPECT_EQ(aSharedCache[1], 0.1) << "Wrong value for key 1 after using in threads";
  EXPECT_EQ(aSharedCache[2], 0.2) << "Wrong value for key 2 after using in threads";
}

TEST(G4Cache, SingleValue)
{
  G4Cache<double> vc(3.431);
  EXPECT_EQ(vc.Get(), 3.431) << "Cache not correctly constructed";

  double theV = vc.Get();
  theV += 1.5;
  vc.Put(theV);
  EXPECT_EQ(vc.Get(), 4.931) << "Cache value incorrectly modified by Put";

  ++vc.Get();
  EXPECT_EQ(vc.Get(), 5.931) << "Cache value incorrectly modified by direct access";
}

TEST(G4Cache, ObjectInstance)
{
  G4Cache<B> bb;
  EXPECT_EQ(bb.Get().b, 1) << "Wrong Cache content for default initialization";

  B ab;
  ab.b = 1234;
  G4Cache<B> abCache;
  abCache.Put(ab);
  G4Cache<B> aCopy(abCache);

  EXPECT_EQ(abCache.Get().b, aCopy.Get().b) << "Copied cache has different value";
  EXPECT_NE(&abCache.Get(), &aCopy.Get()) << "Copied object has same address as original";

  B anewb;
  anewb.b = 4321;
  aCopy.Put(anewb);
  EXPECT_EQ(aCopy.Get().b, 4321) << "Copy does not have new value";
  EXPECT_EQ(abCache.Get().b, 1234) << "After copy Original does not have correct value";

  // Copy assignment
  G4Cache<B> assB = abCache;
  EXPECT_EQ(assB.Get().b, abCache.Get().b) << "Assigned cache has different value";
  EXPECT_NE(&assB.Get(), &abCache.Get()) << "Copied object has same address as original";

  assB.Put(anewb);
  EXPECT_EQ(assB.Get().b, 4321) << "Assigned cache does not have new value";
  EXPECT_EQ(abCache.Get().b, 1234)
    << "After assignement Original does not have anymore correct value";
}

TEST(G4Cache, ObjectPointerValue)
{
  G4Cache<A*> ac;
  EXPECT_EQ(ac.Get(), nullptr) << "Held pointer is not nullptr on default construction";

  A* mya = nullptr;
  ac.Put(mya = new A(2));
  EXPECT_EQ(ac.Get(), mya) << "Wrong pointer returned by Get()";
  EXPECT_EQ(ac.Get()->a, 2) << "Wrong value for pointed-to object";

  mya->a = 1;
  EXPECT_EQ(ac.Get()->a, 1) << "Wrong value for pointed to object";
}

TEST(G4Cache, VectorValue)
{
  G4VectorCache<double> aV;
  EXPECT_EQ(aV.Size(), 0) << "Default constructed G4VectorCache not empty";

  aV.Push_back(1.01);
  EXPECT_EQ(aV.Size(), 1) << "Wrong size after adding an element";
  EXPECT_EQ(aV[0], 1.01) << "First element has wrong value";
  EXPECT_EQ(aV.Pop_back(), 1.01) << "Wrong value of popped element";
  EXPECT_EQ(aV.Size(), 0) << "Not empty after popping only element";

  aV.Push_back(3);
  aV[0] = 2;
  EXPECT_EQ(aV[0], 2) << "[] operator did not modify element correctly";
}

TEST(G4Cache, VectorValueArrayConstruction)
{
  double array[3] = {1.1, 2.2, 3.3};
  G4VectorCache<double> aV(3, array);
  // NB: Use of Matchers:
  // - https://google.github.io/googletest/reference/matchers.html
  // is better here, e.g.
  // EXPECT_THAT(aV, ::testing::Eq(array));
  // but needs inclusion/linking to GMock, and support in tested class
  EXPECT_EQ(aV[0], array[0]);
  EXPECT_EQ(aV[1], array[1]);
  EXPECT_EQ(aV[2], array[2]);

  // Check "stack" functionality
  aV.Push_back(4.4);
  EXPECT_EQ(aV.Pop_back(), 4.4) << "Wrong last element";
  EXPECT_EQ(aV.Pop_back(), 3.3) << "Wrong last element";
  EXPECT_EQ(aV.Pop_back(), 2.2) << "Wrong last element";
  EXPECT_EQ(aV.Pop_back(), 1.1) << "Wrong last element";
  EXPECT_EQ(aV.Size(), 0) << "Cache not empty after last element popped";
}

TEST(G4Cache, CollectionOfCache)
{
  G4Cache<A> anE;
  anE.Put(A(1234));

  std::vector<G4Cache<A>> aVec(10, anE);
  aVec.resize(20);
  for (int i = 0; i < 20; ++i)
  {
    const int val = i > 9 ? 0 : 1234;
    EXPECT_EQ(aVec[i].Get().a, val) << "Value in vector<G4Cache> wrong for element at index " << i;
  }

  G4Cache<A> anE2(5678);
  for (int i = 10; i < 20; ++i)
  {
    aVec[i] = anE2;
  }

  for (int i = 0; i < 20; ++i)
  {
    const int val = i > 9 ? 5678 : 1234;
    EXPECT_EQ(aVec[i].Get().a, val) << "Value in vector<G4Cache> wrong for element at index " << i;
  }
}

TEST(G4Cache, MapValue)
{
  G4MapCache<int, double> aM1;
  aM1[1] = 10.1;
  aM1[10] = 100.1;
  aM1[5] = 50.1;
  EXPECT_EQ(aM1.Size(), 3) << "Wrong map size";

  // Another use-case for Matchers....
  auto e = aM1.Insert(20, 200.1);
  EXPECT_TRUE(e.second) << "Did not insert new value correctly";
  EXPECT_EQ(aM1.Size(), 4) << "Wrong map size after insert";
  EXPECT_TRUE(aM1.Has(20)) << "New value not present in map";
  EXPECT_EQ(aM1[20], 200.1) << "Wrong value for key";

  EXPECT_TRUE((aM1.Begin()->first == 1) && (aM1.Begin()->second == 10.1)) << "Wrong head";
  EXPECT_NE(aM1.Find(5), aM1.End()) << "Held value not found";
  EXPECT_EQ(aM1.Find(111), aM1.End()) << "Bad value found";
  EXPECT_TRUE(aM1.Has(5)) << "Held value not found";
  EXPECT_EQ(aM1.Get(10), 100.1) << "Wrong value returned from Get";

  // Modifiers
  EXPECT_EQ(aM1.Erase(10), 1) << "Value not erased at expected position";
  EXPECT_FALSE(aM1.Has(10)) << "Erased value still found";
  aM1[20] = 199.9;
  EXPECT_EQ(aM1[20], 199.9) << "Value not modified correctly by operator[]";
}

TEST(G4Cache, MapOfPointersToObjects)
{
  G4MapCache<int, A*> mm;
  mm[0] = new A(10);
  mm[1] = new A(11);
  mm[2] = new A(12);
  EXPECT_EQ(mm[0]->a, 10) << "Wrong first element";
  EXPECT_EQ(mm[1]->a, 11) << "Wrong second element";
  EXPECT_EQ(mm[2]->a, 12) << "Wrong third element";
  // NB: MEMORY LEAK HERE
  // Cache does not own pointers...
}

TEST(G4Cache, UnusedCachesDoNotFail)
{
  // This tests that a regression identified as part of Issue #281:
  // 1. We construct N instances of G4Cache.
  //    - These have internal ids 0,...,N-1
  // 2. On construction, the backing store has size S=0
  // 3. We access the G4Cache with id=M
  //    - Backing store size is now S=M+1
  // 4. On destruction of a G4Cache instance, it calls
  //    G4CacheReference<V>::Destroy(unsigned int id, G4bool last)
  //    with id being its identity
  //    - If S < id a fatal G4exception is thrown
  // 5. Thus we can trigger this if:
  //    - S < N-1 => M < N-2
  //    - i.e. if we only accessed the G4Cache elements with id M < N-2
  //
  // This was resolved by always initializing the backing store for
  // a G4Cache instance on construction. The following is a regression test.

  // Unique struct to fully control cache size
  struct Fail
  {};
  // We put this in an assertion just to mark that we don't expect a failure
  // TODO: Need a dedicated G4ExceptionHandler to check for G4Exceptions...
  EXPECT_NO_THROW({
    // Create N cache elements, ids 0, ..., N-1
    constexpr size_t N = 8;
    G4Cache<Fail> vars[N];
    // Access cache value with id = N-3 (max element that can trigger error)
    vars[N - 3].Get();
  });
}

TEST(G4Cache, StaticCachesDoNotFailOnExit)
{
  // Another regression/death test.
  // Generic biasing uses _static_ G4Cache instances, and these may be destructed
  // on program exit _after_ the G4CacheReference static backing store is. Current
  // implementation is partially safe given manual raw memory management. This test
  // reproduces errors seen when trying to move to RAII methods.
  //
  // Here, the backing store will already be empty/cleared when the G4Cache destructor
  // (which calls Destroy) is run.
  // It will cause an fatal G4Exception on _program exit_ if present as the id (0...3).
  // will be less than the store size of 0 at this point.
  static G4Cache<bool> a;
  static G4Cache<bool> b;
  static G4Cache<bool> c;
  static G4Cache<bool> d;
}
