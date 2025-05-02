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
// Unit test for G4ReferenceCountedHandle
//-------------------------------------------------------------------//
#include "G4ReferenceCountedHandle.hh"

#include "G4Threading.hh"

#include <gtest/gtest.h>

#include <string>

//---------------------------------------------------------------------------//
// "Fixtures" for tests

class TesterBase
{
  public:
    TesterBase() = default;
    virtual ~TesterBase() = default;
    virtual int Value() const = 0;
    virtual std::string AsString() const = 0;
};

class TesterString : public TesterBase
{
  public:
    TesterString(const std::string& str) : fData(str) {}

    int Value() const override { return (int)fData.size(); }
    std::string AsString() const override { return fData; }

  private:
    std::string fData;
};

class TesterInt : public TesterBase
{
  public:
    TesterInt(int i = 0) : fData(i) {}
    int Value() const override { return fData; }
    std::string AsString() const override { return std::to_string(fData); }

  private:
    int fData;
};

using Counted = G4ReferenceCountedHandle<TesterBase>;
using CountedString = G4ReferenceCountedHandle<TesterString>;
using CountedInt = G4ReferenceCountedHandle<TesterInt>;

void PassByValueCheckConst(Counted c, int expected_count, TesterBase* expected_handle)
{
  EXPECT_EQ(c.Count(), expected_count) << "Handle has wrong count after passing by value";
  EXPECT_EQ(c(), expected_handle) << "Handle does not point to the right object";
}

void PassByValueCheckModify(Counted c, const std::string& modify)
{
  // Modify locally
  c = new TesterString(modify);
  EXPECT_EQ(c.Count(), 1) << "Handle has wrong count modifying";
}

void PassByRefCheckConst(const Counted& c, int expected_count, TesterBase* expected_handle)
{
  EXPECT_EQ(c.Count(), expected_count) << "Handle has wrong count after passing by reference";
  EXPECT_EQ(c(), expected_handle) << "Handle does not point to the right object";
}

void PassByRefCheckModify(Counted& c, const std::string& modify)
{
  // Modify locally
  c = new TesterString(modify);
  EXPECT_EQ(c.Count(), 1) << "Handle has wrong count modifying";
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(G4ReferenceCountedHandle, Construction)
{
  // Default
  Counted t0;
  EXPECT_EQ(t0.Count(), 0) << "Default constructed handle has non-zero count";
  EXPECT_FALSE(t0) << "Default constructed handle is true";
  EXPECT_EQ(t0(), nullptr) << "Default constructed handle does not hold nullptr";

  // Assignment
  Counted t1 = new TesterString("default");
  EXPECT_EQ(t1.Count(), 1) << "Assigned handle has incorrect count";
  EXPECT_TRUE(t1) << "Assigned handle is false";
  EXPECT_NE(t1(), nullptr) << "Assigned handle holds nullptr";

  // Copy assignment
  t0 = t1;
  EXPECT_EQ(t0.Count(), t1.Count()) << "Handles pointing to the same object have different counts";
  EXPECT_EQ(t1.Count(), 2) << "Handle does not have correct reference count";
  EXPECT_EQ(t1(), t0()) << "Handles do not point to same object";

  // Construction/Destruction
  {
    Counted t2;
    EXPECT_EQ(t2.Count(), 0);
    // Also check that no modifaction has happened behind the scenes
    EXPECT_EQ(t1.Count(), 2);

    t2 = t1;
    EXPECT_EQ(t2.Count(), t1.Count())
      << "Handles pointing to the same object have different counts";
    EXPECT_EQ(t2.Count(), 3) << "Handle does not have correct reference count";
    EXPECT_EQ(t2(), t1()) << "Handles do not point to same object";
  }

  // t2 should have destructed so...
  EXPECT_EQ(t0.Count(), t1.Count()) << "Handles pointing to the same object have different counts";
  EXPECT_EQ(t1.Count(), 2) << "Handle does not have correct reference count";
  EXPECT_EQ(t1(), t0()) << "Handles do not point to same object";
}

TEST(G4ReferenceCountedHandle, Copying)
{
  Counted t1 = new TesterString("SomeString");
  Counted t2 = t1;

  // Repoint...
  t1 = new TesterInt(314);
  EXPECT_EQ(t1.Count(), 1) << "Repointed handle should only have one ref";
  EXPECT_EQ(t2.Count(), 1) << "Copied handle should now have only one ref";

  // Access should also be o.k.
  EXPECT_EQ(t2->Value(), 10);
  EXPECT_EQ(t1->Value(), 314);
}

TEST(G4ReferenceCountedHandle, PassByValue)
{
  auto* handle = new TesterString("PassByValue");
  Counted t1 = handle;
  // Passing by value to a subroutine should result in that having a ref count of 2
  PassByValueCheckConst(t1, 2, handle);
  EXPECT_EQ(t1.Count(), 1);

  // Passing by value to something that modifies it not change things back here
  PassByValueCheckModify(t1, "NewValue");
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t1(), handle);
  EXPECT_EQ(t1->AsString(), "PassByValue");
}

TEST(G4ReferenceCountedHandle, NaiveMTPassByValue)
{
#ifndef G4MULTITHREADED
  GTEST_SKIP() << "Test only valid in multi-threaded mode";
#endif

  auto* handle = new TesterString("PassByValue");
  Counted t1 = handle;

  // Passing by value to a thread should behave the same.
  // NB: G4ReferenceCountedHandle is **NOT THREADSAFE IN TOTAL**:
  // 1. The ref count is not mutexed
  // 2. If a handle is moved to a thread, then deletion there may use the wrong TLS allocator
  //
  // We are o.k. here because we immediately join, and we don't move (or at least, the
  // main thread retains at least one copy of the handle.
  // NB: We have three refs in a subthread because of pass-by-value into the thread first
  G4Thread worker1(PassByValueCheckConst, t1, 3, handle);
  worker1.join();
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t1(), handle);

  G4Thread worker2(PassByValueCheckModify, t1, "ThreadNewValue");
  worker2.join();
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t1(), handle);
}

TEST(G4ReferenceCountedHandle, PassByReference)
{
  auto* handle = new TesterString("PassByRef");
  Counted t1 = handle;

  // Passing by ref to a subroutine should result in that having a ref count of 1
  PassByRefCheckConst(t1, 1, handle);
  EXPECT_EQ(t1.Count(), 1);

  // Passing by non-const ref can repoint things
  Counted t2 = t1;
  PassByRefCheckModify(t1, "RefNewValue");
  EXPECT_EQ(t2.Count(), 1);
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t2->AsString(), "PassByRef");
  EXPECT_EQ(t1->AsString(), "RefNewValue");
}

TEST(G4ReferenceCountedHandle, NaiveMTPassByReference)
{
#ifndef G4MULTITHREADED
  GTEST_SKIP() << "Test only valid in multi-threaded mode";
#endif

  auto* handle = new TesterString("PassByRefMT");
  Counted t1 = handle;

  // Passing to a thread by reference is o.k.
  // As above only thread-safe here because we fork-join immediately and retain
  // lifetime of one ref in the parent thread.
  // Only ever have one ref count because of passing ref around.
  G4Thread worker1(PassByRefCheckConst, std::ref(t1), 1, handle);
  worker1.join();
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t1(), handle);

  // Simarly, passing to thread by ref should also behave the same.
  // Passing by non-const ref can repoint things
  Counted t2 = t1;
  G4Thread worker2(PassByRefCheckModify, std::ref(t1), "MTRefNewValue");
  worker2.join();
  EXPECT_EQ(t2.Count(), 1);
  EXPECT_EQ(t1.Count(), 1);
  EXPECT_EQ(t2->AsString(), "PassByRefMT");
  EXPECT_EQ(t1->AsString(), "MTRefNewValue");

  // - TODO: test also for Allocator case (nominally a death/fail test because
  // new/deletion across threads with thread_local allocator should not work)?
  // - TODO: determine is custom allocator is needed here. There is not much
  // being assigned (one int, two pointers). Held thing is new'd elsewhere,
  // so user responsible for checking that isn't moved across threads?
}
