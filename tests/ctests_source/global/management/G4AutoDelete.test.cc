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
// ---------------------------------------------------------------
// Unit test for G4AutoDelete
// ---------------------------------------------------------------
#include "G4AutoDelete.hh"

#include <gtest/gtest.h>

#include <cstdlib>
#include <functional>
#include <map>

//---------------------------------------------------------------------------//
// "Fixtures" for tests
G4Mutex aMutex = G4MUTEX_INITIALIZER;

struct A
{
    A() = default;
    A(size_t _a) : a(_a) {}
    ~A() = default;
    size_t a{std::numeric_limits<size_t>::max()};
};

struct B
{
    B() = default;
    B(size_t _a) : a(_a) {}
    ~B() = default;
    size_t a{std::numeric_limits<size_t>::max()};
};

struct C
{
    C() = default;
    C(int _a) : a(_a) {}
    ~C() = default;
    size_t a{std::numeric_limits<size_t>::max()};
};

G4Mutex mapMutex = G4MUTEX_INITIALIZER;
std::map<std::size_t, A*> amap;

G4ThreadFunReturnType myfunc(G4ThreadFunArgType /*val*/)
{
  // Something here...
  A* a = new A(std::hash<G4Pid_t>{}(G4Threading::G4GetPidId()));
  G4AutoDelete::Register(a);
  G4AutoDelete::Register(new B(a->a));
  G4AutoDelete::Register(new C(a->a));
  G4AutoLock l(&mapMutex);
  amap[std::hash<G4Pid_t>{}(G4Threading::G4GetPidId())] = a;
  return nullptr;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(G4AutoDelete, Registration)
{
  int nthreads = 2;  // num threads
  auto* tid = new G4Thread[nthreads];

  for (int idx = 0; idx < nthreads; ++idx) {
    G4Thread* tr = &tid[idx];
    G4THREADCREATE(tr, myfunc, static_cast<void*>(nullptr));
  }

  for (int idx = 0; idx < nthreads; ++idx) {
    G4THREADJOIN((tid[idx]));
  }

  auto it = amap.begin();
  for (; it != amap.end(); ++it) {
    EXPECT_EQ(it->first, it->second->a) << "Wrong content of AutoDelete object";
  }

  // NOTE: We don't yet have a way to test static destruction behaves correctly...
}
