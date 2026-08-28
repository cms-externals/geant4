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
//
// ---------------------------------------------------------------
// Unit test for G4ThreadLocalSingleton
// ---------------------------------------------------------------
#include "G4ThreadLocalSingleton.hh"

#include "G4AutoDelete.hh"
#include "G4Threading.hh"

#include <gtest/gtest.h>

//---------------------------------------------------------------------------//
// "Fixtures" for tests

G4Mutex aMutex = G4MUTEX_INITIALIZER;

struct A
{
    A() = default;
    A(int _a) : a(_a) {};
    ~A() = default;
    int a = -1;
};

class G4SingletonExample
{
    friend class G4ThreadLocalSingleton<G4SingletonExample>;

  private:

    G4SingletonExample() = default;
    static G4Mutex ctrm;
    static int ctr;

  public:

    ~G4SingletonExample() = default;
    static G4SingletonExample* GetInstance()
    {
      static G4ThreadLocalSingleton<G4SingletonExample> inst;
      return inst.Instance();
    }
};

G4ThreadFunReturnType myfunc(G4ThreadFunArgType /*val*/)
{
  G4SingletonExample* inst = G4SingletonExample::GetInstance();
  EXPECT_EQ(inst, G4SingletonExample::GetInstance())
    << "Second call to G4SingletonExample::GetInstance() in thread returns different address";
  // Something here...
  A* a = new A(std::hash<G4Pid_t>{}(G4Threading::G4GetPidId()));
  G4AutoDelete::Register(a);
  return nullptr;
}

void foo(G4SingletonExample* test)
{
  EXPECT_EQ(test, G4SingletonExample::GetInstance())
    << "Call to G4SingletonExample::GetInstance() in function returns different address in "
       "function";
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST(G4ThreadLocalSingleton, Main)
{
  static G4ThreadLocalSingleton<A> aSing;
  A* theS = aSing.Instance();
  EXPECT_EQ(theS, aSing.Instance()) << "Second call to singleton returns different address";
  G4SingletonExample* inst = G4SingletonExample::GetInstance();
  EXPECT_EQ(inst, G4SingletonExample::GetInstance())
    << "Third call to G4SingletonExample::GetInstance() returns different address";
  foo(inst);

  int nthreads = 2;  // num threads
  auto* tid = new G4Thread[nthreads];

  for (int idx = 0; idx < nthreads; ++idx)
  {
    G4Thread* tr = &tid[idx];
    G4THREADCREATE(tr, myfunc, static_cast<void*>(inst));
  }

  for (int idx = 0; idx < nthreads; ++idx)
  {
    G4THREADJOIN(tid[idx]);
  }
}
