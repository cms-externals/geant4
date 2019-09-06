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
/// \file RngUnitTest.hh
/// \brief Templates for testing the different random number generators
//
//

#ifndef RngUnitTest_hh_
#define RngUnitTest_hh_

#include <future>
#include <atomic>
#include <vector>
#include <iostream>
#include <cstdint>

#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/DualRand.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RanshiEngine.h"

//----------------------------------------------------------------------------//

namespace G4Test
{

//----------------------------------------------------------------------------//

typedef std::atomic<uintmax_t> atomic_size_t;

//----------------------------------------------------------------------------//

template <typename _RNG_t>
void RngTest(G4SharedFuture<void> fut, uintmax_t n, atomic_size_t* gfail)
{
    // defer the lock
    G4AutoLock l(G4TypeMutex<decltype(std::cout)>(), std::defer_lock);

    // thread id
    static atomic_size_t tid(0);
    uintmax_t _tid = ++tid;
    // inform
    l.lock(); std::cout << "Thread #" << _tid << " is initialized...\n"; l.unlock();

    // barrier until signaled
    fut.wait();
    l.lock(); std::cout << "Thread #" << _tid << " starting...\n"; l.unlock();

    // set the thread-local RNG
    G4Random::setTheEngine(new _RNG_t());

    // count the failures
    atomic_size_t lfail(0);

    for(uintmax_t i = 0; i < n; ++i)
    {
        try {
            G4double rand = G4UniformRand();
            if(!std::isfinite(rand) || rand < 0.0 || !(rand < 1.0))
                ++lfail;
        } catch(...) { ++lfail; }
    }

    l.lock(); std::cout << "--> Total invalid random numbers: "
                        << lfail << " out of " << n << std::endl;

    (*gfail) += lfail.load();
}

//----------------------------------------------------------------------------//

template <typename _RNG_t>
uintmax_t
RngRun(uintmax_t nthread = G4Threading::G4GetNumberOfCores(),
       uintmax_t nitr = 50000)
{
#if defined(G4MULTITHREADED)
    // 16 threads is plenty
    uintmax_t max_threads = 16;
    nthread = (nthread > max_threads) ? max_threads : nthread;
#else
    nthread = 1;
#endif

    // all non-finite RNG failures
    atomic_size_t gfail(0);
    // list of threads
    std::vector<G4Thread*> threads(nthread, nullptr);
    // synchronize the threads with a promise
    G4Promise<void> prom;
    // the threads will wait on this
    G4SharedFuture<void> fut = prom.get_future().share();
    // start threads
    for(auto& itr : threads)
        itr = new G4Thread(RngTest<_RNG_t>, fut, nitr, &gfail);

    // wait a bit until all threads have initialized
    G4ThisThread::sleep_for(std::chrono::milliseconds(500));
    // indicate the threads to start
    prom.set_value();

    for(auto& itr : threads)
    {
        itr->join();
        delete itr;
    }

    // don't return larger than 255
    const uint_least8_t max_ret =
            std::numeric_limits<uint_least8_t>::max();
    return (gfail.load() > max_ret) ? max_ret : gfail.load();
}

//----------------------------------------------------------------------------//

} // namespace G4Test

//----------------------------------------------------------------------------//

#endif // RngUnitTest_hh_

