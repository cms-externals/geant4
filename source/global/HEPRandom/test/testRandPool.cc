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
// Test for RNG algorithms for Geant4 multi-threaded
//
// Author: A. Dotti
//
//
// History: 
//
// 5 March 2015 - First implementation
// ----------------------------------------------------------------------

#include <iostream>
#include <map>
#include <vector>

#include "Randomize.hh"
#include "G4Threading.hh"
#include "G4Timer.hh"
#include "G4AutoLock.hh"

#define SUCCESS 1
#define FAIL 0

// Just to give some time, on my MacBook Pro Intel i7 2.6GHz:
// 2 min for 3M TRIALS, set a reasonable number. Scales linearly

long TRIALS=100000000;
int GROUPSIZE=10;
int NTHREADS=1;
#define STATFILE "/tmp/status.rand"

#include "G4UniformRandPool.hh"

G4Mutex mutex = G4MUTEX_INITIALIZER;

//G4UniformRandPool::PoolSize_t poolsize = G4UniformRandPool::DEFAULT_POOL_SIZE;
int poolsize = G4UNIFORMRANDPOOL_DEFAULT_POOLSIZE;

int engineType = 0;

// per-thread engine
static G4ThreadLocal CLHEP::HepRandomEngine *pdefaultEngine = 0 ;

void InitMe( int engineSwitch = 0 )
{
  if ( ! pdefaultEngine )
  {
      switch (engineSwitch) {
          case 0:
              pdefaultEngine = new CLHEP::Ranlux64Engine();
              break;
          case 1:
              pdefaultEngine = new CLHEP::RanecuEngine();
              break;
          case 2:
              pdefaultEngine = new CLHEP::RanluxEngine();
              break;
          case 3:
              pdefaultEngine = new CLHEP::MTwistEngine();
              break;
          case 4:
              pdefaultEngine = new CLHEP::RanshiEngine();
              break;
          case 5:
              pdefaultEngine = new CLHEP::HepJamesRandom();
              break;
          case 6:
              pdefaultEngine = new CLHEP::DualRand();
              break;
          case 7:
              pdefaultEngine = new CLHEP::MixMaxRng();
              break;
          //case 8:
          //     pdefaultEngine = new G4MKL_engine();
          //    break;
          default:
              std::cerr<<"Non valid selection for RandomEngine, valid values from 0 (default) to 4"<<std::endl;
              abort();
              break;
      }
  } 
  G4Random::setTheEngine( pdefaultEngine );
}

std::map<G4Pid_t,double> finalRandomValue;

//======================================
// Functions executed in threads
//======================================

// Argument to thread function
struct arg_t
{
    long seed;
    int size;
    long histories;
    long group;
};


typedef int Ret_t;
typedef Ret_t (*Test_t)(void);

G4ThreadFunReturnType mythreadfuncSTD(G4ThreadFunArgType arg)
{
    InitMe( engineType );
    arg_t* conf = (arg_t*)arg;
    long seed = conf->seed;
    G4Random::setTheSeed( seed );
    double val=0;
    for ( int i = 0 ; i < conf->histories ; ++i )
    {
        val=G4UniformRand();
    }
    G4AutoLock l(&mutex);
    finalRandomValue[ G4Threading::G4GetPidId() ] = val;
    return (G4ThreadFunReturnType)NULL;
}

G4ThreadFunReturnType mythreadfuncSTD2(G4ThreadFunArgType arg)
{
    InitMe( engineType );
    arg_t* conf = (arg_t*)arg;
    long seed = conf->seed;
    G4Random::setTheSeed( seed );
    double* val = new double[conf->group];
    for ( int i = 0 ; i < conf->histories/conf->group ; ++i )
    {
        G4Random::getTheEngine()->flatArray(conf->group,val);
    }
    G4AutoLock l(&mutex);
    finalRandomValue[ G4Threading::G4GetPidId() ] = val[conf->group-1];
    delete[] val;
    return (G4ThreadFunReturnType)NULL;
}


G4ThreadFunReturnType mythreadfunc(G4ThreadFunArgType arg)
{
    InitMe( engineType );
    arg_t* conf = (arg_t*)arg;
    long seed = conf->seed;
    G4Random::setTheSeed( seed );
    G4UniformRandPool pool( conf->size );
    double val=0;
    for ( int i = 0 ; i < conf->histories ; ++i )
    {
        val=pool.GetOne();
    }
    G4AutoLock l(&mutex);
    finalRandomValue[ G4Threading::G4GetPidId() ] = val;
    return (G4ThreadFunReturnType)NULL;
}

// Simply run and store final value
G4ThreadFunReturnType mythreadfunc2(G4ThreadFunArgType arg)
{
    InitMe( engineType );
    arg_t* conf = (arg_t*)arg;
    long seed = conf->seed;
    G4Random::setTheSeed( seed );
    G4UniformRandPool pool( conf->size );
    double * val = new G4double[conf->group];
    for ( int i = 0 ; i < conf->histories/conf->group ; ++i )
    {
        pool.GetMany(val,conf->group);
    }
    G4AutoLock l(&mutex);
    finalRandomValue[ G4Threading::G4GetPidId() ] = val[conf->group-1];
    delete[] val;
    return (G4ThreadFunReturnType)NULL;
}

Ret_t threadsSTD0()
{
    finalRandomValue.clear();
    G4Thread * threads = new G4Thread[NTHREADS];
    long seed = 123456789;
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        arg_t myarg = { seed , poolsize, TRIALS , GROUPSIZE };
#if defined(G4MULTITHREADED)
        G4Thread* tr = &threads[i];
        G4THREADCREATE( tr ,  mythreadfuncSTD , &myarg );
#else
	mythreadfuncSTD((G4ThreadFunArgType)&myarg);
#endif
    }
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        G4THREADJOIN( threads[i] );
    }
    std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
    double val = (*it).second;
    ++it;
    for ( ; it != finalRandomValue.end() ; ++it )
    {
        if ( (*it).second != val ) return FAIL;
    }
    delete[] threads;
    return SUCCESS;
}

Ret_t threadsSTD1()
{
    finalRandomValue.clear();
    G4Thread * threads = new G4Thread[NTHREADS];
    long seed = 123456789;
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        arg_t myarg = { seed , poolsize, TRIALS , GROUPSIZE };
#if defined(G4MULTITHREADED)
        G4Thread* tr = &threads[i];
        G4THREADCREATE( tr ,  mythreadfuncSTD2 , &myarg );
#else
	mythreadfuncSTD2((G4ThreadFunArgType)&myarg);
#endif
    }
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        G4THREADJOIN( threads[i] );
    }
    std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
    double val = (*it).second;
    ++it;
    for ( ; it != finalRandomValue.end() ; ++it )
    {
        if ( (*it).second != val ) return FAIL;
    }
    delete[] threads;
    return SUCCESS;
}


Ret_t threads1()
{
  finalRandomValue.clear();
  G4Thread * threads = new G4Thread[NTHREADS];
  long seed = 123456789;
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { seed , poolsize, TRIALS , GROUPSIZE };
#if defined(G4MULTITHREADED)
      G4Thread* tr = &threads[i];
      G4THREADCREATE( tr ,  mythreadfunc , &myarg );
#else
      mythreadfunc((G4ThreadFunArgType)&myarg);
#endif
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  double val = (*it).second;
  ++it;
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( (*it).second != val ) return FAIL;
    }
  delete[] threads;
  return SUCCESS;
}

Ret_t threads2()
{
    finalRandomValue.clear();
    G4Thread* threads = new G4Thread[NTHREADS];
    long seed = 123456789;
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        arg_t myarg = { seed , poolsize, TRIALS , GROUPSIZE };
#if defined(G4MULTITHREADED)
        G4Thread* tr = &threads[i];
        G4THREADCREATE( tr ,  mythreadfunc2 , &myarg );
#else
	mythreadfunc2((G4ThreadFunArgType)&myarg);
#endif
    }
    for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        G4THREADJOIN( threads[i] );
    }
    std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
    double val = (*it).second;
    ++it;
    for ( ; it != finalRandomValue.end() ; ++it )
    {
        if ( (*it).second != val ) return FAIL;
    }
    delete[] threads;
    return SUCCESS;
}
//===================================
// End Tests
//===================================

struct tt {
  Test_t aTest;
  const char* name;
  const char* mess;
};

//The list of tests to be performed
tt tests[] = {
    {&threadsSTD0,"Default","Retrieve single random number from G4Random."},
    {&threadsSTD1,"DefaultMany","Retrieve set random number from G4Random via array interface."},
    {&threads1 ,"GetOne","Retrieve single random number from pool."} ,
    {&threads2,"GetMany","Retrieve set of random numbers from pool."},
  {NULL,NULL,NULL}
 };

Ret_t testme( const tt& t) 
{
  G4Timer timer;
  timer.Start();
    std::cout<<std::left<<"\t  Sub-Test: ";
    std::cout<<std::left<<std::setfill(' ')<<std::setw(20)<<t.name;
  Ret_t result = t.aTest();
  if ( result == FAIL ) 
    {
      std::cout<<std::right<<" FAILED"<<std::flush;
      std::cerr<<"Sub-Test: "<<t.name<<":"<<t.mess<<" FAILED"<<std::endl;
    }
  else 
    {
        std::cout<<std::right<<" SUCCESS";
    }
  timer.Stop();
  std::cout<<"  (t="<<timer<<")"<<std::endl;
  return result;
}


int main(int argc,char** argv)
{
    if ( argc > 1 ) TRIALS = atol(argv[1]);
    if ( argc > 2 ) GROUPSIZE = atoi(argv[2]);
    if ( argc > 3 ) NTHREADS = atoi(argv[3]);
    bool result = SUCCESS;
#if not defined(G4MULTITHREADED)
    if ( NTHREADS > 1) {
        std::cout<<"Testing with sequential build, using only 1 thread"<<std::endl;
        NTHREADS = 1;
    }
#endif
    std::cout<<"Testing Geant4 RNG"<<std::endl;
    std::cout<<"Running with "<<NTHREADS<<" threads, number histories: "<<TRIALS<<std::endl;
    std::cout<<"For array interface, request: "<<GROUPSIZE<<" numbers at once"<<std::endl;
    std::cout<<"Tests: "<<std::endl;
    int i = 0;
    while (true) {
        if ( tests[i].aTest == NULL ) break;
        std::cout<<std::setw(17)<<std::setfill(' ')<<std::left<<tests[i].name;
        std::cout<<std::setw(3)<<" : ";
        std::cout<<std::right<<tests[i].mess<<std::endl;
        ++i;
    }
//    G4UniformRandPool::PoolSize_t ps[] = {
    int ps[] = {
        G4UNIFORMRANDPOOL_DEFAULT_POOLSIZE , G4UNIFORMRANDPOOL_TINY_POOLSIZE,
        G4UNIFORMRANDPOOL_SMALL_POOLSIZE , G4UNIFORMRANDPOOL_MEDIUM_POOLSIZE ,
        G4UNIFORMRANDPOOL_LARGE_POOLSIZE , G4UNIFORMRANDPOOL_HUGE_POOLSIZE };
    for ( int pc = 0 ; pc < 1 ; ++pc ) {
  //for ( int pc = 0 ; pc < 6 ; ++pc ) {
        poolsize = ps[pc];
        std::cout<<"Testing with pool size: "<<poolsize<<std::endl;
   	const int max = 8;//9
        for ( int t = 0 ; t < max ; ++t ) //Loop on engine types
        {
           //t = 7;
            if ( pdefaultEngine )
            {
                delete pdefaultEngine;
                pdefaultEngine = NULL;
            }
            engineType = t;
            InitMe( t );
            std::cout<<"\t================================================="<<std::endl;
            std::cout<<"\tTesting Random Engine: "<<pdefaultEngine->name()<<std::endl;
            std::cout<<"\t================================================="<<std::endl;
            G4Timer someT;
            someT.Start();
            i = 0;
            while (true) {
                if ( tests[i].aTest == NULL ) break;
                result &= testme(tests[i]);
                ++i;
            }
            someT.Stop();
            std::cout<<"\tDone, user elapsed time: "<<someT.GetUserElapsed()<<" s."<<std::endl;
        }
    }
    //std::cout<<"To Conclude some timing on G4MKL_Engine::flat() vs GetOne():"<<std::endl;
    //G4MKL_engine* en = new G4MKL_engine();
    //G4Random::setTheEngine(en);
    //G4UniformRandPool rp;
    //G4Timer t;
    //double v=0;
    //t.Start();
    //for ( int i = 0 ; i < TRIALS ; ++i )
    //  {
    //    v+=en->novirt();
    //  }
    //t.Stop();
    //std::cout<<"Via direct non-virtual interface: "<<t<<std::endl;

    //t.Start();
    //for ( int i = 0 ; i < TRIALS ; ++i )
    //  {
    //    v+=en->flat();
    //  }
    //t.Stop();
    //std::cout<<"Via virtual flat() interface: "<<t<<std::endl;
    //t.Start();
    //for ( int i = 0 ; i < TRIALS ; ++i )
    //  {
    //    v+=rp.GetOne();
    //  }
    //t.Stop();
    //std::cout<<"Via GetOne interface: "<<t<<std::endl;
    //delete en;

    //std::cout<<G4UniformRandPool::flat()<<std::endl;
    //G4double a[10];
    //G4UniformRandPool::flatArray(10,a);
    if ( result == SUCCESS ) return 0;
    else return 1;
}
