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
// Unit test for Threading
// ---------------------------------------------------------------
#include <functional>
#include <iostream>

#include "G4Threading.hh"

G4Mutex g4autolockdebug = G4MUTEX_INITIALIZER;

#define TID std::hash<std::thread::id>{}(G4Threading::G4GetPidId())

#define MESSAGE( msg ) {		  \
        G4MUTEXLOCK(&g4autolockdebug); \
        std::cout << "ThreadID: " << TID << " " << msg<<std::endl; \
        G4MUTEXUNLOCK(&g4autolockdebug); }

// Example of usage of threading classes of G4

// Create a global mutex
G4Mutex mutex = G4MUTEX_INITIALIZER;

// Define a thread-function using G4 types
G4ThreadFunReturnType myfunc( G4ThreadFunArgType val )
{
  double value = *(double*)val;
  MESSAGE( "value is:"<<value );
  return /*(G4ThreadFunReturnType)*/NULL;
}

// Example: spawn 10 threads that execute myfunc
//
int main(int,char**)
{
  MESSAGE( "Starting program ");
  int nthreads = 10;
  G4Thread* tid = new G4Thread[nthreads];
  double *valss = new double[nthreads];
  for ( int idx = 0 ; idx < nthreads ; ++idx )
  {
    valss[idx] = (double)idx;
    G4Thread* tr = &tid[idx];
    G4THREADCREATE( tr , myfunc, &(valss[idx]) );
  }
  for ( int idx = 0 ; idx < nthreads ; ++idx )
  {
    G4THREADJOIN( tid[idx] );
  } 
  MESSAGE( "Program ended ");
  return 0;
}
