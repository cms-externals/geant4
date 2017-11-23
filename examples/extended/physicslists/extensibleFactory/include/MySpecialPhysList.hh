//
// $Id: MySpecialPhysList.hh 66241 2017-11-02 18:34:42Z rhatcher $
//
//---------------------------------------------------------------------------
//
// ClassName:  MySpecialPhysList
//
// Author: 2017-11-02 R. Hatcher
//   Example "alternative" physics list by typedef'ing QBBC
//
//----------------------------------------------------------------------------
//
#ifndef TMySpecialPhysList_h
#define TMySpecialPhysList_h 1

#include "MySpecialPhysList.icc"

// Users would define their own physics list in
//   MySpecialPhysList.icc and MySpecialPhysList.hh
//
// The only requirement for registering with the extensible factory
// is that the physics list constructor must accept a single G4int argument
// which is the verbosity, i.e. :-
//
//  class MySpecialPhysicsList : public G4ModularPhysList
//  {
//     public:
//       MySpecialPhysicsList(G4int ver = 1 [, any defaulted args ] );
//       virtual ~MySpecialPhysicsList();
//       virtual void SetCuts();
//     ....

#include "QBBC.hh"
typedef   QBBC MySpecialPhysList;

#endif

