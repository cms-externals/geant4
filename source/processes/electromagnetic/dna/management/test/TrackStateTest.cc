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
/*
 * TrackStateTest.cc
 *
 *  Created on: 7 avr. 2014
 *      Author: kara
 */

#include "G4TrackState.hh"
#include "globals.hh"
#include <cassert>

////////////////////////:
// Define class and State
class Class : public G4TrackStateDependent<Class>
{
public:
  Class();
  void NewTrackState();
  class State
  {
  public:
    State()
    {
      test = 0;
    }
    State(int i);
    int test;
  };
  static int ID;
  int fInstanceID;
};

int Class::ID = 0;

Class::Class()
{
  fInstanceID = Class::ID;
  Class::ID++;
}

Class::State::State(int i)
{
  test = i;
}

//LINKSTATE(Class, Class::State)

template<>
class G4TrackState<Class> : public G4TrackStateBase<Class>,
                            public Class::State
{
  typedef Class::State State;
  friend class G4TrackStateDependent<Class> ;
public:
  G4TrackState() :
      G4TrackStateBase<Class>(), Class::State()
  {
  }
  G4TrackState(int i) :
      G4TrackStateBase<Class>(), Class::State(i)
  {
  }
  virtual ~G4TrackState()
  {
  }
  virtual int GetID()
  {
    return G4TrackStateID<Class>::GetID();
  }
  static int ID()
  {
    return G4TrackStateID<Class>::GetID();
  }
protected:
};

// NewTrackState, LoadTrackState, etc... can be redefined if needed
void Class::NewTrackState()
{
  fpTrackState = StateTypeHandle(new StateType(fInstanceID));
}

// End of class and state definition
//////

void CreateState(Class& c, G4TrackStateManager& man)
{
  c.NewTrackState();
  c.SaveTrackState(man);
  c.PopTrackState();
}

//----------

struct Class2{
  class State
  {
  public:
    State()
    {
      test = 0;
    }
    State(int i);
    int test;
  };
};

template<>
class G4TrackState<Class2> : public G4TrackStateBase<Class2>,
public Class2::State
{
  typedef Class2::State State;
  friend class G4TrackStateDependent<Class2> ;
public:
  G4TrackState() :
  G4TrackStateBase<Class2>(), Class2::State()
  {
  }
  G4TrackState(int i) :
  G4TrackStateBase<Class2>(), Class2::State(i)
  {
  }
  virtual ~G4TrackState()
  {
  }
  virtual int GetID()
  {
    return G4TrackStateID<Class2>::GetID();
  }
  static int ID()
  {
    return G4TrackStateID<Class2>::GetID();
  }
protected:
};

//#include "G4ITNavigator.hh"

int main()
{
  G4TrackStateManager man;
  Class c[10];

  for (int i = 0; i < 10; i++)
  {
    CreateState(c[i], man);
    assert(c[i].GetConcreteTrackState().get() == 0);
  }

  for (int i = 0; i < 10; i++)
  {
    G4VTrackStateHandle state1 = man.GetTrackState(&c[i]);
    G4shared_ptr<G4TrackState<Class> > state_converted = ConvertToConcreteTrackState<Class>(state1);
    assert(state_converted->test == c[i].fInstanceID);
    G4cout << state_converted->test  << G4endl;
  }
  for (int i = 0; i < 10; i++)
  {
    c[i].LoadTrackState(man);
    Class::StateTypeHandle state = c[i].GetConcreteTrackState();
    assert(state->test == c[i].fInstanceID);
  }

//  constexpr int testConstExpr1 = G4TrackState<Class>::ID();
//  constexpr int testConstExpr2 = G4TrackState<Class2>::ID();

//  G4cout << testConstExpr1 << G4endl;
//  G4cout << testConstExpr2 << G4endl;
  
  G4cout << G4TrackState<Class>::ID() << G4endl;
  G4cout << G4TrackState<Class2>::ID() << G4endl;

  //assert(G4TrackState<Class>::ID() == 0);
  //assert(G4TrackState<Class2>::ID() == 1);
  
  assert(G4TrackState<Class2>::ID() == G4TrackState<Class>::ID()+1);
  
//  G4cout << G4TrackState<G4ITNavigator>::ID() << G4endl;
  
  G4cout << "test OK" << G4endl;

  return 0;
}
