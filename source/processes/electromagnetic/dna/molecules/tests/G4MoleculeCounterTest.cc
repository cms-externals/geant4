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
 * G4MoleculeCounterTest.cc
 *
 *  Created on: 5 mars 2014
 *      Author: kara
 */

#define MOLECULE_COUNTER_TESTING
#include "G4MoleculeCounter.hh"
#include "G4Electron_aq.hh"
#include "G4OH.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Track.hh"

void Test1()
{
  G4MoleculeCounter* counter = G4MoleculeCounter::Instance();
  counter->Use(true);

  G4MolecularConfiguration* e_aq =
      G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4Electron_aq::Definition());
  G4MolecularConfiguration* e_aq2 =
      G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4Electron_aq::Definition());

  assert((e_aq == e_aq2) == true);
  assert((e_aq < e_aq2) == false);
  assert((e_aq2 < e_aq) == false);

  /*
   * Add both e_aq & OH
   */

  const size_t N = 5;
  double time[N] =
  { 1, 5, 10, 12, 20 };

  for(size_t i = 0; i < N; i++)
  {
    counter->AddAMoleculeAtTime(e_aq, time[i]);
  }

  G4MolecularConfiguration* OH =
        G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4OH::Definition());

  const size_t N2 = 5;
  double time2[N2] =
  { 0, 2, 10, 11, 25 };

  for(size_t i = 0; i < N2; i++)
  {
    counter->AddAMoleculeAtTime(OH, time2[i]);
  }

  const NbMoleculeAgainstTime& timeMap = counter->GetNbMoleculeAgainstTime(OH);

  NbMoleculeAgainstTime::const_iterator it = timeMap.begin();
  NbMoleculeAgainstTime::const_iterator endTime_it = timeMap.end();

  {
    int i = 1;

    for(; it != endTime_it; it++)
    {
      G4cout << "Testing : " << it->second << " == " << i << G4endl;
      assert(it->second == i);
      i++;
    }
  }

  /*
   * Test OH
   */

  G4cout << "OH" << G4endl;
//  time2[N2] = {0, 2, 10, 11, 25};
  const size_t Ntests2 = 10;
  double testTimes2[Ntests2] =
  { -1, 0, 1, 4, 9.9, 11, 12, 24, 26, 27 };
  double resultTests2[Ntests2] =
  { 0, 1, 1, 2, 2, 4, 4, 4, 5, 5 };

  for(size_t i = 0; i < Ntests2; i++)
  {
    G4cout << testTimes2[i] << " : "
    << counter->GetNMoleculesAtTime(OH, testTimes2[i]) << G4endl;
//    assert(counter->GetNMoleculesAtTime(OH, testTimes2[i]) == resultTests2[i]);
  }

  assert(counter->GetNMoleculesAtTime(OH, -1) == 0);
  assert(counter->GetNMoleculesAtTime(OH, 0) == 1);
  assert(counter->GetNMoleculesAtTime(OH, 1) == 1);
  assert(counter->GetNMoleculesAtTime(OH, 4) == 2);
  assert(counter->GetNMoleculesAtTime(OH, 9.9) == 2);
  assert(counter->GetNMoleculesAtTime(OH, 11) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 12) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 24) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 26) == 5);
  assert(counter->GetNMoleculesAtTime(OH, 27) == 5);

  /*
   * Test e_aq
   */

  const size_t Ntests = 10;
  double testTimes[Ntests] =
  { 0, 1, 3, 5, 6, 10, 12, 15, 20, 25 };
  double resultTests[Ntests] =
  { 0, 1, 1, 2, 2, 3, 4, 4, 5, 5 };

  for(size_t i = 0; i < Ntests; i++)
  {
    G4cout << testTimes[i] << " : "
           << counter->GetNMoleculesAtTime(e_aq, testTimes[i]) << G4endl;
  }

  assert(counter->GetNMoleculesAtTime(e_aq, 0) == 0);
  assert(counter->GetNMoleculesAtTime(e_aq, 1) == 1);
  assert(counter->GetNMoleculesAtTime(e_aq, 3) == 1);
  assert(counter->GetNMoleculesAtTime(e_aq, 5) == 2);
  assert(counter->GetNMoleculesAtTime(e_aq, 6) == 2);
  assert(counter->GetNMoleculesAtTime(e_aq, 10) == 3);
  assert(counter->GetNMoleculesAtTime(e_aq, 12) == 4);
  assert(counter->GetNMoleculesAtTime(e_aq, 15) == 4);
  assert(counter->GetNMoleculesAtTime(e_aq, 20) == 5);
  assert(counter->GetNMoleculesAtTime(e_aq, 25) == 5);
}

void Test2()
{
  G4MoleculeCounter* counter = G4MoleculeCounter::Instance();
  counter->Use(true);

  G4MolecularConfiguration* e_aq =
      G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4Electron_aq::Definition());
  G4MolecularConfiguration* e_aq2 =
      G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4Electron_aq::Definition());

  assert((e_aq == e_aq2) == true);
  assert((e_aq < e_aq2) == false);
  assert((e_aq2 < e_aq) == false);

  /*
   * Add both e_aq & OH
   */

  const size_t N = 5;
  double time[N] =
  { 1, 5, 10, 12, 20 };

  double remove_times[N] =
  { 24, 25, 27, 27.5, 29};
  /* 4, 3, 2, 1, 0*/

  std::vector<G4Track*> tracks_eq(N*2);

  G4ThreeVector position;

  for(size_t i = 0; i < N; i++)
  {
//    counter->AddAMoleculeAtTime(e_aq, time[i]);
    tracks_eq[i] = (new G4Molecule(e_aq))->BuildTrack(time[i], position);
  }

  for(size_t i = 0; i < N; i++)
  {
    tracks_eq[i]->SetGlobalTime(remove_times[i]);
    delete tracks_eq[i];
  }

  //-------------------------------------------------------------

  G4MolecularConfiguration* OH =
        G4MolecularConfiguration::GetOrCreateMolecularConfiguration(G4OH::Definition());

  const size_t N2 = 5;
  double time2[N2] =
  { 0, 2, 10, 11, 25 };

  std::vector<G4Track*> tracks_OH(N2);

  for(size_t i = 0; i < N2; i++)
  {
//    counter->AddAMoleculeAtTime(OH, time2[i]);
    tracks_OH[i] = (new G4Molecule(OH))->BuildTrack(time2[i], position);
  }

  const NbMoleculeAgainstTime& timeMap = counter->GetNbMoleculeAgainstTime(OH);

  NbMoleculeAgainstTime::const_iterator it = timeMap.begin();
  NbMoleculeAgainstTime::const_iterator endTime_it = timeMap.end();

  {
    int i = 1;

    for(; it != endTime_it; it++)
    {
      G4cout << "Testing : " << it->second << " == " << i << G4endl;
      assert(it->second == i);
      i++;
    }
  }
  //-------------------------------------------------------------

  /*
   * Test OH
   */

  G4cout << "OH" << G4endl;
//  time2[N2] = {0, 2, 10, 11, 25};
  const size_t Ntests2 = 10;
  double testTimes2[Ntests2] =
  { -1, 0, 1, 4, 9.9, 11, 12, 24, 26, 27 };
  double resultTests2[Ntests2] =
  { 0, 1, 1, 2, 2, 4, 4, 4, 5, 5 };

  for(size_t i = 0; i < Ntests2; i++)
  {
    G4cout << testTimes2[i] << " : "
    << counter->GetNMoleculesAtTime(OH, testTimes2[i]) << G4endl;
//    assert(counter->GetNMoleculesAtTime(OH, testTimes2[i]) == resultTests2[i]);
  }

  assert(counter->GetNMoleculesAtTime(OH, -1) == 0);
  assert(counter->GetNMoleculesAtTime(OH, 0) == 1);
  assert(counter->GetNMoleculesAtTime(OH, 1) == 1);
  assert(counter->GetNMoleculesAtTime(OH, 4) == 2);
  assert(counter->GetNMoleculesAtTime(OH, 9.9) == 2);
  assert(counter->GetNMoleculesAtTime(OH, 11) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 12) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 24) == 4);
  assert(counter->GetNMoleculesAtTime(OH, 26) == 5);
  assert(counter->GetNMoleculesAtTime(OH, 27) == 5);

  //-------------------------------------------------------------

  /*
   * Test e_aq
   */

//  const size_t N = 5;
//  double time[N] =
//  { 1, 5, 10, 12, 20 };
//
//  double remove_times[N] =
//  { 24, 25, 27, 27.5, 29};
//  /* 4, 3, 2, 1, 0*/

  const size_t Ntests = 14;
  double testTimes[Ntests] =
  { 0, 1, 3, 5, 6, 10, 12, 15, 20, 24, 26, 27, 27.5, 29 };
  double resultTests[Ntests] =
  { 0, 1, 1, 2, 2,  3,  4,  4,  5,  4,  3,  2,   1,  0 };


  for(size_t i = 0; i < Ntests; i++)
  {
    G4cout << i << " "
           << testTimes[i] << " : "
           << counter->GetNMoleculesAtTime(e_aq, testTimes[i]) << G4endl;

    assert(counter->GetNMoleculesAtTime(e_aq, testTimes[i]) == resultTests[i]);
  }

//  assert(true == false);

  assert(counter->GetNMoleculesAtTime(e_aq, 0) == 0);
  assert(counter->GetNMoleculesAtTime(e_aq, 1) == 1);
  assert(counter->GetNMoleculesAtTime(e_aq, 3) == 1);
  assert(counter->GetNMoleculesAtTime(e_aq, 5) == 2);
  assert(counter->GetNMoleculesAtTime(e_aq, 6) == 2);
  assert(counter->GetNMoleculesAtTime(e_aq, 10) == 3);
  assert(counter->GetNMoleculesAtTime(e_aq, 12) == 4);
  assert(counter->GetNMoleculesAtTime(e_aq, 15) == 4);
  assert(counter->GetNMoleculesAtTime(e_aq, 20) == 5);
//  assert(counter->GetNMoleculesAtTime(e_aq, 25) == 5);
}


int main()
{
  G4MoleculeCounter::Instance()->CheckTimeForConsistency(false);
  Test1();
  G4MoleculeCounter::Instance()->ResetCounter();
  Test1();
  G4MoleculeCounter::Instance()->ResetCounter();
  Test2();
  G4MoleculeCounter::DeleteInstance();

  G4cout << "test OK" << G4endl;
  return 0;
}

