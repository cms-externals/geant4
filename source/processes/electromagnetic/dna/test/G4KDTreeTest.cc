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
 * G4KDTreeTest.cc
 *
 *  Created on: 9 avr. 2015
 *      Author: matkara
 */

#include "G4KDTree.hh"
#include <cassert>
#include <vector>
#include "Randomize.hh"
#include "G4RandomDirection.hh"

using namespace std;

// global variables
G4KDTree gTree;
std::vector<G4ThreeVector*> generatedPositions;
const int nPositions = 100000;
const double boxSize = 100000;
const int nTrials = 1000;

void GenerateRandomPositions()
{
  //generatedPositions.reserve(nPositions);
  for(int i = 0 ; i < nPositions ; ++i)
  {
    G4ThreeVector* pos = new G4ThreeVector(G4RandomDirection()*G4UniformRand()*boxSize);
   // generatedPositions[i] = pos;
    generatedPositions.push_back(pos);
    gTree.Insert(pos);
  }
}

G4ThreeVector* FindExpectedAnser(G4ThreeVector& trialPoint)
{
//    size_t _out_i = 0;
    double dis = 3e8;
    G4ThreeVector* answer(0);

    for(size_t i = 0 ; i != generatedPositions.size() ; i++)
    {
        double tmp= (*(generatedPositions[i])-trialPoint).mag();
        if(tmp < dis)
        {
//          _out_i = i;
          dis = tmp;
          answer = generatedPositions[i];
        }
    }

    return answer;
}

G4ThreeVector* FindExpectedAnser(size_t index)
{
//    size_t _out_i = 0;
    double dis = 3e8;
    G4ThreeVector* answer(0);

    G4ThreeVector& trialPoint = *(generatedPositions[index]);

    for(size_t i = 0 ; i != generatedPositions.size() ; i++)
    {
        if(i == index) continue;
        double tmp= (*(generatedPositions[i])-trialPoint).mag();
        if(tmp < dis)
        {
//          _out_i = i;
          dis = tmp;
          answer = generatedPositions[i];
        }
    }

    return answer;
}


void TestReply(G4ThreeVector& testPoint,
               G4ThreeVector* expectedAnswer)
{
  G4KDTreeResultHandle result = gTree.Nearest(testPoint);
  G4ThreeVector* reply = result->GetItem<G4ThreeVector>();
  assert(*reply == *expectedAnswer);

  if(*reply != *expectedAnswer)
  {
    G4Exception("TestReply", "NOT_MATCH_1", FatalException, "Not match");
  }
}

void TestReply(G4ThreeVector* testPoint,
               G4ThreeVector* expectedAnswer)
{
  G4KDTreeResultHandle result = gTree.Nearest(*testPoint);
  //G4ThreeVector* reply = result->GetItem<G4ThreeVector>();
  G4KDNode_Base* node = result->GetNode();
  assert(node != 0);
  if(node == 0)
  {
    G4Exception("TestReply", "NODE_NULL", FatalException, "Node ptr null");
  }

  result = gTree.Nearest(node);
  G4ThreeVector* reply = result->GetItem<G4ThreeVector>();
  assert(*reply == *expectedAnswer);
  if(*reply != *expectedAnswer)
  {
    G4Exception("TestReply", "NOT_MATCH_2", FatalException, "Not match");
  }
}

void RandomTrials()
{
  G4cout << "Generate " << nTrials << " random trials" << G4endl;
  for(int i = 0 ; i < nTrials ; ++i)
  {
    G4cout << "trial: " << i << "\r";
    G4ThreeVector trialPoint (G4RandomDirection()*G4UniformRand()*boxSize);
    G4ThreeVector* expectedAnswer = FindExpectedAnser(trialPoint);
    assert(expectedAnswer != 0);

    if(expectedAnswer == 0)
    {
      G4Exception("RandomTrials", "EXPECTED_ANSW", FatalException, "no expected answer!");
    }

    TestReply(trialPoint, expectedAnswer);
  }
}

void TrialsFromGeneratedPoints()
{
  srand (time(NULL));
  size_t index;
  G4cout << "Generate " << nTrials << " trials from already inserted points" << G4endl;
  for(int i = 0 ; i < nTrials ; ++i)
  {
    G4cout << "trial: " << i << "\r";
    index = rand() % nPositions;
    G4ThreeVector* expectedAnswer = FindExpectedAnser(index);
    assert(expectedAnswer != 0);
    if(expectedAnswer == 0)
    {
      G4Exception("TrialsFromGeneratedPoints", "EXPECTED_ANSW", FatalException, "no expected answer!");
    }
    TestReply(generatedPositions[index], expectedAnswer);
  }
}


int main()
{

/////////////////////////////////////
// QUICK TESTS
//  G4ThreeVector center(0, 0, 0);
//
//  G4ThreeVector* pos1 = new G4ThreeVector(-2, -2, 0);
//  tree.Insert(pos1);
//  TestReply(center, pos1);
//
//  G4ThreeVector* pos2 = new G4ThreeVector(-1, -1, 0);
//  tree.Insert(pos2);
//  TestReply(center, pos2);
//
//  G4ThreeVector* pos3 = new G4ThreeVector(-1, 0.5, 0);
//  tree.Insert(pos3);
//  TestReply(center, pos3);
/////////////////////////////////////

  GenerateRandomPositions();
  RandomTrials();
  TrialsFromGeneratedPoints();

}
