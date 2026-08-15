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
#ifndef DATAMANAGER_HH
#define DATAMANAGER_HH

#include "globals.hh"

#include <fstream>
#include <utility>
#include <vector>

class DataManager;

typedef std::pair<G4double, DataManager*> coll;
typedef std::vector<coll> collector;

class zCompare
{
  public:

    bool operator()(const coll& l, const coll& r) const { return keyLess(l.first, r.first); }

    bool operator()(const coll& l, const coll::first_type& k) const { return keyLess(l.first, k); }

    bool operator()(const coll::first_type& k, const coll& r) const { return keyLess(k, r.first); }

  private:

    bool keyLess(const coll::first_type& k1, const coll::first_type& k2) const { return k1 < k2; }
};

class DataManager
{
  public:

    DataManager(G4double zLow, G4double zUp);
    virtual ~DataManager();

    virtual void ScoreEnergyDeposit(G4double en, G4double x, G4double y, G4double z,
                                    const G4String& type = "all");
    virtual void ScoreParticleEnergy(G4double en, G4double x, G4double y, G4double z,
                                     const G4String& type);
    virtual void PrintResults();

    void AddDataCollector(DataManager* comp);
    bool HasDataCollector(G4double z);
    bool OverlapsWithOtherChild(G4double zLow, G4double zUp);
    bool IsContained(G4double zLow, G4double zUp);
    DataManager* MatchingChildDataCollector(G4double zLow, G4double zUp);

    G4double GetUpperBound() { return zAxisUpperBound; }
    G4double GetLowerBound() { return zAxisLowerBound; }

  private:

    collector dataCollectors;
    collector garbage;

    G4double zAxisLowerBound;
    G4double zAxisUpperBound;
    G4double eps;
};

#endif  // DATAMANAGER_HH
