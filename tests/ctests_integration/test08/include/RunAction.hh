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
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#  define RunAction_h 1

#  include "G4AccArray.hh"
#  include "G4AccMap.hh"
#  include "G4AccUnorderedMap.hh"
#  include "G4AccValue.hh"
#  include "G4AccVector.hh"
#  include "G4UserRunAction.hh"
#  include "globals.hh"

// user defined accumulable
#  include "ProcCounterAccumulable.hh"

#  include <limits>
#  include <vector>

class G4Run;
class G4LogicalVolume;

/// @brief Define std::hash template specialization for G4String
/// so that it can be used  as a key in unordered associative containers
namespace std
{
template<>
struct hash<G4String>
{
    std::size_t operator()(const G4String& str) const { return std::hash<std::string>()(str); }
};
}  // namespace std

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:

    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void AddNsec(G4int nsec);
    void AddEdep(G4double edep);
    void SetPassed(G4bool passed);
    void CountEvent();
    void CountProcess(G4String procName);

  private:

    // methods
    void DefineUnits();
    G4bool Check(G4VAccumulable* accumulable, const G4String& message) const;
    G4int GetId(G4VAccumulable* accumulable, G4String& message) const;
    void TestGetById();
    G4String GetName(G4VAccumulable* accumulable, G4String& message) const;
    void TestGetByName();
    void PrintAccRange(G4int firstId, G4int count, const G4String& message) const;
    void PrintAllAccumulables() const;

    // Accummulables of double Value type
    // B1 example use case
    G4AccValue<G4double> fEdep{0.};
    // B1 example use case
    G4AccValue<G4double> fEdep2{"Edep2", 0.};
    // Merge mode: kMaximum
    G4AccValue<G4double> fMaxEdep{0., G4MergeMode::kMaximum};
    // Merge mode: kMinimum
    G4AccValue<G4double> fMinEdep{std::numeric_limits<double>::max(), G4MergeMode::kMinimum};
    // Test prefix increment operator
    G4AccValue<G4double> fEventCounter{0., G4MergeMode::kAddition};
    // Test postfix increment operator
    G4AccValue<G4double> fEventCounter2{0., G4MergeMode::kAddition};
    // Save of all double value accumulables to simplify testing
    std::vector<G4AccValue<G4double>*> fDValues;

    // Accummulables of other than double Value type
    // G4int, default merge mode: kAddition
    G4AccValue<G4int> fNsec{0};
    // G4bool, default merge mode: kAddition
    G4AccValue<G4bool> fPassed{false};

    // User defined accumulable
    ProcCounterAccumulable* fProcCounter{nullptr};
    ProcCounterAccumulable* fProcCounterTest1{nullptr};

    // Vectors
    // Vector accumulable ctor 1
    G4AccVector<G4double> fEdepVectorCtor1;
    // Vector accumulable ctor 3 - initialized with count, value
    G4AccVector<G4double> fEdepVectorCtor3;
    // Vector accumulable ctor 3n - initialized with count, value, name
    G4AccVector<G4double> fEdepVectorCtor3n;
    // Vector accumulable ctor 4 - initialized with count
    G4AccVector<G4double> fEdepVectorCtor4;
    // Vector accumulable ctor 4n - initialized with count, value, name
    G4AccVector<G4double> fEdepVectorCtor4n;
    // Vector accumulable ctor 10 - constructed with std::initializer_list
    G4AccVector<G4double> fEdepVectorCtor10{{0., 0.}};
    // Vector accumulable ctor 10n - constructed with std::initializer_list, name
    G4AccVector<G4double> fEdepVectorCtor10n{"EdepVector10n", {0., 0.}};
    // Save of all vectors to simplify testing
    std::vector<G4AccVector<G4double>*> fVectors;

    // Arrays
    // Array intialized with default ctor
    G4AccArray<G4double, 2> fEdepArrayCtor1;
    // Array intialized with default ctor with name
    G4AccArray<G4double, 2> fEdepArrayCtor1n;
    // Array intialized with ctor 2 with the initializer list
    G4AccArray<G4double, 2> fEdepArrayCtor2{0., 0.};
    // Array intialized with ctor 2 with the initializer list and name
    G4AccArray<G4double, 2> fEdepArrayCtor3{"EdepArray3", 0., 0.};
    // Save of all arrays to simplify testing
    std::vector<G4AccArray<G4double, 2>*> fArrays;

    // Maps
    // Map intialized with default ctor
    G4AccMap<G4String, G4int> fProcCounterMapCtor1;
    // Map intialized with default ctor with name
    G4AccMap<G4String, G4int> fProcCounterMapCtor1n;
    // Map intialized with ctor 10
    G4AccMap<G4String, G4int> fProcCounterMapCtor10{
      {G4String("CoulombScat"), 0}, {G4String("Rayl"), 0}, {G4String("FictiveProc"), 500}};
    // Map intialized with  ctor 10 with name
    G4AccMap<G4String, G4int> fProcCounterMapCtor10n{
      "ProcCounterMap10n",
      {{G4String("CoulombScat"), 0}, {G4String("Rayl"), 0}, {G4String("FictiveProc"), 500}}};
    // Save of all maps to simplify testing
    std::vector<G4AccMap<G4String, G4int>*> fMaps;

    // Unordered maps
    // std::string as key
    // Unordered map intialized with default ctor
    G4AccUnorderedMap<G4String, G4int> fProcCounterUMapCtor1;
    // Map intialized with default ctor with name
    G4AccUnorderedMap<G4String, G4int> fProcCounterUMapCtor1n;
    // Map intialized with ctor 10
    G4AccUnorderedMap<G4String, G4int> fProcCounterUMapCtor10{
      {G4String("CoulombScat"), 0}, {G4String("Rayl"), 0}, {G4String("FictiveProc"), 500}};
    // Map intialized with  ctor 10 with name
    G4AccUnorderedMap<G4String, G4int> fProcCounterUMapCtor10n{
      "ProcCounterUMap10n",
      {{G4String("CoulombScat"), 0}, {G4String("Rayl"), 0}, {G4String("FictiveProc"), 500}}};
    // Save of all maps to simplify testing
    std::vector<G4AccUnorderedMap<G4String, G4int>*> fUnorderedMaps;
};

#  include "RunAction.icc"

#endif
