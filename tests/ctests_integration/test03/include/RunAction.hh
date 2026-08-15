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

#  include "G4AnalysisUtilities.hh"
#  include "G4UImanager.hh"
#  include "G4UserRunAction.hh"
#  include "globals.hh"

#  include <sstream>

class G4Run;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following
/// physics quantities:
/// - Edep in absorber
/// - Track length in absorber
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed
/// dispersion is printed.
///

class EventAction;

class RunAction : public G4UserRunAction
{
  public:

    RunAction(EventAction* eventAction);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

  private:

    template<typename HT>
    void TestGetHt(const G4int id) const;

    void TestWriting() const;
    void TestReading() const;
    void TestGetCommands() const;
    void PrintStatistics() const;

    EventAction* fEventAction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename HT>
void RunAction::TestGetHt(G4int id) const
{
  // apply /analysis/hn|pn/get id
  G4UImanager* uiManager = G4UImanager::GetUIpointer();
  auto hnType = G4Analysis::GetHnType<HT>();
  auto command = "/analysis/" + hnType + "/get";
  uiManager->ApplyCommand(command + " " + std::to_string(id));

  // get command value
  auto htAdress = uiManager->GetCurrentValues(command.c_str());
  if (htAdress.empty())
  {
    G4cerr << "Get " + hnType + " id = " << id << " address failed." << G4endl;
    return;
  }

  // get histogram/profile
  void* htPtr = nullptr;
  std::istringstream is(htAdress);
  is >> htPtr;
  auto ht = static_cast<HT*>(htPtr);

  // print histogram/profile title or error message if failure
  if (ht != nullptr)
  {
    G4cout << "Got " + hnType + " id = " << id << ": " << ht->title() << G4endl;
  }
  else
  {
    G4cerr << "Failed to get " + hnType + " id = " << id << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
