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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#  define RunAction_h 1

#  include "G4Run.hh"
#  include "G4UserRunAction.hh"
#  include "globals.hh"

class Run;
class RunActionMessenger;
class DetectorConstruction;
class HistoManager;
class TestSeries;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:

    explicit RunAction(DetectorConstruction*);
    ~RunAction() override;

    G4Run* GenerateRun() override;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void CreateRangeTest(G4double refRange, G4double relError);
    void CreateEnergyLossTest(G4double refLoss, G4double relError);

    void SetDNATest(G4bool val) { testDNA = val; }

  private:

    void TestDNAStopping();

    DetectorConstruction* fDetector;
    RunActionMessenger* fMessenger;
    HistoManager* fHistoManager;
    Run* fRun{nullptr};

    TestSeries* testRange;
    TestSeries* testEnergyLoss;
    G4bool testDone{false};
    G4bool testDNA{false};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
