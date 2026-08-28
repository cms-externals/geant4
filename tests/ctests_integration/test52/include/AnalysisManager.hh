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
#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "globals.hh"

#include "DataManager.hh"

class AnalysisManager : public DataManager
{
  public:

    static AnalysisManager* Instance(G4double zLow = 0.0, G4double zUp = 10.0);
    static void Destroy();

    void CreateCalorimeter(G4double pos, G4double thickn, G4double rad);
    void ScoreEnergyDeposit(G4double en, G4double x, G4double y, G4double z,
                            const G4String& type = "");
    void ScoreParticleEnergy(G4double en, G4double x, G4double y, G4double z, const G4String& type);
    void ScoreEnteringParticles(G4double en, G4double x, G4double y, G4double z,
                                const G4String& type);
    void ScoreExitingParticles(G4double en, G4double x, G4double y, G4double z,
                               const G4String& type, G4int id);
    void PrintResults();

  protected:

    AnalysisManager(G4double zLow, G4double zUp);
    ~AnalysisManager();

  private:

    static AnalysisManager* instance;

    std::vector<G4double> binBoundaries;

    G4double totEnergyDeposit;
    G4double totEnergyDepositSumSquares;

    G4double elecEnergyEnterTarget;
    G4double photEnergyExitTarget;
    G4double photEnergyExitTargetSumSquares;
    G4double elecEnergyExitTarget;
    G4double elecEnergyExitTargetSumSquares;
    G4double primElecEnergyExitTarget;

    G4int nmbEnergyDeposits;
    G4int nmbElecEnterTarget;
    G4int nmbElecExitTarget;
    G4int nmbPrimElecExitTarget;
    G4int nmbPhotExitTarget;
};

#endif  // ANALYSISMANAGER_HH
