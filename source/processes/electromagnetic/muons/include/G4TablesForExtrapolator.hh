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
//---------------------------------------------------------------------------
//
// ClassName:    G4TablesForExtrapolator
//
// Description:  This class keep dedx, range, inverse range tables
//               for extrapolator
//
// Author:       24.10.14 V.Ivanchenko
//
// Modification:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4TABLESFOREXTRAPOLATOR_HH
#define G4TABLESFOREXTRAPOLATOR_HH

#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"

#include <vector>

class G4ParticleDefinition;
class G4ProductionCuts;
class G4ProductionCutsTable;
class G4MaterialCutsCouple;
class G4Material;
class G4EmDataHandler;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

enum ExtTableType
{
  fDedxElectron = 0,
  fDedxPositron,
  fDedxProton,
  fDedxMuon,
  fRangeElectron,
  fRangePositron,
  fRangeProton,
  fRangeMuon,
  fInvRangeElectron,
  fInvRangePositron,
  fInvRangeProton,
  fInvRangeMuon,
  fMscElectron
};

class G4TablesForExtrapolator
{
  public:

    explicit G4TablesForExtrapolator(G4int verb, G4int bins, G4double e1, G4double e2);

    ~G4TablesForExtrapolator() = default;

    const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const;

    void Initialisation();

    // hide assignment operator
    G4TablesForExtrapolator& operator=(const G4TablesForExtrapolator& right) = delete;
    G4TablesForExtrapolator(const G4TablesForExtrapolator&) = delete;

  private:

    void PrepareTable(std::size_t idx);

    void ComputeElectronDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

    void ComputeMuonDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

    void ComputeProtonDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

    void ComputeTrasportXS(const G4ParticleDefinition* part, G4PhysicsTable* table);

    G4DataVector cuts;

    const G4ParticleDefinition* electron{nullptr};
    const G4ParticleDefinition* positron{nullptr};
    const G4ParticleDefinition* muonPlus{nullptr};
    const G4ParticleDefinition* muonMinus{nullptr};
    const G4ParticleDefinition* proton{nullptr};

    G4EmDataHandler* theData;
    G4ProductionCutsTable* fCoupleTable;

    G4PhysicsTable* dedxElectron{nullptr};
    G4PhysicsTable* dedxPositron{nullptr};
    G4PhysicsTable* dedxMuon{nullptr};
    G4PhysicsTable* dedxProton{nullptr};
    G4PhysicsTable* rangeElectron{nullptr};
    G4PhysicsTable* rangePositron{nullptr};
    G4PhysicsTable* rangeMuon{nullptr};
    G4PhysicsTable* rangeProton{nullptr};
    G4PhysicsTable* invRangeElectron{nullptr};
    G4PhysicsTable* invRangePositron{nullptr};
    G4PhysicsTable* invRangeMuon{nullptr};
    G4PhysicsTable* invRangeProton{nullptr};
    G4PhysicsTable* mscElectron{nullptr};

    G4double emin;
    G4double emax;

    G4int verbose;
    G4int nbins;
    std::size_t nCouples{0};
    G4bool splineFlag{false};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
