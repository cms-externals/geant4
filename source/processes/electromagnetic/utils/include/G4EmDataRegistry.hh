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
// 18 April 2024   V.Ivanchenko
//
// This is a singleton class to store shared EM physics data.
// The access to data is possible only by data name.
// Once created data pointers does not changed or deleted until the
// end of the job, this class is responsible for destruction.
//

#ifndef G4EMDATAREGISTRY_HH
#define G4EMDATAREGISTRY_HH

#include "G4EmDataHandler.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"

#include <vector>

class G4EmParameters;
class G4LossTableBuilder;
class G4LossTableManager;
class G4ParticleDefinition;

class G4EmDataRegistry
{
  public:

    static G4EmDataRegistry* Instance();

    ~G4EmDataRegistry();

    // create or access data handler
    // nTable - number of tables in a new handler
    G4EmDataHandler* GetHandlerByName(const G4String&, std::size_t nTable);

    // register a new handler
    void Register(G4LossTableManager*, G4int threadID);

    // register a new handler
    void Register(G4EmDataHandler*);

    // register a new handler
    void DeRegister(G4EmDataHandler*);

    // initialise base materials
    void InitialiseBaseMaterials(const G4PhysicsTable* table = nullptr);

    // access methods
    const std::vector<G4int>* GetCoupleIndexes() const;

    const std::vector<G4double>* GetDensityFactors() const;

    const std::vector<G4bool>* GetFluctuationFlags() const;

    const std::vector<G4bool>* GetFlags() const;

    G4bool GetFlag(std::size_t idx) const;

    G4bool GetBaseMaterialFlag() const;

    G4LossTableBuilder* GetLossTableBuilder() const;

    // method called once per run before initialisation
    void CheckBaseMaterials();

    // build tables in parallel
    void BuildTablesInParallel(const G4int verbose);

    const std::vector<const G4ParticleDefinition*>& ListForParallelBuild() const;

    G4bool TablesAreBuilt(const G4ParticleDefinition*, const G4int id);

    void ResetFlag(const G4ParticleDefinition*);

    // hide assignment operator
    G4EmDataRegistry& operator=(const G4EmDataRegistry&) = delete;
    G4EmDataRegistry(const G4EmDataRegistry&) = delete;

  private:

    G4EmDataRegistry();

    G4EmDataHandler* EmDataHandler(const G4String&);

    static G4EmDataRegistry* instance;
    G4EmParameters* theParameters;
    G4LossTableBuilder* theTableBuilder;

    G4bool isInitialized{false};
    G4bool baseMatFlag{false};
    G4bool parallelIni{false};
    G4int fNumberOfManagers{0};
    G4int verbose{0};

    std::vector<G4EmDataHandler*> fDataHandlers;
    std::vector<G4LossTableManager*> fManagers;
    std::vector<const G4ParticleDefinition*> fParticle;
    std::vector<G4int> isBuilt;
    std::vector<G4int> fID;
    std::vector<G4double>* theDensityFactor;
    std::vector<G4int>* theDensityIdx;
    std::vector<G4bool>* theFlag;
    std::vector<G4bool>* theFluct;
};

#endif
