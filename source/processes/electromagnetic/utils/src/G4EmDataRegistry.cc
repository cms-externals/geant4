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

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Author:  V.Ivanchenko 18.04.2024
//
//----------------------------------------------------------------------------

#include "G4EmDataRegistry.hh"

#include "G4AutoLock.hh"
#include "G4EmParameters.hh"
#include "G4LossTableBuilder.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProductionCutsTable.hh"

G4EmDataRegistry* G4EmDataRegistry::instance = nullptr;

namespace
{
G4Mutex theEmData = G4MUTEX_INITIALIZER;
const G4int nManagerMax = 8;
const G4int nHandlerMax = 50;
const G4int pdg[8] = {11, -11, -13, 13, 211, -211, 321, -321};
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmDataRegistry* G4EmDataRegistry::Instance()
{
  if (instance == nullptr)
  {
    static G4EmDataRegistry manager;
    instance = &manager;
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataRegistry::G4EmDataRegistry()
{
  fManagers.reserve(nManagerMax);
  fID.reserve(nManagerMax);
  fParticle.reserve(nManagerMax);
  isBuilt.reserve(nManagerMax);
  fDataHandlers.reserve(nHandlerMax);
  theDensityFactor = new std::vector<G4double>;
  theDensityIdx = new std::vector<G4int>;
  theFlag = new std::vector<G4bool>;
  theFluct = new std::vector<G4bool>;
  theParameters = G4EmParameters::Instance();
  theTableBuilder = new G4LossTableBuilder();
  theTableBuilder->SetRegistry(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataRegistry::~G4EmDataRegistry()
{
  for (auto const& p : fDataHandlers)
  {
    delete p;
  }
  fDataHandlers.clear();
  delete theDensityFactor;
  delete theDensityIdx;
  delete theFlag;
  delete theFluct;
  delete theTableBuilder;
  theDensityFactor = nullptr;
  theDensityIdx = nullptr;
  theFlag = nullptr;
  theFluct = nullptr;
  theTableBuilder = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler* G4EmDataRegistry::GetHandlerByName(const G4String& nam, std::size_t n)
{
  // handler already exist
  G4EmDataHandler* ptr = EmDataHandler(nam);
  if (nullptr != ptr)
  {
    return ptr;
  }

  // create a new handler
  G4AutoLock l(&theEmData);
  ptr = EmDataHandler(nam);
  if (nullptr == ptr)
  {
    ptr = new G4EmDataHandler(n, nam);
  }
  l.unlock();

  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::Register(G4EmDataHandler* ptr)
{
  if (nullptr == ptr)
  {
    return;
  }
  for (auto const& p : fDataHandlers)
  {
    if (p == ptr)
    {
      return;
    }
  }
  fDataHandlers.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::DeRegister(G4EmDataHandler* ptr)
{
  if (nullptr == ptr)
  {
    return;
  }
  std::size_t n = fDataHandlers.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    if (fDataHandlers[i] == ptr)
    {
      fDataHandlers[i] = nullptr;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::Register(G4LossTableManager* ptr, G4int id)
{
  if (nullptr == ptr || fNumberOfManagers >= nManagerMax)
  {
    return;
  }
  for (auto const& p : fManagers)
  {
    if (p == ptr)
    {
      return;
    }
  }
  fManagers.push_back(ptr);
  fID.push_back(id);
  ++fNumberOfManagers;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler* G4EmDataRegistry::EmDataHandler(const G4String& nam)
{
  G4EmDataHandler* ptr = nullptr;
  if (fDataHandlers.empty())
  {
    return ptr;
  }
  for (auto const& p : fDataHandlers)
  {
    if (p->GetName() == nam)
    {
      ptr = p;
      break;
    }
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4int>* G4EmDataRegistry::GetCoupleIndexes() const
{
  return theDensityIdx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4double>* G4EmDataRegistry::GetDensityFactors() const
{
  return theDensityFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4bool>* G4EmDataRegistry::GetFluctuationFlags() const
{
  return theFluct;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4bool>* G4EmDataRegistry::GetFlags() const
{
  return theFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataRegistry::GetFlag(std::size_t idx) const
{
  return (idx < theFlag->size()) ? (*theFlag)[idx] : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmDataRegistry::GetBaseMaterialFlag() const
{
  return baseMatFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LossTableBuilder* G4EmDataRegistry::GetLossTableBuilder() const
{
  return theTableBuilder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmDataRegistry::CheckBaseMaterials()
{
  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t nCouples = theCoupleTable->GetTableSize();
  if (0 == nCouples)
  {
    return;
  }

  G4AutoLock l(&theEmData);
  std::size_t nFlags = theFlag->size();

  // extend vectors
  if (nFlags < nCouples)
  {
    theFlag->resize(nCouples, true);
    theFluct->resize(nCouples, theParameters->LossFluctuation());
    theParameters->DefineFluctuationFlags(theFluct);
    theDensityFactor->resize(nCouples, 1.0);
    theDensityIdx->resize(nCouples, 0);
  }

  if (!baseMatFlag)
  {
    if (nFlags < nCouples)
    {
      for (std::size_t i = nFlags; i < nCouples; ++i)
      {
        if (nullptr
            != theCoupleTable->GetMaterialCutsCouple((G4int)i)->GetMaterial()->GetBaseMaterial())
        {
          baseMatFlag = true;
          break;
        }
      }
    }
  }
  isInitialized = false;
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::InitialiseBaseMaterials(const G4PhysicsTable* table)
{
  if (nullptr == table || isInitialized)
  {
    return;
  }

  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t nCouples = theCoupleTable->GetTableSize();
  std::size_t nFlags = theFlag->size();
  if (0 == nCouples)
  {
    return;
  }
  else if (nFlags < nCouples)
  {
    CheckBaseMaterials();
  }

  G4AutoLock l(&theEmData);

  // define default flag and index of used material cut couple
  for (std::size_t i = 0; i < nCouples; ++i)
  {
    (*theFlag)[i] = table->GetFlag(i);
    (*theDensityIdx)[i] = (G4int)i;
  }

  // use base materials
  if (baseMatFlag)
  {
    for (std::size_t i = 0; i < nCouples; ++i)
    {
      // base material is needed only for a couple which is not
      // initialised and for which tables will be computed
      auto couple = theCoupleTable->GetMaterialCutsCouple((G4int)i);
      auto pcuts = couple->GetProductionCuts();
      auto mat = couple->GetMaterial();
      auto bmat = mat->GetBaseMaterial();

      // base material exists - find it and check if it can be reused
      if (nullptr != bmat)
      {
        for (std::size_t j = 0; j < nCouples; ++j)
        {
          if (j == i)
          {
            continue;
          }
          auto bcouple = theCoupleTable->GetMaterialCutsCouple((G4int)j);

          if (bcouple->GetMaterial() == bmat && bcouple->GetProductionCuts() == pcuts)
          {
            // based couple exist in the same region
            (*theDensityFactor)[i] = mat->GetDensity() / bmat->GetDensity();
            (*theDensityIdx)[i] = (G4int)j;
            (*theFlag)[i] = false;

            // ensure that there will no double initialisation
            (*theDensityFactor)[j] = 1.0;
            (*theDensityIdx)[j] = (G4int)j;
            (*theFlag)[j] = true;
            break;
          }
        }
      }
    }
  }
  isInitialized = true;
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::BuildTablesInParallel(const G4int verb)
{
  if (2 >= fNumberOfManagers)
  {
    return;
  }
  auto ptable = G4ParticleTable::GetParticleTable();
  for (G4int i = 0; i < nManagerMax; ++i)
  {
    fParticle.push_back(ptable->FindParticle(pdg[i]));
    isBuilt.push_back(0);
  }
  parallelIni = true;
  verbose = verb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataRegistry::TablesAreBuilt(const G4ParticleDefinition* part, const G4int id)
{
  if (!parallelIni || id >= nManagerMax)
  {
    return true;
  }

  G4AutoLock l(&theEmData);
  G4bool res = true;
  G4int i = 0;
  for (; i < nManagerMax; ++i)
  {
    if (part == fParticle[i] && 0 == isBuilt[i])
    {
      res = false;
      break;
    }
  }
  l.unlock();
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::ResetFlag(const G4ParticleDefinition* part)
{
  G4AutoLock l(&theEmData);
  for (G4int i = 0; i < nManagerMax; ++i)
  {
    if (part == fParticle[i])
    {
      isBuilt[i] = 1;
      for (G4int j = i + 1; j < nManagerMax; ++j)
      {
        if (0 == isBuilt[j])
        {
          return;
        }
      }
      break;
    }
  }
  parallelIni = false;

  if (1 < verbose)
  {
    G4cout << "### G4EmDataRegistry::BuildTablesInParallel done" << G4endl;
  }
}

const std::vector<const G4ParticleDefinition*>& G4EmDataRegistry::ListForParallelBuild() const
{
  return fParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
