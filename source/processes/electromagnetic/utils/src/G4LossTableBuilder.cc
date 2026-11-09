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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4LossTableBuilder
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 23-01-03 V.Ivanchenko Cut per region
// 21-07-04 V.Ivanchenko Fix problem of range for dedx=0
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 07-12-04 Fix of BuildDEDX table (V.Ivanchenko)
// 27-03-06 Add bool options isIonisation (V.Ivanchenko)
// 16-01-07 Fill new (not old) DEDX table (V.Ivanchenko)
// 12-02-07 Use G4LPhysicsFreeVector for the inverse range table (V.Ivanchenko)
// 24-06-09 Removed hidden bin in G4PhysicsVector (V.Ivanchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableBuilder.hh"

#include "G4EmDataRegistry.hh"
#include "G4EmParameters.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::G4LossTableBuilder(G4bool)
{
  theParameters = G4EmParameters::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::SetRegistry(G4EmDataRegistry* ptr)
{
  theRegistry = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4int>* G4LossTableBuilder::GetCoupleIndexes() const
{
  return theRegistry->GetCoupleIndexes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4double>* G4LossTableBuilder::GetDensityFactors() const
{
  return theRegistry->GetDensityFactors();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4bool>* G4LossTableBuilder::GetFluctuationFlags() const
{
  return theRegistry->GetFluctuationFlags();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LossTableBuilder::GetFlag(std::size_t idx) const
{
  return theRegistry->GetFlag(idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LossTableBuilder::GetBaseMaterialFlag() const
{
  return theRegistry->GetBaseMaterialFlag();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildDEDXTable(G4PhysicsTable* dedxTable,
                                        const std::vector<G4PhysicsTable*>& list) const
{
  std::size_t n_processes = list.size();
  if (1 >= n_processes)
  {
    return;
  }

  std::size_t nCouples = dedxTable->size();
  G4int verb = theParameters->Verbose();
  if (1 < verb)
  {
    G4cout << "G4LossTableBuilder::BuildDEDXTable: Nproc= " << n_processes
           << " nCouples=" << nCouples << " Nv= " << dedxTable->size() << G4endl;
  }
  if (0 >= nCouples)
  {
    return;
  }

  // flags are not checked for the sum of dedx tables
  for (std::size_t i = 0; i < nCouples; ++i)
  {
    auto pv0 = static_cast<G4PhysicsLogVector*>((*(list[0]))[i]);
    if (pv0 == nullptr)
    {
      continue;
    }
    std::size_t npoints = pv0->GetVectorLength();
    auto pv = new G4PhysicsLogVector(*pv0);
    for (std::size_t j = 0; j < npoints; ++j)
    {
      G4double dedx = 0.0;
      for (std::size_t k = 0; k < n_processes; ++k)
      {
        const G4PhysicsVector* pv1 = (*(list[k]))[i];
        dedx += (*pv1)[j];
      }
      pv->PutValue(j, dedx);
    }
    if (splineFlag)
    {
      pv->FillSecondDerivatives();
    }
    G4PhysicsTableHelper::SetPhysicsVector(dedxTable, i, pv);
  }
  if (2 < verb)
  {
    G4cout << "### G4LossTableBuilder::BuildDEDXTable " << G4endl;
    G4cout << *dedxTable << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildRangeTable(const G4PhysicsTable* dedxTable,
                                         G4PhysicsTable* rangeTable, G4bool useBaseMaterial) const
// Build range table from the energy loss table
{
  const std::size_t nCouples = dedxTable->size();
  G4int verb = theParameters->Verbose();
  if (1 < verb)
  {
    G4cout << "### G4LossTableBuilder::BuildRangeTable: nCouples=" << nCouples
           << "  ThreadID=" << G4Threading::G4GetThreadId() << G4endl;
  }
  if (0 >= nCouples)
  {
    return;
  }

  const std::size_t n = 100;
  const G4double del = 1.0 / (G4double)n;

  auto const& theFlag = theRegistry->GetFlags();

  for (std::size_t i = 0; i < nCouples; ++i)
  {
    auto pv = static_cast<G4PhysicsLogVector*>((*dedxTable)[i]);
    if (nullptr == pv || (useBaseMaterial && !(*theFlag)[i]))
    {
      continue;
    }
    std::size_t npoints = pv->GetVectorLength();
    std::size_t bin0 = 0;
    G4double elow = pv->Energy(0);
    G4double ehigh = pv->Energy(npoints - 1);
    G4double dedx1 = (*pv)[0];

    // protection for specific cases dedx=0
    if (dedx1 <= 0.0)
    {
      for (std::size_t k = 1; k < npoints; ++k)
      {
        ++bin0;
        elow = pv->Energy(k);
        dedx1 = (*pv)[k];
        if (dedx1 > 0.0)
        {
          break;
        }
      }
      npoints -= bin0;
    }

    // initialisation of a new vector
    if (npoints < 3)
    {
      npoints = 3;
      bin0 = 1;
    }

    delete (*rangeTable)[i];
    G4PhysicsLogVector* v = (0 == bin0)
                              ? new G4PhysicsLogVector(*pv)
                              : new G4PhysicsLogVector(elow, ehigh, npoints - 1, splineFlag);

    // assumed dedx proportional to beta for the low edge range
    G4double energy1 = v->Energy(0);
    G4double range = 2. * energy1 / dedx1;
    v->PutValue(0, range);

    if (2 < verb)
    {
      G4cout << "New Range vector Npoints=" << v->GetVectorLength() << " coupleIdx=" << i
             << " spline=" << v->GetSpline() << " Elow=" << v->GetMinEnergy()
             << " Ehigh=" << v->GetMaxEnergy() << " DEDX(Elow)=" << dedx1 << " R(Elow)=" << range
             << G4endl;
    }

    for (std::size_t j = 1; j < npoints; ++j)
    {
      G4double energy2 = v->Energy(j);
      G4double de = (energy2 - energy1) * del;
      G4double energy = energy2 + de * 0.5;
      G4double sum = 0.0;
      std::size_t idx = j - 1;
      for (std::size_t k = 0; k < n; ++k)
      {
        energy -= de;
        dedx1 = pv->Value(energy, idx);
        if (dedx1 > 0.0)
        {
          sum += de / dedx1;
        }
      }
      range += sum;
      v->PutValue(j, range);
      energy1 = energy2;
    }
    if (splineFlag)
    {
      v->FillSecondDerivatives();
    }
    G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
  }
  if (2 < verb)
  {
    G4cout << "### Range table" << G4endl;
    G4cout << *rangeTable << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
                                                G4PhysicsTable* invRangeTable,
                                                G4bool useBaseMaterial) const
// Build inverse range table from the energy loss table
{
  std::size_t nCouples = rangeTable->size();
  if (0 >= nCouples)
  {
    return;
  }

  auto const& theFlag = theRegistry->GetFlags();

  for (std::size_t i = 0; i < nCouples; ++i)
  {
    G4PhysicsVector* pv = (*rangeTable)[i];
    if (nullptr == pv || (useBaseMaterial && !(*theFlag)[i]))
    {
      continue;
    }
    std::size_t npoints = pv->GetVectorLength();

    delete (*invRangeTable)[i];
    auto v = new G4PhysicsFreeVector(npoints, splineFlag);

    for (std::size_t j = 0; j < npoints; ++j)
    {
      G4double e = pv->Energy(j);
      G4double r = (*pv)[j];
      v->PutValues(j, r, e);
    }
    if (splineFlag)
    {
      v->FillSecondDerivatives();
    }
    v->EnableLogBinSearch(theParameters->NumberForFreeVector());

    G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
  }
  if (2 < theParameters->Verbose())
  {
    G4cout << "### Inverse range table" << G4endl;
    G4cout << *invRangeTable << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::InitialiseBaseMaterials(const G4PhysicsTable* table) const
{
  theRegistry->InitialiseBaseMaterials(table);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4LossTableBuilder::BuildTableForModel(G4PhysicsTable* aTable, G4VEmModel* model,
                                                       const G4ParticleDefinition* part,
                                                       G4double emin, G4double emax,
                                                       G4bool spline) const
{
  G4PhysicsTable* table = G4PhysicsTableHelper::PreparePhysicsTable(aTable);
  if (aTable != nullptr && aTable != table)
  {
    aTable->clearAndDestroy();
    delete aTable;
  }
  if (nullptr == table)
  {
    return table;
  }

  //  auto const & theFlag = theRegistry->GetFlags();
  G4int nbins = theParameters->NumberOfBinsPerDecade();

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();
  G4int verb = theParameters->Verbose();
  if (2 < verb)
  {
    G4cout << "G4LossTableBuilder::BuildTableForModel Ncouple=" << numOfCouples
           << " model: " << model->GetName() << "  " << part->GetParticleName()
           << " size=" << table->size() << G4endl;
  }
  G4PhysicsLogVector* aVector = nullptr;

  for (G4int i = 0; i < (G4int)numOfCouples; ++i)
  {
    if (table->GetFlag(i))
    {
      // create physics vector and fill it
      auto couple = theCoupleTable->GetMaterialCutsCouple(i);
      delete (*table)[i];

      // if start from zero then change the scale
      const G4Material* mat = couple->GetMaterial();

      G4double tmin = std::max(std::max(emin, model->MinPrimaryEnergy(mat, part)), CLHEP::eV);
      G4int n = nbins;

      if (tmin >= emax)
      {
        aVector = nullptr;
      }
      else
      {
        n *= G4lrint(std::log10(emax / tmin));
        n = std::max(n, 3);
        aVector = new G4PhysicsLogVector(tmin, emax, n, spline);
      }

      if (nullptr != aVector)
      {
        for (G4int j = 0; j <= n; ++j)
        {
          G4double e = aVector->Energy(j);
          G4double y = model->Value(couple, part, e);
          aVector->PutValue(j, y);
        }
        if (spline)
        {
          aVector->FillSecondDerivatives();
        }
      }
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }
  if (2 < verb)
  {
    G4cout << "G4LossTableBuilder::BuildTableForModel done for " << part->GetParticleName()
           << " and " << model->GetName() << "  " << table << G4endl;
    G4cout << *table << G4endl;
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
