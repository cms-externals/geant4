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
// File name:     G4LossTableManager
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications: by V.Ivanchenko
//
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableManager.hh"

#include "G4Electron.hh"
#include "G4ElectronIonPair.hh"
#include "G4EmConfigurator.hh"
#include "G4EmCorrections.hh"
#include "G4EmParameters.hh"
#include "G4EmSaturation.hh"
#include "G4EmTableType.hh"
#include "G4EmTableUtil.hh"
#include "G4EmUtility.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NIELCalculator.hh"
#include "G4Neutron.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VEmProcess.hh"
#include "G4VMultipleScattering.hh"
#include "G4VSubCutProducer.hh"
#include "G4VXRayModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

static std::once_flag applyOnce;
G4ThreadLocal G4LossTableManager* G4LossTableManager::instance = nullptr;

G4LossTableManager* G4LossTableManager::Instance()
{
  if (nullptr == instance)
  {
    static G4ThreadLocalSingleton<G4LossTableManager> inst;
    instance = inst.Instance();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::~G4LossTableManager()
{
  for (auto const& p : loss_vector)
  {
    delete p;
  }
  for (auto const& p : msc_vector)
  {
    delete p;
  }
  for (auto const& p : emp_vector)
  {
    delete p;
  }
  for (auto const& p : p_vector)
  {
    delete p;
  }
  for (auto const& p : xray_vector)
  {
    delete p;
  }

  std::size_t mod = mod_vector.size();
  std::size_t fmod = fmod_vector.size();
  for (std::size_t a = 0; a < mod; ++a)
  {
    if (nullptr != mod_vector[a])
    {
      for (std::size_t b = 0; b < fmod; ++b)
      {
        if ((G4VEmModel*)(fmod_vector[b]) == mod_vector[a])
        {
          fmod_vector[b] = nullptr;
        }
      }
      delete mod_vector[a];
      mod_vector[a] = nullptr;
    }
  }
  for (auto const& p : fmod_vector)
  {
    delete p;
  }

  Clear();
  delete emCorrections;
  delete emConfigurator;
  delete emElectronIonPair;
  delete nielCalculator;
  delete atomDeexcitation;
  delete subcutProducer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::G4LossTableManager()
{
  theParameters = G4EmParameters::Instance();
  theRegistry = G4EmDataRegistry::Instance();
  emCorrections = new G4EmCorrections(verbose);

  // only one thread is the master
  std::call_once(applyOnce, [this]() { isMaster = true; });
  std::size_t n = 70;
  loss_vector.reserve(n);
  part_vector.reserve(n);
  base_part_vector.reserve(n);
  msc_vector.reserve(10);
  emp_vector.reserve(16);
  mod_vector.reserve(150);
  fmod_vector.reserve(60);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Clear()
{
  all_tables_are_built = false;
  currentLoss = nullptr;
  currentParticle = nullptr;
  if (n_loss > 0)
  {
    loss_map.clear();
    loss_vector.clear();
    part_vector.clear();
    base_part_vector.clear();
    n_loss = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEnergyLossProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : loss_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VEnergyLossProcess : " << p->GetProcessName()
           << "  idx= " << n_loss << G4endl;
  }
  ++n_loss;
  loss_vector.push_back(p);
  part_vector.push_back(nullptr);
  base_part_vector.push_back(nullptr);
  if (nullptr == theElectron)
  {
    theElectron = G4Electron::Electron();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ResetParameters()
{
  // initialisation once per run
  if (!resetParam)
  {
    return;
  }
  resetParam = false;
  startInitialisation = true;
  if (0 <= run)
  {
    theParameters->SetIsPrintedFlag(true);
  }
  if (isMaster)
  {
    verbose = theParameters->Verbose();

    // defined base material flag
    theRegistry->CheckBaseMaterials();

    // dump EM parameters only once
    if (verbose > 0 && -1 == run)
    {
      theParameters->Dump();
    }
  }
  else
  {
    verbose = theParameters->WorkerVerbose();
  }
  emCorrections->SetVerbose(verbose);

  if (nullptr != nielCalculator)
  {
    nielCalculator->Initialise();
  }
  if (nullptr != emConfigurator)
  {
    emConfigurator->SetVerbose(verbose);
  }
  if (nullptr != emElectronIonPair)
  {
    emElectronIonPair->SetVerbose(verbose);
  }
  if (nullptr != atomDeexcitation)
  {
    atomDeexcitation->SetVerboseLevel(verbose);
    atomDeexcitation->InitialiseAtomicDeexcitation();
  }
  if (1 < verbose)
  {
    G4cout << "====== G4LossTableManager::ResetParameters "
           << " Nloss=" << loss_vector.size() << " run=" << run << " master=" << isMaster << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEnergyLossProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (G4int i = 0; i < n_loss; ++i)
  {
    if (loss_vector[i] == p)
    {
      loss_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VMultipleScattering* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : msc_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VMultipleScattering : " << p->GetProcessName()
           << "  idx= " << msc_vector.size() << G4endl;
  }
  msc_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VMultipleScattering* p)
{
  if (nullptr == p)
  {
    return;
  }
  std::size_t msc = msc_vector.size();
  for (std::size_t i = 0; i < msc; ++i)
  {
    if (msc_vector[i] == p)
    {
      msc_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : emp_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VEmProcess : " << p->GetProcessName()
           << "  idx= " << emp_vector.size() << G4endl;
  }
  emp_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  std::size_t emp = emp_vector.size();
  for (std::size_t i = 0; i < emp; ++i)
  {
    if (emp_vector[i] == p)
    {
      emp_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : p_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VProcess : " << p->GetProcessName()
           << "  idx= " << p_vector.size() << G4endl;
  }
  p_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VProcess* p)
{
  if (nullptr == p)
  {
    return;
  }
  std::size_t emp = p_vector.size();
  for (std::size_t i = 0; i < emp; ++i)
  {
    if (p_vector[i] == p)
    {
      p_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmModel* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : mod_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  mod_vector.push_back(p);
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VEmModel : " << p->GetName() << "  " << p << "  "
           << mod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmModel* p)
{
  std::size_t n = mod_vector.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    if (mod_vector[i] == p)
    {
      mod_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmFluctuationModel* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : fmod_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  fmod_vector.push_back(p);
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VEmFluctuationModel : " << p->GetName() << "  "
           << fmod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmFluctuationModel* p)
{
  std::size_t n = fmod_vector.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    if (fmod_vector[i] == p)
    {
      fmod_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VXRayModel* p)
{
  if (nullptr == p)
  {
    return;
  }
  for (auto& ptr : xray_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::Register G4VXRayModel : " << p->GetName() << "  "
           << xray_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VXRayModel* p)
{
  std::size_t n = xray_vector.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    if (xray_vector[i] == p)
    {
      xray_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::RegisterExtraParticle(const G4ParticleDefinition* part,
                                               G4VEnergyLossProcess* p)
{
  if (nullptr == p || nullptr == part)
  {
    return;
  }
  for (auto& ptr : loss_vector)
  {
    if (ptr == p)
    {
      return;
    }
  }
  if (verbose > 1)
  {
    G4cout << "G4LossTableManager::RegisterExtraParticle " << part->GetParticleName()
           << "  G4VEnergyLossProcess : " << p->GetProcessName() << "  idx= " << n_loss << G4endl;
  }
  ++n_loss;
  loss_vector.push_back(p);
  part_vector.push_back(part);
  base_part_vector.push_back(p->BaseParticle());
  all_tables_are_built = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEnergyLossProcess*
G4LossTableManager::GetEnergyLossProcess(const G4ParticleDefinition* aParticle)
{
  if (aParticle != currentParticle)
  {
    currentParticle = aParticle;
    if (theElectron == aParticle)
    {
      currentLoss = electronLoss;
    }
    else if (G4EmTableUtil::IsIon(aParticle))
    {
      currentLoss = ionLoss;
    }
    else
    {
      auto const& pos = loss_map.find(aParticle);
      if (pos != loss_map.end())
      {
        currentLoss = (*pos).second;
      }
      else
      {
        currentLoss = nullptr;
      }
    }
  }
  return currentLoss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle,
                                             G4VEnergyLossProcess* p)
{
  if (1 < verbose)
  {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " << particle->GetParticleName()
           << " and " << p->GetProcessName() << " run= " << run << "   loss_vector "
           << loss_vector.size() << " run=" << run << " master=" << isMaster << G4endl;
  }

  // start initialisation for the first run
  if (-1 == run)
  {
    if (nullptr != emConfigurator)
    {
      emConfigurator->PrepareModels(particle, p);
    }
    // define particles for given process only once
    for (G4int j = 0; j < n_loss; ++j)
    {
      if (p == loss_vector[j])
      {
        if (nullptr == part_vector[j])
        {
          part_vector[j] = particle;
          if (particle == theElectron)
          {
            if (p->IsIonisationProcess())
            {
              electronLoss = p;
            }
          }
          else if (particle->GetParticleName() == "GenericIon")
          {
            theGenericIon = particle;
            ionLoss = p;
          }
        }
      }
    }
  }
  ResetParameters();
  if (isMaster)
  {
    for (auto& ptr : loss_vector)
    {
      ptr->SetTablesAreBuilt(false);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle, G4VEmProcess* p)
{
  if (1 < verbose)
  {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " << particle->GetParticleName()
           << " and " << p->GetProcessName() << " run=" << run << " master=" << isMaster << G4endl;
  }

  // start initialisation for the first run
  if (-1 == run && nullptr != emConfigurator)
  {
    emConfigurator->PrepareModels(particle, p);
  }

  ResetParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* part,
                                             G4VMultipleScattering* p)
{
  if (1 < verbose)
  {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " << part->GetParticleName() << " and "
           << p->GetProcessName() << " run=" << run << " master=" << isMaster << G4endl;
  }

  // start initialisation for the first run
  if (-1 == run && nullptr != emConfigurator)
  {
    emConfigurator->PrepareModels(part, p);
  }

  ResetParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::BuildPhysicsTable(const G4ParticleDefinition*)
{
  if (startInitialisation && nullptr != emConfigurator)
  {
    emConfigurator->Clear();
  }
  if (startInitialisation)
  {
    resetParam = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::LocalPhysicsTables(const G4ParticleDefinition* aParticle,
                                            G4VEnergyLossProcess* p)
{
  if (1 < verbose)
  {
    G4cout << "### G4LossTableManager::LocalPhysicsTable() for " << aParticle->GetParticleName()
           << " and process " << p->GetProcessName() << " worker thread" << G4endl;
  }
  currentParticle = nullptr;
  if (startInitialisation)
  {
    if (-1 == run)
    {
      if (nullptr != emConfigurator)
      {
        emConfigurator->Clear();
      }
      firstParticle = aParticle;
    }
    ++run;
    if (1 < verbose)
    {
      G4cout << "===== G4LossTableManager::LocalPhysicsTable() for run " << run
             << " =====" << G4endl;
    }
    startInitialisation = false;
    resetParam = true;
  }

  for (G4int i = 0; i < n_loss; ++i)
  {
    if (p == loss_vector[i])
    {
      p->SetTablesAreBuilt(true);
      if (0 == run)
      {
        part_vector[i] = p->Particle();
        base_part_vector[i] = p->BaseParticle();
        if (p->IsIonisationProcess())
        {
          loss_map[part_vector[i]] = p;
        }
      }
      if (1 < verbose)
      {
        G4cout << i << ".   " << p->GetProcessName() << "  for "
               << part_vector[i]->GetParticleName()
               << "  isIonisation= " << p->IsIonisationProcess() << G4endl;
      }
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::BuildPhysicsTable(const G4ParticleDefinition* aParticle,
                                           G4VEnergyLossProcess* p)
{
  if (1 < verbose)
  {
    G4cout << "### G4LossTableManager::BuildPhysicsTable() for " << aParticle->GetParticleName()
           << " and process " << p->GetProcessName() << G4endl;
  }
  currentParticle = nullptr;

  // initialisation of tables is performed once for all dedx processes
  if (startInitialisation)
  {
    if (-1 == run)
    {
      if (nullptr != emConfigurator)
      {
        emConfigurator->Clear();
      }
      firstParticle = aParticle;
    }
    ++run;
    startInitialisation = false;
    resetParam = true;
    if (1 < verbose)
    {
      G4cout << "    run=" << run << " atomDeexcitation: " << atomDeexcitation << G4endl;
    }

    for (G4int i = 0; i < n_loss; ++i)
    {
      auto el = loss_vector[i];
      el->SetTablesAreBuilt(false);
      if (0 < run)
      {
        continue;
      }
      part_vector[i] = el->Particle();
      base_part_vector[i] = el->BaseParticle();
      if (1 < verbose)
      {
        G4cout << i << ".  " << el->GetProcessName() << " for " << el->Particle()->GetParticleName()
               << "  isIonisation= " << el->IsIonisationProcess();
        if (nullptr != base_part_vector[i])
        {
          G4cout << " baseParticle " << base_part_vector[i]->GetParticleName();
        }
        G4cout << G4endl;
      }
    }

    G4bool buildCSDA = theParameters->BuildCSDARange();

    // parallel initialisation only for the 1st run
    if (0 == run && theParameters->UseParallelInitialisation())
    {
      if (1 < verbose)
      {
        G4cout << "### G4LossTableManager::BuildPhysicsTable() "
               << "parallel initialisation is enabled." << G4endl;
      }
    }

    // build all dedx, range, inverse range tables
    for (G4int i = 0; i < n_loss; ++i)
    {
      if (!loss_vector[i]->TablesAreBuilt() && nullptr == base_part_vector[i])
      {
        G4int idx =
          G4EmUtility::BuildTables(part_vector[i], part_vector, loss_vector, verbose, buildCSDA);
        if (1 < verbose)
        {
          G4cout << "### Build Table for i=" << i << "  idx=" << idx << G4endl;
        }
        if (idx >= 0)
        {
          if (0 == run)
          {
            loss_map[part_vector[idx]] = loss_vector[idx];
          }
          CopyTables(part_vector[idx], loss_vector[idx]);
        }
      }
    }
  }
  if (1 < verbose)
  {
    G4cout << "### G4LossTableManager::BuildPhysicsTable done." << G4endl;
    if (nullptr != subcutProducer)
    {
      G4cout << "     SubCutProducer <" << subcutProducer->GetName() << ">" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::CopyTables(const G4ParticleDefinition* part,
                                    G4VEnergyLossProcess* base_proc)
{
  for (G4int j = 0; j < n_loss; ++j)
  {
    G4VEnergyLossProcess* proc = loss_vector[j];

    // base particle approach
    if (part == base_part_vector[j])
    {
      proc->SetTablesAreBuilt(true);
      proc->SetLambdaTable(base_proc->LambdaTable());
      if (proc->IsIonisationProcess())
      {
        proc->SetDEDXTable(base_proc->IonisationTable(), fRestricted);
        proc->SetRangeTableForLoss(base_proc->RangeTableForLoss());
        proc->SetInverseRangeTable(base_proc->InverseRangeTable());
        if (0 == run)
        {
          loss_map[part_vector[j]] = proc;
        }
        proc->SetDEDXTable(base_proc->DEDXunRestrictedTable(), fTotal);
        proc->SetCSDARangeTable(base_proc->CSDARangeTable());
      }
      if (1 < verbose)
      {
        G4cout << "   CopyTables for " << proc->GetProcessName() << " for "
               << part_vector[j]->GetParticleName() << " base_part= " << part->GetParticleName()
               << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ParticleHaveNoLoss(const G4ParticleDefinition* aParticle)
{
  G4ExceptionDescription ed;
  ed << "Energy loss process not found for " << aParticle->GetParticleName() << " !";
  G4Exception("G4LossTableManager::ParticleHaveNoLoss", "em0001", FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEnergyLossProcess*>& G4LossTableManager::GetEnergyLossProcessVector()
{
  return loss_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEmProcess*>& G4LossTableManager::GetEmProcessVector()
{
  return emp_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VMultipleScattering*>& G4LossTableManager::GetMultipleScatteringVector()
{
  return msc_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmSaturation* G4LossTableManager::EmSaturation() const
{
  return theParameters->GetEmSaturation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmConfigurator* G4LossTableManager::EmConfigurator()
{
  if (nullptr == emConfigurator)
  {
    emConfigurator = new G4EmConfigurator(verbose);
  }
  return emConfigurator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectronIonPair* G4LossTableManager::ElectronIonPair()
{
  if (nullptr == emElectronIonPair)
  {
    emElectronIonPair = new G4ElectronIonPair(verbose);
  }
  return emElectronIonPair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetNIELCalculator(G4NIELCalculator* ptr)
{
  if (nullptr != ptr && ptr != nielCalculator)
  {
    delete nielCalculator;
    nielCalculator = ptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NIELCalculator* G4LossTableManager::NIELCalculator()
{
  if (nullptr == nielCalculator)
  {
    nielCalculator = new G4NIELCalculator(nullptr, verbose);
  }
  return nielCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableManager::SetAtomDeexcitation(G4VAtomDeexcitation* p)
{
  if (atomDeexcitation != p)
  {
    delete atomDeexcitation;
    atomDeexcitation = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetSubCutProducer(G4VSubCutProducer* p)
{
  if (subcutProducer != p)
  {
    delete subcutProducer;
    subcutProducer = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableManager::PrintEWarning(G4String tit, G4double /*val*/)
{
  G4String ss = "G4LossTableManager::" + tit;
  G4ExceptionDescription ed;
  /*
  ed << "Parameter is out of range: " << val
     << " it will have no effect!\n" << " ## "
     << " nbins= " << nbinsLambda
     << " nbinsPerDecade= " << nbinsPerDecade
     << " Emin(keV)= " << minKinEnergy/keV
     << " Emax(GeV)= " << maxKinEnergy/GeV;
  */
  G4Exception(ss, "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableManager::DumpHtml()
{
  // Automatic generation of html documentation page for physics lists
  // List processes and models for the most important
  // particles in descending order of importance
  // NB. for model names with length > 18 characters the .rst file needs
  // to be edited by hand. Or modify G4EmModelManager::DumpModelList

  char* dirName = std::getenv("G4PhysListDocDir");
  char* physList = std::getenv("G4PhysListName");
  if (dirName && physList)
  {
    G4String physListName = G4String(physList);
    G4String pathName = G4String(dirName) + "/" + physListName + ".rst";

    std::ofstream outFile;
    outFile.open(pathName);

    outFile << physListName << G4endl;
    outFile << std::string(physListName.length(), '=') << G4endl;

    std::vector<G4ParticleDefinition*> particles{
      G4Gamma::Gamma(),
      G4Electron::Electron(),
      G4Positron::Positron(),
      G4Proton::ProtonDefinition(),
      G4MuonPlus::MuonPlusDefinition(),
      G4MuonMinus::MuonMinusDefinition(),
    };

    std::vector<G4VEmProcess*> emproc_vector = GetEmProcessVector();
    std::vector<G4VEnergyLossProcess*> enloss_vector = GetEnergyLossProcessVector();
    std::vector<G4VMultipleScattering*> mscat_vector = GetMultipleScatteringVector();

    for (auto theParticle : particles)
    {
      outFile << G4endl << "**" << theParticle->GetParticleName() << "**" << G4endl << G4endl
              << " .. code-block:: none" << G4endl;

      G4ProcessManager* pm = theParticle->GetProcessManager();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int plen = pm->GetProcessListLength();

      for (auto emproc : emproc_vector)
      {
        for (G4int i = 0; i < plen; ++i)
        {
          G4VProcess* proc = (*pv)[i];
          if (proc == emproc)
          {
            outFile << G4endl;
            proc->ProcessDescription(outFile);
            break;
          }
        }
      }

      for (auto mscproc : mscat_vector)
      {
        for (G4int i = 0; i < plen; ++i)
        {
          G4VProcess* proc = (*pv)[i];
          if (proc == mscproc)
          {
            outFile << G4endl;
            proc->ProcessDescription(outFile);
            break;
          }
        }
      }

      for (auto enlossproc : enloss_vector)
      {
        for (G4int i = 0; i < plen; ++i)
        {
          G4VProcess* proc = (*pv)[i];
          if (proc == enlossproc)
          {
            outFile << G4endl;
            proc->ProcessDescription(outFile);
            break;
          }
        }
      }
    }
    outFile.close();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
