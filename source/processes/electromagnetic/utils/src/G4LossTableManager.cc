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
#include "G4SystemOfUnits.hh"

#include "G4VMultipleScattering.hh"
#include "G4VEmProcess.hh"
#include "G4EmUtility.hh"

#include "G4EmParameters.hh"
#include "G4EmSaturation.hh"
#include "G4EmConfigurator.hh"
#include "G4ElectronIonPair.hh"
#include "G4NIELCalculator.hh"
#include "G4EmCorrections.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VSubCutProducer.hh"
#include "G4VXRayModel.hh"

#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4EmTableUtil.hh"
#include "G4EmTableType.hh"
#include "G4Region.hh"
#include "G4PhysicalConstants.hh"

#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4OpticalPhoton.hh"
#include "G4Neutron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Threading.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

static std::once_flag applyOnce;
std::vector<const G4ParticleDefinition*> G4LossTableManager::part_vector;
std::vector<const G4ParticleDefinition*> G4LossTableManager::base_part_vector;

G4ThreadLocal G4LossTableManager* G4LossTableManager::instance = nullptr;

G4LossTableManager* G4LossTableManager::Instance()
{
  if(nullptr == instance) {
    static G4ThreadLocalSingleton<G4LossTableManager> inst;
    instance = inst.Instance();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::~G4LossTableManager()
{
  for (auto const & p : loss_vector) { delete p; }
  for (auto const & p : msc_vector) { delete p; }
  for (auto const & p : emp_vector) { delete p; }
  for (auto const & p : p_vector) { delete p; }
  for (auto const & p : xray_vector) { delete p; }

  std::size_t mod = mod_vector.size();
  std::size_t fmod = fmod_vector.size();
  for (std::size_t a=0; a<mod; ++a) {
    if( nullptr != mod_vector[a] ) { 
      for (std::size_t b=0; b<fmod; ++b) {
        if((G4VEmModel*)(fmod_vector[b]) == mod_vector[a]) {
          fmod_vector[b] = nullptr;
        }
      }
      delete mod_vector[a]; 
      mod_vector[a] = nullptr;
    }
  }
  for (auto const & p : fmod_vector) { delete p; }

  delete emCorrections;
  delete emConfigurator;
  delete emElectronIonPair;
  delete nielCalculator;
  delete atomDeexcitation;
  delete subcutProducer;

  all_tables_are_built = false;
  currentLoss = nullptr;
  currentParticle = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::G4LossTableManager()
{
  // only one thread is the master 
  std::call_once(applyOnce, [this]() { isMaster = true; });
  std::size_t n = 70;
  loss_vector.reserve(n);
  if (isMaster) {
    part_vector.reserve(n);
    base_part_vector.reserve(n);
  }
  msc_vector.reserve(10);
  emp_vector.reserve(16);
  mod_vector.reserve(150);
  fmod_vector.reserve(60);

  theParameters = G4EmParameters::Instance();
  theRegistry = G4EmDataRegistry::Instance();
  threadID = G4Threading::G4GetThreadId();
  theRegistry->Register(this, threadID);
  emCorrections = new G4EmCorrections(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEnergyLossProcess* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : loss_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEnergyLossProcess : " 
           << p->GetProcessName() << "  idx= " << n_loss << G4endl;
  }
  ++n_loss;
  loss_vector.push_back(p);
  if (isMaster) {
    part_vector.push_back(nullptr);
    base_part_vector.push_back(nullptr);
  }
  if (nullptr == theElectron) {
    theElectron = G4Electron::Electron();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ResetParameters()
{
  // initialisation once per run
  if (!resetParam) { return; }
  resetParam = false;
  startInitialisation = true;
  if (0 <= run) {
    theParameters->SetIsPrintedFlag(true);
  }
  if (isMaster) {
    verbose = theParameters->Verbose();

    // defined base material flag
    theRegistry->CheckBaseMaterials();

    // dump EM parameters only once
    if (verbose > 0 && -1 == run) {
      theParameters->Dump();
    }
  }
  else {
    verbose = theParameters->WorkerVerbose();
  }
  emCorrections->SetVerbose(verbose);

  for (G4int i=0; i<n_loss; ++i) {
    loss_vector[i]->ResetFlagPrepared();
  }

  if (nullptr != nielCalculator) { nielCalculator->Initialise(); }
  if (nullptr != emConfigurator) { emConfigurator->SetVerbose(verbose); }
  if (nullptr != emElectronIonPair) { emElectronIonPair->SetVerbose(verbose); }
  if (nullptr != atomDeexcitation) {
    atomDeexcitation->SetVerboseLevel(verbose);
    atomDeexcitation->InitialiseAtomicDeexcitation();
  }
  if (1 < verbose) {
    G4cout << "====== G4LossTableManager::ResetParameters " 
           << " Nloss=" << loss_vector.size()
           << " run=" << run << " master=" << isMaster
	   << " id=" << G4Threading::G4GetThreadId()      
           << G4endl;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEnergyLossProcess* p)
{
  if (nullptr == p) { return; }
  for (G4int i=0; i<n_loss; ++i) {
    if (loss_vector[i] == p) { 
      loss_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VMultipleScattering* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : msc_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VMultipleScattering : " 
           << p->GetProcessName() << "  idx= " << msc_vector.size() << G4endl;
  }
  msc_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VMultipleScattering* p)
{
  if (nullptr == p) { return; }
  std::size_t msc = msc_vector.size();
  for (std::size_t i=0; i<msc; ++i) {
    if(msc_vector[i] == p) { 
      msc_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmProcess* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : emp_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmProcess : " 
           << p->GetProcessName() << "  idx= " << emp_vector.size() << G4endl;
  }
  emp_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmProcess* p)
{
  if (nullptr == p) { return; }
  std::size_t emp = emp_vector.size();
  for (std::size_t i=0; i<emp; ++i) {
    if(emp_vector[i] == p) { 
      emp_vector[i] = nullptr; 
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VProcess* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : p_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VProcess : " 
           << p->GetProcessName() << "  idx= " << p_vector.size() << G4endl;
  }
  p_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VProcess* p)
{
  if (nullptr == p) { return; }
  std::size_t emp = p_vector.size();
  for (std::size_t i=0; i<emp; ++i) {
    if(p_vector[i] == p) { 
      p_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmModel* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : mod_vector) { if (ptr == p) { return; } }
  mod_vector.push_back(p);
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmModel : " 
           << p->GetName() << "  " << p << "  " << mod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmModel* p)
{
  std::size_t n = mod_vector.size();
  for (std::size_t i=0; i<n; ++i) {
    if(mod_vector[i] == p) { 
      mod_vector[i] = nullptr; 
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmFluctuationModel* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : fmod_vector) { if (ptr == p) { return; } }
  fmod_vector.push_back(p);
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmFluctuationModel : " 
           << p->GetName() << "  " << fmod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmFluctuationModel* p)
{
  std::size_t n = fmod_vector.size();
  for (std::size_t i=0; i<n; ++i) {
    if(fmod_vector[i] == p) {
      fmod_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VXRayModel* p)
{
  if (nullptr == p) { return; }
  for (auto & ptr : xray_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VXRayModel : " 
           << p->GetName() << "  " << xray_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VXRayModel* p)
{
  std::size_t n = xray_vector.size();
  for (std::size_t i=0; i<n; ++i) {
    if (xray_vector[i] == p) {
      xray_vector[i] = nullptr;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::RegisterExtraParticle(
     const G4ParticleDefinition* part,
     G4VEnergyLossProcess* p)
{ 
  if (nullptr == p || nullptr == part) { return; }
  for (auto & ptr : loss_vector) { if (ptr == p) { return; } }
  if (verbose > 1) {
    G4cout << "G4LossTableManager::RegisterExtraParticle "
           << part->GetParticleName() << "  G4VEnergyLossProcess : " 
           << p->GetProcessName() << "  idx= " << n_loss << G4endl;
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
  if (aParticle != currentParticle) {
    currentParticle = aParticle;
    if (theElectron == aParticle) {
      currentLoss = electronLoss;
    }
    else if (G4EmTableUtil::IsIon(aParticle)) {
      currentLoss = ionLoss;
    }
    else {
      auto const & pos = loss_map.find(aParticle);
      if (pos != loss_map.end()) {
        currentLoss = (*pos).second;
      }
      else {
        currentLoss = nullptr;
      }
    }
  }
  return currentLoss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* part,
                                        G4VEnergyLossProcess* p)
{
  if (1 < verbose) {
    G4cout << "### G4LossTableManager::PreparePhysicsTable for " 
           << part->GetParticleName() 
           << " and " << p->GetProcessName() << " run=" << run 
           << " loss_vector=" << n_loss
           << " run=" << run << " master=" << isMaster
	   << " id=" << G4Threading::G4GetThreadId()
           << G4endl;
  }
  
  // Start initialisation for the first run
  // In the following runs the list of particles, processes, and models
  // are unchanged. Materials and cuts may be changed.
  
  if ( -1 == run ) {
    if (nullptr != emConfigurator) {
      emConfigurator->PrepareModels(part, p);
    }

    // define particles for given process only once
    for (G4int i=0; i<n_loss; ++i) {
      if (p != loss_vector[i]) { continue; }
      if (isMaster) {
        if (nullptr == part_vector[i]) { 
          part_vector[i] = part;
        }
        base_part_vector[i] = p->BaseParticle();
      }
      if (p->IsIonisationProcess()) {
        loss_map[part] = p;
        if (part == theElectron) {
          electronLoss = p;
        }
        else if (part->GetParticleName() == "GenericIon") {
          theGenericIon = part;
          ionLoss = p;
        }
      }
    }
  }
  ResetParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* part)
{
  for (G4int i=0; i<n_loss; ++i) {
    if (part != part_vector[i]) { continue; }
    loss_vector[i]->PreparePhysicsTable(*part);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle,
                                        G4VEmProcess* p)
{
  if (1 < verbose) {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " 
           << particle->GetParticleName() 
           << " and " << p->GetProcessName()
           << " run=" << run << " master=" << isMaster
           << G4endl;
  }

  // start initialisation for the first run
  if (-1 == run && nullptr != emConfigurator) {
    emConfigurator->PrepareModels(particle, p);
  }

  ResetParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* part,
                                             G4VMultipleScattering* p)
{
  if (1 < verbose) {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " 
           << part->GetParticleName() 
           << " and " << p->GetProcessName()
           << " run=" << run << " master=" << isMaster
           << G4endl;
  }

  // start initialisation for the first run
  if (-1 == run && nullptr != emConfigurator) {
    emConfigurator->PrepareModels(part, p);
  } 
  
  ResetParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::BuildPhysicsTable(const G4ParticleDefinition*)
{
  if (startInitialisation && nullptr != emConfigurator) {
    emConfigurator->Clear();
  }
  if (startInitialisation) { resetParam = true; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::WorkerPhysicsTables(
     const G4ParticleDefinition* aParticle,
     G4VEnergyLossProcess* p)
{
  if (1 < verbose) {
    G4cout << "### G4LossTableManager::WorkerPhysicsTable() for "
           << aParticle->GetParticleName()
           << " and process " << p->GetProcessName()
           << "  threadID=" << threadID << G4endl;
  }
  currentParticle = nullptr;
  if (startInitialisation) {
    if (-1 == run) {
      if (nullptr != emConfigurator) { emConfigurator->Clear(); }
      firstParticle = aParticle;
    }
    ++run;
    startInitialisation = false;
    
    resetParam = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::BuildPhysicsTable(
     const G4ParticleDefinition* aParticle,
     G4VEnergyLossProcess* p)
{
  if (1 < verbose) {
    G4cout << "### G4LossTableManager::BuildPhysicsTable() for "
           << aParticle->GetParticleName()
           << " and process " << p->GetProcessName()
           << " startInitialisation=" << startInitialisation 
	         << " id=" << G4Threading::G4GetThreadId()
           << G4endl;
  }
  currentParticle = nullptr;

  // initialisation of tables is performed once for all dedx processes
  if (startInitialisation) {
    if (-1 == run) {
      if ( nullptr != emConfigurator) { emConfigurator->Clear(); }
      firstParticle = aParticle;
    }
    ++run;
    startInitialisation = false;
    resetParam = true;
    if (1 < verbose) {
      G4cout << "    run=" << run << " atomDeexcitation: "
             << atomDeexcitation << G4endl;
    }

    // parallel initialisation only for the 1st run
    if (0 == run && theParameters->UseParallelInitialisation()) {
      if (1 < verbose) {
        G4cout << "### G4LossTableManager::BuildPhysicsTable() "
               << "parallel initialisation is enabled." << G4endl;
      }
      theRegistry->BuildTablesInParallel(verbose);
    }
    auto plist = theRegistry->ListForParallelBuild();

    // build all dedx, range, inverse range tables in the master thread
    for (G4int i=0; i<n_loss; ++i) {
      if (!loss_vector[i]->TablesAreBuilt() && nullptr == base_part_vector[i]) {
	      if (!plist.empty()) {
	        for (auto const & ptr : plist) {
	          if (ptr == part_vector[i]) {
	            continue;
	          }
	        }
	      }
        BuildTables(part_vector[i], true);
      }
    }
  }
  if (1 < verbose) {
    G4cout << "### G4LossTableManager::BuildPhysicsTable done for all particles."
           << G4endl;
    if (nullptr != subcutProducer) {
      G4cout << "     SubCutProducer <" << subcutProducer->GetName()
             << ">" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4LossTableManager::BuildTables(const G4ParticleDefinition* part, G4bool ok)
{
  if (1 < verbose) {
    G4cout << "### G4LossTableManager::BuildTables for "
           << part->GetParticleName() << G4endl;
  }

  std::vector<G4PhysicsTable*> t_list;
  t_list.reserve(5);
  std::vector<G4VEnergyLossProcess*> loss_list;
  loss_list.reserve(5);
  G4VEnergyLossProcess* em = nullptr;
  G4VEnergyLossProcess* p = nullptr;
  G4int iem = -1;
  G4PhysicsTable* dedx = nullptr;

  G4ProcessVector* pvec = part->GetProcessManager()->GetProcessList();
  G4int nvec = (G4int)pvec->size();

  for (G4int i = 0; i < n_loss; ++i) {
    p = loss_vector[i];
    G4bool yes = (part == part_vector[i]);

    // possible case of process sharing between particle/anti-particle
    if (!yes) {
      auto ptr = static_cast<G4VProcess*>(p);
      for (G4int j=0; j<nvec; ++j) {
        if (ptr == (*pvec)[j]) {
          yes = true;
          break;
        }
      }
    }
    // the process does not used for the particle
    if (!yes) { continue; }

    // ionisation process for this particle
    if (p->IsIonisationProcess() || nullptr == em) {
      em = p;
      iem = i;
    }

    // tables may be shared between particle/anti-particle
    if (!p->TablesAreBuilt()) {
      // build new table
      dedx = p->BuildDEDXTable(fRestricted);
      p->SetDEDXTable(dedx, fRestricted, ok);
      p->SetTablesAreBuilt(true);
    } else {
      // use existing table
      dedx = p->DEDXTable();
    }
    t_list.push_back(dedx);
    loss_list.push_back(p);
  }

  G4int n_dedx = (G4int)t_list.size();
  if (0 == n_dedx || nullptr == em) {
    G4cout << "### G4LossTableManager::BuildTables WARNING: no DEDX processes for " 
           << part->GetParticleName() << G4endl;
    return;
  }

  G4bool buildCSDA = theParameters->BuildCSDARange();

  if (2 < verbose) {
    G4cout << "     Start to build the sum of " << n_dedx << " processes"
           << " iem=" << iem << " em: " << em->GetProcessName()
           << " buildCSDARange=" << buildCSDA << G4endl;
  }

  dedx = em->DEDXTable();
  em->SetDEDXTable(dedx, fIsIonisation, true);
  auto tableBuilder = G4EmDataRegistry::Instance()->GetLossTableBuilder();
  tableBuilder->SetSplineFlag(em->Spline());

  if (1 < n_dedx) {
    dedx = nullptr;
    dedx = G4PhysicsTableHelper::PreparePhysicsTable(dedx);
    tableBuilder->BuildDEDXTable(dedx, t_list);
    em->SetDEDXTable(dedx, fRestricted, true);
  }

  G4PhysicsTable* range = em->RangeTableForLoss();
  if (nullptr == range) {
    range = G4PhysicsTableHelper::PreparePhysicsTable(range);
  }

  G4PhysicsTable* invrange = em->InverseRangeTable();
  if (nullptr == invrange) {
    invrange = G4PhysicsTableHelper::PreparePhysicsTable(invrange);
  }

  tableBuilder->BuildRangeTable(dedx, range);
  tableBuilder->BuildInverseRangeTable(range, invrange);

  em->SetRangeTableForLoss(range);
  em->SetInverseRangeTable(invrange);

  // build lambda table for all processes from the list
  for (auto & ptr : loss_list) {
    ptr->SetLambdaTable(ptr->BuildLambdaTable(fRestricted), ok);
  }

  // CSDA dedx and range tables
  if (buildCSDA) {
    std::vector<G4PhysicsTable*> listCSDA;
    listCSDA.reserve(n_dedx);
    for (auto & ptr : loss_list) {
      dedx = ptr->BuildDEDXTable(fTotal);
      ptr->SetDEDXTable(dedx, fTotal, true);
      listCSDA.push_back(dedx);
    }
    G4PhysicsTable* dedxCSDA = em->DEDXunRestrictedTable();
    if (1 < n_dedx) {
      dedxCSDA = nullptr;
      dedxCSDA = G4PhysicsTableHelper::PreparePhysicsTable(dedxCSDA);
      tableBuilder->BuildDEDXTable(dedxCSDA, listCSDA);
      em->SetDEDXTable(dedxCSDA, fTotal, true);
    }
    G4PhysicsTable* rCSDA = em->CSDARangeTable();
    if (nullptr == rCSDA) {
      rCSDA = G4PhysicsTableHelper::PreparePhysicsTable(rCSDA);
    }
    tableBuilder->BuildRangeTable(dedxCSDA, rCSDA);
    em->SetCSDARangeTable(rCSDA);
  }
  // fill derived processes
  CopyTables(part, em);

  if (1 < verbose) {
    G4cout << "### G4LossTableManager::BuildTables: Tables are built for "
           << part->GetParticleName() << " Nproc=" << n_dedx
           << " ionisation process: " << em->GetProcessName()
           << " idx=" << iem << " threadID=" << threadID << G4endl;
  }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::CopyTables(const G4ParticleDefinition* part,
                                    G4VEnergyLossProcess* base_proc)
{
  // base particle approach
  for (G4int j=0; j<n_loss; ++j) {
    if (part == base_part_vector[j]) {
      G4VEnergyLossProcess* proc = loss_vector[j];
      proc->SetTablesAreBuilt(true);
      proc->SetLambdaTable(base_proc->LambdaTable(), true);
      if (proc->IsIonisationProcess()) { 
        proc->SetDEDXTable(base_proc->IonisationTable(), fRestricted, true);
        proc->SetRangeTableForLoss(base_proc->RangeTableForLoss());
        proc->SetInverseRangeTable(base_proc->InverseRangeTable());
        proc->SetDEDXTable(base_proc->DEDXunRestrictedTable(), fTotal, true);
        proc->SetCSDARangeTable(base_proc->CSDARangeTable());
      }
      if (1 < verbose) {
         G4cout << "   CopyTables for " << proc->GetProcessName()
                << " for " << part_vector[j]->GetParticleName()
                << " base_part= " << part->GetParticleName()
                << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ParticleHaveNoLoss(
     const G4ParticleDefinition* aParticle)
{
  G4ExceptionDescription ed;
  ed << "Energy loss process not found for " << aParticle->GetParticleName() 
     << " !";
  G4Exception("G4LossTableManager::ParticleHaveNoLoss", "em0001",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEnergyLossProcess*>& 
G4LossTableManager::GetEnergyLossProcessVector() const
{
  return loss_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEmProcess*>&
G4LossTableManager::GetEmProcessVector() const
{
  return emp_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VMultipleScattering*>& 
G4LossTableManager::GetMultipleScatteringVector() const
{
  return msc_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<const G4ParticleDefinition*>& 
G4LossTableManager::GetParticleVector() const
{
  return part_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmSaturation* G4LossTableManager::EmSaturation() const
{
  return theParameters->GetEmSaturation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmConfigurator* G4LossTableManager::EmConfigurator()
{
  if (nullptr == emConfigurator) {
    emConfigurator = new G4EmConfigurator(verbose); 
  }
  return emConfigurator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectronIonPair* G4LossTableManager::ElectronIonPair()
{
  if (nullptr == emElectronIonPair) { 
    emElectronIonPair = new G4ElectronIonPair(verbose);
  }
  return emElectronIonPair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetNIELCalculator(G4NIELCalculator* ptr)
{
  if(nullptr != ptr && ptr != nielCalculator) {
    delete nielCalculator;
    nielCalculator = ptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NIELCalculator* G4LossTableManager::NIELCalculator()
{
  if (nullptr == nielCalculator) { 
    nielCalculator = new G4NIELCalculator(nullptr, verbose); 
  }
  return nielCalculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4LossTableManager::SetAtomDeexcitation(G4VAtomDeexcitation* p)
{
  if(atomDeexcitation != p) {
    delete atomDeexcitation;
    atomDeexcitation = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetSubCutProducer(G4VSubCutProducer* p) 
{
  if(subcutProducer != p) {
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
  if (dirName && physList) {
    G4String physListName = G4String(physList);
    G4String pathName = G4String(dirName) + "/" + physListName + ".rst";

    std::ofstream outFile;
    outFile.open(pathName);
   
    outFile << physListName << G4endl;
    outFile << std::string(physListName.length(), '=') << G4endl;

    std::vector<const G4ParticleDefinition*> particles {
        G4Gamma::Gamma(),
        G4Electron::Electron(),
        G4Positron::Positron(),
        G4Proton::Proton(),
        G4MuonPlus::MuonPlus(),
        G4MuonMinus::MuonMinus(),
      };
   
    std::vector<G4VEmProcess*> emproc_vector = GetEmProcessVector();
    std::vector<G4VEnergyLossProcess*> enloss_vector = 
      GetEnergyLossProcessVector();
    std::vector<G4VMultipleScattering*> mscat_vector =
      GetMultipleScatteringVector();
    
    for (auto theParticle : particles) {
      outFile << G4endl << "**" << theParticle->GetParticleName()
              << "**" << G4endl << G4endl << " .. code-block:: none" << G4endl;

      G4ProcessManager* pm = theParticle->GetProcessManager();
      G4ProcessVector*  pv = pm->GetProcessList();
      G4int plen = pm->GetProcessListLength();

      for (auto emproc : emproc_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == emproc) {
            outFile << G4endl;
            proc->ProcessDescription(outFile);
            break;
          }
        }
      }

      for (auto mscproc : mscat_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == mscproc) {
            outFile << G4endl;
            proc->ProcessDescription(outFile);
            break;
          }
        }
      }

      for (auto enlossproc : enloss_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == enlossproc) {
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

