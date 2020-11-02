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

#include "G4RunManagerFactory.hh"
#include "G4EnvironmentUtils.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4TaskRunManager.hh"
#include "G4Threading.hh"
#include "templates.hh"

//============================================================================//

namespace
{
  // failure message
  static void fail(const std::string& _prefix, const std::string& _name,
                   const std::set<std::string>& _opts, G4int _num)
  {
    G4ExceptionDescription msg;
    msg << _prefix << ": \"" << _name << "\". "
        << "Must be one of: ";
    std::stringstream ss;
    for(const auto& itr : _opts)
      ss << ", \"" << itr << "\"";
    msg << ss.str().substr(2);
    auto mnum = std::string("RunManagerCreator000") + std::to_string(_num);
    G4Exception("G4RunManagerFactory::CreateRunManager", mnum.c_str(),
                FatalException, msg);
  }

  static G4RunManager* master_run_manager              = nullptr;
  static G4MTRunManager* mt_master_run_manager         = nullptr;
  static G4RunManagerKernel* master_run_manager_kernel = nullptr;
}  // namespace

//============================================================================//

G4RunManager* G4RunManagerFactory::CreateRunManager(G4RunManagerType _type,
                                                    G4VUserTaskQueue* _queue,
                                                    G4bool fail_if_unavail,
                                                    G4int nthreads)
{
  G4RunManager* rm = nullptr;

  auto rm_type  = G4GetEnv<std::string>("G4RUN_MANAGER_TYPE", GetName(_type));
  auto force_rm = G4GetEnv<std::string>("G4FORCE_RUN_MANAGER_TYPE", "",
                                        "Forcing G4RunManager type...");

  // if forced, set rm_type to specification, enable failure if not available,
  // and set enum value to default so that
  if(force_rm.length() > 0)
  {
    rm_type         = force_rm;
    _type           = GetType(force_rm);
    fail_if_unavail = true;
  }
  else if(rm_type.empty())
  {
    rm_type = GetDefault();
    _type   = GetType(rm_type);
  }

  // if string is empty
  if(rm_type.empty())
  {
    if(_type == G4RunManagerType::Default)
      rm_type = GetDefault();
    else
      rm_type = GetName(_type);
  }

  // convert default to a type
  if(_type == G4RunManagerType::Default)
    _type = GetType(rm_type);

  // if type is still default and fail
  if(fail_if_unavail && _type == G4RunManagerType::Default)
    fail("Run manager type is not available", rm_type, GetOptions(), 0);

  auto opts = GetOptions();
  if(fail_if_unavail && opts.find(GetName(_type)) == opts.end())
  {
    // if type selected is not available and fail_if_unavail is true, fail
    fail("Run manager type is not available", GetName(_type), opts, 1);
  }
  else if(!fail_if_unavail && opts.find(GetName(_type)) == opts.end())
  {
    // if type selected is not available and fail_if_unavail is false, set to
    // default
    rm_type = GetDefault();
    _type   = GetType(rm_type);
  }

  switch(_type)
  {
    case G4RunManagerType::Serial:
      rm = new G4RunManager();
      break;
    case G4RunManagerType::MT:
#if defined(G4MULTITHREADED)
      rm = new G4MTRunManager();
#endif
      break;
    case G4RunManagerType::Tasking:
#if defined(G4MULTITHREADED)
      rm = new G4TaskRunManager(_queue, false);
#endif
      break;
    case G4RunManagerType::TBB:
#if defined(G4MULTITHREADED) && defined(GEANT4_USE_TBB)
      rm = new G4TaskRunManager(_queue, true);
#endif
      break;
    case G4RunManagerType::Default:
      break;
  }

  if(!rm)
    fail("Failure creating run manager", GetName(_type), GetOptions(), 3);

  auto mtrm = dynamic_cast<G4MTRunManager*>(rm);
  if(nthreads > 0 && mtrm)
    mtrm->SetNumberOfThreads(nthreads);

  master_run_manager        = rm;
  mt_master_run_manager     = mtrm;
  master_run_manager_kernel = rm->kernel;

  G4ConsumeParameters(_queue);
  return rm;
}

//============================================================================//

std::string G4RunManagerFactory::GetDefault()
{
#if defined(G4MULTITHREADED)
  // For version 10.7, default is set to MT
  //  return "Tasking";
  return "MT";
#else
  return "Serial";
#endif
}

//============================================================================//

std::set<std::string> G4RunManagerFactory::GetOptions()
{
  static auto _instance = []() {
    std::set<std::string> options = { "Serial" };
#if defined(G4MULTITHREADED)
    options.insert({ "MT", "Tasking" });
#  if defined(GEANT4_USE_TBB)
    options.insert("TBB");
#  endif
#endif
    return options;
  }();
  return _instance;
}

//============================================================================//

G4RunManagerType G4RunManagerFactory::GetType(const std::string& key)
{
  // IGNORES CASE!
  static const auto opts = std::regex::icase;

  if(std::regex_match(key, std::regex("^(Serial).*", opts)))
    return G4RunManagerType::Serial;
  else if(std::regex_match(key, std::regex("^(MT).*", opts)))
    return G4RunManagerType::MT;
  else if(std::regex_match(key, std::regex("^(Task).*", opts)))
    return G4RunManagerType::Tasking;
  else if(std::regex_match(key, std::regex("^(TBB).*", opts)))
    return G4RunManagerType::TBB;

  return G4RunManagerType::Default;
}

//============================================================================//

std::string G4RunManagerFactory::GetName(G4RunManagerType _type)
{
  switch(_type)
  {
    case G4RunManagerType::Serial:
      return "Serial";
    case G4RunManagerType::MT:
      return "MT";
    case G4RunManagerType::Tasking:
      return "Tasking";
    case G4RunManagerType::TBB:
      return "TBB";
    default:
      break;
  };
  return "";
}

//============================================================================//

G4RunManager* G4RunManagerFactory::GetMasterRunManager()
{
#if !defined(G4MULTITHREADED)
  // if serial build just return G4RunManager
  return G4RunManager::GetRunManager();
#else
  // if the application used G4RunManagerFactory to create the run-manager
  if(master_run_manager)
    return master_run_manager;

  // if the application did not use G4RunManagerFactory and is MT
  if(G4Threading::IsMultithreadedApplication())
  {
    auto mt_rm = GetMTMasterRunManager();
    if(mt_rm)
      return mt_rm;
  }

  // if the application did not use G4RunManagerFactory and is serial
  return G4RunManager::GetRunManager();
#endif
}

//============================================================================//

G4MTRunManager* G4RunManagerFactory::GetMTMasterRunManager()
{
#if defined(G4MULTITHREADED)
  // if the application used G4RunManagerFactory to create the run-manager
  if(mt_master_run_manager)
    return mt_master_run_manager;

  // if the application did not use G4RunManagerFactory
  if(G4Threading::IsMultithreadedApplication())
  {
    auto task_rm = G4TaskRunManager::GetMasterRunManager();
    if(task_rm)
      return task_rm;
    return G4MTRunManager::GetMasterRunManager();
  }
#endif

  return nullptr;
}

//============================================================================//

G4RunManagerKernel* G4RunManagerFactory::GetMasterRunManagerKernel()
{
#if !defined(G4MULTITHREADED)
  // if serial build just return G4RunManager
  return G4RunManager::GetRunManager()->kernel;
#else
  // if the application used G4RunManagerFactory to create the run-manager
  if(master_run_manager_kernel)
    return master_run_manager_kernel;

  // if the application did not use G4RunManagerFactory and is MT
  if(G4Threading::IsMultithreadedApplication())
  {
    auto mt_rm = GetMTMasterRunManager();
    if(mt_rm)
      return mt_rm->kernel;
  }

  // if the application did not use G4RunManagerFactory and is serial
  return G4RunManager::GetRunManager()->kernel;
#endif
}

//============================================================================//
