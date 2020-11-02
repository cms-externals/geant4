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

// Author: Ivana Hrivnacova, 15/09/2020 (ivana@ipno.in2p3.fr)

#include "G4RootNtupleFileManager.hh"
#include "G4RootFileManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootPNtupleManager.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <cstdio>

using namespace G4Analysis;
using std::make_shared;

G4RootNtupleFileManager* G4RootNtupleFileManager::fgMasterInstance = nullptr;

//_____________________________________________________________________________
G4RootNtupleFileManager::G4RootNtupleFileManager(const G4AnalysisManagerState& state)
 : G4VNtupleFileManager(state, "root"),
   fIsInitialized(false),
   fNofNtupleFiles(0),
   fNtupleRowWise(false),
   fNtupleRowMode(true),
   fNtupleMergeMode(G4NtupleMergeMode::kNone),
   fNtupleManager(nullptr),
   fSlaveNtupleManager(nullptr),
   fFileManager(nullptr)
{
  if ( G4Threading::IsMasterThread() && fgMasterInstance ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4RootNtupleFileManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4RootNtupleFileManager::G4RootNtupleFileManager()",
                "Analysis_F001", FatalException, description);
  }
  if ( G4Threading::IsMasterThread() ) fgMasterInstance = this;

  // Do not merge ntuples by default
  // Merging may require user code migration as analysis manager
  // must be created both on master and workers.
  auto mergeNtuples = false;
  SetNtupleMergingMode(mergeNtuples, fNofNtupleFiles);
}

//_____________________________________________________________________________
G4RootNtupleFileManager::~G4RootNtupleFileManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;  
}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleMergingMode(G4bool mergeNtuples, 
                                                 G4int nofNtupleFiles)

{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("set", "ntuple merging mode", "");
#endif

  auto canMerge = true;

  // Illegal situations
  if ( mergeNtuples && ( ! G4Threading::IsMultithreadedApplication() ) ) {
    if ( nofNtupleFiles > 0 ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Merging ntuples is not applicable in sequential application." 
        << G4endl 
        << "      " << "Setting was ignored.";
        G4Exception("G4RootNtupleFileManager::SetNtupleMergingMode()",
                  "Analysis_W013", JustWarning, description);
    }
    canMerge = false;      
  }

  // Illegal situations
  if ( mergeNtuples && G4Threading::IsMultithreadedApplication() &&
       ( ! fgMasterInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "Merging ntuples requires G4AnalysisManager instance on master." 
      << G4endl 
      << "      " << "Setting was ignored.";
      G4Exception("G4RootNtupleFileManager::SetNtupleMergingMode()",
                "Analysis_W013", JustWarning, description);
    canMerge = false;      
  }

  G4String mergingMode;
  if ( ( ! mergeNtuples ) || ( ! canMerge ) ) {
    fNtupleMergeMode = G4NtupleMergeMode::kNone;
    mergingMode = "G4NtupleMergeMode::kNone";
  }
  else {
    // Set the number of reduced ntuple files
    fNofNtupleFiles = nofNtupleFiles;
  
    // Check the number of reduced ntuple files
    // if ( fNofNtupleFiles < 0 || fNofNtupleFiles > nofThreads ) {
    if ( fNofNtupleFiles < 0  ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Number of reduced files must be [0, nofThreads]."
        << G4endl 
        << "      " << "Cannot set  " <<  nofNtupleFiles
        // << " files when nofThreads is " << nofThreads << G4endl   
        << " files" << G4endl   
        << "      " << "Ntuples will be merged in a single file.";
        G4Exception("G4RootNtupleFileManager::SetNtupleMergingMode()",
                  "Analysis_W013", JustWarning, description);
      fNofNtupleFiles = 0;
    }
  
    // Forced merging mode
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    if ( isMaster ) {
      fNtupleMergeMode = G4NtupleMergeMode::kMain;
      mergingMode = "G4NtupleMergeMode::kMain";
    } else {
      fNtupleMergeMode = G4NtupleMergeMode::kSlave;    
      mergingMode = "G4NtupleMergeMode::kSlave";
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()
      ->Message("set", "ntuple merging mode", mergingMode);
#endif
}

//_____________________________________________________________________________
G4int G4RootNtupleFileManager::GetNtupleFileNumber()
{
  if ( ! fNofNtupleFiles ) return 0;

  G4int nofMainManagers = fNofNtupleFiles;
  if ( ! nofMainManagers ) nofMainManagers = 1;

  // Debug - check G4Threading::GetNumberOfRunningWorkerThreads()
  G4cout << "In GetNtupleFileNumber: "
         << G4Threading::GetNumberOfRunningWorkerThreads() << G4endl;

  auto fileNumber = G4Threading::G4GetThreadId() % nofMainManagers;
  return fileNumber;
}

// 
// public methods
//

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleMerging(G4bool mergeNtuples, 
                                          G4int  nofNtupleFiles)

{
  if ( fIsInitialized ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Cannot change merging mode." << G4endl
        << "      " << "The function must be called before OpenFile().";
      G4Exception("G4RootNtupleFileManager::SetNtupleMerging",
                "Analysis_W013", JustWarning, description);
      return;
  }

  // Set ntuple merging mode 
  SetNtupleMergingMode(mergeNtuples, nofNtupleFiles);
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleRowWise(G4bool rowWise, G4bool rowMode) 
{
#ifdef G4VERBOSE
  // Print info even when setting makes no effect
  // (as we do not get printed the default setting in the output)
  G4String rowWiseMode;
  if ( rowWise ) {
    rowWiseMode = "row-wise with extra branch";
  } 
  else if ( rowMode ) {
    rowWiseMode = "row-wise";
  }
  else {
    rowWiseMode = "column-wise"; 
  }

  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("set", "ntuple merging row mode", rowWiseMode);
#endif

  // Do nothing if the mode is not changed
  if ( fNtupleRowWise == rowWise && fNtupleRowMode == rowMode ) return;

  fNtupleRowWise = rowWise;
  fNtupleRowMode = rowMode;

  if ( fNtupleManager ) {
    fNtupleManager->SetNtupleRowWise(rowWise, rowMode);
  }

  if ( fSlaveNtupleManager ) {
    fSlaveNtupleManager->SetNtupleRowWise(rowWise, rowMode);
  }
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetBasketSize(unsigned int basketSize) 
{
  fFileManager->SetBasketSize(basketSize);
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetBasketEntries(unsigned int basketEntries)
{
  fFileManager->SetBasketEntries(basketEntries);
}

//_____________________________________________________________________________
std::shared_ptr<G4VNtupleManager> G4RootNtupleFileManager::CreateNtupleManager()
{
// Create and return ntuple manager.
// If ntuple marging is activated, managers of different types are created
// on master/worker.

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple managers", "");
#endif

  // Check that file manager and anaysis manager are set !

  std::shared_ptr<G4VNtupleManager> activeNtupleManager = nullptr;
  switch ( fNtupleMergeMode )
  {
    case G4NtupleMergeMode::kNone:
      fNtupleManager 
        = make_shared<G4RootNtupleManager>(
            fState, fBookingManager, 0, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;

    case G4NtupleMergeMode::kMain: {
      G4int nofMainManagers = fNofNtupleFiles;
      if ( ! nofMainManagers ) nofMainManagers = 1;
             // create one manager if merging required into the histos & profiles files
      fNtupleManager 
        = make_shared<G4RootNtupleManager>(
            fState, fBookingManager, nofMainManagers, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;
    }

    case G4NtupleMergeMode::kSlave:
      fNtupleManager = fgMasterInstance->fNtupleManager;
        // The master class is used only in Get* functions
      auto mainNtupleManager 
        = fNtupleManager->GetMainNtupleManager(GetNtupleFileNumber()); 
      fSlaveNtupleManager 
        = make_shared<G4RootPNtupleManager>(
            fState, fBookingManager, mainNtupleManager, fNtupleRowWise, fNtupleRowMode); 
      activeNtupleManager = fSlaveNtupleManager;
      break;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("create", "ntuple managers", "");
#endif

  fIsInitialized = true; 

  return activeNtupleManager;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtOpenFile(const G4String& fileName)
{
  // Check if fNtupleBookingManager is set

  auto finalResult = true;

  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone )  {

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "analysis file", fileName);
#endif

    // Create file if it does not exist
    auto file = fFileManager->GetTFile(fileName, false);
    if (! file ) {
      file = fFileManager->CreateTFile(fileName);
    }
    fNtupleManager->SetNtupleFile(file);

    fNtupleManager->CreateNtuplesFromBooking(
      fBookingManager->GetNtupleBookingVector());

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "analysis file", fileName, finalResult);
#endif
  
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {

#ifdef G4VERBOSE
    // G4String fileName = fFileManager->GetFullFileName();
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "main analysis file", fileName);
#endif

    fFileManager->SetNofNtupleFiles(fNofNtupleFiles);
    // Create file if it does not exist
    auto file = fFileManager->GetTFile(fileName, false);
    if (! file ) {
      file = fFileManager->CreateTFile(fileName);
    }
    fNtupleManager->SetNtupleFile(file);
    auto result = fFileManager->OpenNtupleFiles();
    finalResult = finalResult && result;

    fNtupleManager->CreateNtuplesFromBooking(
      fBookingManager->GetNtupleBookingVector());

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "main analysis file", fileName, finalResult);
#endif  
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave )  {
    G4cout << "Going to create slave ntuples from main" << G4endl;
    // No file is open by Slave manager
    fSlaveNtupleManager->CreateNtuplesFromMain();
  }

  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::WriteNtuple()
{
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone ) {
    return true;
  }
  
  auto finalResult = true;

  G4String ntupleType;
  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) ntupleType = "main ntuples";
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) ntupleType = "slave ntuples";

#ifdef G4VERBOSE 
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("merge", ntupleType, "");
#endif

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Merge();
    finalResult = result && finalResult;
  }  
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) {
    auto result = fSlaveNtupleManager->Merge();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("merge", ntupleType, "");
#endif
  
  return finalResult;
}
    
//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtWrite() 
{

  auto finalResult = true;

  // Ntuples
  auto result = WriteNtuple();
  finalResult = finalResult && result;

  // File
  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave )  {
    // notify that file is not empty
    G4String fileName = fFileManager->GetFullFileName();
    result = fFileManager->SetIsEmpty(fileName, false);
    finalResult = finalResult && result;
  }

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtCloseFile(G4bool reset)
{
  auto finalResult = true;

  if ( ! reset ) {
    // The ntuples must be always reset when closing file) 
    auto result = Reset();
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4RootNtupleFileManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
    }
    finalResult = finalResult && result;
  }

  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave )  {
    // close file
    auto result = fFileManager->CloseFile(); 
    finalResult = finalResult && result;
  }

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone || 
       fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  

  return finalResult;
}

