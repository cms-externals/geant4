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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4XmlFileManager.hh"
#include "G4XmlHnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/waxml/begend"

using namespace tools;

//_____________________________________________________________________________
G4XmlFileManager::G4XmlFileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager<std::ofstream>(state),
   fHnFile(nullptr)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4XmlHnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4XmlHnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4XmlHnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4XmlHnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4XmlHnFileManager<histo::p2d>>(this);
}

//_____________________________________________________________________________
G4XmlFileManager::~G4XmlFileManager()
{}

// 
// protected methods
//

//_____________________________________________________________________________
std::shared_ptr<std::ofstream> G4XmlFileManager::CreateFileImpl(const G4String& fileName)
{
  std::shared_ptr<std::ofstream> file = std::make_shared<std::ofstream>(fileName);
  if ( file->fail() ) {
    file = nullptr;
    G4ExceptionDescription description;
    description << "      " << "Cannot create file " << fileName;
    G4Exception("G4XmlFileManager::CreateFileImpl()",
              "Analysis_W001", JustWarning, description);
    return nullptr;
  }

  waxml::begin(*file);
  return file;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::WriteFileImpl(std::shared_ptr<std::ofstream> /*file*/)
{
  // Nothing to be done here
  return true;  
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseFileImpl(std::shared_ptr<std::ofstream> file)    
{
  if ( ! file ) return false;

  // close file
  waxml::end(*file);
  file->close(); 

  return true;
}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlFileManager::OpenFile(const G4String& fileName)
{
  // Keep and locks file name
  fFileName =  fileName;
  auto name = GetFullFileName(fFileName);
  G4cout << "G4XmlFileManager::OpenFile " << name << G4endl;

  // Create histograms file (on master)
  if ( fState.GetIsMaster() ) {
    // Create file (and save in in the file map (on master only)
    auto file = CreateTFile(name);
    if ( file && (! fHnFile) ) {
      // Backward compatibility
      // If a default file is not yet defined
      fHnFile = file;
    }
  }

  fIsOpenFile = true;

  return true;
}  
  
//_____________________________________________________________________________
G4bool G4XmlFileManager::WriteFile() 
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseFile()
{
  
  auto fileName = GetFullFileName(fFileName);
  auto result = G4VTFileManager<std::ofstream>::CloseFile(fileName);

  // fLockFileName = false;
  fIsOpenFile = false;

  return result;  
} 

//_____________________________________________________________________________
G4bool G4XmlFileManager::CreateNtupleFile(
  G4TNtupleDescription<waxml::ntuple>* ntupleDescription)
{
  auto ntupleFileName = ntupleDescription->fFileName;
  if ( ntupleFileName.empty() ) {
    auto ntupleName = ntupleDescription->fNtupleBooking.name();
    ntupleFileName = GetNtupleFileName(ntupleName);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "ntuple file", ntupleFileName);
#endif

  auto ntupleFile = new std::ofstream(ntupleFileName);
  if ( ntupleFile->fail() ) {
    delete ntupleFile;
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << ntupleFileName;
    G4Exception("G4XmlFileManager::CreateNtupleFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  
  waxml::begin(*ntupleFile);
  ntupleDescription->fFile = ntupleFile;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("create", "ntuple file", ntupleFileName);
#endif

  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseNtupleFile(
  G4TNtupleDescription<waxml::ntuple>* ntupleDescription)
{
  // Do nothing if there is no file
  if ( ! ntupleDescription->fFile ) return true;

  auto ntupleFileName = ntupleDescription->fFileName;
  if ( ntupleFileName.empty() ) {
    auto ntupleName = ntupleDescription->fNtupleBooking.name();
    ntupleFileName = GetNtupleFileName(ntupleName);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("close", "ntuple file", ntupleFileName);
#endif

  // close file
  waxml::end(*(ntupleDescription->fFile));
  ntupleDescription->fFile->close(); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("close", "ntuple file", ntupleFileName);
#endif

  return true; 
}
