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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4RootFileManager.hh"
#include "G4RootHnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/to"
#include "tools/zlib"

#include <iostream>
#include <cstdio>

using namespace tools;

//_____________________________________________________________________________
G4RootFileManager::G4RootFileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager<G4RootFile>(state),
   fFile(nullptr),
   fNofNtupleFiles(0),
   // fNtupleFiles(),
   // fMainNtupleDirectories(),
   fBasketSize(0),
   fBasketEntries(0)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4RootHnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4RootHnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4RootHnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4RootHnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4RootHnFileManager<histo::p2d>>(this);
}

//_____________________________________________________________________________
G4RootFileManager::~G4RootFileManager()
{}

//
// private methods
//

//_____________________________________________________________________________
tools::wroot::directory*  G4RootFileManager::CreateDirectory(
  std::shared_ptr<tools::wroot::file> rfile,
#ifdef G4VERBOSE
  const G4String& directoryName, const G4String& objectType) const
#else
  const G4String& directoryName, const G4String& /*objectType*/) const
#endif
{
  if ( ! rfile ) return nullptr;

  if ( directoryName == "" ) {
    // Do not create a new directory if its name is not set
    return &(rfile->dir());
  }  
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "directory for " + objectType, directoryName);
#endif
  
  auto directory = rfile->dir().mkdir(directoryName);
  if ( ! directory ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot create directory " << directoryName;
    G4Exception("G4RootFileManager::CreateDirectory()",
              "Analysis_W001", JustWarning, description);
    return nullptr;       
  }       
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("create", "directory for " + objectType, directoryName);
  }    
#endif
  return directory;
}

//_____________________________________________________________________________
std::shared_ptr<G4RootFile> 
G4RootFileManager::CreateFileImpl(const G4String& fileName)
{
  // create file
  std::shared_ptr<wroot::file> file = std::make_shared<wroot::file>(G4cout, fileName);
  file->add_ziper('Z',compress_buffer);
  file->set_compression(fState.GetCompressionLevel());
  
  if ( ! file->is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot create file " << fileName;
    G4Exception("G4RootAnalysisManager::CreateFileImpl()",
                "Analysis_W001", JustWarning, description);
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  }

  // create histo directory
  tools::wroot::directory* hdirectory 
    = CreateDirectory(file, fHistoDirectoryName, "histograms");
  if ( ! hdirectory ) {
    // Warning is issued in CreateDirectory
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  } 

  // create ntuple directory
  tools::wroot::directory* ndirectory 
    = CreateDirectory(file, fNtupleDirectoryName, "ntuples");
  if ( ! ndirectory ) {
    // Warning is issued in CreateDirectory
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  } 

  return std::make_shared<G4RootFile>(file, hdirectory, ndirectory);
}

//_____________________________________________________________________________
G4bool G4RootFileManager::WriteFileImpl(std::shared_ptr<G4RootFile> file)
{
// New prototype: called by G4TFileManager base classe

  if ( ! file ) return false;

  unsigned int n;
  return std::get<0>(*file)->write(n);
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFileImpl(std::shared_ptr<G4RootFile> file)    
{
// New prototype: called by G4TFileManager base classe

  // G4cout << "G4RootFileManager::CloseFileImpl " << file.get() << G4endl;

  if ( ! file ) return false;

  // close file
  std::get<0>(*file)->close();

  return true;
}

//
// public methods
//
//_____________________________________________________________________________
G4bool G4RootFileManager::OpenFile(const G4String& fileName)
{
// Open default file

  // Keep file name
  fFileName =  fileName;
  auto name = GetFullFileName();
  
  if ( fFile ) {    
    G4ExceptionDescription description;
    description 
      << "File " << fileName << " already exists."; 
    G4Exception("G4RootFileManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
  }

  // Create histograms file (on master)
  if ( fState.GetIsMaster() ) {
    // Create file (and save in in the file map)
    fFile = CreateTFile(name);
    if ( ! fFile ) {
      G4ExceptionDescription description;
      description 
        << "Failed to create file " << fileName; 
      G4Exception("G4RootFileManager::OpenFile()",
                  "Analysis_W001", JustWarning, description);
      return false;
    }
  }

  fLockDirectoryNames = true;
  fIsOpenFile = true;

  return true;
}

//_____________________________________________________________________________
G4bool G4RootFileManager:: WriteFile()
{
// Write default file

  auto fileName = GetFullFileName(fFileName);

  return G4VTFileManager<G4RootFile>::WriteFile(fileName);
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFile()
{
// Close default file

  auto fileName = GetFullFileName(fFileName);
  auto result = G4VTFileManager<G4RootFile>::CloseFile(fileName);

  fIsOpenFile = false;
  fFile.reset();

  return result;
}

//_____________________________________________________________________________
G4bool G4RootFileManager::OpenNtupleFiles()
{
  auto finalResult = true;

  for ( auto i = 0; i < fNofNtupleFiles; i++ ) {
    
    auto name = GetNtupleFileName(i);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "main ntuple file", name);
#endif

    // create new file
    // auto rfile = std::make_shared<wroot::file>(G4cout, name);
    auto rfile = CreateTFile(name);
    if ( ! rfile ) {
      G4ExceptionDescription description;
      description 
        << "Failed to create ntuple file " << name; 
      G4Exception("G4RootFileManager::OpenNtupleFiles()",
                  "Analysis_W001", JustWarning, description);
      return false;
    }

    // fNtupleFiles.push_back(rfile);
    // fMainNtupleDirectories.push_back(directory);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("create", "main ntuple file", name);
#endif

  }

  return finalResult;
}  

//_____________________________________________________________________________
std::shared_ptr<G4RootFile> 
G4RootFileManager::GetNtupleFile(G4int index) const
{ 
  if ( index==0 && ( ! fNofNtupleFiles ) ) return fFile;    

  if ( index < 0 || index >= fNofNtupleFiles ) {
    G4String inFunction = "G4RootFileManager::GetNtupleFile()";
    G4ExceptionDescription description;
    description << "      " << "ntuple file " << index << " does not exist.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return nullptr;         
  }

  return GetTFile(GetNtupleFileName(index));
}
