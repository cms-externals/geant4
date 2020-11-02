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

// The manager for Hdf5 file output operations.

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5FileManager_h
#define G4Hdf5FileManager_h 1

#include "G4VTFileManager.hh"
#include "globals.hh"

#include "tools/hdf5/ntuple" // for hid_t

#include <memory>
#include <tuple>

using G4Hdf5File = std::tuple<hid_t, hid_t, hid_t>;

class G4Hdf5FileManager : public G4VTFileManager<G4Hdf5File>
{
  public:
    explicit G4Hdf5FileManager(const G4AnalysisManagerState& state);
    ~G4Hdf5FileManager();

    using G4VTFileManager<G4Hdf5File>::WriteFile;
    using G4VTFileManager<G4Hdf5File>::CloseFile;

    // Methods to manipulate output files
    virtual G4bool OpenFile(const G4String& fileName) final;
    virtual G4bool WriteFile() final;
    virtual G4bool CloseFile() final; 

    virtual G4String GetType() const final { return "Hdf5"; }
    virtual G4String GetFileType() const final { return "hdf5"; }

    // Set methods
    void  SetBasketSize(unsigned int basketSize);

    // Get methods
    std::shared_ptr<G4Hdf5File> GetFile() const;
    hid_t GetHistoDirectory() const;
    hid_t GetNtupleDirectory() const;
    unsigned int GetBasketSize() const; 
    
  protected:
    // // Methods derived from templated base class
    virtual std::shared_ptr<G4Hdf5File> CreateFileImpl(const G4String& fileName) final;
    virtual G4bool WriteFileImpl(std::shared_ptr<G4Hdf5File> file) final;
    virtual G4bool CloseFileImpl(std::shared_ptr<G4Hdf5File> file) final;    

  private:
    hid_t CreateDirectory(hid_t& file, const G4String& directoryName, 
             const G4String& objectType);
    hid_t GetFileHid() const;

    // constants
    static const G4String fgkDefaultDirectoryName;
    
    // data members
    std::shared_ptr<G4Hdf5File>  fFile;
    unsigned int fBasketSize;
};

// inline functions

//_____________________________________________________________________________
inline void G4Hdf5FileManager::SetBasketSize(unsigned int basketSize)  
{ fBasketSize = basketSize; }

//_____________________________________________________________________________
inline std::shared_ptr<G4Hdf5File> G4Hdf5FileManager::GetFile() const
{ return fFile; }

//_____________________________________________________________________________
inline unsigned int G4Hdf5FileManager::GetBasketSize() const
{ return fBasketSize; }

#endif
