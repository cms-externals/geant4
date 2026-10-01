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
/// \file PSFReader.hh
/// \brief Header file for PSFReader for Molecular DNA simulation

//Created: Ji Won Choi, Adrien Paillet

#ifndef PSFREADER_HH
#define PSFREADER_HH

#include "G4RootAnalysisReader.hh"
#include "G4PrimaryVertex.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Structure to gather data about particles red in the PSF 
struct PhaseSpaceParticle {
  G4int pdg;
  G4ThreeVector pos;
  G4ThreeVector dir;
  G4double eng;
  G4double weight;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// singleton pattern
// the class must be available from the main
// and each thread must share the same vector
class PSFReader {

public:
  PSFReader() = default;
  ~PSFReader() = default;
  static PSFReader &Instance() {
    static PSFReader instance;
    return instance;
  }

  // Actually reads the PSF and stores the data in g_particles
  void InitializePhaseSpace();

  // get/set public functions
  inline G4bool IsLoaded() const { return g_isLoaded; }
  inline size_t GetNumberOfEntries() const { return g_particles.size(); }
  inline G4String GetPhaseSpacePattern() const { return fpattern; }
  inline G4String GetPSFsamplingMode() const { return fPSFsamplingMode; }
  inline G4long GetNofOriginal() const { return fnumOfOriginal; }
  inline G4String GetPSFNTupleName() const     { return fPSFNTupleName;     }

  inline void SetPhaseSpacePattern(const G4String &pattern) { fpattern = pattern; }
  inline void SetPSFsamplingMode(const G4String &samplingMode) { fPSFsamplingMode = samplingMode; }
  inline void SetNofOriginal(const G4long &nOfOriginal) { fnumOfOriginal = nOfOriginal; }
  inline void SetPSFNTupleName(const G4String& PSFNTupleName)  {  fPSFNTupleName = PSFNTupleName;  }

  // must be called by the GeneratePrimaries method
  G4PrimaryVertex *GetNextPrimary();

private:
  // PSF-reading functions
  void ReadRootFiles(const G4String &pattern);

  void ReadTextFiles(const G4String &pattern);

  void Read_GRAS_CSV_Files(const G4String &pattern);

  // index-modifying functions
  void Shuffle() const;

  void SampleWithReplacement() const;

  void OversampleWithReplacement(const size_t &nOfEvents) const;

  // position modifying function
  void mapOnASphere(G4double offset_x, G4double offset_y, G4double offset_z, G4double radius_x,
                    G4double radius_y, G4double radius_z);

  // Phase-space data shared across all worker threads.
  std::vector<PhaseSpaceParticle> g_particles;
  G4bool g_isLoaded = false;
  std::atomic<size_t> g_nextIndex{0};
  G4String fpattern = "";
  G4long fnumOfOriginal;
  G4String fPSFsamplingMode = "Default";
  G4String fPSFNTupleName = "PSF";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
