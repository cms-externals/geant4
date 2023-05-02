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
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4LevelReader
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
// Helper class to read Geant4 nuclear level data  
// 

#ifndef G4LEVELREADER_HH
#define G4LEVELREADER_HH 1

#include "globals.hh"
#include "G4LevelManager.hh"
#include <iosfwd>
#include <vector>

class G4NuclearLevelData;

class G4LevelReader 
{

public:

  explicit G4LevelReader(G4NuclearLevelData*);

  // create run manager using G4LEVELGAMMADATA data for Z and A
  const G4LevelManager* CreateLevelManager(G4int Z, G4int A);

  // create run manager using whatever data
  const G4LevelManager* MakeLevelManager(G4int Z, G4int A,
					 const G4String& filename);

  inline void SetVerbose(G4int val);

  G4LevelReader(const G4LevelReader & right) = delete;  
  const G4LevelReader& operator=(const G4LevelReader &right) = delete;
  G4bool operator==(const G4LevelReader &right) const = delete;
  G4bool operator!=(const G4LevelReader &right) const = delete;
  
private:

  G4bool ReadData(std::istringstream& dataFile, G4double& x);

  G4bool ReadDataItem(std::istream& dataFile, G4double& x);

  G4bool ReadDataItem(std::istream& dataFile, G4float& x);

  G4bool ReadDataItem(std::istream& dataFile, G4int& x);

  G4bool ReadDataItem(std::istream& dataFile, G4String& x);
  
  const std::vector<G4float>* NormalizedICCProbability(G4int Z);

  const G4LevelManager* LevelManager(G4int Z, G4int A,
				     std::ifstream& infile);

  G4NuclearLevelData* fData;

  G4double fEnergy = 0.;
  G4double fTransEnergy = 0.;
  G4double fTime = 0.;
  G4double fTimeFactor;

  G4float fProb = 0.f;
  G4float fSpin = 0.f;
  G4float fAlpha = 0.f;
  G4float fAlphaMax;
  G4float fRatio = 0.f;
  G4float fNorm1 = 0.f;
  G4float fICC[10] = {0.f};

  G4int nbufmax = 20;
  G4int nbuf1 = 14;
  G4int nbuf2 = 8;

  G4int fVerbose = 1;
  G4int fLevelMax = 632;
  G4int fTransMax = 145;
  G4int ntrans = 0;
  G4int i1 = 0;
  G4int i2 = 0;
  G4int k = 0;
  G4int kk = 0;
  G4int tnum = 0;

  char buffer[20] = {' '};
  char buff1[14] = {' '};
  char buff2[8] = {' '};
  char bufp[3] = {' '};

  std::vector<G4double> vEnergy;
  std::vector<G4int> vSpin;
  std::vector<const G4NucLevel*> vLevel;

  std::vector<G4int> vTrans;
  std::vector<G4float> vRatio;
  std::vector<G4float> vGammaCumProbability;
  std::vector<G4float> vGammaProbability;
  std::vector<const std::vector<G4float>*> vShellProbability;

  G4String fPol = "  ";
  G4String fDirectory;
};

inline void G4LevelReader::SetVerbose(G4int val)
{
  fVerbose = val;
}

#endif
