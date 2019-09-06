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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------

//for the documentation of this class see also G4EMDataSet class

#ifndef FluoDataSet_hh
#define FluoDataSet_hh 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4VEMDataSet.hh"

class G4VDataSetAlgorithm;

class XrayFluoDataSet : public G4VEMDataSet {

public:

  XrayFluoDataSet(G4int Z,
		  G4DataVector* points, 
	      G4DataVector* values,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);
  
  XrayFluoDataSet(G4int Z,
		  const G4String& dataFile,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  ~XrayFluoDataSet();
 
  //find the value corresponding to the energy e in the set
  //identified by id
  G4double FindValue(G4double e, G4int id = 0) const;
  
  void PrintData() const;

  const G4DataVector& GetEnergies(G4int i) const { return *energies; }
  const G4DataVector& GetData(G4int i) const { return *data; }

private:

  void LoadData(const G4String& dataFile);
  G4int z;
  G4int FindBinLocation(G4double energy) const;

  G4DataVector* energies; // Owned pointer
  G4DataVector* data;     // Owned pointer

  const G4VDataSetAlgorithm* algorithm; // Not owned pointer 
  
  G4double unit1;
  G4double unit2;

  size_t numberOfBins;

};
 
#endif
