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
// ------------------------------------------------------------
// Geant4 class header file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#ifndef AXPETMaterial_h
#define AXPETMaterial_h 1

#include "globals.hh"

class G4Material;

class AXPETMaterial
{
  public:
  AXPETMaterial();
  ~AXPETMaterial();
  
  public:
  void  DefineMaterials();
  G4Material* GetMat(const G4String&); // Returns the material

  private:
  G4Material* aluminum;
  G4Material* matAir;
  G4Material* LYSO;
  G4Material* EJ280;
  G4Material* Vacuum;
  
};
#endif
