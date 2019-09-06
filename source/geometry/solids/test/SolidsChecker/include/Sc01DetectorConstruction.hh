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
//  GEANT 4 class header file 
//
//      This class is a class derived from G4VUserDetectorConstruction
//      for constructing all particles and processes.
//
//  History
//        first version              09 Sept. 1998 by S.Magni
//        modified for geometry test  11.02.04 V. Grichine 
// ------------------------------------------------------------

#ifndef Sc01DetectorConstruction_h
#define Sc01DetectorConstruction_h 1

#include "Sc01DetectorMessenger.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VisAttributes;

class Sc01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Sc01DetectorConstruction();
    ~Sc01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     G4VPhysicalVolume* SelectDetector (const G4String& val);
     G4VPhysicalVolume* SelectTubeSector();
     void               SwitchDetector();
     void               CleanGeometry();
     void               SetupGeometry();
     void               SetMaterial();
     G4double           GetHallSize() { return fHallSize; }

     G4LogicalVolume* GetConePolycone();

  private:

     Sc01DetectorMessenger* detectorMessenger;

     G4VSolid* Solid;
     G4LogicalVolume* LogicalVolume;
     G4VPhysicalVolume *WorldVolume, *PlacedVolume;
     G4RotationMatrix* rot;
     G4Box *b1, *b2;

     G4Material* Water;
     G4Material* Water1;

     G4OpticalSurface* aSurface;
     G4LogicalBorderSurface *bSurface1, *bSurface2;
     G4VisAttributes* visAttr;
     G4double fHallSize;
};

#endif
