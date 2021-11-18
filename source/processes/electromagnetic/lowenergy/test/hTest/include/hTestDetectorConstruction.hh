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
#ifndef hTestDetectorConstruction_h
#define hTestDetectorConstruction_h 1

// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestDetectorConstruction -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
#include "globals.hh"
#include "G4ios.hh"

class hTestCalorimeterSD;
class hTestEventAction;
class hTestDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    hTestDetectorConstruction();
   ~hTestDetectorConstruction();

  public:
     
    inline void SetAbsorberMaterial(const G4String& v) {nameMatAbsorber = v;};
    inline G4Material* GetAbsorberMaterial()  const {return AbsorberMaterial;};
    void SetNumberOfAbsorbers(G4int);     
    inline G4int GetNumberOfAbsorbers() const {return NumberOfAbsorbers;};
    void SetAbsorberThickness(G4double);     
    inline G4double  GetAbsorberThickness() const {return AbsorberThickness;};
    void SetAbsorberSizeXY   (G4double);            
    inline G4double    GetAbsorberSizeXY() const {return SizeXY;};
    inline void SetWorldMaterial(const G4String& v) {nameMatWorld = v;};
    inline G4Material* GetWorldMaterial() const {return WorldMaterial;};
    void SetWorldSizeZ   (G4double);
    inline G4double    GetWorldSizeZ() const {return WorldSizeZ;}; 
    void SetGap(G4double val);
    inline G4double    GetGap() const {return gap;}; 
    inline void SetNumberOfEvents(G4int val) {nEvents = val;};
    inline G4int GetNumberOfEvents() const {return nEvents;};
    inline void SetEventAction(hTestEventAction* p) {theEvent = p;};
    inline hTestEventAction* GetEventAction() const {return theEvent;};
    inline void SetVerbose(G4int val) {myVerbose = val;};
    inline G4int GetVerbose() const {return myVerbose;};
    inline void SetNumAbsorbersSaved(G4int val) {nAbsSaved = val;};
    inline G4int GetNumAbsorbersSaved() const {return nAbsSaved;};
    inline void SetFirstEventToDebug(G4int val) {nFirstEvtToDebug = val;};
    inline G4int GetFirstEventToDebug() const {return nFirstEvtToDebug;};
    inline void SetLastEventToDebug(G4int val) {nLastEvtToDebug = val;};
    inline G4int GetLastEventToDebug() const {return nLastEvtToDebug;};
    inline G4double  GetMaxDeltaEnergy() const {return maxDelta;};
    inline void SetMaxDeltaEnergy(G4double val) {maxDelta = val;};            

    inline const G4VPhysicalVolume* GetPhysWorld() const {return physWorld;};
    inline const G4LogicalVolume*   GetAbsorber()  const {return logicAbs;};

    void SetMagField(G4double,G4int);
     
    G4VPhysicalVolume* Construct();
    void UpdateGeometry();
    void PrintGeomParameters(); 
                                     
  private:

    // Methods

     void DefineMaterials();
     G4Material* GetMaterial(const G4String&);
     void ComputeGeomParameters();
     G4VPhysicalVolume* ConstructGeometry();     
     void GeometryIsChanged();     
     void MaterialIsChanged();     

    // Members
     
     G4String           nameMatAbsorber;
     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;  // 
     G4double           SizeXY;
     G4double           gap;

     G4int              NumberOfAbsorbers;

     G4double           ZposAbs;
     
     G4String           nameMatWorld;
     G4Material*        WorldMaterial;
     G4double           WorldSizeZ;     
                
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physWorld;     //pointer to the physical World

     G4Box*             solidAbs;      //pointer to the solid Absorber
     G4LogicalVolume*   logicAbs;      //pointer to the logical Absorber
     G4VPhysicalVolume* physAbs;       //pointer to the physical Absorber
     
     G4UniformMagField* magField;      //pointer to the magnetic field
     
     hTestDetectorMessenger* detectorMessenger;  
     hTestCalorimeterSD* calorimeterSD;  
     hTestEventAction* theEvent;

     G4int myVerbose;      
     G4int nEvents;
     G4bool detIsConstructed;
     G4int nAbsSaved;
     G4int nFirstEvtToDebug;
     G4int nLastEvtToDebug;
     G4double maxDelta;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

