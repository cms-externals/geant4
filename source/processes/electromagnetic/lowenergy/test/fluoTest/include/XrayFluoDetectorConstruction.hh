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

#ifndef XrayFluoDetectorConstruction_h
#define XrayFluoDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class XrayFluoDetectorMessenger;
class XrayFluoHPGeSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  XrayFluoDetectorConstruction();
  ~XrayFluoDetectorConstruction();
  
public:
  
  G4VPhysicalVolume* Construct();
  
  void UpdateGeometry();
  
public:
  
  void PrintApparateParameters(); 
  
  G4double GetWorldSizeZ()           {return WorldSizeZ;}; 
  G4double GetWorldSizeXY()          {return WorldSizeXY;};
  
  G4double GetDeviceThickness()      {return DeviceThickness;}; 
  G4double GetDeviceSizeX()          {return DeviceSizeX;};
  G4double GetDeviceSizeY()          {return DeviceSizeY;};
  G4double GetPixelSizeXY()          {return PixelSizeXY;};
  G4double GetContactSizeXY()        {return ContactSizeXY;};
  
  G4int GetNbOfPixels()              {return NbOfPixels;}; 
  G4int GetNbOfPixelRows()           {return NbOfPixelRows;}; 
  G4int GetNbOfPixelColumns()        {return NbOfPixelColumns;}; 
  
  G4Material* GetOhmicPosMaterial()  {return OhmicPosMaterial;};
  G4double    GetOhmicPosThickness() {return OhmicPosThickness;};      
  
  G4Material* GetOhmicNegMaterial()  {return OhmicNegMaterial;};
  G4double    GetOhmicNegThickness() {return OhmicNegThickness;};      
  
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};  
  const G4VPhysicalVolume* GetHPGe()        {return physiHPGe;};
  const G4VPhysicalVolume* GetSample()     {return physiSample;};
  const G4VPhysicalVolume* GetDia1()        {return physiDia1;};
  const G4VPhysicalVolume* GetDia3()        {return physiDia3;};
  
  const G4VPhysicalVolume* GetphysiPixel()  {return physiPixel;};           
  const G4VPhysicalVolume* GetOhmicPos()    {return physiOhmicPos;};
  const G4VPhysicalVolume* GetOhmicNeg()    {return physiOhmicNeg;};
  
  
private:
  
  G4double           DeviceSizeX;
  G4double           DeviceSizeY;
  G4double           DeviceThickness;
  
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidHPGe; //pointer to the solid Sensor
  G4LogicalVolume*   logicHPGe; //pointer to the logical Sensor
  G4VPhysicalVolume* physiHPGe; //pointer to the physical Sensor
  
  G4Box*             solidSample;    //pointer to the solid Sample
  G4LogicalVolume*   logicSample;    //pointer to the logical Sample
  G4VPhysicalVolume* physiSample;    //pointer to the physical Sample
  
  G4Tubs*             solidDia1; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia1; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia1; //pointer to the physical Diaphragm 
  
  G4Tubs*             solidDia3; //pointer to the solid  Diaphragm
  G4LogicalVolume*   logicDia3; //pointer to the logical  Diaphragm
  G4VPhysicalVolume* physiDia3; //pointer to the physical Diaphragm  
  
  G4Box*             solidOhmicPos;
  G4LogicalVolume*   logicOhmicPos; 
  G4VPhysicalVolume* physiOhmicPos; 
  
  G4Box*             solidOhmicNeg;
  G4LogicalVolume*   logicOhmicNeg; 
  G4VPhysicalVolume* physiOhmicNeg;     
  
  G4Box*             solidPixel;   
  G4LogicalVolume*   logicPixel;  
  G4VPhysicalVolume* physiPixel;    
 
  //pointers to the materials used 
  G4Material*        OhmicPosMaterial;
  G4Material*        OhmicNegMaterial; 
  G4Material*        pixelMaterial;
  G4Material*        sampleMaterial;
  G4Material*        Dia1Material;
  G4Material*        Dia3Material;
  G4Material*        defaultMaterial;
  G4double           OhmicPosThickness;
  
  G4double           OhmicNegThickness;
  
  G4int              PixelCopyNb;
  G4int              NbOfPixels;
  G4int              NbOfPixelRows;
  G4int              NbOfPixelColumns;
  G4double           PixelThickness;
  
  G4double           PixelSizeXY;
  G4double           ContactSizeXY;
  
public:

  G4Material* GetSampleMaterial()  {return sampleMaterial;};
  G4Material* GetPixelMaterial()  {return pixelMaterial;}; 
  G4Material* GetDia1Material()  {return Dia1Material;}; 
  G4Material* GetDia3Material()  {return Dia3Material;}; 
  
private:

  G4double           SampleThickness;
  G4double           SampleSizeXY;
  G4double           SiSizeXY; 
  G4double           Dia1Thickness;
  G4double           Dia1SizeXY;
  G4double           Dia3Thickness;
  G4double           Dia3SizeXY;
  G4double           DiaInnerSize;
  G4double           Dia3InnerSize;
  G4double           DistSi;
public: 
  
  
  G4double GetSampleThickness()         {return SampleThickness;};
  G4double GetSampleSizeXY()              {return SampleSizeXY;};
  
  G4double GetDia1Thickness()         {return Dia1Thickness;};
  G4double GetDia1SizeXY()              {return Dia1SizeXY;};
  
  G4double GetDia3Thickness()         {return Dia3Thickness;};
  G4double GetDia3SizeXY()              {return Dia3SizeXY;};
  
  
private:

  
  G4double           ThetaHPGe;
  G4double           ThetaDia1;
  G4double           ThetaDia3;
  
  G4double           DistDe;
  G4double           DistDia;
  G4double           Dia3Dist;
  G4double           PhiHPGe;
  G4double           PhiDia1;
  G4double           PhiDia3;
  G4double AlphaDia1;
  G4double AlphaDia3;
  
  
  G4RotationMatrix   zRotPhiHPGe;
  G4RotationMatrix   zRotPhiDia1;
  G4RotationMatrix   zRotPhiDia3;
  G4double           WorldSizeXY;
  G4double           WorldSizeZ;
  
  
  XrayFluoDetectorMessenger* detectorMessenger; //pointer to the Messenger

  XrayFluoHPGeSD* HPGeSD;  //pointer to the sensitive detector
  
private:
  
  void DefineMaterials();
  G4VPhysicalVolume* ConstructApparate();

  //calculates some quantities used to construct geometry
  void ComputeApparateParameters();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void XrayFluoDetectorConstruction::ComputeApparateParameters()
{     
  // Compute derived parameters of the apparate
  
  DeviceThickness = PixelThickness+OhmicNegThickness+OhmicPosThickness;
  DeviceSizeY =NbOfPixelRows*max(ContactSizeXY,PixelSizeXY);
  DeviceSizeX =NbOfPixelColumns*max(ContactSizeXY,PixelSizeXY);
  
  WorldSizeZ = (2 * (DistDe  +1.4142 *(max(max(DeviceThickness,DeviceSizeY), DeviceSizeX))))+5*m; 
  WorldSizeXY = 2 * (DistDe +1.4142 *Dia1SizeXY)+5*m;
  
}

#endif






