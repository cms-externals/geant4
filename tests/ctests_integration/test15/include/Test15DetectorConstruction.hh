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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef Test15DetectorConstruction_h
#  define Test15DetectorConstruction_h 1

#  include "G4VUserDetectorConstruction.hh"
#  include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

// class Test15LeadSD;
// class Test15SampleSD;

class Test15DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    Test15DetectorConstruction();
    ~Test15DetectorConstruction();

  public:

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4VPhysicalVolume* GetWorldVolume();
    G4VPhysicalVolume& GetWorldVolumeAddress() const;

  private:

    void DefineMaterials();

    void CreateGrichine();
    void CreateLayer1();
    void CreateLayer2();
    void CreateLayer3();
    void CreateLayer4();
    void CreateLayer5();
    void CreateLayer6();
    void CreateLayer7();
    void CreateLayer8();
    void CreateLayer9();
    void CreateLayer10();
    void CreateLayer11();

    G4bool olap_test;
    G4int replica_idx_A;
    G4int replica_idx_B;
    G4int replica_idx_C;
    G4int replica_sample;

    G4double blockWidthA;
    G4double blockLengthA;
    G4double blockHeightA;
    G4double blockWidthB;
    G4double blockLengthB;
    G4double blockHeightB;
    G4double blockWidthC;
    G4double blockLengthC;
    G4double blockHeightC;

    G4Material* world_mat;  // materials used
    G4Material* lab_mat;
    G4Material* lead_mat;
    G4Material* sample_mat;
    G4Material* sampleSphere_mat;
    G4Material* sampleTube_mat;

    // G4double worldRadius;                // sizes
    // G4double worldHeight;

    G4LogicalVolume* world_log;  // pointers
    G4VPhysicalVolume* world_phys;

    G4LogicalVolume* lab_log;
    G4VPhysicalVolume* lab_phys;

    G4LogicalVolume* lead_log;
    G4VPhysicalVolume* lead_phys;

    G4LogicalVolume* blockA_log;
    G4VPhysicalVolume* blockA_phys;
    G4LogicalVolume* blockB_log;
    G4VPhysicalVolume* blockB_phys;
    G4LogicalVolume* beamblockB_log;
    G4VPhysicalVolume* beamblockB_phys;
    G4LogicalVolume* blockC_log;
    G4VPhysicalVolume* blockC_phys;

    G4LogicalVolume* sample_log;
    G4VPhysicalVolume* sample_phys;
    G4VPhysicalVolume* sample_phys2;
    G4LogicalVolume* sampleSphere_log;
    G4LogicalVolume* sampleSphere_log2;

    G4LogicalVolume* sampleTube_log;
    G4VPhysicalVolume* sampleTube_phys;
};

#endif
