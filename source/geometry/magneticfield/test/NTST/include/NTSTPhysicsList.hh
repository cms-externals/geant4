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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTPhysicsList_h
#define NTSTPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class NTSTLooperDeath;

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4eMultipleScattering;
class G4MuMultipleScattering;
class G4hMultipleScattering;

class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4MuIonisation;
class G4MuBremsstrahlung;
class G4MuPairProduction;

class G4hIonisation;

class NTSTPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTPhysicsList: public G4VUserPhysicsList
{
public:
  NTSTPhysicsList();
  ~NTSTPhysicsList();
  //
  // set methods
  //
  inline void SetBgsTran(const G4bool on) { useBgsTran = on; }
  void SetMinimumEnergyCut(const G4double e);
  void SetMaximumEnergyCut(const G4double e);
  inline void SetLengthCut(const G4double e) { Cut = e ; }
  void SetLooperCut(const G4double e);
    
protected:
  void AddBgsTransportation();
  
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
 
  void SetCuts();
    
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBarions();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructLeptHad();
  void ConstructHad();
    
public:
  // this method allow to set on/off the processes
  void SetStatusEmProcess();
    
private:
  G4bool   useBgsTran;
  G4double MinimumEnergyCut;
  G4double MaximumEnergyCut;
  G4double Cut;
  G4double LooperCut;

  // processes

  NTSTLooperDeath*       theLooperDeath;

  G4PhotoElectricEffect* thePhotoElectricEffect;
  G4ComptonScattering*   theComptonScattering;
  G4GammaConversion*     theGammaConversion;
    
  G4eMultipleScattering* theeminusMultipleScattering;
  G4eIonisation*         theeminusIonisation;
  G4eBremsstrahlung*     theeminusBremsstrahlung;
    
  G4eMultipleScattering* theeplusMultipleScattering;
  G4eIonisation*         theeplusIonisation;
  G4eBremsstrahlung*     theeplusBremsstrahlung;
  G4eplusAnnihilation*   theeplusAnnihilation;
    
  //  G4hMultipleScattering* theHadronMultipleScattering;
  NTSTPhysicsListMessenger* physicsListMessenger;
};

#endif



