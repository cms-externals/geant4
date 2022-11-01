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
// ------------------------------------------------------------
//
//  To print cross sections per atom and mean free path for simple material
//
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Timer.hh"

#include "globals.hh"

#include "G4PenelopeRayleighModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4PenelopeOscillatorManager.hh"

#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"

#include "G4SystemOfUnits.hh"

int main() {

  // define materials
  //
  G4double Z, A;

  G4Element* AuEl = new G4Element ("Gold","Au",79.,196.97*g/mole);
  G4Element* HEl = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* OEl = new G4Element("Oxygen","O",8.,15.9994*g/mole);
  G4Element* NEl = new G4Element("Nitrogen","N",7,14.0067*g/mole);

  G4Material* Gold = new G4Material("Gold",19.3*g/cm3,1);
  Gold->AddElement(AuEl,1);

  G4Material* Iodine =
  new G4Material("Iodine", Z=53., A=126.90*g/mole, 4.93*g/cm3);

  G4Material* Water = new G4Material("Water",1.0*g/cm3,2);
  Water->AddElement(HEl,2);
  Water->AddElement(OEl,1);
  
  G4Material* Water2 = new G4Material("WaterFractionMass",1.0*g/cm3,2);
  Water2->AddElement(HEl,11.189*perCent);
  Water2->AddElement(OEl,88.811*perCent);

  G4Material* Air = new G4Material("Air",1.290*mg/cm3,2,
				   kStateGas,300.00*kelvin,1.0*atmosphere);
  Air->AddElement(OEl,20.0*perCent);
  Air->AddElement(NEl,80.0*perCent);
  //Air->AddElement(OEl,1);
  //Air->AddElement(NEl,4);

  G4Material* material = Water;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // initialise gamma processes (models)
  //
  G4ParticleDefinition* gamma = G4Gamma::Gamma();  
  G4PenelopeComptonModel* comp = new 
    G4PenelopeComptonModel();

  // compute CrossSection per atom and MeanFreePath
  //
  G4double Emin = 100*eV, Emax = 100*GeV;
  G4int nstep = 5000;
  G4double logStep = std::log10(Emax/Emin)/((G4double) nstep);

  G4DataVector dummy;
  comp->Initialise(gamma,dummy);
  comp->SetVerbosityLevel(0);

  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4cout << atomDensity/(1./cm3) << " atoms/cm3" << G4endl;

  G4Timer* theTimer = new G4Timer();
  theTimer->Start();

  std::ofstream file("testCo.dat");
  for (G4int i=1;i<=nstep;i++)
    {
      G4double energy = Emin*std::pow(10,i*logStep);
      G4double XSVOLUME = comp->CrossSectionPerVolume(material,gamma,energy);
      G4double meanSigma = XSVOLUME/atomDensity;
      file << energy/eV << " " << XSVOLUME/(1./cm) << " " << 
	meanSigma/barn << G4endl;
	//ray->ComputeMeanFreePath(gamma,energy,material)/cm << 
	//G4endl;
    }	
  file.close();
  theTimer->Stop();

  G4cout << "Time: " << theTimer->GetSystemElapsed() << " " << 
    theTimer->GetRealElapsed() << G4endl;

  G4cout << "Water at 100 keV: " << 
    comp->CrossSectionPerVolume(Water,gamma,100*keV)*cm << " 1/cm" 
	 << G4endl;

  G4PenelopeOscillatorManager* theManager = 
    G4PenelopeOscillatorManager::GetOscillatorManager();
  theManager->Dump(material); 


  /*
  //Now test XSs for Oxygen
  G4int ZZ=79;
  std::ofstream ff("G4GoldPE.dat");
  for (G4int i=1;i<=nstep;i++)
    {
      G4double energy = Emin*std::pow(10,i*logStep);
      G4double XS0 = conv->ComputeCrossSectionPerAtom(gamma,
						      energy,ZZ);
      ff << energy/eV << " " << XS0/barn << G4endl;
    }	
  ff.close();
  */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                                 
return EXIT_SUCCESS;
}
