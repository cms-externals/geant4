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
//
#include "G4Material.hh"
#include "G4Element.hh"

//#include "G4PEEffectModel.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4BetheHeitlerModel.hh"

#include "G4eeToTwoGammaModel.hh"

#include "G4MollerBhabhaModel.hh"
//#include "G4eBremsstrahlungModel.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"

#include "G4MuBetheBlochModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4NistManager.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeIonisationCrossSection.hh"
#include "G4AtomicShellEnumerator.hh"
#include "G4LivermoreIonisationCrossSection.hh"
#include "G4PenelopeIonisationXSHandler.hh"
#include "G4SystemOfUnits.hh"

int main() {

  G4UnitDefinition::BuildUnitsTable();

  // define materials
  //
  G4double Z, A;

  G4Element* AuEl = new G4Element ("Gold","Au",79.,196.97*g/mole);
  G4Element* HEl = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* OEl = new G4Element("Oxygen","O",8.,15.9994*g/mole);
  G4Element* AlEl = new G4Element("Aluminum","Al",13,26.98*g/mole);

  G4Material* Gold = new G4Material("Gold",19.3*g/cm3,1);
  Gold->AddElement(AuEl,1);

  G4Material* Iodine =
  new G4Material("Iodine", Z=53., A=126.90*g/mole, 4.93*g/cm3);

  G4Material* Water = new G4Material("Water",1.0*g/cm3,2);
  Water->AddElement(HEl,2);
  Water->AddElement(OEl,1);
  
  G4Material* Aluminum = new G4Material("Aluminum",2.7*g/cm3,1);
  Aluminum->AddElement(AlEl,100*perCent);

  G4Material* Water2 = new G4Material("WaterFractionMass",1.0*g/cm3,2);
  Water2->AddElement(HEl,11.189*perCent); 
  Water2->AddElement(OEl,88.811*perCent);
  
  //Polyethylene092:

  G4Material* Polyethylene092 = new G4Material("Polyethylene", 0.92*g/cm3,
					       2);
  G4NistManager* man = G4NistManager::Instance();
   
  Polyethylene092->AddElement(man->FindOrBuildElement("C"), 0.1428);
  Polyethylene092->AddElement(man->FindOrBuildElement("H"), 0.8572);

  //Polyethylene092->AddMaterial(man->FindOrBuildMaterial("G4_C"), 0.1428);
  //Polyethylene092->AddMaterial(man->FindOrBuildMaterial("G4_H"), 0.8572);




  //Joe Perl test
  G4Material* TestMaterial = new G4Material("TestMaterial", 1.0*g/cm3,
					    22);
  TestMaterial->AddElement(man->FindOrBuildElement("Ag"), 0.000786);
  TestMaterial->AddElement(man->FindOrBuildElement("As"), 0.002689);
  TestMaterial->AddElement(man->FindOrBuildElement("B"), 0.097731);
  TestMaterial->AddElement(man->FindOrBuildElement("Ba"), 0.033471);
  TestMaterial->AddElement(man->FindOrBuildElement("Be"), 0.00037);
  TestMaterial->AddElement(man->FindOrBuildElement("Cd"), 0.000038);
  TestMaterial->AddElement(man->FindOrBuildElement("Co"), 0.026509);
  TestMaterial->AddElement(man->FindOrBuildElement("Cr"), 0.366652);
  TestMaterial->AddElement(man->FindOrBuildElement("Hg"), 0.000018);
  //G4cout << "UNO" << G4endl;
  TestMaterial->AddElement(man->FindOrBuildElement("Li"), 0.00474);
  TestMaterial->AddElement(man->FindOrBuildElement("Mn"), 0.18652);
  TestMaterial->AddElement(man->FindOrBuildElement("Mo"), 0.0006);
  TestMaterial->AddElement(man->FindOrBuildElement("Ni"), 0.19674);
  TestMaterial->AddElement(man->FindOrBuildElement("Pb"), 0.032769);
  TestMaterial->AddElement(man->FindOrBuildElement("Pd"), 0.003392);
  TestMaterial->AddElement(man->FindOrBuildElement("Sb"), 0.001226);
  //G4cout << "DUE" << G4endl;
  TestMaterial->AddElement(man->FindOrBuildElement("Se"), 0.001437);
  TestMaterial->AddElement(man->FindOrBuildElement("Sn"), 0.000358);
  TestMaterial->AddElement(man->FindOrBuildElement("Sr"), 0.03264);
  TestMaterial->AddElement(man->FindOrBuildElement("Tl"), 0.000805);
  TestMaterial->AddElement(man->FindOrBuildElement("V"), 0.010029);
  TestMaterial->AddElement(man->FindOrBuildElement("W"), 0.000479);

  G4Material* material = Water2;
  G4MaterialCutsCouple* theCouple = new G4MaterialCutsCouple(material);
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  


  G4PenelopeOscillatorManager* theManager = 
    G4PenelopeOscillatorManager::GetOscillatorManager();
  theManager->SetVerbosityLevel(2);  

  //G4PenelopeOscillatorTable* theTable1 = 
  //theManager->GetOscillatorTableIonisation(material); 
  theManager->Dump(material);
  
  //Instantiate PIXE cross sections handlers.
  G4PenelopeIonisationCrossSection* pixeXS = new G4PenelopeIonisationCrossSection();
  pixeXS->SetVerbosityLevel(0);
  G4LivermoreIonisationCrossSection *livXS = new G4LivermoreIonisationCrossSection();

  //Dump shell cross sections
  G4double Emin = 500*eV, Emax = 100*MeV;
  G4int nstep = 5000;
  G4int iZ = 79;
  G4double logStep = std::log10(Emax/Emin)/((G4double) nstep);
  std::ofstream ff("test.dat");
  for (G4int i=1;i<=nstep;i++)
    {
      G4double energy = Emin*std::pow(10,i*logStep);
      G4double p1 = pixeXS->CrossSection(iZ,
					 fKShell,
					 energy,
					 511.0*keV,
					 material);
      G4double p2 = livXS->CrossSection(iZ,
					 fKShell,
					 energy,
					 511.0*keV,
					 material);
    
      if (p1> 1e-3*barn || p2>1e-3*barn)
	ff << energy/keV << " " << p1/barn << " " << p2/barn << G4endl;
    }
  ff.close();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                                 
return EXIT_SUCCESS;
}
