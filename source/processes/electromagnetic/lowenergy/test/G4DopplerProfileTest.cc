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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DataSetManagerTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 8 March 2008
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"
#include "G4VEMDataSet.hh"

int main()
{
  //  G4cout.setf( ios::scientific, ios::floatfield );

  // EADL data 
  // G4String file = "/fluor/binding";
  // G4ShellData* dataSet = new G4ShellData();

  // Doppler binding data
  G4String file = "/doppler/shell-doppler";
  G4DopplerProfile* dataSet = new G4DopplerProfile(1,100);

  //  dataSet->LoadData(file);
  dataSet->PrintData();
  
  // Dump all Doppler profiles
  for (G4int zDoppler=1; zDoppler<101; zDoppler++)
    {
      //    const G4VEMDataSet* profiles = dataSet->Profiles(zDoppler);
      //     profiles->PrintData();
    }

  G4cout << "Enter Z" << G4endl;
  G4int Z;
  G4cin >> Z;

  G4int n = dataSet->NumberOfProfiles(Z);
  G4cout << "Z = " << Z << " has " << n << " shells" << G4endl;
  const G4VEMDataSet* profiles = dataSet->Profiles(Z);
  profiles->PrintData();

  G4cout << "Enter shell index " << G4endl;
  G4int i;
  G4cin >> i;

  const G4VEMDataSet* profile = dataSet->Profile(Z,i);
  profile->PrintData();

  // Test random selection on single and composite profile
  G4double rando = profile->RandomSelect();
  rando = profiles->RandomSelect(i);

  rando = dataSet->RandomSelectMomentum(Z,i);

  /*
  G4int id = dataSet->ShellId(Z,i);
  G4double e = dataSet->BindingEnergy(Z,i) / keV;

  G4cout << "Shell id = " << id 
	 << ", Binding energy = " << e << " keV" <<G4endl;

  std::vector<G4double> idVector = dataSet->ShellIdVector(Z);
  for (G4int ind=0; ind<n; ind++)
    {
      G4int idx = (G4int) idVector[ind];
      G4cout << "Id vector(" << ind << ") = " << idx << G4endl; 
    }

  // Select random shell

  
  for (G4int iter=0; iter<100; iter++)
    {
      G4int shellIndex = dataSet->RandomSelectShell(Z);
      G4cout << "Random selected shell: " << shellIndex << G4endl;
      // G4cout << "random = " << random 
      //        << ", i = " << shellIndex 
      //        << ", prob = " << prob[shellIndex] << G4endl;
    }

  G4cout << "Dump data (1) or stop (2)" << G4endl;
  G4int k;
  G4cin >> k;

  if (k == 1) dataSet->PrintData();  
  */

  delete dataSet;

  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








