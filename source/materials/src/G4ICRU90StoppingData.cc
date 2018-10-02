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
// $Id: G4ICRU90StoppingData.cc
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:    G4ICRU90StoppingData
//
// Description:  Data on electroninc stopping power from ICRU 90
//
// Author:       Lucas Norberto Burigo
//
// Creation date: 03.09.2018
//
// Modifications: 25.09.2018 V.Ivanchenko adopted for material sub-library
// 
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ICRU90StoppingData.hh" 

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU90StoppingData::G4ICRU90StoppingData() : isInitialized(false)
{
  // 1st initialisation 
  for(size_t i=0; i<nvectors; ++i) { 
    materials[i]    = nullptr; 
    sdata_proton[i] = nullptr; 
    sdata_alpha[i]  = nullptr;
  }
  FillData();

  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU90StoppingData::~G4ICRU90StoppingData()
{
  for(size_t i=0; i<nvectors; ++i) { 
    delete sdata_proton[i]; 
    delete sdata_alpha[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4ICRU90StoppingData::Initialise()
{
  if(isInitialized) { return; }
  // this method may be called several times during initialisation
  G4int nmat = G4Material::GetNumberOfMaterials();
  if(nmat == (G4int)nvectors) { return; }

  static const G4String nameNIST_ICRU90[3] = {"G4_WATER","G4_AIR","G4_GRAPHITE"};

  // loop via material list to add extra data
  for(G4int i=0; i<nmat; ++i) {
    const G4Material* mat = (*(G4Material::GetMaterialTable()))[i];

    G4bool isThere = false;  
    for(G4int j=0; j<nvectors; ++j) {
      if(mat == materials[j]) {
	isThere = true;
	break;
      }
    }
    if(!isThere) {
      // check list of NIST materials
      G4String mname = mat->GetName();
      for(G4int j=0; j<nvectors; ++j) {
        if(mname == nameNIST_ICRU90[j]) {
          materials[j] = mat;
          break;
	}
      }
    }
    isInitialized = (materials[0] && materials[1] && materials[2]);
    if(isInitialized) { return; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU90StoppingData::GetElectronicDEDXforProton(
         const G4Material* mat, G4double kinEnergy) const
{
  G4int idx = GetIndex(mat);
  return (idx < 0) ? 0.0 : GetDEDX(sdata_proton[idx], kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU90StoppingData::GetElectronicDEDXforAlpha(
         const G4Material* mat, G4double scaledKinEnergy) const
{
  G4int idx = GetIndex(mat);
  return (idx < 0) ? 0.0 : GetDEDX(sdata_alpha[idx], scaledKinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU90StoppingData::FillData() 
{
  G4double T0_proton[57] = {  0.001000, 0.001500, 0.002000, 0.003000, 0.004000, 0.005000, 0.006000, 0.008000, 0.010000, 0.015000, 0.020000, 0.030000, 0.040000, 0.050000, 0.060000, 0.080000, 0.100000, 0.150000, 0.200000, 0.300000, 0.400000, 0.500000, 0.600000, 0.800000, 1.000000, 1.500000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 8.000000, 10.000000, 15.000000, 20.000000, 30.000000, 40.000000, 50.000000, 60.000000, 80.000000, 100.000000, 150.000000, 200.000000, 300.000000, 400.000000, 500.000000, 600.000000, 800.000000, 1000.000000, 1500.000000, 2000.000000, 3000.000000, 4000.000000, 5000.000000, 6000.000000, 8000.000000, 10000.000000  };

  G4double T0_alpha[49] = {  0.001000, 0.001500, 0.002000, 0.003000, 0.004000, 0.005000, 0.006000, 0.008000, 0.010000, 0.015000, 0.020000, 0.030000, 0.040000, 0.050000, 0.060000, 0.080000, 0.100000, 0.150000, 0.200000, 0.300000, 0.400000, 0.500000, 0.600000, 0.800000, 1.000000, 1.500000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 8.000000, 10.000000, 15.000000, 20.000000, 30.000000, 40.000000, 50.000000, 60.000000, 80.000000, 100.000000, 150.000000, 200.000000, 300.000000, 400.000000, 500.000000, 600.000000, 800.000000, 1000.000000};

  static const G4float e0_proton[57] = {  119.700000f, 146.700000f, 169.300000f, 207.400000f, 239.500000f, 267.800000f, 293.300000f, 338.700000f, 378.700000f, 450.400000f, 506.700000f, 590.500000f, 648.300000f, 687.700000f, 713.200000f, 734.100000f, 729.000000f, 667.200000f, 592.200000f, 476.500000f, 401.200000f, 349.800000f, 312.100000f, 258.700000f, 222.700000f, 168.200000f, 137.000000f, 101.700000f, 81.920000f, 69.050000f, 59.940000f, 47.810000f, 40.040000f, 28.920000f, 22.930000f, 16.520000f, 13.120000f, 10.980000f, 9.514000f, 7.618000f, 6.441000f, 4.815000f, 3.975000f, 3.117000f, 2.686000f, 2.431000f, 2.265000f, 2.069000f, 1.962000f, 1.850000f, 1.820000f, 1.828000f, 1.861000f, 1.898000f, 1.934000f, 1.998000f, 2.052000f  };

  static const G4float e0_alpha[49] = {  87.500000f, 108.600000f, 126.700000f, 157.300000f, 183.500000f, 206.700000f, 227.900000f, 265.900000f, 299.600000f, 372.300000f, 434.300000f, 539.500000f, 629.000000f, 708.100000f, 779.800000f, 906.800000f, 1018.000000f, 1247.000000f, 1429.000000f, 1693.000000f, 1861.000000f, 1961.000000f, 2008.000000f, 2002.000000f, 1922.000000f, 1626.000000f, 1381.000000f, 1071.000000f, 885.800000f, 760.600000f, 669.600000f, 545.300000f, 463.400000f, 342.300000f, 274.700000f, 200.100000f, 159.300000f, 133.200000f, 115.000000f, 91.190000f, 76.140000f, 54.930000f, 43.690000f, 31.840000f, 25.630000f, 21.790000f, 19.170000f, 15.830000f, 13.790000f  };

  static const G4float e1_proton[57] = {  133.700000f, 163.800000f, 189.100000f, 231.600000f, 267.500000f, 299.000000f, 327.600000f, 378.200000f, 422.900000f, 503.600000f, 567.300000f, 662.800000f, 729.000000f, 774.000000f, 802.600000f, 824.100000f, 814.500000f, 736.000000f, 658.500000f, 543.500000f, 464.300000f, 406.500000f, 362.400000f, 299.700000f, 257.400000f, 193.400000f, 156.900000f, 116.000000f, 93.190000f, 78.420000f, 68.010000f, 54.170000f, 45.320000f, 32.690000f, 25.890000f, 18.640000f, 14.790000f, 12.380000f, 10.720000f, 8.578000f, 7.250000f, 5.417000f, 4.470000f, 3.504000f, 3.018000f, 2.731000f, 2.544000f, 2.323000f, 2.203000f, 2.065000f, 2.017000f, 1.998000f, 2.010000f, 2.029000f, 2.050000f, 2.090000f, 2.124000f  };

  static const G4float e1_alpha[49] = {  98.910000f, 122.800000f, 143.100000f, 177.600000f, 206.900000f, 233.000000f, 256.800000f, 299.300000f, 337.000000f, 418.100000f, 487.100000f, 603.900000f, 703.000000f, 790.500000f, 869.600000f, 1009.000000f, 1131.000000f, 1383.000000f, 1582.000000f, 1873.000000f, 2062.000000f, 2178.000000f, 2240.000000f, 2256.000000f, 2190.000000f, 1877.000000f, 1599.000000f, 1239.000000f, 1021.000000f, 874.700000f, 768.700000f, 623.900000f, 529.000000f, 389.400000f, 311.800000f, 226.700000f, 180.200000f, 150.600000f, 130.000000f, 103.000000f, 85.930000f, 61.940000f, 49.230000f, 35.860000f, 28.850000f, 24.520000f, 21.560000f, 17.800000f, 15.500000f  };

  static const G4float e2_proton[57] = {  118.500000f, 145.100000f, 167.600000f, 205.300000f, 237.000000f, 265.000000f, 290.300000f, 335.200000f, 374.800000f, 435.300000f, 481.100000f, 550.300000f, 611.000000f, 663.000000f, 697.600000f, 726.700000f, 726.600000f, 671.200000f, 598.600000f, 483.300000f, 407.100000f, 354.600000f, 315.900000f, 262.400000f, 226.100000f, 171.100000f, 139.400000f, 103.500000f, 83.260000f, 70.100000f, 60.810000f, 48.450000f, 40.550000f, 29.250000f, 23.170000f, 16.680000f, 13.230000f, 11.070000f, 9.591000f, 7.675000f, 6.486000f, 4.843000f, 3.994000f, 3.125000f, 2.687000f, 2.426000f, 2.255000f, 2.050000f, 1.936000f, 1.806000f, 1.760000f, 1.744000f, 1.756000f, 1.775000f, 1.796000f, 1.834000f, 1.868000f  };

  static const G4float e2_alpha[49] = {  192.300000f, 228.900000f, 259.000000f, 308.300000f, 348.900000f, 384.000000f, 415.300000f, 469.900000f, 517.000000f, 615.100000f, 695.500000f, 826.200000f, 932.700000f, 1024.000000f, 1104.000000f, 1240.000000f, 1354.000000f, 1574.000000f, 1731.000000f, 1929.000000f, 2027.000000f, 2063.000000f, 2060.000000f, 1993.000000f, 1891.000000f, 1615.000000f, 1387.000000f, 1085.000000f, 898.200000f, 772.000000f, 680.400000f, 554.400000f, 471.200000f, 347.800000f, 278.800000f, 202.800000f, 161.200000f, 134.800000f, 116.300000f, 92.140000f, 76.890000f, 55.430000f, 44.050000f, 32.090000f, 25.810000f, 21.930000f, 19.280000f, 15.910000f, 13.840000f  };

  AddData(57, T0_proton, e0_proton);
  AddData(57, T0_proton, e1_proton);
  AddData(57, T0_proton, e2_proton);

  AddData(49, T0_alpha, e0_alpha);
  AddData(49, T0_alpha, e1_alpha);
  AddData(49, T0_alpha, e2_alpha);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LPhysicsFreeVector* G4ICRU90StoppingData::AddData(G4int n, const G4double* e, 
						    const G4float* dedx)
{
  static const G4double fac = CLHEP::MeV*CLHEP::cm2/CLHEP::g;

  G4LPhysicsFreeVector* data = new G4LPhysicsFreeVector(n, e[0], e[n-1]);
  for(G4int i=0; i<n; ++i) { 
    data->PutValues(i, e[i]*CLHEP::MeV, ((G4double)dedx[i])*fac); 
  }
  data->SetSpline(true);
  data->FillSecondDerivatives();
  return data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
