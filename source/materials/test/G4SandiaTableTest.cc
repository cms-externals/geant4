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
// 
//
////////////////////////////////////////////////////////////////////////
//
//  This program illustrates the different ways to define photoabsorption 
//  cross section according G4Sandiatable
//
// History:
//
// 15.09.99 V.Grichine, start from G4MaterialTest.cc
//

#include "G4ios.hh"
#include <iomanip>
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4Material.hh"
#include "G4SandiaTable.hh"

int main() 
{
  // set output format

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  G4String name, symbol;             // a=mass of a mole;
  G4double a, z, density;            // z=mean number of protons;  
  G4int iz, n;                       // iz=nb of protons  in an isotope; 
                                     // n=nb of nucleons in an isotope;

  G4int ncomponents, natoms, nel ;
  G4double abundance, fractionmass ;
  G4double temperature, pressure ;

  G4UnitDefinition::BuildUnitsTable();

//
// define Elements
//

  a = 1.01*g/mole;
  G4Element* elHold  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 1.01*g/mole;
  G4Isotope* ih1 = new G4Isotope("Hydrogen",iz=1,n=1,a);

  a = 2.01*g/mole;
  G4Isotope* ih2 = new G4Isotope("Deuterium",iz=1,n=2,a);

  G4Element* elH = new G4Element(name="Hydrogen",symbol="H",2);
  elH->AddIsotope(ih1,.999);
  elH->AddIsotope(ih2,.001);


  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

  a = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron"    ,symbol="Fe", z=26., a);

  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);
  
  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);


  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);

  a = 69.723*g/mole;
  G4Element* elGa  = new G4Element(name="Ga", symbol="Ga", z=31., a);

  a = 74.9216*g/mole;
  G4Element* elAs  = new G4Element(name="As", symbol="As", z=33., a);


//
// define an Element from isotopes, by relative abundance 
//

  G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* elU  = new G4Element(name="enriched Uranium", 
                                symbol="U", ncomponents=2);
  elU->AddIsotope(U5, abundance= 90.*perCent);
  elU->AddIsotope(U8, abundance= 10.*perCent);


// G4cout << *(G4Isotope::GetIsotopeTable()) << G4endl;

// G4cout << *(G4Element::GetElementTable()) << G4endl;

//
// define simple materials
//

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead "     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//
 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(elC, natoms=9);
  Sci->AddElement(elH, natoms=10);

  density = 2.200*g/cm3;
  G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

  density = 1.290*mg/cm3;
  G4Material* oldAir = new G4Material(name="Air  "  , density, ncomponents=2);
  oldAir->AddElement(elN, fractionmass=0.7);
  oldAir->AddElement(elO, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

  density = 0.200*g/cm3;
  G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=0.625);
  Aerog->AddMaterial(H2O , fractionmass=0.374);
  Aerog->AddElement (elC , fractionmass=0.1*perCent);

//
// examples of gas in non STP conditions
//

  density     = 27.*mg/cm3;
  pressure    = 50.*atmosphere;
  temperature = 325.*kelvin;
  G4Material* CO2 = new G4Material(name="Carbonic gas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
  CO2->AddElement(elC, natoms=1);
  CO2->AddElement(elO, natoms=2);
 
  density     = 0.3*mg/cm3;
  pressure    = 2.*atmosphere;
  temperature = 500.*kelvin;
  G4Material* steam = new G4Material(name="Water steam ", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;                      //from PhysicalConstants.h

  G4Material* beam = new G4Material(name="Beam ", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  beam->AddMaterial(oldAir, fractionmass=1.);

  // maylar

  density = 1.39*g/cm3;
  G4Material* Maylar = new G4Material(name="Maylar", density, nel=3);
  Maylar->AddElement(elO,2);
  Maylar->AddElement(elC,5);
  Maylar->AddElement(elH,4);

  // Kapton Dupont de Nemur (density: 1.396-1.430, get middle )

  density = 1.413*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

  // Germanium as detector material

  density = 5.323*g/cm3;
  a = 72.59*g/mole;
  G4Material* Ge = new G4Material(name="Ge", z=32., a, density);

  // GaAs detectors

  density = 5.32*g/cm3;
  G4Material* GaAs = new G4Material(name="GaAs",density, nel=2);
  GaAs->AddElement(elGa,1);
  GaAs->AddElement(elAs,1);

  // Diamond detectors

  density = 3.5*g/cm3;
  G4Material* Diamond = new G4Material(name="Diamond",density, nel=1);
  Diamond->AddElement(elC,1);


  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
				      kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
				       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);

  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = new G4Material(name="XeCO2CF4", XeCO2CF4_density,
					ncomponents=3,
					kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);
      
  density = 0.935*g/cm3;
  G4Material* TRT_CH2 = new G4Material(name="TRT_CH2",density, nel=2);
  TRT_CH2->AddElement(elC,1);
  TRT_CH2->AddElement(elH,2);

  density = 0.059*g/cm3;
  G4Material* Radiator = new G4Material(name="Radiator",density, nel=2);
  Radiator->AddElement(elC,1);
  Radiator->AddElement(elH,2);

  density = 0.145*g/cm3;
  G4Material* CarbonFiber = new G4Material(name="CarbonFiber",density, nel=1);
  CarbonFiber->AddElement(elC,1);

  // Dry air (average composition)


  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

  // Helium as detector gas, STP

  density = 0.178*mg/cm3 ;
  a = 4.0026*g/mole ;
  G4Material* He  = new G4Material(name="He",z=2., a, density );

  // Neon as detector gas, STP

  density = 0.900*mg/cm3 ;
  a = 20.179*g/mole ;
  G4Material* Ne  = new G4Material(name="Ne",z=10., a, density );

  // Krypton as detector gas, STP

  density = 3.700*mg/cm3 ;
  a = 83.80*g/mole ;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  // Carbone dioxide, CO2 STP

  density = 1.977*mg/cm3 ;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2) ;
  CarbonDioxide->AddElement(elC,1) ;
  CarbonDioxide->AddElement(elO,2) ;

  // Metane, STP

  density = 0.7174*mg/cm3 ;
  G4Material* metane = new G4Material(name="CH4",density,nel=2) ;
  metane->AddElement(elC,1) ;
  metane->AddElement(elH,4) ;

  // Propane, STP

  density = 2.005*mg/cm3 ;
  G4Material* propane = new G4Material(name="C3H8",density,nel=2) ;
  propane->AddElement(elC,3) ;
  propane->AddElement(elH,8) ;

  // iso-Butane (methylpropane), STP

  density = 2.67*mg/cm3 ;
  G4Material* isobutane = new G4Material(name="isoC4H10",density,nel=2) ;
  isobutane->AddElement(elC,4) ;
  isobutane->AddElement(elH,10) ;

// 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm

  density = 4.9196*mg/cm3 ;

  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  , density,
                                                             ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,       fractionmass = 0.971 ) ;
  XeCH4C3H8->AddMaterial( metane,   fractionmass = 0.010 ) ;
  XeCH4C3H8->AddMaterial( propane,  fractionmass = 0.019 ) ;

  // Propane in MWPC, 2 atm, 20 C

  //  density = 3.758*mg/cm3 ;
  density = 3.736*mg/cm3 ;
  G4Material* propaneDet = new G4Material(name="detC3H8",density,nel=2) ;
  propaneDet->AddElement(elC,3) ;
  propaneDet->AddElement(elH,8) ;

  // 80% Ar + 20% CO2, STP

  density = 1.8223*mg/cm3 ;
  G4Material* Ar20CO2 = new G4Material(name="Ar20CO2"  , density,
                                                             ncomponents=2);
  Ar20CO2->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 93% Ar + 7% CH4, STP

  density = 1.709*mg/cm3 ;
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density,
                                                             ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 ) ;

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3 ;
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , density,
                                                             ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3 ;
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density,
                                                             ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 ) ;

  // 80% He + 20% CO2, STP

  density = 0.5378*mg/cm3 ;
  G4Material* He20CO2 = new G4Material(name="He20CO2"  , density,
                                                             ncomponents=2);
  He20CO2->AddMaterial( He,              fractionmass = 0.265 ) ;
  He20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.735 ) ;



//
// Print the table of materials
//

// G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//
////////////////////////////////////////////////////////////////////////
//
//
// Checking Sandia table coefficients
//
  G4int numberOfMat, iMat, matIndex, nbOfElements, sanIndex, row, iSan;
  G4double unit;
  G4String materialName = "Air";
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  numberOfMat = theMaterialTable->size() ;

  for(iMat=0;iMat<numberOfMat;iMat++)
  {
    if(materialName == (*theMaterialTable)[iMat]->GetName() )
    {
      matIndex = (*theMaterialTable)[iMat]->GetIndex() ;
      break ;
    }
  }
  
//
////////////////////////////////////////////////////////////////////////
//
// Sandia cof according old PAI stuff
//
  for(iMat=0;iMat<numberOfMat;iMat++)
  {
     G4String matName = (*theMaterialTable)[iMat]->GetName();
     matIndex = (*theMaterialTable)[iMat]->GetIndex();
     nbOfElements = (*theMaterialTable)[iMat]->GetNumberOfElements();
     density = (*theMaterialTable)[iMat]->GetDensity();
     
     G4cout<<matIndex<<"\t"<<matName<<G4endl<<G4endl;
     
     G4cout<<"Sandia cof according old PAI stuff"<<G4endl<<G4endl;

     G4int* thisMaterialZ = new G4int[nbOfElements];
     for(iSan=0;iSan<nbOfElements;iSan++)
     {
        thisMaterialZ[iSan] = (G4int)(*theMaterialTable)[iMat]->
                                      GetElement(iSan)->GetZ();
     }     
     G4SandiaTable sandia(matIndex) ;

     //  sanIndex = sandia.SandiaIntervals(thisMaterialZ,nbOfElements);    
     //  sanIndex = sandia.SandiaMixing( thisMaterialZ ,
     //                             (*theMaterialTable)[iMat]->GetFractionVector() ,
     //			     nbOfElements,sanIndex) ;
     sanIndex = sandia.GetMaxInterval() ;
     G4cout<<"fMaxInterval = "<<sanIndex<<G4endl<<G4endl;

     for(row = 0; row < sanIndex - 1 ; row++)
     {
       G4cout<<row+1<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,0)/keV;
       
       unit = cm2/g;
       for(iSan = 1; iSan < 5; iSan++)
       {
         unit *= keV;         
         G4cout<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,iSan)/unit;
       }
       G4cout<<G4endl ;
     }
     G4cout<<G4endl ;
     
//
////////////////////////////////////////////////////////////////////////
//
// Sandia cof according ComputeMatSandiaMatrix()
//
   G4SandiaTable* sanMatrix = G4Material::GetMaterial(matName)->
                              GetSandiaTable();
   sanIndex = sanMatrix->GetMatNbOfIntervals();
      
   G4cout<<"Sandia cof according ComputeMatSandiaMatrix()"<<G4endl<<G4endl;

     for (row=0; row<sanIndex; row++) {
       G4cout<<row+1<<"\t"<<sanMatrix->GetSandiaCofForMaterial(row,0)/keV;
	       
       unit = cm2/g;
       for (iSan=1; iSan<5; iSan++) {
         unit *= keV; 
         G4cout<<"\t"<<(sanMatrix->GetSandiaCofForMaterial(row,iSan))
	                                                    /(density*unit);
       }
       G4cout<<G4endl;
     }      
     G4cout<<G4endl;     
  }
  return EXIT_SUCCESS;
}

//
//
///////////////////////  end of G4SandiaTableTest.cc  /////////////////////////
