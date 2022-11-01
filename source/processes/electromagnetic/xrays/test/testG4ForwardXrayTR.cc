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



#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "g4templates.hh"
#include "globals.hh"
#include "G4Timer.hh"

#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsVector.hh"
     
#include "G4eEnergyLoss.hh"
#include "G4eIonisation.hh"
#include "G4hEnergyLoss.hh"
#include "G4PAIenergyLoss.hh"
#include "G4hIonisation.hh"
#include "G4PAIonisation.hh"
#include "G4ForwardXrayTR.hh"
#include "G4TransitionRadiation.hh"

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4GPILSelection.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
                  
//    It tests the G4PAIenergyLoss,G4PAIonisation and G4ForwardXrayTR processes 
//    created by L.Urban on 06/06/97 and modified by V.Grichine on 17.12.97
//
//    Modifications:
//    23-09-97, geometry adapted for the touchable
//    14.11.97  PAIonisation class is responsible for <dE/dx> calculation
//    16.12.97  G4ForwardXrayTR class for generation of X-ray TR
//

G4VPhysicalVolume* BuildVolume(G4Material* matworld)
//  it builds a simple box filled with material matword .......
{
  G4Box *myWorldBox= new G4Box ("WBox",10000.*cm,10000.*cm,10000.*cm);

  G4LogicalVolume *myWorldLog = new G4LogicalVolume(myWorldBox,matworld,
                                                    "WLog",0,0,0) ;

  G4PVPlacement *myWorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                                 "WPhys",
                                                 myWorldLog,
                                                 0,false,0) ;

  
  return myWorldPhys ;

}
                                              

int main()
{
  //-------- set output format-------
   G4cout.setf( std::ios::scientific, std::ios::floatfield );
  //---write results to the file  hloss -----
   std::ofstream outFile("XrayTR.cc", std::ios::out ) ;
   outFile.setf( std::ios::scientific, std::ios::floatfield );

  //--------- Material definition ---------

  G4double a, z, ez, density ,temperature,pressure;
  G4State state ;
  G4String name, symbol;
  G4int nel;

  G4Element*   elH = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);

  a = 6.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", ez=6., a);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", ez=8., a);

  G4Material* C2H2 = new G4Material ("Polypropelene" , 0.91*g/cm3, 2);
  C2H2->AddElement(elC,2);
  C2H2->AddElement(elH,2);


  density = 1.29e-03*g/cm3;
  state = kStateGas ;
  temperature = 273.*kelvin ;
  pressure = 1.*atmosphere ;
  G4Material* Air = new G4Material(name="Air", density, nel=2 ,
                                   state ,temperature , pressure ) ;
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);


  /*  ***************************************************************

  density = 5.85e-3*g/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon", ez=54., a, density);

  density = 1.782e-3*g/cm3 ;
  a = 39.948*g/mole ;
  G4Material* Ar  = new G4Material(name="Argon", ez=18., a, density);

  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);


  G4Material* H2O = new G4Material ("Water" , 1.*g/cm3, 2);
  H2O->AddElement(elH,2);
  H2O->AddElement(elO,1);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);

  a = 196.97*g/mole;
  density = 19.32*g/cm3;
  G4Material* Au = new G4Material(name="Gold", z=79., a, density);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  a = 0. ;          
  density = 0. ;         
  G4Material* Vac= new G4Material(name="Vacuum",z=0., a, density,kVacuum);

  ******************************************************************** */

  const G4MaterialTable* theMaterialTable ;
  G4Material* apttoMaterial ;
  G4String MaterialName ; 
  G4Timer theTimer ;




//--------- Particle definition ---------

  G4ParticleDefinition* theGamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* theElectron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* thePositron = G4Positron::PositronDefinition();
  G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* theAntiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* thePionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* thePionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* theKaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* theKaonMinus = G4KaonMinus::KaonMinusDefinition();

  G4double* GammaKineticEnergyCuts ;
  G4double* ElectronKineticEnergyCuts ;
  G4double* PositronKineticEnergyCuts ;
  G4double* ParticleKineticEnergyCuts ;

  theMaterialTable = G4Material::GetMaterialTable() ;

  G4double cutinrange,CutInRangeele,CutInRangepos ;

  G4eIonisation theElectronIonisation , thePositronIonisation ;
  G4ProcessManager* theElectronProcessManager = theElectron->GetProcessManager();
  theElectronProcessManager->AddProcess(&theElectronIonisation,-1,0,0) ;
  G4ProcessManager* thePositronProcessManager = thePositron->GetProcessManager();
  thePositronProcessManager->AddProcess(&thePositronIonisation,-1,0,0) ;

  G4ParticleWithCuts* theParticle ;
    
  G4double energy, momentum, mass;
  G4ProcessVector* palongget ;
  G4ProcessVector* palongdo ;
  G4ProcessVector* ppostget ;
  G4ProcessVector* ppostdo ;

  G4String confirm ;

  G4cout << " Do you want the proton as particle (yes/no)? " << std::flush;
  //G4cin >> confirm ;
  confirm = "yes"  ;
  if(confirm == "yes")
  {
    mass=theProton->GetPDGMass();
    theParticle = theProton;
  }
  else
  {
     G4cout << " Do you want the antiproton as particle (yes/no)? " << std::flush;
     G4cin >> confirm ;
     if(confirm == "yes")
     {
        mass=theAntiProton->GetPDGMass();
        theParticle = theAntiProton;
     }
     else
     {
      G4cout << " Do you want the pi+ as particle (yes/no)? " << std::flush;
      G4cin >> confirm ;
      if(confirm == "yes")
      {
      mass=thePionPlus->GetPDGMass();
      theParticle = thePionPlus;
      }
      else
      {
        G4cout << " Do you want the pi- as particle (yes/no)? " << std::flush;
        G4cin >> confirm ;
        if(confirm == "yes")
        {
        mass=thePionMinus->GetPDGMass();
        theParticle = thePionMinus;
        } 
        else
        {
          G4cout << " Do you want the K+ as particle (yes/no)? " << std::flush;
          G4cin >> confirm ;
          if(confirm == "yes")
          {
          mass=theKaonPlus->GetPDGMass();
          theParticle = theKaonPlus;
          }
          else
          {
            G4cout << " Do you want the K- as particle (yes/no)? " << std::flush;
            G4cin >> confirm ;
            if(confirm == "yes")
            {
            mass=theKaonMinus->GetPDGMass();
            theParticle = theKaonMinus;
            }
            else
            {
             G4cout << " There is no other particle in the test." << G4endl;
             return EXIT_FAILURE;
            }
          }
        }
       }
      }
     }
     
  
    energy = 1000.0*GeV + mass ;  // was 1.0*GeV now 1.0*TeV
    momentum=std::sqrt(energy*energy-mass*mass) ;
  
    G4ParticleMomentum theMomentum(momentum,0.,0.);
  
    G4double pModule = theMomentum.mag();

    G4DynamicParticle aParticle(theParticle,energy,theMomentum);

    aParticle.SetKineticEnergy(energy-mass);


    // G4hIonisation theParticleIonisation  ;

  G4PAIonisation theParticleIonisation  ;

  G4ProcessManager* theParticleProcessManager = theParticle->GetProcessManager();
  theParticleProcessManager->AddProcess(&theParticleIonisation,-1,0,0) ;

  G4ForceCondition cond ;
  G4ForceCondition* condition = &cond ;

  G4double currentSafety ;
  G4double& refsafety=currentSafety;


  G4cout << "cut for GAMMA in mm =" ;
  // G4cin >> cutinrange ;
  cutinrange = 0.1 ;
  theGamma->SetCuts(cutinrange) ;
    G4cout << "gamma,cut in range(mm)=" << theGamma->GetCuts() << G4endl ;
    outFile << "  ---------------------------------------" << G4endl ;
    outFile << "  gamma,cut in range(mm)=" << theGamma->GetCuts() << G4endl ;

  GammaKineticEnergyCuts = theGamma->GetCutsInEnergy() ;
  for (G4int icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut] << G4endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut] << G4endl ;
  }


  G4cout << "cut for ELECTRON in mm =" ;
  // G4cin >> cutinrange ; 
  cutinrange = 1.0 ;
  CutInRangeele = cutinrange ;
  theElectron->SetCuts(cutinrange) ;
    G4cout << "electron,cut in range(mm)=" << theElectron->GetCuts() << G4endl ;
    outFile << "  ---------------------------------------" << G4endl ;
    outFile << "  electron,cut in range(mm)=" << theElectron->GetCuts() << G4endl ;

  ElectronKineticEnergyCuts = theElectron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut] << G4endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut] << G4endl ;
  }

  G4cout << "cut for POSITRON in mm =" ;
  //  G4cin >> cutinrange ; 
  cutinrange = 1.0 ;
  CutInRangepos = cutinrange ;
  thePositron->SetCuts(cutinrange) ;
    G4cout << "positron,cut in range(mm)=" << thePositron->GetCuts() << G4endl ;
    outFile << "  ---------------------------------------" << G4endl ;
    outFile << "  positron,cut in range(mm)=" << thePositron->GetCuts() << G4endl ;

  PositronKineticEnergyCuts = thePositron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut] << G4endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut] << G4endl ;
  }

  G4cout << "cut for hadrons in mm =" ;
  //  G4cin >> cutinrange ; 
  cutinrange = 1.0 ;
  theParticle->SetCuts(cutinrange) ;
  G4cout << "after particle setcuts " << G4endl;


    G4cout << "cut in range(mm)=" << theParticle->GetLengthCuts() << G4endl ;
    outFile << "  ---------------------------------------" << G4endl ;
    outFile << "  cut in range(mm)=" << theParticle->GetLengthCuts() << G4endl ;

  ParticleKineticEnergyCuts = theParticle->GetEnergyCuts() ;

  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         ParticleKineticEnergyCuts[icut] << G4endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         ParticleKineticEnergyCuts[icut] << G4endl ;
  }

    G4cout << "  ------         ----- " << G4endl ;
    outFile << "  " << G4endl;

    // Creation of X-ray transition radiation object
    G4cout<<"Just before TR object creation"<<G4endl ;

  G4ForwardXrayTR theXrayTR  ;

    G4cout<<"Just after TR object creation"<<G4endl ;
  theParticleProcessManager->AddProcess(&theXrayTR,-1,0,0) ;




    /* *********************************************************************
    // Output to file of fPAItransferBank data

    G4bool isOutRange ;
    G4PhysicsLogVector* 
    aLogVector = new G4PhysicsLogVector(G4PAIonisation::GetMinKineticEnergy(),
                                        G4PAIonisation::GetMaxKineticEnergy(),
                                        G4PAIonisation::GetBinNumber()     ) ;

    for(G4int materialIndex=0 ;
              materialIndex < theMaterialTable->length() ; materialIndex++)
    {
      apttoMaterial = (*theMaterialTable)[materialIndex] ;
      MaterialName = apttoMaterial->GetName() ;
      outFile<<"PAI transfer data for the material = "
             <<MaterialName<<G4endl<<G4endl ;

      outFile<<"Particle Tkin"<<"\t\t"
             <<"Transfer energy, keV"<<"\t"
             <<"Primary ionisation, 1/mm"<<G4endl<<G4endl ;

  //  G4cout<<"Read of TotBin ? TotBin = "<<G4PAIonisation::GetBinNumber()<<G4endl;

      for(G4int iTkin=48;iTkin<G4PAIonisation::GetBinNumber();iTkin++)
      {

        outFile<<aLogVector->GetLowEdgeEnergy(iTkin)<<G4endl<<G4endl ;

        G4PhysicsVector* trVec = (*G4PAIenergyLoss::GetPAItransferBank())
                       (materialIndex*G4PAIonisation::GetBinNumber() + iTkin) ;
	G4cout<<materialIndex*G4PAIonisation::GetBinNumber() + iTkin
	       <<"\t\t"<<trVec->GetVectorLength()<<G4endl ;
        for(G4int iTr=0;iTr<trVec->GetVectorLength();iTr++)
        {
          outFile<<"\t\t"<<trVec->GetLowEdgeEnergy(iTr)*1000.
                 <<"\t\t"<<(*trVec)(iTr)<<G4endl ;
        }
               
      }
      outFile<<G4endl ;
    } 


    G4cout<<"Start dE/dx distribution"<<G4endl ;
    G4int iLoss, iStat, iStatMax = 10000;
    G4double energyLoss[200], delta ;
    G4int spectrum[200] ;
    const G4DynamicParticle* particlePtr = &aParticle ;

    for(iLoss=0;iLoss<200;iLoss++)
    {
        energyLoss[iLoss] = 0.5*iLoss*keV ;
      spectrum[iLoss] = 0 ;
    }
    for(iStat=0;iStat<iStatMax;iStat++)
    {
      delta = G4PAIenergyLoss::GetLossWithFluct(10.0*mm,particlePtr,Xe) ;
      for(iLoss=0;iLoss<200;iLoss++)
      {
        if(delta <= energyLoss[iLoss]) break ;
      }
      spectrum[iLoss-1]++ ;
    }
    G4double meanLoss = 0.0 ;
    outFile<<"Energy loss, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl ;
    for(iLoss=0;iLoss<200;iLoss++)
    {
      outFile<<energyLoss[iLoss]*1000.0<<"\t\t"<<spectrum[iLoss]<<G4endl ;
      meanLoss +=energyLoss[iLoss]*spectrum[iLoss] ;
    }
    G4cout<<"Mean loss over spectrum = "<<meanLoss*1000./iStatMax<<" keV"<<G4endl ;


    ****************************************************************  */

    // Output of TR information

    G4int iMat, jMat, iPlace, iTkin, iTR ;
    G4PhysicsLogVector* vectorTkin = new 
    G4PhysicsLogVector(G4ForwardXrayTR::GetMinProtonTkin(),
                       G4ForwardXrayTR::GetMaxProtonTkin(),
                       G4ForwardXrayTR::GetTotBin()          ) ;

     for(iMat=0;iMat<theMaterialTable->length();iMat++)
     {
       for(jMat=0;jMat<theMaterialTable->length();jMat++)
       {
         if(jMat == iMat) continue ;
        
         if(jMat < iMat)
	 {
           outFile<<(*theMaterialTable)[iMat]->GetName()<<" -> "<<"\t\t"
		  <<(*theMaterialTable)[jMat]->GetName()<<G4endl<<G4endl ;
           outFile<<"Lorentz Factor"<<"\t\t"
               <<"energy TR number"<<"\t\t"  
		  <<"theta TR number"<<"\t\t"<<G4endl<<G4endl ;  
          for(iTkin=0;iTkin<G4ForwardXrayTR::GetTotBin();iTkin++)
	   {
           iPlace =(iMat*(theMaterialTable->length()-1)+jMat)*
                     G4ForwardXrayTR::GetTotBin()+iTkin           ;

          outFile<<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(0)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(0)<<G4endl ;
          G4cout<<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(0)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(0)<<G4endl ;
	   }
	 }
         else
	 {
           outFile<<(*theMaterialTable)[iMat]->GetName()<<" -> "<<"\t\t"
		  <<(*theMaterialTable)[jMat]->GetName()<<G4endl<<G4endl ;
           outFile<<"Lorentz Factor"<<"\t\t"
               <<"energy TR number"<<"\t\t"  
		  <<"theta TR number"<<"\t\t"<<G4endl<<G4endl ;  
           for(iTkin=0;iTkin<G4ForwardXrayTR::GetTotBin();iTkin++)
	   {
             iPlace =(iMat*(theMaterialTable->length()-1)+jMat-1)*
                     G4ForwardXrayTR::GetTotBin()+iTkin          ;

          outFile<<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(0)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(0)<<G4endl ;
          G4cout<<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(0)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(0)<<G4endl ;
	   }
	 }
         outFile<<G4endl ; 
         G4cout<<G4endl ; 
       } 
     } 

     for(iMat=0;iMat<theMaterialTable->length();iMat++)
     {
       for(jMat=0;jMat<theMaterialTable->length();jMat++)
       {
         if(jMat == iMat) continue ;
        
         if(jMat < iMat)
	 {
           outFile<<(*theMaterialTable)[iMat]->GetName()<<" -> "<<"\t\t"
		  <<(*theMaterialTable)[jMat]->GetName()<<G4endl<<G4endl ;
          for(iTkin=0;iTkin<G4ForwardXrayTR::GetTotBin();iTkin++)
	   {
             if(iTkin == 9 || iTkin == 19 || iTkin == 29
                           || iTkin == 39 || iTkin == 49 )
	     {
	       outFile<<"iTkin = "<<iTkin<<"\t\t"
                      <<"Lorentz factor = "
                      <<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)
                      <<G4endl<<G4endl ;
           outFile<<"TR Energy"<<"\t\t"
                  <<" > energy TR number"<<"\t\t"
                  <<"TR ThetaSq"<<"\t\t"
		  <<" > theta TR number"<<"\t\t"<<G4endl<<G4endl ;  
               iPlace =(iMat*(theMaterialTable->length()-1)+jMat)*
                     G4ForwardXrayTR::GetTotBin()+iTkin           ;

               for(iTR=0;iTR<G4ForwardXrayTR::GetBinTR();iTR++)
	       {
          outFile<<(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace)->GetLowEdgeEnergy(iTR)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(iTR)<<"\t\t"
                 <<(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace)->GetLowEdgeEnergy(iTR)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(iTR)<<G4endl ;
	       }
	     }
	   }
	 }
         else
	 {
           outFile<<(*theMaterialTable)[iMat]->GetName()<<" -> "<<"\t\t"
		  <<(*theMaterialTable)[jMat]->GetName()<<G4endl<<G4endl ;
           for(iTkin=0;iTkin<G4ForwardXrayTR::GetTotBin();iTkin++)
	   {
             if(iTkin == 9 || iTkin == 19 || iTkin == 29
                           || iTkin == 39 || iTkin == 49 )
	     {
	       outFile<<"iTkin = "<<iTkin<<"\t\t"
                      <<"Lorentz factor = "
                      <<1+(vectorTkin->GetLowEdgeEnergy(iTkin)/proton_mass_c2)
                      <<G4endl<<G4endl ;
           outFile<<"TR Energy"<<"\t\t"
                  <<" > energy TR number"<<"\t\t"
                  <<"TR ThetaSq"<<"\t\t"
		  <<" > theta TR number"<<"\t\t"<<G4endl<<G4endl ;  
               iPlace =(iMat*(theMaterialTable->length()-1)+jMat-1)*
                     G4ForwardXrayTR::GetTotBin()+iTkin               ;

               for(iTR=0;iTR<G4ForwardXrayTR::GetBinTR();iTR++)
	       {
          outFile<<(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace)->GetLowEdgeEnergy(iTR)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetEnergyDistrTable())
                   (iPlace))(iTR)<<"\t\t"
                 <<(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace)->GetLowEdgeEnergy(iTR)<<"\t\t"
                 <<(*(*G4ForwardXrayTR::GetAngleDistrTable())
                   (iPlace))(iTR)<<G4endl ;
	       }
	     }
	   }
	 }
         outFile<<G4endl ; 
         G4cout<<G4endl ; 
       } 
     } 

    G4cout<<"Start TR energy distribution"<<G4endl ;
    G4int iLossTR, iStatTR, iStatMaxTR = 1000;
    G4double energyTR[200], deltaTR ;
    G4int spectrumTR[200] ;
    // const G4DynamicParticle* particlePtr = &aParticle ;

    for(iLossTR=0;iLossTR<200;iLossTR++)
    {
        energyTR[iLossTR] = 0.5*iLossTR*keV ;
      spectrumTR[iLossTR] = 0 ;
    }
    for(iStatTR=0;iStatTR<iStatMaxTR;iStatTR++)
    {
      deltaTR = theXrayTR.GetEnergyTR(0,1,39) ;
      for(iLossTR=0;iLossTR<200-1;iLossTR++)
      {
        if(deltaTR <= energyTR[iLossTR]) break ;
      }
      spectrumTR[iLossTR]++ ;
    }
    G4double meanLossTR = 0.0 ;
    outFile<<"TR energy, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl ;
    for(iLossTR=0;iLossTR<200;iLossTR++)
    {
      outFile<<energyTR[iLossTR]*1000.0<<"\t\t"<<spectrumTR[iLossTR]<<G4endl ;
      meanLossTR +=energyTR[iLossTR]*spectrumTR[iLossTR] ;
    }
    G4cout<<"Mean TR energy over spectrum = "
        <<meanLossTR*1000./iStatMaxTR<<" keV"<<G4endl ;


    G4int message ;
    G4cout<<"Continue ?(enter int>0)  Break ? (enter int<=0)"<<G4endl ;
    G4cin>>message ;
    if(message<= 0)      // break program at this point
    {
      return 1 ;    
    }
 


    outFile << " ionisation test  **************************************" << G4endl ;
    outFile << "  " << G4endl;
    outFile << "   particle = " << 
             aParticle.GetDefinition()->GetParticleName() << G4endl ;
    outFile << G4endl;
    
    palongget = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeGPIL);
    ppostget = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeGPIL);
    palongdo = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeDoIt);
    ppostdo = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeDoIt);

//---------------------------------- Physics --------------------------------

  G4int itry=1, Ntry=1, Nstart, ir;
  G4double r ;

//**************************************************************************
  const G4int Nbin=97 ;
  G4double TkinMeV[Nbin]  =
              {0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,
               0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,
               0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,
               1.,1.5,2.,3.,4.,5.,6.,8.,
               10.,15.,20.,30.,40.,50.,60.,80.,
               100.,150.,200.,300.,400.,500.,600.,800.,
               1.0e3,1.5e3,2.0e3,3.0e3,4.0e3,5.0e3,6.0e3,8.0e3,
               1.0e4,1.5e4,2.0e4,3.0e4,4.0e4,5.0e4,6.0e4,8.0e4,
               1.0e5,1.5e5,2.0e5,3.0e5,4.0e5,5.0e5,6.0e5,8.0e5,
               1.0e6,1.5e6,2.0e6,3.0e6,4.0e6,5.0e6,6.0e6,8.0e6,
               1.0e7,1.5e7,2.0e7,3.0e7,4.0e7,5.0e7,6.0e7,8.0e7,
               1.0e8,1.5e8,2.0e8,3.0e8,4.0e8,5.0e8,6.0e8,8.0e8,
               1.0e9} ; 
         
    G4int J=-1 ;

    G4double lambda,trueStep,geomStep,stepLimit,
             previousStepSize,currentMinimumStep ;
    G4ParticleChange* aParticleChange ;
    
    G4double T,dEdx,range ;

    NEXTMATERIAL: ;
    J = J+1 ;
    if ( J >= theMaterialTable->length() )
      { G4cout << "that was the last material in the table --> STOP" << G4endl;
        return EXIT_FAILURE ; }  

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName = apttoMaterial->GetName() ; 
    G4cout << "material=" << MaterialName << G4endl ;
    G4cout << "Do you want the Energyloss test 1. for this material?" << G4endl ;
    G4cout << "type a positive number if the answer is YES" << G4endl ;
    G4cout << "type a negative number if the answer is NO " << G4endl ;
    G4int icont ;
    G4cin >> icont ;
    if ( icont < 0 )
        goto NEXTMATERIAL ;

//---------- Volume definition ---------------------

    G4VPhysicalVolume* myVolume ;

    myVolume = BuildVolume(apttoMaterial) ;

//--------- track and Step definition (for this test ONLY!)------------
    G4ThreeVector aPosition(0.,0.,0.);
    const G4ThreeVector aDirection(0.,0.,1.) ;
    G4double aTime = 0. ;

    G4Track* tracke = new G4Track(&aParticle,aTime,aPosition) ;
    G4Track& trackele = (*tracke) ;
    //(*tracke).SetVolume(myVolume) ;
    G4GRSVolume* touche = new G4GRSVolume(myVolume, NULL, aPosition);   
    (*tracke).SetTouchable(touche);        
    (*tracke).SetMomentumDirection(aDirection) ;


    G4Step* Step = new G4Step() ;
    G4Step& Step = (*Step) ;

    G4StepPoint* aPoint = new G4StepPoint();
    (*aPoint).SetPosition(aPosition) ;
    G4double safety = 10000.*cm ;
    (*aPoint).SetSafety(safety) ;

    (*Step).SetPostStepPoint(aPoint) ;
   
//**************************************************************************

    G4cout <<  G4endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 1." << G4endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    G4cout << G4endl ;
    G4cout << "kin.en.(MeV)    dE/dx(MeV/mm)    range(mm)    Step(mm)" << G4endl ;
    G4cout << G4endl ;
 
    outFile <<  G4endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 1." << G4endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;
    outFile << "kin.en.(MeV)    dE/dx(MeV/mm)    range(mm)    Step(mm)" << G4endl ;
    outFile << G4endl ;
 

    G4GPILSelection* selection ;

    for ( G4int i=0 ; i<Nbin ; i++)
    {
      trueStep = cutinrange ; 
      previousStepSize = cutinrange ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = (*palongget)(0)->AlongStepGetPhysicalInteractionLength( 
                                                         trackele,            
                                                         previousStepSize,
                                                         currentMinimumStep,
                                                         refsafety,
                                                         selection) ;
 
      dEdx = theParticleIonisation.GetdEdx() ;
      range = theParticleIonisation.GetRangeNow() ;
      T = TkinMeV[i] ;

       G4cout <<" " <<  T << "  " << dEdx/mm << "  " ;
       G4cout << range/mm << "  " << stepLimit/mm << G4endl ; 

       outFile <<" " <<  T << "  " << dEdx/mm << "  " ;
       outFile << range/mm << "  " << stepLimit/mm << G4endl ; 

    }

    G4cout <<  G4endl;
    outFile << G4endl;

    ENERGYLOSS2: ;

    G4cout << "material=" << MaterialName << G4endl ;
    G4cout << "Do you want the Energyloss test 2. for this material?" << G4endl ;
    G4cout << "type a positive number if the answer is YES" << G4endl ;
    G4cout << "type a negative number if the answer is NO " << G4endl ;
    G4cin >> icont ;
    if ( icont < 0 )
        goto ENERGYLOSS3 ;

    G4double TMeV,stepmm,stepmx,meanloss,lossnow ;
  

    G4cout << "give an energy value in MeV " ;
    G4cin >> TMeV ;

      trueStep = cutinrange ; 
      previousStepSize = cutinrange ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackele,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              selection);

    G4cout << " give a steplength in mm , the max. meaningful Step is " << stepmx << " mm" <<G4endl;
    G4cout << "Step:" ;
    G4cin >> stepmm ;
 
   (*Step).SetTrack(tracke) ;
   (*Step).SetStepLength(stepmm);

 
      aParticleChange = (G4ParticleChange*)((*palongdo)(0)->
                                  AlongStepDoIt(trackele,Step));
      meanloss = theParticleIonisation.GetMeanLoss() ;
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();

    G4cout <<  G4endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 2." << G4endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    G4cout << G4endl ;
    G4cout << "kin.en.(MeV)    Step(mm)   meanloss(MeV)  act.loss(MeV)" << G4endl ;
    G4cout << TMeV << "   " << stepmm << "  " << meanloss << "  " << lossnow << G4endl ;
    G4cout << " status change:" << (*aParticleChange).GetStatusChange() << G4endl ;
    G4cout << G4endl ;
 
    outFile <<  G4endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 2." << G4endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;
    outFile << "kin.en.(MeV)    Step(mm)   meanloss(MeV)  act.loss(MeV)" << G4endl ;
    outFile << TMeV << "   " << stepmm << "  " << meanloss << "  " << lossnow << G4endl ;
    outFile << " status change:" << (*aParticleChange).GetStatusChange() << G4endl ;
    outFile << G4endl ;
 
    goto ENERGYLOSS2 ;

    ENERGYLOSS3: ;

    G4cout << "material=" << MaterialName << G4endl ;
    G4cout << "Do you want the Energyloss test 3. for this material?" << G4endl ;
    G4cout << "type a positive number if the answer is YES" << G4endl ;
    G4cout << "type a negative number if the answer is NO " << G4endl ;
    G4cin >> icont ;
    if ( icont < 0 )
        goto DELTARAY1 ;

    G4cout << "give an energy value in MeV " ;
    G4cin >> TMeV ;

      trueStep = cutinrange ; 
      previousStepSize = cutinrange ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackele,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              selection);

    G4cout << " give a steplength in mm , the max. meaningful Step is " << stepmx << " mm" <<G4endl;
    G4cout << "Step:" ;
    G4cin >> stepmm ;
 
   (*Step).SetTrack(tracke) ;
   (*Step).SetStepLength(stepmm);


    G4cout << " give number of events you want " ;
    G4int nbev,ibev ;
    G4cin >> nbev ;

    meanloss=0.;
    theTimer.Start();

    for ( ibev=0; ibev<nbev; ibev++)
    { 
      aParticleChange =(G4ParticleChange*)((*palongdo)(0)->
                           AlongStepDoIt(trackele,Step));
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();

    meanloss += lossnow ;
   }

    theTimer.Stop();
    meanloss /= nbev ;
    G4cout <<  G4endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 3." << G4endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    G4cout << G4endl ;
    G4cout << "kin.en.(MeV)    Step(mm)   meanloss(MeV) time/event(sec) " << G4endl ;
    G4cout << TMeV << "   " << stepmm << "  " << meanloss << "  " << 
    theTimer.GetUserElapsed()/nbev << G4endl ;
    G4cout << G4endl ;
 
    outFile <<  G4endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 3." << G4endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;
    outFile << "kin.en.(MeV)    Step(mm)   meanloss(MeV) time/event(sec) " << G4endl ;
    outFile << TMeV << "   " << stepmm << "  " << meanloss << "  " << 
    theTimer.GetUserElapsed()/nbev << G4endl ;
    outFile << G4endl ;
 
    goto ENERGYLOSS3 ;

    DELTARAY1: ;
    G4cout << "material=" << MaterialName << G4endl ;
    G4cout << "Do you want the delta ray test 1. for this material?" << G4endl ;
    G4cout << "type a positive number if the answer is YES" << G4endl ;
    G4cout << "type a negative number if the answer is NO " << G4endl ;
    G4cin >> icont ;
    if ( icont < 0 )
        goto DELTARAY2 ;


    G4cout <<  G4endl;
    G4cout <<"  " << MaterialName  << "  delta ray test 1." << G4endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    G4cout << G4endl ;
    G4cout << "kin.en.(MeV)      mean free path(mm)" << G4endl ;
    G4cout << G4endl ;
 
    outFile <<  G4endl;
    outFile <<"  " << MaterialName  << "  delta ray test 1." << G4endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;
    outFile << "kin.en.(MeV)      mean free path(mm)" << G4endl ;
    outFile << G4endl ;

    for ( i=0 ; i<Nbin ; i++)
    {

      previousStepSize = cutinrange ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = theParticleIonisation.GetMeanFreePath(                                           
                                                         trackele,            
                                                         previousStepSize,
                                                         condition) ;                                           

      T = TkinMeV[i] ;

      G4cout <<" " <<  T << "      " <<  stepLimit/mm << G4endl ;

      outFile <<" " <<  T << "      " << stepLimit/mm << G4endl ;
 
    }

    G4cout <<  G4endl;
    outFile << G4endl;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    DELTARAY2: ;

    G4cout << "material=" << MaterialName << G4endl ;
    G4cout << "Do you want the deltaray test 2. for this material?" << G4endl ;
    G4cout << "type a positive number if the answer is YES" << G4endl ;
    G4cout << "type a negative number if the answer is NO " << G4endl ;
    G4cin >> icont ;

    G4double newenergy,dx,dy,dz,Tdelta,ddx,ddy,ddz ;
    G4int nd ;
    const G4ThreeVector* momdir ;
    G4ParticleMomentum ddir ;

    if ( icont < 0 )
        goto BREMS1 ;  

    G4cout << "give an energy value in MeV " ;
    G4cin >> TMeV ;
                           
     stepmm = 1. ;

   (*Step).SetTrack(tracke) ;
   (*Step).SetStepLength(stepmm);


      (*tracke).SetKineticEnergy(TMeV) ;
      aParticleChange = (G4ParticleChange*)((*ppostdo)(0)->
                          PostStepDoIt(trackele,Step));

     newenergy=(*aParticleChange).GetEnergyChange() ;
     momdir=(*aParticleChange).GetMomentumChange();
     dx = (*momdir).x();
     dy = (*momdir).y();
     dz = (*momdir).z();
     nd=aParticleChange->GetNumberOfSecondaries();
 
    if(nd>0)
    {
     Tdelta=aParticleChange->GetSecondary(0)->GetKineticEnergy();
     ddir=aParticleChange->GetSecondary(0)->
                              GetMomentumDirection();
     ddx = (ddir).x();
     ddy = (ddir).y();
     ddz = (ddir).z();
    }

    G4cout <<  G4endl;
    G4cout <<"  " << MaterialName  << "  delta ray test 2." << G4endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    G4cout << G4endl ;
    G4cout << "T=" << TMeV << "   newT=" << newenergy << "  (MeV)" << G4endl ;
    G4cout << " status change:" << (*aParticleChange).GetStatusChange() << G4endl ;
    if(nd>0)
    G4cout << "Tdelta=" << Tdelta << G4endl ;
    G4cout << "new direction:" << dx << "  " << dy << "  " << dz << G4endl;
    if(nd>0)
    G4cout << "delta direction:" << ddx << "  " << ddy << "  " << ddz << G4endl ;
    G4cout << G4endl ;
 
    outFile <<  G4endl;
    outFile <<"  " << MaterialName  << "  delta ray test 2." << G4endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;
    outFile << "T=" << TMeV << "   newT=" << newenergy << "   (MeV)" << G4endl;
    outFile << " status change:" << (*aParticleChange).GetStatusChange() << G4endl ;
    if(nd>0)
    outFile << "Tdelta=" << Tdelta << G4endl ;
    outFile << "new direction:" << dx << "  " << dy << "  " << dz << G4endl;
    if(nd>0)
    outFile << "delta direction:" << ddx << "  " << ddy << "  " << ddz << G4endl ;
    outFile << G4endl ;

    (*aParticleChange).Clear();

    goto DELTARAY2 ;

    BREMS1: ;

    if( J < theMaterialTable->length()-1 )
       goto NEXTMATERIAL ;
 
  return EXIT_SUCCESS;
}
