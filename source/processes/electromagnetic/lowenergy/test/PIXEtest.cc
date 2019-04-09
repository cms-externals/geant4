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
// --------------------------------------------------------------------
//
// GEANT4     Test file
//
// File name: PIXEtest
//
// Author:    Alfonso Mantero  (alfonso.mantero@ge.infn.it)
// 
// History:
// --------
// 
// 22 Apr 2009 ALF 1st implementation based on work of Simona Saliceti
// 24 Apr 2009 ALF revision
// --------------------------------------------------------------------
//
// Test Description: 
// ------------------
// Test of second implementation of the Empiric Model for shell cross sections in proton ionisation
// --------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "G4VhShellCrossSection.hh"
#include "G4teoCrossSection.hh"
#include "G4empCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
//#include "AIDA/AIDA.h"

int main()
{ 
   G4int Z; 
   G4double incidentEnergy;
   G4double mass;
   G4double deltaEnergy;
   size_t shellNumber;
   G4String fileName;
   G4String partType;
   G4String shellName;
   G4String model;
   G4String nameId;
   G4VhShellCrossSection* shellCS;
   G4ParticleDefinition* particle;

   //    G4VhShellCrossSection* shellExp = new G4hShellCrossSectionDoubleExp();
   //  here you decide the implementation: G4teoCrossSection is ECPSSR theory: 
   //  you can choose the "analytical" implementetation of the theory or the "interpolated" 
   //  implementation from REIS et al.
   //  G4hShellCrossSectionDoubleExp is previous work with ownmade fitting functions to Paul for Alpha.

   G4cout << "Enter model (Analytical / ECPSSR_FormFactor / empirical): " << G4endl;
   G4cin >> model;

   if (model == "Analytical" || model == "ECPSSR_FormFactor") {

     shellCS = new G4teoCrossSection(model);
   }
   else if (model == "empirical") {
     shellCS = new G4empCrossSection();
   }

   G4cout << "Enter particle (p/a): " << G4endl;
   G4cin >> partType;
  
   if (partType == "p") {
   //  G4Proton* aProtone = G4Proton::Proton();
     particle = G4Proton::Proton();
   
   }
   else if (partType == "a") {
     particle  = G4Alpha::Alpha();
   }

   mass = particle->GetPDGMass();

   std::vector<G4double> energies;
   
   for (G4double i=0; i<401; i=i+1) 
     {
       energies.push_back(std::pow(10,(-2+i/100)) *MeV);

       //       energies.push_back(std::pow(10,(0.05*i+1)) *keV);
     } 
  

   G4cout << "Enter shell Index (0=K 3=L3 4-8=M1-M5): " << G4endl;
   G4cin >> shellNumber;

   if (shellNumber == 0) {shellName = "K";}
   else if (shellNumber == 1) {shellName = "L1";}
   else if (shellNumber == 2) {shellName = "L2";}
   else if (shellNumber == 3) {shellName = "L3";}
   else if (shellNumber == 4) {shellName = "M1";}
   else if (shellNumber == 5) {shellName = "M2";}
   else if (shellNumber == 6) {shellName = "M3";}
   else if (shellNumber == 7) {shellName = "M4";}
   else if (shellNumber == 8) {shellName = "M5";}


   if (model  == "Analytical" ) { nameId = "A";}
   else if (model  == "ECPSSR_FormFactor" ) { nameId = "I";}
   else if (model  == "empirical" ) { nameId = "E";}
   G4String fileNameTxt = fileName;
   char buffer[3];
   std::ofstream myfile;

   //Z is the atomic number
   for (Z = 6; Z<=92; Z++)
     { 
       G4cout << "Z = " << Z << G4endl;
       snprintf(buffer, 3, "%d", Z);
       
       fileNameTxt = nameId + "-" + shellName + "-" + partType + "-" + buffer + ".dat";
       
       myfile.open (fileNameTxt);
       G4cout << "************ here! ***************" << G4endl;       
       
       //Cross section for each incident energy
       for (size_t k=0; k<energies.size();k++)
	 {
	   incidentEnergy = energies[k];//*MeV;
	   deltaEnergy = 0.0;

	   /*
	   
	   // put zero to mimick behavior of REIS + Analytical

	   if ((incidentEnergy > 0.1 * MeV && incidentEnergy < 10*MeV) && (model  == "Analytical" )) {
	     
	     G4cout << incidentEnergy/MeV << "Energy "<<  model << " model" << G4endl;
	     myfile << incidentEnergy/MeV << "\t\t" << 0 << G4endl;
	   }

	   else { */
	     std::vector<G4double> CS = shellCS->GetCrossSection(Z,incidentEnergy,mass,deltaEnergy,false);
	     myfile << incidentEnergy/MeV << "\t\t" << CS[shellNumber]/barn << G4endl;  

	     //	   }
	   	   
	   
	 }
       
       myfile.close();
       fileNameTxt = fileName;     
       
     } 
   delete shellCS;
   
   G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}
