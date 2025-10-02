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
#ifndef G4MODELS_IE_H
#define G4MODELS_IE_H

const int NModels_IE = 2;
// std::string ModelName_IE[3]  = { "inclxx", "bertini", "ftfp" };
// std::string ModelName_IE[3]  = { "bertini", "ftfp", "inclxx" };
// std::string ModelName_IE[4]  = { "bertini", "ftfp", "ftfp_tune3", "fluka.cern.4.4.0" };
std::string ModelName_IE[2]  = { "ftfp", "fluka.cern.4.4.0" };
// std::string ModelName_IE[1] = {"ftfp"};
//
/* ---> for VMP business
std::string ModelName_IE[9]  = { "ftfp-default",
                                 "ftfp-pr0-0.1-pr1-0.1", 
				 "ftfp-pr0-0.5-pr1-0.5",
				 "ftfp-pr0-0.1-pr1-0.5",
				 "ftfp-pr0-0.4-pr1-0.25",
                                 "ftfp-pr0-a1-6.5", "ftfp-pr0-a1-19.5", 
				 "ftfp-pr0-b1-0.85", "ftfp-pr0-b1-2.65"  };
*/
// int         ColorModel_IE[5] = { kMagenta, 7, kRed, kBlue /* kBlack */ , 14 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel_IE[6] = { kMagenta, kRed, kGreen, kBlue, 14, 7 }; // 14 = grey, 7 = light "sky"-blue
// ---> for VMP --> int         ColorModel_IE[5] = { kRed, kGreen, kBlue /* kBlack */ , 7, kMagenta }; 
int         SymbModel_IE[4]     = { 8, 21, 23, 25 };

#endif
