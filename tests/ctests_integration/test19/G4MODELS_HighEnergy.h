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
#ifndef G4MODELS_HE_H
#define G4MODELS_HE_H

const int NModels_HE = 4;
//std::string ModelName_HE[3]  = { "qgsp", "ftfp", "qgsp-g4lund-str-fragm" };
std::string ModelName_HE[4]  = { "qgsp", "ftfp", "ftfp_tune3", "fluka4.4.0" };
//int         ColorModel_HE[5] = { kMagenta, 7, kRed, kBlue /* kBlack */ , 14 }; // 14 = grey, 7 = light "sky"-blue
// kMagenta is for BERTINI !!!
// int         ColorModel_HE[5] = { kBlue /* kBlack */ , 7, kRed, 14, kGreen }; // 14 = grey, 7 = light "sky"-blue
// int         SymbModel_HE[4]     = { 8, 21, 25, 23 };
// int         ColorModel_HE[5] = { 7, kRed, kGreen, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel_HE[6] = { 7, kRed, kGreen, kBlue, 14, kMagenta }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel_HE[4]     = { 8, 25, 21, 23 };

#endif
