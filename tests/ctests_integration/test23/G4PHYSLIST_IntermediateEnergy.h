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
#ifndef G4PHYSLIST_IE_H
#define G4PHYSLIST_IE_H

const int NModels_IE = 4;
// --> std::string ModelName_IE[2]  = { "ftfp_bert", "NuBeam" };
std::string ModelName_IE[4]  = { "ftfp_bert", "qgsp_bert", "Shielding", "ShieldingM" };
// --> !!! std::string ModelName_IE[3]  = { "ftfp_bert", "Shielding", "ShieldingM" };
// int         ColorModel_IE[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
// int         ColorModel_IE[5] = { 7, kMagenta, kGreen, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel_IE[5] = { 7, kRed, kGreen, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel_IE[4]     = { 8, 21, 23, 25 };

#endif
