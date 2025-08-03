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
#ifndef G4PHYSLIST_HE_H
#define G4PHYSLIST_HE_H

const int NModels_HE = 3;
std::string ModelName_HE[4]  = { "ftfp_bert", "qgsp_bert", "NuBeam", "NuBeam-3-3.5" };
//
//int         ColorModel_HE[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
///
// magenta is reserved for Bertini !
// ftfp_bert is light blue (7) like ftfp in model-level
// NuBeam is red is in qgsp-g4lund-str-fragm
// NuBeam-3-3.5 is Green as a "green light to improved NuBeam"
// qgsp_bert is black like qgsp in model-level
//
int         ColorModel_HE[5] = { 7, kBlack, kRed, kGreen, 14 }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel_HE[4]     = { 8, 25, 21, 23 };

/*
   enum EMarkerStyle {kDot=1, kPlus, kStar, kCircle=4, kMultiply=5,
                      kFullDotSmall=6, kFullDotMedium=7, kFullDotLarge=8,
                      kFullCircle=20, kFullSquare=21, kFullTriangleUp=22,
                      kFullTriangleDown=23, kOpenCircle=24, kOpenSquare=25,
                      kOpenTriangleUp=26, kOpenDiamond=27, kOpenCross=28,
                      kFullStar=29, kOpenStar=30, kOpenTriangleDown=32,
                      kFullDiamond=33, kFullCross=34};
*/

#endif
