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
#ifndef G4VAL_NA49_H
#define G4VAL_NA49_H


// read exp.data
//

int    NPointsNA49         = 0;

float* xF              = 0;

float* dXSecdxF     = 0;
float* err_dXSecdxF = 0;
float* dNdxF           = 0;
float* err_dNdxF       = 0;
float* pT              = 0;
float* err_pT          = 0;
float* pT2             = 0;
float* err_pT2         = 0;

float* XSec            = 0;
float* err_XSec        = 0;

// for pions (pi+ & pi-)
//
float* Y =0;
float* dNdY[2]     = { 0, 0 };
float* err_dNdY[2] = { 0, 0 };

int NSubSetsNA49 = 0;
double*     DDiff_xF = 0;
TObjArray*  DDiffDataHolder = 0;

static int isNA49UtilLoaded = 0;

#endif
