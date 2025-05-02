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
#ifndef G4VAL_HARP_H
#define G4VAL_HARP_H


// read exp.data
//
// FW data (forward = small angle with respect to projectile)
//
const int NSetsFW = 4;
//
// LA data (lateral = larger angle with respect to projectile)
//
const int NSetsLA = 9;

// general purpose counter
//
static int NSetsHARP = 0;

float** AngleBinHARP = 0;
int*    NPointsHARP = 0; 
float** XMinHARP = 0;
float** XMaxHARP = 0;
float** YHARP = 0;
float** EYHARP = 0;
float   NORM_UNCERTANTY = 0.02;

// NOTE: Based on the info from the ref. below, errors on the HARP data are published
//       as a total of stat&sys (except 2% overal normalization uncertanty which is not included)
//       See here:
//       http://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2008-26.pdf
//       Journal publications are not very explicit on this matter

#endif
