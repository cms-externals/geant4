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
// 20100804  M. Kelsey -- Add name string to ctor
// 20110719  M. Kelsey -- Add initial state code to ctor
// 20110720  M. Kelsey -- Bugs in the very last 9-body final state for pp and nn
// 20110725  M. Kelsey -- Instantiate cross-section object for self-registration
// 20110916  M. Kelsey -- Drop self-registration due to platform inconsistencies
// 20120831  M. Kelsey -- Replace T1NNChannel with separate pp and nn files.
// 20120907  M. Kelsey -- Subclass and overload findCrossSection() function.

#include "G4CascadePPChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // p p : Outgoing particle types of a given multiplicity
  static const G4int pp2bfs[1][2] =
  {{pro,pro}};

  static const G4int pp3bfs[6][3] =
  {{pro,pro,pi0}, {pro,neu,pip}, {pro,lam,kpl},
   {pro,s0,kpl},  {pro,sp,k0},   {neu,sp,kpl}};

  static const G4int pp4bfs[18][4] =
  {{pro,pro,pip,pim}, {pro,neu,pip,pi0}, {pro,pro,pi0,pi0},
   {neu,neu,pip,pip}, {pro,lam,kpl,pi0}, {pro,lam,k0,pip},
   {neu,lam,kpl,pip}, {neu,s0,kpl,pip},  {pro,s0,kpl,pi0},
   {pro,s0,k0,pip},   {pro,sm,kpl,pip},  {pro,sp,k0,pi0},
   {neu,sp,k0,pip},   {pro,sp,kpl,pim},  {neu,sp,kpl,pi0},
   {pro,pro,k0,k0b},  {pro,pro,kpl,kmi}, {pro,neu,kpl,k0b}};

  static const G4int pp5bfs[32][5] =
  {{pro,pro,pip,pim,pi0}, {pro,pro,pi0,pi0,pi0}, {pro,neu,pip,pip,pim},
   {pro,neu,pip,pi0,pi0}, {neu,neu,pip,pip,pi0}, {pro,lam,kpl,pip,pim},
   {pro,lam,kpl,pi0,pi0}, {pro,lam,k0,pip,pi0},  {pro,s0,kpl,pip,pim},
   {pro,s0,kpl,pi0,pi0},  {pro,s0,k0,pip,pi0},   {pro,sp,k0,pip,pim},
   {pro,sp,k0,pi0,pi0},   {pro,sp,kpl,pim,pi0},  {pro,sm,kpl,pip,pi0},
   {pro,sm,k0,pip,pip},   {neu,lam,kpl,pip,pi0}, {neu,lam,k0,pip,pip},
   {neu,s0,kpl,pip,pi0},  {neu,s0,k0,pip,pip},   {neu,sp,k0,pip,pi0},
   {neu,sp,kpl,pip,pim},  {neu,sp,kpl,pi0,pi0},  {neu,sm,kpl,pip,pip},
   {pro,pro,pip,k0,kmi},  {pro,pro,pim,kpl,k0b}, {pro,pro,pi0,k0,k0b},
   {pro,pro,pi0,kpl,kmi}, {pro,neu,pip,k0,k0b},  {pro,neu,pip,kpl,kmi},
   {pro,neu,pi0,kpl,k0b}, {neu,neu,pip,kpl,k0b}};

  static const G4int pp6bfs[7][6] =
  {{pro,pro,pip,pip,pim,pim}, {pro,pro,pip,pim,pi0,pi0},
   {pro,pro,pi0,pi0,pi0,pi0}, {pro,neu,pip,pip,pim,pi0},
   {pro,neu,pip,pi0,pi0,pi0}, {neu,neu,pip,pip,pip,pim},
   {neu,neu,pip,pip,pi0,pi0}};

  static const G4int pp7bfs[8][7] =
  {{pro,pro,pip,pip,pim,pim,pi0}, {pro,pro,pip,pim,pi0,pi0,pi0},
   {pro,pro,pi0,pi0,pi0,pi0,pi0}, {pro,neu,pip,pip,pip,pim,pim},
   {pro,neu,pip,pip,pim,pi0,pi0}, {pro,neu,pip,pi0,pi0,pi0,pi0},
   {neu,neu,pip,pip,pip,pim,pi0}, {neu,neu,pip,pip,pi0,pi0,pi0}};

  static const G4int pp8bfs[10][8] =
  {{pro,pro,pip,pip,pip,pim,pim,pim}, {pro,pro,pip,pip,pim,pim,pi0,pi0},
   {pro,pro,pip,pim,pi0,pi0,pi0,pi0}, {pro,pro,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,neu,pip,pip,pip,pim,pim,pi0}, {pro,neu,pip,pip,pim,pi0,pi0,pi0},
   {pro,neu,pip,pi0,pi0,pi0,pi0,pi0}, {neu,neu,pip,pip,pip,pip,pim,pim},
   {neu,neu,pip,pip,pip,pim,pi0,pi0}, {neu,neu,pip,pip,pi0,pi0,pi0,pi0}};

  static const G4int pp9bfs[11][9] =
  {{pro,pro,pip,pip,pip,pim,pim,pim,pi0}, {pro,pro,pip,pip,pim,pim,pi0,pi0,pi0},
   {pro,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pro,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
   {pro,neu,pip,pip,pip,pip,pim,pim,pim}, {pro,neu,pip,pip,pip,pim,pim,pi0,pi0},
   {pro,neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {pro,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
   {neu,neu,pip,pip,pip,pip,pim,pim,pi0}, {neu,neu,pip,pip,pip,pim,pi0,pi0,pi0},
   {neu,neu,pip,pip,pi0,pi0,pi0,pi0,pi0}};
}

namespace {
  // Total p p cross sections as a function of kinetic energy
  static const G4double ppTotXSec[30] = 
  // Stepanov cross sections below 400 MeV
  {17613.0, 863.3, 674.6, 495.2, 376.0, 285.4, 205.8, 135.7,  93.7,  69.1,
      55.2,  44.5,  38.8,  35.1,  33.0,  32.0,  44.0,  47.04, 44.86, 46.03,
      44.09, 41.81, 41.17, 40.65, 40.15, 40.18, 39.26, 38.36, 38.39, 38.41};

  static const G4double ppCrossSections[93][30] = {
  //
  // multiplicity 2 (1 channel)
  //
  //  p p
  // Stepanov cross sections below 400 MeV 
   {17613.0, 863.3, 674.6, 495.2, 376.0, 285.4, 205.8, 135.7, 93.7, 69.1,
       55.2,  44.5,  38.8,  35.1,  32.3,  26.1,  25.0,  23.5, 21.0, 18.0,
       16.0,  14.3,  12.5,  11.2,  10.3,   9.6,   9.0,   8.5,  8.0,  7.7 },
  //
  // multiplicity 3 (6 channels)
  //
  //  p p pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  1.4,  4.0,  4.3,  4.0,  4.0,
     3.6,  3.0,  2.8,  2.5,  1.7,  1.3,  1.1,  1.0,  0.9,  0.85 },

  //  p n pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.7,  4.5, 15.0, 19.1, 18.0, 16.0,
    13.0, 10.0,  8.2,  6.0,  4.3,  3.3,  2.6,  2.0,  1.65, 1.4 },

  //  p L K+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.03, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04 ,0.04, 0.04, 0.03 },

  //  p S0 K+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S+ K0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02 },

  //  n S+ K+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.02, 0.06, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02 },
  //
  // multiplicity 4 (18 channels)
  //
  //  p p pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  1.9,
     2.8,  3.0,  3.0,  2.8,  2.5,  2.1,  1.9,  1.6,  1.4,  1.2 },

  //  p n pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.05, 0.6,  3.5,
     4.0,  3.9,  3.5,  3.1,  2.8,  2.4,  2.2,  1.9,  1.7,  1.5 },

  //  p p pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 0.76,
     1.1,  1.2,  1.2,  1.1,  1.0,  0.84, 0.76, 0.64, 0.56, 0.48 },

  //  n n pi+ pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.24, 1.4,
     1.6,  1.6,  1.4,  1.2,  1.1,  1.0,  0.88, 0.76, 0.68, 0.6 },

  //  L K+ p pi0
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.05, 0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  L K0 p pi+
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.02, 0.06, 0.09, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 },

  //  L K+ n pi+
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.01, 0.04, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03 },

  //  S0 K+ n pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01 },

  //  S0 K+ p pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01,0.01 },

  //  S0 K0 p pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02 },

  //  S- K+ p pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01 },

  //  S+ K0 p pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  S+ K0 n pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  S+ K+ p pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01 },

  //  S+ K+ n pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01 },

  //  p p K0 K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p p K+ K-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p n K+ K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },
  //
  // multiplicity 5 (32 channels)
  //
  //  p p pi+ pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.06,
     0.4,  1.1,  1.8,  2.4,  2.4,  2.2,  2.0,  1.7,  1.5,  1.3 },

  //  p p pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.02,
     0.12, 0.33, 0.54, 0.72, 0.72, 0.66, 0.6,  0.51, 0.45, 0.39 },

  //  p n pi+ pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.12, 0.26,
     0.7,  1.6,  2.4,  2.6,  2.3,  2.0,  1.8,  1.6,  1.4,  1.2 },

  //  p n pi+ pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.04, 0.08,
     0.21, 0.48, 0.72, 0.78, 0.69, 0.6,  0.54, 0.48, 0.42, 0.36 },

  //  n n pi+ pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.01, 0.04,
     0.24, 0.66, 1.08, 1.44, 1.44, 1.32, 1.2,  1.0,  0.9,  0.78 },

  //  p L K+ pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.05, 0.04, 0.04, 0.03, 0.02 },

  //  p L K+ pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01 },

  //  p L K0 pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02 },

  //  p S0 K+ pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S0 K+ pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p S0 K0 pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S+ K0 pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S+ K0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.02, 0.01, 0.01, 0.01, 0.01 },

  //  p S+ K+ pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S- K+ pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p S- K0 pi+ pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  n L K+ pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n L K0 pi+ pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n S0 K+ pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S0 K0 pi+ pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S+ K0 pi+ pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S+ K+ pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.04, 0.03, 0.02, 0.02, 0.01 },

  //  n S+ K+ pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  n S- K+ pi+ pi+
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01 },

  //  p p pi+ K0 K-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.04, 0.04, 0.03 },

  //  p p pi- K+ K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  p p pi0 K0 K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.04, 0.03 },

  //  p p pi0 K+ K-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.04, 0.03 },

  //  p n pi+ K0 K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.03, 0.02, 0.02 },

  //  p n pi+ K+ K-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.05, 0.03, 0.02, 0.02 },

  //  p n pi0 K+ K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },

  //  n n pi+ K+ K0bar
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.03, 0.03, 0.02, 0.02 },
  //
  // multiplicity 6 (7 channels)
  //
  //  p p pi+ pi+ pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p p pi+ pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  p p pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.01, 0.02, 0.05, 0.1,  0.13, 0.12, 0.11, 0.1,  0.1,  0.09 },

  //  p n pi+ pi+ pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.06, 0.1,  0.18, 0.38, 0.49, 0.46, 0.43, 0.40, 0.38, 0.36 },

  //  p n pi+ pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  n n pi+ pi+ pi+ pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },

  //  n n pi+ pi+ pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.03, 0.05, 0.09, 0.19, 0.25, 0.23, 0.22, 0.2,  0.19, 0.18 },
  //
  // multiplicity 7 (8 channels)
  //
  //  p p pi+ pi+ pi- pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  p p pi+ pi- pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.1, 0.30, 0.42, 0.42, 0.42, 0.40, 0.37 },

  //  p p pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.05, 0.14, 0.20, 0.22, 0.20, 0.19, 0.18 },

  //  p n pi+ pi+ pi+ pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.19, 0.31, 0.41, 0.44, 0.47, 0.45, 0.45 },

  //  p n pi+ pi+ pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.12, 0.18, 0.24, 0.26, 0.23, 0.28, 0.26 },

  //  p n pi+ pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.02, 0.06, 0.08, 0.12, 0.13, 0.14, 0.13, 0.13 },

  //  n n pi+ pi+ pi+ pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.06, 0.17, 0.5,  0.7,  0.7,  0.69, 0.66, 0.62 },

  //  n n pi+ pi+ pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.04, 0.1, 0.30, 0.42, 0.42, 0.41, 0.40, 0.37 },
  //
  // multiplicity 8 (10 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p p pi+ pi+ pi- pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p p pi+ pi- pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p p pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.02, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  p n pi+ pi+ pi+ pi- pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  p n pi+ pi+ pi- pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  p n pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },

  //  n n pi+ pi+ pi+ pi+ pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.08, 0.18, 0.27, 0.30, 0.27, 0.24 },

  //  n n pi+ pi+ pi+ pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.04, 0.12, 0.15, 0.18, 0.15, 0.15 },

  //  n n pi+ pi+ pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.06, 0.09, 0.12, 0.09, 0.09 },
  //
  // multiplicity 9 (11 channels)
  //
  //  p p pi+ pi+ pi+ pi- pi- pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.07, 0.11, 0.14, 0.15, 0.15, 0.15 },

  //  p p pi+ pi+ pi- pi- pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.09, 0.11, 0.12, 0.12, 0.12 },

  //  p p pi+ pi- pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.04, 0.06, 0.07, 0.07, 0.07, 0.07 },

  //  p p pi0 pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04 },

  //  p n pi+ pi+ pi+ pi+ pi- pi- pi-
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  p n pi+ pi+ pi+ pi- pi- pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.03, 0.08, 0.20, 0.25, 0.29, 0.29, 0.29 },

  //  p n pi+ pi+ pi- pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.05, 0.12, 0.15, 0.17, 0.17, 0.17 },

  //  p n pi+ pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.07, 0.09, 0.10, 0.10, 0.10 },

  //  n n pi+ pi+ pi+ pi+ pi- pi- pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.06, 0.15, 0.19, 0.22, 0.22, 0.22 },

  //  n n pi+ pi+ pi+ pi- pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.02, 0.05, 0.12, 0.15, 0.17, 0.17, 0.17 },

  //  n n pi+ pi- pi0 pi0 pi0 pi0 pi0
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.01, 0.03, 0.07, 0.09 ,0.10, 0.10, 0.10 }};
}

// Initialize p-p cross-section table

const G4CascadePPChannelData::data_t
G4CascadePPChannelData::data(pp2bfs, pp3bfs, pp4bfs, pp5bfs, pp6bfs, pp7bfs,
			     pp8bfs, pp9bfs, ppCrossSections, ppTotXSec,
			     pro*pro, "ProtonProton");


// Overload base class interpolator to use function for 0-10 MeV total, elastic

G4double 
G4CascadePPChannel::findCrossSection(G4double ke,
                                     const G4double (&xsec)[30]) const {
  if (ke < 0.01 && (xsec == ppTotXSec || xsec == ppCrossSections[0])) {
    // Stepanov's function for ke < 10 MeV, up to zero-energy value
    const G4double kemin = 4.0/ppTotXSec[0];
    return (ke>0.001 ? (9.0692 - 0.0050574/ke)/ke + 6.9466 :
	    ke>kemin ? 4.0/ke : ppTotXSec[0]);
  }
  return G4PionNucSampler::findCrossSection(ke, xsec);	// Call through to base
}
