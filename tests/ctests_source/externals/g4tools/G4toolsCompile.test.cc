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

// This is a totally dumb, but simple way, to check that the tools headers compile
// correctly in case they are not included in a source file.
//
// It also provides a ways to run clang-tidy on individual tools headers as they must be
// #included in a source file for clang-tidy to see them. To run clang-tidy on an individual
// tools header, setup for Geant4 and clang-tidy as normal (see CODING_GUIDELINES.rst), then run
//
// $ run-clang-tidy tests/ctests_source/externals -header-filter="tools/HEADERTOTEST"
//
// We do leave out some .icc files and those which pull in externals like zlib, xml...

// TODO: Expand this list as needed
#include "tools/aida_ntuple"

// This is a test suite as such, but a test to check that the tools headers compile.
// We put an empty case in here just so we have something to compile and link.
#include <gtest/gtest.h>

TEST(tools_clang_tidy, Basic) {}
