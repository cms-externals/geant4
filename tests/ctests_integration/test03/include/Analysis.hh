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
/// \file Analysis.hh
/// \brief Selection of the analysis technology

#ifndef Analysis_h
#  define Analysis_h 1

#  ifdef TEST_ANALYSIS_CSV
#    include "G4CsvAnalysisManager.hh"
#    include "G4CsvAnalysisReader.hh"
using G4AnalysisManager = G4CsvAnalysisManager;
using G4AnalysisReader = G4CsvAnalysisReader;
#  endif

#  ifdef TEST_ANALYSIS_HDF5
#    include "G4Hdf5AnalysisManager.hh"
#    include "G4Hdf5AnalysisReader.hh"
using G4AnalysisManager = G4Hdf5AnalysisManager;
using G4AnalysisReader = G4Hdf5AnalysisReader;
#  endif

#  ifdef TEST_ANALYSIS_ROOT
#    include "G4RootAnalysisManager.hh"
#    include "G4RootAnalysisReader.hh"
using G4AnalysisManager = G4RootAnalysisManager;
using G4AnalysisReader = G4RootAnalysisReader;
#  endif

#  ifdef TEST_ANALYSIS_XML
#    include "G4XmlAnalysisManager.hh"
#    include "G4XmlAnalysisReader.hh"
using G4AnalysisManager = G4XmlAnalysisManager;
using G4AnalysisReader = G4XmlAnalysisReader;
#  endif

#  ifdef TEST_ANALYSIS_GENERIC
#    include "G4GenericAnalysisManager.hh"
#    include "G4RootAnalysisReader.hh"
using G4AnalysisManager = G4GenericAnalysisManager;
using G4AnalysisReader = G4RootAnalysisReader;
#  endif

#endif
