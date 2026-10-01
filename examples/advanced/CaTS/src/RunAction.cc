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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//   March 28th, 2026: Modified by Ilker Parmaksiz for latest Opticks
// ********************************************************************
//
/// \file RunAction.cc
/// \brief Implementation of the CaTS::RunAction class

#include <G4UserRunAction.hh>
#include "G4Run.hh"
#ifdef WITH_G4OPTICKS
#include "G4TransportationManager.hh"
#include "G4CXOpticks.hh"
#endif
#ifdef WITH_ROOT
#include "G4AnalysisManager.hh"
#include "RootIO.hh"
#endif
// project headers
#include "RunAction.hh"
#include "ConfigurationManager.hh"

RunAction::RunAction()
: G4UserRunAction() {
}

void RunAction::BeginOfRunAction(const G4Run*) {
#ifdef WITH_ROOT
    if (ConfigurationManager::getInstance()->isdoAnalysis()) {
        // Create the generic analysis manager
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->SetDefaultFileType("root");
        G4cout << "Using " << analysisManager->GetType() << G4endl;
        analysisManager->SetVerboseLevel(1);
        G4String HistoFileName =
                ConfigurationManager::getInstance()->getHistoFileName();
        G4cout << "Opening Analysis output File: " << HistoFileName << G4endl;
        analysisManager->SetFileName(HistoFileName);
        analysisManager->OpenFile();
        //
        // Book histograms, ntuple
        //
        // Creating 1D histograms
        analysisManager->CreateH1("ENeutron", "Energy of created Neutrons", 200, 0,
                100);
        analysisManager->CreateH1("EProton", "Energy of created Protons", 200, 0,
                100);
    }
#endif
#ifdef WITH_G4OPTICKS
    if (ConfigurationManager::getInstance()->isEnable_opticks()) {
        if (!geo_initialized) {
            G4cout << "\n\n###[ RunAction::BeginOfRunAction G4CXOpticks.setGeometry\n\n"
                    << G4endl;
            G4VPhysicalVolume* world =
                    G4TransportationManager::GetTransportationManager()
                    ->GetNavigatorForTracking()
                    ->GetWorldVolume();
            assert(world);
            G4CXOpticks* g4cx = G4CXOpticks::Get();
            g4cx->SetGeometry(world);

            G4cout << "\n\n###] RunAction::BeginOfRunAction G4CXOpticks.setGeometry\n\n"
                    << G4endl;
            geo_initialized = true;
        }
    }
#endif
}

void RunAction::EndOfRunAction(const G4Run*) {
#ifdef WITH_G4OPTICKS
    if (ConfigurationManager::getInstance()->isEnable_opticks()) {
        if (ConfigurationManager::getInstance()->isEnable_verbose()) {
            G4cout << "\n\n###[ RunAction::EndOfRunAction G4CXOpticks.Finalize\n\n"
                    << G4endl;
        }
        G4CXOpticks::Finalize();
        if (ConfigurationManager::getInstance()->isEnable_verbose()) {
            G4cout << "\n\n###] RunAction::EndOfRunAction G4CXOpticks.Finalize\n\n"
                    << G4endl;
        }
    }
#endif
#ifdef WITH_ROOT
    if (ConfigurationManager::getInstance()->isEnable_verbose()) {
        G4cout << "##############RunAction::EndOfRunAction" << G4endl;
    }
    if (ConfigurationManager::getInstance()->isWriteHits()) {
        if (G4Threading::IsMultithreadedApplication()) {
            RootIO::GetInstance()->Merge();
        } else {
            RootIO::GetInstance()->Close();
        }
    }
    if (ConfigurationManager::getInstance()->isdoAnalysis()) {
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->Write();
        analysisManager->CloseFile();
    }
#endif
}
