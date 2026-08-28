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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UnitsTable.hh"

#include "Analysis.hh"
#include "ApplicationParameters.hh"
#include "DetectorParameters.hh"
#include "EventAction.hh"

using namespace ApplicationParameters;
using namespace DetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction) : G4UserRunAction(), fEventAction(eventAction)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Test writing
  if (TestWrite) TestWriting();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestWriting() const
{
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
#if (defined(TEST_ANALYSIS_GENERIC))
  analysisManager->SetDefaultFileType("root");
#endif
  // Create directories
  if (TestDirectory)
  {
    analysisManager->SetHistoDirectoryName(HistoDirName);
    analysisManager->SetNtupleDirectoryName(NtupleDirName);
  }
  analysisManager->SetVerboseLevel(VerboseLevel);
  if (MergeNtuple)
  {
    analysisManager->SetNtupleMerging(true, NofNtupleFiles);
    analysisManager->SetNtupleRowWise(NtupleRowWise, NtupleRowMode);
  }
  G4cout << "Using " << analysisManager->GetType() << " verboseLevel = " << VerboseLevel;
  if (MergeNtuple)
  {
    G4cout << " mergeNtuple = " << MergeNtuple << " nofNtupleFiles = " << NofNtupleFiles << G4endl;
  }
  else
  {
    G4cout << G4endl;
  }

  // Book histograms, ntuple
  //

  // H1 properties

  if (TestH1)
  {
    // ID = 0
    analysisManager->CreateH1("Edep", "Edep in absorber", 100, 0., 100 * MeV);
    // ID = 1
    analysisManager->CreateH1("1_TrackL", "TrackL", 100, 0., NofZLayers * LayerZSize);
    // ID = 2  unit
    analysisManager->CreateH1("2_TrackL_Unit", "TrackL unit=cm", 100, 0., NofZLayers * LayerZSize,
                              "cm");
    // ID = 3  fcn (log10)
    analysisManager->CreateH1("3_TrackL_Fcn", "TrackL fcn=log10", 100, 1., NofZLayers * LayerZSize,
                              "none", "log10");
    // ID = 4  binScheme
    analysisManager->CreateH1("4_TrackL_BinScheme", "TrackL binScheme=log", 100, 1.,
                              NofZLayers * LayerZSize, "none", "none", "log");

    // ID = 5  unit & fcn (log10)
    analysisManager->CreateH1("5_TrackL_Unit_Fcn", "TrackL unit=cm & fcn=log10", 100, 1.,
                              NofZLayers * LayerZSize, "cm", "log10");
    // ID = 6  unit & binScheme
    analysisManager->CreateH1("6_TrackLength_Unit_BinScheme", "TrackL unit=cm & binScheme=log", 100,
                              1., NofZLayers * LayerZSize, "cm", "none", "log");

    // ID = 7  as ID =1
    analysisManager->CreateH1("7_TrackL", "Not in test", 100, 0., NofZLayers * LayerZSize);
    // ID = 8  as ID =1
    analysisManager->CreateH1("8_TrackL", "Not in test", 100, 0., NofZLayers * LayerZSize);
    // ID = 9  as ID =1
    analysisManager->CreateH1("9_TrackL", "Not in test", 100, 0., NofZLayers * LayerZSize);

    // Histogram dimensions

    // 1D histograms
    //
    analysisManager->CreateH1("PX", "PX in absorber layer #2",  // ID = 10
                              100, -10. * MeV, 10 * MeV);
    analysisManager->CreateH1("PX_set", "PX in absorber layer #2",  // ID = 11
                              100, -10. * MeV, 10 * MeV);
    analysisManager->CreateH1("PX_copy", "PX in absorber layer #2",  // ID = 12
                              100, -10. * MeV, 10 * MeV);
    analysisManager->CreateH1("PX_set_copy", "PX in absorber layer #2",  // ID = 13
                              100, -10. * MeV, 10 * MeV);

    // Histogram defined via edges
    // same as ID = 0
    std::vector<G4double> edges;
    G4double value = 0.;
    for (G4int i = 0; i < 100; ++i)
    {
      edges.push_back(value);
      value += 1. * MeV;
    }
    analysisManager->CreateH1("EdepEdges", "Edep in absorber - edges",  // ID 14
                              edges);

    // Non-equidistant edges
    std::vector<G4double> edges2 = {0.,        10. * MeV, 30. * MeV, 40. * MeV,
                                    60. * MeV, 70. * MeV, 90. * MeV, 100. * MeV};
    analysisManager->CreateH1("EdepEdges2", "Edep in absorber - edges 2",  // ID 15
                              edges2);

    // Ascii
    analysisManager->SetH1Ascii(10, TestAscii);

    // Plotting
    if (TestPlot)
    {
      analysisManager->SetH1Plotting(10, true);
      analysisManager->SetH1Plotting(11, true);
      analysisManager->SetH1Plotting(12, true);
      analysisManager->SetH1Plotting(13, true);
    }
  }

  if (TestH2)
  {
    // 2D histograms
    //
    analysisManager->CreateH2("PX_PY", "PX PY in absorber layer #2",  // ID = 0
                              100, -10. * MeV, 10 * MeV, 100, -10. * MeV, 10 * MeV);

    analysisManager->CreateH2("PX_PY_set", "PX PY in absorber layer #2",  // ID = 1
                              100, -10. * MeV, 10 * MeV, 100, -10. * MeV, 10 * MeV);

    analysisManager->CreateH2("PX_PY_set_dim", "PX PY in absorber layer #3",  // ID = 2
                              100, -10. * MeV, 10 * MeV, 100, -10. * MeV, 10 * MeV);

    // Ascii
    analysisManager->SetH2Ascii(0, TestAscii);

    // Plotting
    if (TestPlot)
    {
      analysisManager->SetH2Plotting(0, true);
      analysisManager->SetH2Plotting(1, true);
    }
  }

  if (TestH3)
  {
    // 3D histograms
    //
    analysisManager->CreateH3("PX_PY_PZ", "PX PY PZ in absorber layer #2",  // ID = 0
                              10, -10. * MeV, 10 * MeV, 10, -10. * MeV, 10 * MeV, 10, -10. * MeV,
                              10 * MeV);
    analysisManager->CreateH3("PX_PY_PZ_set", "PX PY PZ in absorber layer #2",  // ID = 1
                              10, -10. * MeV, 10 * MeV, 10, -10. * MeV, 10 * MeV, 10, -10. * MeV,
                              10 * MeV);
    analysisManager->CreateH3("PX_PY_PZ_set_dim", "PX PY PZ in absorber layer #3",  // ID = 2
                              10, -10. * MeV, 10 * MeV, 10, -10. * MeV, 10 * MeV, 10, -10. * MeV,
                              10 * MeV);

    // Ascii
    analysisManager->SetH3Ascii(0, TestAscii);
  }

  if (TestP1)
  {
    // 1D profiles
    //
    analysisManager->CreateP1("Edep_Vs_Z", "Longit Edep Z pseudo-profile",  // ID = 0
                              NofZLayers, 0, NofZLayers * LayerZSize, 0., 10. * GeV);
    analysisManager->CreateP1("Edep_Vs_Z_set", "Longit Edep Z pseudo-profile #2",  // ID = 1
                              NofZLayers, 0, NofZLayers * LayerZSize, 0., 10. * GeV);
    analysisManager->CreateP1("Edep_Vs_Z_set_dim", "Longit Edep Z pseudo-profile #3",  // ID = 2
                              NofZLayers, 0, NofZLayers * LayerZSize, 0., 10. * GeV);

    // Ascii
    analysisManager->SetP1Ascii(0, TestAscii);

    // Plotting
    if (TestPlot)
    {
      analysisManager->SetP1Plotting(0, true);
      analysisManager->SetP1Plotting(1, true);
    }
  }

  if (TestP2)
  {
    // 2D profiles
    analysisManager->CreateP2(
      "Edep_Vs_XY", "Longit Edep XY pseudo-profile",  // ID = 0
      NofXYLayers, -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, NofXYLayers,
      -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, 0, 10 * GeV);
    analysisManager->CreateP2(
      "Edep_Vs_XY_set", "Longit Edep XY pseudo-profile #2",  // ID = 1
      NofXYLayers, -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, NofXYLayers,
      -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, 0, 10 * GeV);
    analysisManager->CreateP2(
      "Edep_Vs_XY_set_dim", "Longit Edep XY pseudo-profile #3",  // ID = 2
      NofXYLayers, -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, NofXYLayers,
      -NofXYLayers * LayerXYSize / 2, NofXYLayers * LayerXYSize / 2, 0, 10 * GeV);
    // Ascii
    analysisManager->SetP2Ascii(0, TestAscii);
  }

  // if ( TestH1 ) {
  //   // Ilegal definitions
  //   // (objects are not created)
  //   //
  //   // log10 and xmin = 0
  //   analysisManager->CreateH1("7_TrackL_Fcn",
  //                             "TrackL fcn=log10",
  //                             100, 0., NofZLayers*LayerZSize, "none", "log10");

  //   // fcn (log10) & binScheme
  //   analysisManager->CreateH1("8_TrackLength_Unit_Fcn_BinScheme",
  //                             "TrackL fcn=log10 & binScheme=log",
  //                             100, 1., NofZLayers*LayerZSize, "none", "log10", "log");
  // }

  if (TestNtuple)
  {
    // Creating ntuples
    // ntupleId = 0
    analysisManager->CreateNtuple("EDep", "Edep and Track Length");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleSColumn("Label");
    analysisManager->CreateNtupleDColumn("EdepVector", fEventAction->GetEdepVector());
    analysisManager->CreateNtupleDColumn("TrackLengthVector", fEventAction->GetTrackLengthVector());
    analysisManager->CreateNtupleSColumn("LabelVector", fEventAction->GetLabelVector());
    analysisManager->FinishNtuple();
    // ntupleId = 1
    analysisManager->CreateNtuple("TrackL", "Track Length");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->FinishNtuple();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestReading() const
{
  // Read something from another analysis file
  G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();
  analysisReader->SetVerboseLevel(ReadVerboseLevel);

  G4String histoDir = "";
  G4String ntupleDir = "";
  if (TestDirectory)
  {
    histoDir = HistoDirName;
    ntupleDir = NtupleDirName;
  }

  // Define base file name
  analysisReader->SetFileName(FileName);

  if (TestH1)
  {
    G4int h1Id = analysisReader->ReadH1("Edep", "", histoDir);
    if (h1Id >= 0)
    {
      auto h1 = analysisReader->GetH1(h1Id);
      if (h1)
      {
        G4cout << "   H1: "
               << "   mean: " << h1->mean() << " rms: " << h1->rms() << G4endl;
      }
    }
  }

  if (TestH2)
  {
    G4int h2Id = analysisReader->ReadH2("PX_PY", "", histoDir);
    if (h2Id >= 0)
    {
      auto h2 = analysisReader->GetH2(h2Id);
      if (h2)
      {
        G4cout << "   H2: "
               << "   mean_x: " << h2->mean_x() << " rms_x: " << h2->rms_x() << G4endl << "       "
               << "   mean_y: " << h2->mean_y() << " rms_y: " << h2->rms_y() << G4endl;
      }
    }
  }

  if (TestH3)
  {
    G4int h3Id = analysisReader->ReadH3("PX_PY_PZ", "", histoDir);
    if (h3Id >= 0)
    {
      auto h3 = analysisReader->GetH3(h3Id);
      if (h3)
      {
        G4cout << "   H3: "
               << "   mean_x: " << h3->mean_x() << " rms_x: " << h3->rms_x() << G4endl << "       "
               << "   mean_y: " << h3->mean_y() << " rms_y: " << h3->rms_y() << G4endl << "       "
               << "   mean_z: " << h3->mean_z() << " rms_z: " << h3->rms_z() << G4endl;
      }
    }
  }

  if (TestP1)
  {
    G4int p1Id = analysisReader->ReadP1("Edep_Vs_Z", "", histoDir);
    if (p1Id >= 0)
    {
      auto p1 = analysisReader->GetP1(p1Id);
      if (p1)
      {
        G4cout << "   P1: "
               << "   mean: " << p1->mean() << " rms: " << p1->rms() << G4endl;
      }
    }
  }

  if (TestP2)
  {
    G4int p2Id = analysisReader->ReadP2("Edep_Vs_XY", "", histoDir);
    if (p2Id >= 0)
    {
      auto p2 = analysisReader->GetP2(p2Id);
      if (p2)
      {
        G4cout << "   P2: "
               << "   mean_x: " << p2->mean_x() << " rms_x: " << p2->rms_x() << G4endl << "       "
               << "   mean_y: " << p2->mean_y() << " rms_y: " << p2->rms_y() << G4endl;
      }
    }
  }

  if (TestNtuple)
  {
    // TrackL
    G4int ntupleId;
#if (defined(TEST_ANALYSIS_ROOT) || defined(TEST_ANALYSIS_GENERIC))
#  ifdef TEST_ANALYSIS_ROOT
    // Root output filename depends on NofNtupleFiles,
    // lets always take testAnalysis.root
    ntupleId = analysisReader->GetNtuple("TrackL", "testAnalysis.root", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_GENERIC
    // Root output filename depends on NofNtupleFiles,
    // lets always take testAnalysis.root
    ntupleId = analysisReader->GetNtuple("TrackL", "testAnalysisG.root", ntupleDir);
#  endif
#else
    if (!G4Threading::IsMultithreadedApplication()
        || (G4Threading::IsMultithreadedApplication() && !isMaster))
    {
      // In MT application only ntuples are written only on workers
      // (for other than Root output types)
      ntupleId = analysisReader->GetNtuple("TrackL", "", ntupleDir);
      // file name can be omitted then the file name will be
      // defined from the base file name like in writing phase
    }
    else
    {
      // In MT application on master thread
      // test reading from a file with explicitly given name
      // (specific to the output format)
#  ifdef TEST_ANALYSIS_HDF5
      ntupleId = analysisReader->GetNtuple("TrackL", "testAnalysis_t1.hdf5", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_XML
      ntupleId = analysisReader->GetNtuple("TrackL", "testAnalysis_nt_TrackL_t1.xml", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_CSV
      ntupleId = analysisReader->GetNtuple("TrackL", "testAnalysis_nt_TrackL_t1.csv", ntupleDir);
#  endif
    }
#endif
    if (ntupleId >= 0)
    {
      G4double trackL;
      analysisReader->SetNtupleDColumn("Labs", trackL);
      G4int counter = 0;
      G4cout << "Ntuple TrackL, reading selected column Labs" << G4endl;
      while (analysisReader->GetNtupleRow() && counter < 10)
      {
        G4cout << counter++ << "th entry: "
               << "  TrackL: " << trackL << std::endl;
      }
    }
    else
    {
      G4cout << "Failed to read TrackL ntuple." << G4endl;
    }

    // EDep
#if (defined(TEST_ANALYSIS_ROOT) || defined(TEST_ANALYSIS_GENERIC))
#  ifdef TEST_ANALYSIS_ROOT
    // Root output filename depends on NofNtupleFiles,
    // lets always take testAnalysis.root
    ntupleId = analysisReader->GetNtuple("EDep", "testAnalysis.root", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_GENERIC
    // Root output filename depends on NofNtupleFiles,
    // lets always take testAnalysis.root
    ntupleId = analysisReader->GetNtuple("EDep", "testAnalysisG.root", ntupleDir);
#  endif
#else
    if (!G4Threading::IsMultithreadedApplication()
        || (G4Threading::IsMultithreadedApplication() && !isMaster))
    {
      // In MT application only ntuples are written only on workers
      ntupleId = analysisReader->GetNtuple("EDep", "", ntupleDir);
      // file name can be omitted then the file name will be
      // define from the base file name like in writing phase
    }
    else
    {
      // In MT application on master thread
      // test reading from a file with explicitly given name
      // (specific to the output format)
#  ifdef TEST_ANALYSIS_HDF5
      ntupleId = analysisReader->GetNtuple("EDep", "testAnalysis_t1.hdf5", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_XML
      ntupleId = analysisReader->GetNtuple("EDep", "testAnalysis_nt_EDep_t1.xml", ntupleDir);
#  endif
#  ifdef TEST_ANALYSIS_CSV
      ntupleId = analysisReader->GetNtuple("EDep", "testAnalysis_nt_EDep_t1.csv", ntupleDir);
#  endif
    }
#endif

    if (ntupleId >= 0)
    {
      G4String label;
      analysisReader->SetNtupleSColumn(ntupleId, "Label", label);
      std::vector<double> edepVector;
      analysisReader->SetNtupleDColumn(ntupleId, "EdepVector", edepVector);
      std::vector<std::string> labelVector;
      analysisReader->SetNtupleSColumn(ntupleId, "LabelVector", labelVector);
      G4int counter = 0;
      G4cout << "Ntuple EDep, reading selected columns Label, EdepVector" << G4endl;
      while (analysisReader->GetNtupleRow(ntupleId) && counter < 10)
      {
        G4cout << counter++ << "th entry: "
               << "  Label: " << label << G4endl;
        G4cout << "             EdepVector: ";
        for (G4int i = 0; i < G4int(edepVector.size()); ++i)
        {
          G4cout << edepVector[i] << "  ";
        }
        G4cout << G4endl;

        G4cout << "             LabelVector: ";
        for (G4int i = 0; i < G4int(labelVector.size()); ++i)
        {
          G4cout << labelVector[i] << "  ";
        }
        G4cout << G4endl;
      }
    }
    else
    {
      G4cout << "Failed to read EDep ntuple." << G4endl;
    }
  }

  // check get methods
  for (G4int i = 0; i < analysisReader->GetNofNtuples(); ++i)
  {
    auto ntuple = analysisReader->GetNtuple(i);
    if (ntuple)
    {
      G4cout << "Ntuple " << i << " title: " << ntuple->title() << G4endl;
    }
    else
    {
      G4cout << "Failed to get ntuple " << i << G4endl;
    }
  }

  // close files
  analysisReader->CloseFiles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestGetCommands() const
{
  // test /analysis/hn|pn/get command

  if (TestH1)
  {
    auto id = 0;
    TestGetHt<tools::histo::h1d>(id);
  }
  if (TestH2)
  {
    auto id = 0;
    TestGetHt<tools::histo::h2d>(id);
  }
  if (TestH3)
  {
    auto id = 0;
    TestGetHt<tools::histo::h3d>(id);
  }
  if (TestP1)
  {
    auto id = 0;
    TestGetHt<tools::histo::p1d>(id);
  }
  if (TestP2)
  {
    auto id = 0;
    TestGetHt<tools::histo::p2d>(id);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintStatistics() const
{
  // print histogram statistics

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (isMaster)
  {
    if (TestH1)
    {
      G4int id = 0;
      auto h1 = analysisManager->GetH1(id);
      if (h1)
      {
        G4cout << "   H1: "
               << "   mean: " << h1->mean() << " rms: " << h1->rms() << G4endl << "       "
               << "       "
               << "   id: " << id
               << " id by name:: " << analysisManager->GetH1Id(analysisManager->GetH1Name(id))
               << G4endl;
      }
    }
    if (TestH2)
    {
      G4int id = 0;
      auto h2 = analysisManager->GetH2(id);
      if (h2)
      {
        G4cout << "   H2: "
               << "   mean_x: " << h2->mean_x() << " rms_x: " << h2->rms_x() << G4endl << "       "
               << "   mean_y: " << h2->mean_y() << " rms_y: " << h2->rms_y() << G4endl << "       "
               << "   id: " << id
               << " id by name:: " << analysisManager->GetH2Id(analysisManager->GetH2Name(id))
               << G4endl;
      }
    }

    if (TestH3)
    {
      G4int id = 0;
      auto h3 = analysisManager->GetH3(id);
      if (h3)
      {
        G4cout << "   H3: "
               << "   mean_x: " << h3->mean_x() << " rms_x: " << h3->rms_x() << G4endl << "       "
               << "   mean_y: " << h3->mean_y() << " rms_y: " << h3->rms_y() << G4endl << "       "
               << "   mean_z: " << h3->mean_z() << " rms_z: " << h3->rms_z() << G4endl << "       "
               << "   id: " << id
               << " id by name:: " << analysisManager->GetH3Id(analysisManager->GetH3Name(id))
               << G4endl;
      }
    }

    if (TestP1)
    {
      G4int id = 0;
      auto p1 = analysisManager->GetP1(id);
      if (p1)
      {
        G4cout << "   P1: "
               << "   mean: " << p1->mean() << " rms: " << p1->rms() << G4endl << "       "
               << "   id: " << id
               << " id by name:: " << analysisManager->GetP1Id(analysisManager->GetP1Name(id))
               << G4endl;
      }
    }

    if (TestP2)
    {
      G4int id = 0;
      auto p2 = analysisManager->GetP2(id);
      if (p2)
      {
        G4cout << "   P2: "
               << "   mean_x: " << p2->mean_x() << " rms_x: " << p2->rms_x() << G4endl << "       "
               << "   mean_y: " << p2->mean_y() << " rms_y: " << p2->rms_y() << G4endl << "       "
               << "   id: " << id
               << " id by name:: " << analysisManager->GetP2Id(analysisManager->GetP2Name(id))
               << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Test reading
  if (TestRead) TestReading();

  if (TestWrite)
  {
    G4String fileName = FileName;
    // Append filename with "W" if both Read and Write tests are selected
    if (TestRead) fileName.append("W");

    G4AnalysisManager::Instance()->OpenFile(fileName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  if (TestWrite)
  {
    // print some information
    PrintStatistics();

    // save histograms & ntuple
    G4AnalysisManager::Instance()->Write();

    // test /analysis/hn|pn/get command
    TestGetCommands();

    // close file
    G4AnalysisManager::Instance()->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
