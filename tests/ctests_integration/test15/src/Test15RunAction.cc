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
// RunAction program
// --------------------------------------------------------------

#include "Test15RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "Test15Run.hh"

#include <fstream>
#include <iomanip>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test15RunAction::Test15RunAction() : G4UserRunAction(), analysisManager(nullptr)
{
  number_shells = 26;
  //  G4double shell_thickness = 10.0*mm;
  G4double shell_thickness = 2.0 * mm;

  shell_outer_radius = 457.0 * mm;
  shell_inner_radius = shell_outer_radius - shell_thickness;

  G4double radii_start[] = {200.0 * cm, 190.0 * cm, 185.0 * cm, 175.0 * cm, 165.0 * cm, 150.0 * cm,
                            140.0 * cm, 130.0 * cm, 120.0 * cm, 110.0 * cm, 100.0 * cm, 90.0 * cm,
                            80.0 * cm,  70.0 * cm,  60.0 * cm,  50.0 * cm,  45.7 * cm,  40.0 * cm,
                            30.0 * cm,  25.0 * cm,  20.0 * cm,  15.0 * cm,  10.0 * cm,  8.0 * cm,
                            5.0 * cm,   3.0 * cm};

  for (G4int i = 0; i < number_shells; ++i)
  {
    outer_radius[i] = radii_start[i];
    inner_radius[i] = radii_start[i] - shell_thickness;
  }

  for (G4int i = 0; i < 4; ++i)
  {
    local_energy_integral[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test15RunAction::~Test15RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* Test15RunAction::GenerateRun()
{
  return new Test15Run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // Get analysis manager
  analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Open an output file
  //
  G4String fileName = "Test15_output";
  analysisManager->OpenFile(fileName);

  // Create directories
  // analysisManager->SetHistoDirectoryName("histograms");
  // analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //

  // Creating histograms - primary beam investigation:
  analysisManager->CreateH1("1", "Gamma Edep /keV", 1000, 0., 1000 * keV);
  analysisManager->CreateH1("2", "Neutron ener vs. 1/mom /eV", 100000, 0., 1000000.);
  analysisManager->CreateH1("3", "Electron Edep /keV", 1000, 0., 1000 * keV);
  analysisManager->CreateH1("4", "Positron Edep /keV", 1000, 0., 1000 * keV);
  analysisManager->CreateH1("5", "Other Edep /keV", 1000, 0., 1000 * keV);

  analysisManager->CreateH1("6", "Particle Stack", 12, 0.5, 12.5);
  analysisManager->CreateH1("7", "Neutrons/event", 30, 0., 30.);
  analysisManager->CreateH1("8", "Protons/event", 30, 0., 30.);

  analysisManager->CreateH2("11", "Neutron Energy vs. Time", 100, -1.0, 4.0, 100, -2.0, 7.0);
  analysisManager->CreateH2("12", "OTHER particle Energy vs. Time", 100, -1.0, 4.0, 100, -2.0, 7.0);

  // Creating ntuples
  //
  analysisManager->CreateNtuple("Test15_Secondaries", "Secondary Particle Info");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleIColumn("particle");
  analysisManager->CreateNtupleDColumn("momentum");
  analysisManager->CreateNtupleIColumn("parentid");
  analysisManager->CreateNtupleDColumn("e_prim");
  analysisManager->CreateNtupleIColumn("parent");
  analysisManager->CreateNtupleDColumn("e_parent");
  analysisManager->CreateNtupleIColumn("numgen");
  analysisManager->CreateNtupleIColumn("event");
  analysisManager->FinishNtuple();  // ntupleID: 0 - filled

  analysisManager->CreateNtuple("Test15 Energy Time", "Neutron Time");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("primary");
  analysisManager->FinishNtuple();  // ntupleID: 1

  analysisManager->CreateNtuple("Test15 Exiting", "Neutrons Exiting");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->FinishNtuple();  // ntupleID: 2 - filled

  analysisManager->CreateNtuple("Test15 Flux 4002", "Neutrons Test15 flux");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("tarcflux");
  analysisManager->CreateNtupleDColumn("errstat");
  analysisManager->CreateNtupleDColumn("errsyst");
  analysisManager->CreateNtupleDColumn("g4flux");
  analysisManager->CreateNtupleDColumn("g4perp");
  analysisManager->CreateNtupleDColumn("gfluence");
  analysisManager->CreateNtupleDColumn("g4err");
  analysisManager->CreateNtupleDColumn("rawflux");
  analysisManager->CreateNtupleDColumn("trceflux");
  analysisManager->CreateNtupleDColumn("g4eflux");
  analysisManager->CreateNtupleDColumn("gstep");
  analysisManager->CreateNtupleDColumn("gfl_cyl");
  analysisManager->CreateNtupleDColumn("g4front");
  analysisManager->CreateNtupleDColumn("g4_shell");
  analysisManager->CreateNtupleDColumn("g4_shell_err");
  analysisManager->FinishNtuple();  // ntupleID: 3

  analysisManager->CreateNtuple("Test15 Flux 4004", "Neutrons Test15 flux");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("tarcflux");
  analysisManager->CreateNtupleDColumn("errstat");
  analysisManager->CreateNtupleDColumn("errsyst");
  analysisManager->CreateNtupleDColumn("g4flux");
  analysisManager->CreateNtupleDColumn("g4perp");
  analysisManager->CreateNtupleDColumn("gfluence");
  analysisManager->CreateNtupleDColumn("g4err");
  analysisManager->CreateNtupleDColumn("rawflux");
  analysisManager->CreateNtupleDColumn("gstep");
  analysisManager->CreateNtupleDColumn("gfl_cyl");
  analysisManager->CreateNtupleDColumn("g4front");
  analysisManager->CreateNtupleDColumn("g4_shell");
  analysisManager->CreateNtupleDColumn("g4_shell_err");
  analysisManager->FinishNtuple();  // ntupleID: 4

  analysisManager->CreateNtuple("Test15 Flux 4005", "Neutrons Test15 flux");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("tarcflux");
  analysisManager->CreateNtupleDColumn("errstat");
  analysisManager->CreateNtupleDColumn("errsyst");
  analysisManager->CreateNtupleDColumn("g4flux");
  analysisManager->CreateNtupleDColumn("g4perp");
  analysisManager->CreateNtupleDColumn("gfluence");
  analysisManager->CreateNtupleDColumn("g4zflux");
  analysisManager->CreateNtupleDColumn("g4err");
  analysisManager->CreateNtupleDColumn("flux5cm");
  analysisManager->CreateNtupleDColumn("err5cm");
  analysisManager->CreateNtupleDColumn("rawflux");
  analysisManager->CreateNtupleDColumn("gstep");
  analysisManager->CreateNtupleDColumn("gfl_cyl");
  analysisManager->CreateNtupleDColumn("g4front");
  analysisManager->CreateNtupleDColumn("g4_shell");
  analysisManager->CreateNtupleDColumn("g4_shell_err");
  analysisManager->FinishNtuple();  // ntupleID: 5

  analysisManager->CreateNtuple("Test15 Created Neutrons", "Created Neutrons");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleIColumn("particle");
  analysisManager->CreateNtupleDColumn("momentum");
  analysisManager->CreateNtupleDColumn("zmom");
  analysisManager->FinishNtuple();  // ntupleID: 6

  analysisManager->CreateNtuple("Test15 Created Neutrons", "Created Neutrons");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("starte");
  analysisManager->CreateNtupleIColumn("trackid");
  analysisManager->CreateNtupleIColumn("parentid");
  analysisManager->CreateNtupleDColumn("fluxe");
  analysisManager->CreateNtupleDColumn("fluxidx");
  analysisManager->CreateNtupleDColumn("zmom");
  analysisManager->CreateNtupleDColumn("startt");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("e_parent");
  analysisManager->CreateNtupleIColumn("parent");
  analysisManager->CreateNtupleIColumn("step");
  analysisManager->CreateNtupleIColumn("dupli");
  analysisManager->FinishNtuple();  // ntupleID: 7

  analysisManager->CreateNtuple("Test15 Radial Shell Fluence", "Radial Shell Fluence");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("fluence");
  analysisManager->CreateNtupleDColumn("true_e");
  analysisManager->CreateNtupleDColumn("true_f");
  analysisManager->FinishNtuple();  // ntupleID: 8

  analysisManager->CreateNtuple("Test15 Radial Fluence Data", "Radial Fluence Data");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("data");
  analysisManager->CreateNtupleDColumn("error");
  analysisManager->FinishNtuple();  // ntupleID: 9

  analysisManager->CreateNtuple("Test15 Radial Fluence He3", "Radial Fluence He3");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("data");
  analysisManager->CreateNtupleDColumn("error");
  analysisManager->FinishNtuple();  // ntupleID: 10

  analysisManager->CreateNtuple("Test15 3.5GeV He3 experimental data", "Radial He3 Exp Data");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("data");
  analysisManager->CreateNtupleDColumn("stat_err");
  analysisManager->CreateNtupleDColumn("syst_err");
  analysisManager->FinishNtuple();  // ntupleID: 11

  analysisManager->CreateNtuple("Test15 Radial Fluence Li", "Radial Fluence Li");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("data");
  analysisManager->CreateNtupleDColumn("stat_err");
  analysisManager->CreateNtupleDColumn("syst_err");
  analysisManager->FinishNtuple();  // ntupleID: 12

  analysisManager->CreateNtuple("Test15 Radial Fluence", "Radial Fluence");
  analysisManager->CreateNtupleDColumn("radius");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("fluence");
  analysisManager->CreateNtupleDColumn("he_data");
  analysisManager->FinishNtuple();  // ntupleID: 13

  analysisManager->CreateNtuple("Test15 Energy Time", "OTHER Time");
  analysisManager->CreateNtupleDColumn("energy");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("primary");
  analysisManager->FinishNtuple();  // ntupleID: 14
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::EndOfRunAction(const G4Run* aRun)
{
  FillRadialExperimentalData();

  G4int nofEvents = aRun->GetNumberOfEvent();
  G4cout << " Number of events is: " << nofEvents << G4endl;
  if (nofEvents == 0) return;

  const Test15Run* tarcRun = static_cast<const Test15Run*>(aRun);
  G4double ExitFlux = tarcRun->GetExitingFlux();
  G4double ExitGrichineFlux = tarcRun->GetExitingGrichineFlux();
  G4double ExitEnergy = tarcRun->GetExitingEnergy();
  G4double ExitCheckFlux = tarcRun->GetExitingCheckFlux();

  // print
  //
  if (isMaster)
  {
    createNeutronFluxHisto(aRun->GetNumberOfEvent(), tarcRun);
    createRadialFluxHisto(aRun->GetNumberOfEvent(), tarcRun);
    secondarySummary(aRun->GetNumberOfEvent(), tarcRun);
    G4cout << G4endl << "--------------------End of Global Run-----------------------" << G4endl
           << "  The run was " << nofEvents << " events " << G4endl << G4endl
           << " Integral Neutron Flux @ 46cm: " << tarcRun->GetIntegralFlux_46cm() << G4endl
           << " Integral EFLUX @ 46cm: " << tarcRun->GetIntegralEFlux_46cm() << G4endl;
  }
  else
  {
    G4cout << G4endl << "--------------------End of Local Run------------------------" << G4endl
           << "  The run was " << nofEvents << " ";
  }
  G4cout << " Exiting Grichine Flux : " << ExitGrichineFlux << G4endl
         << " Exiting Flux : " << ExitFlux << G4endl << " Total Exiting Energy : " << ExitEnergy
         << G4endl << " Exiting Check Flux : " << ExitCheckFlux << G4endl
         << "------------------------------------------------------------" << G4endl << G4endl;

  G4cout << " Gamma Edep : mean = " << analysisManager->GetH1(1)->mean()
         << " rms = " << analysisManager->GetH1(1)->rms() << G4endl;

  G4cout << " Neutron Lethargy : mean = " << analysisManager->GetH1(2)->mean()
         << " rms = " << analysisManager->GetH1(2)->rms() << G4endl;

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::FillRadialExperimentalData()
{
  G4double radial_mean_energy[10] = {0.1,   1.5,   5.0,    10.0,    18.0,
                                     100.0, 480.0, 1000.0, 10000.0, 50000.0};

  G4double radial_distance_1[76] = {
    -187.5, -164.5, -154.6, -132.7, -124.8, -121.2, -111.5, -107.4, -105.3, -97.5, -95.2,
    -94.0,  -90.3,  -76.9,  -70.8,  -70.8,  -64.1,  -60.5,  -58.6,  -54.6,  -50.3, -45.6,
    -40.4,  -27.0,  27.0,   40.4,   45.6,   46.0,   50.3,   54.6,   58.6,   60.5,  60.5,
    67.5,   69.1,   70.8,   70.8,   76.9,   79.7,   83.9,   85.2,   90.3,   90.3,  90.3,
    94.0,   95.2,   102.0,  105.3,  105.3,  111.5,  112.5,  117.4,  120.2,  124.8, 127.5,
    127.5,  127.5,  131.0,  133.5,  137.7,  143.3,  144.1,  144.1,  145.6,  150.2, 153.9,
    154.6,  154.6,  154.6,  164.5,  164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_1[76] = {
    178127,  208387,  252267,  652345,  721648,  695630,  861361,  951128,  992298,  1098173,
    1166686, 1075235, 1214726, 1341157, 1143483, 1401114, 1389402, 1479789, 1591117, 1793428,
    1734488, 1660183, 1906804, 1998677, 1946504, 1799165, 1965170, 1404543, 1709736, 1479024,
    1612297, 1743914, 1638378, 1465282, 1534629, 1367310, 1484603, 1318944, 1307921, 1048934,
    1221003, 1291466, 1157944, 1106473, 1070874, 1101744, 992297,  844377,  925665,  891084,
    859081,  815832,  725850,  671149,  730157,  603684,  532314,  556694,  609798,  540992,
    480779,  445007,  350657,  469112,  242026,  348872,  251507,  269597,  243829,  198383,
    242238,  279758,  216177,  184382,  106216,  159332};

  G4double radial_error_1[76] = {
    20916,  15236,  19480,  26868,  31621,  37039,  37286, 43864, 41004, 43472,  50820,
    52536,  59741,  80500,  106866, 60049,  75086,  75960, 69629, 94551, 72263,  71642,
    100561, 105661, 108645, 97110,  101310, 214244, 71958, 93239, 68516, 117159, 79586,
    64304,  62931,  55589,  73533,  64449,  52070,  68685, 59340, 92312, 46655,  61307,
    48136,  53235,  53476,  130788, 39518,  36844,  34146, 32958, 28554, 30435,  28978,
    25423,  35862,  32606,  26898,  21415,  23472,  19050, 16957, 18525, 18278,  16337,
    11483,  12693,  9988,   9593,   10798,  11411,  12280, 7422,  20669, 21266};

  std::ofstream fLowRadialFile;
  G4String name1("LowEnergyRadialData.dat");
  fLowRadialFile.open(name1);

  for (G4int i = 0; i < 76; ++i)
  {
    copy_radial_fluence_1[i] = radial_fluence_1[i];
    copy_radial_error_1[i] = radial_error_1[i];

    analysisManager->FillNtupleDColumn(9, 0, radial_distance_1[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[0]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_1[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_1[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[0] << " radius: " << radial_distance_1[i]
                   << " Data: " << radial_fluence_1[i] << G4endl;
  }

  G4double radial_distance_2[102] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,  -64.1,  -60.5,
    -58.6,  -54.6,  -54.6,  -50.3,  -50.3,  -45.6,  -40.4,  -27.0,  -27.0,  27.0,   40.4,   45.6,
    46.0,   50.3,   54.6,   58.6,   60.5,   60.5,   60.5,   67.5,   67.5,   69.1,   70.8,   70.8,
    70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   85.2,   90.3,   90.3,   90.3,   90.3,   92.8,
    94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,  112.5,  117.4,  117.4,
    117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,  127.5,  131.0,  133.5,  137.7,
    143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,  154.6,  154.6,  164.5,
    164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_2[102] = {
    52037,  60333,  73270,  75772,  172600, 208772, 208477, 262848, 251775, 285943, 269209, 291771,
    352464, 340796, 362096, 366357, 443692, 374540, 402242, 481558, 477824, 445496, 420413, 475073,
    489428, 479350, 540276, 524329, 515078, 538414, 573929, 563951, 610220, 597918, 546907, 536203,
    435567, 498337, 517116, 455746, 491059, 474549, 484480, 415905, 433361, 422865, 388025, 434043,
    412295, 364314, 395237, 354786, 378500, 351309, 374363, 324092, 308537, 341238, 351774, 368103,
    304353, 321900, 257844, 283505, 271252, 268307, 252442, 242475, 246286, 218185, 216957, 213167,
    227283, 186689, 191253, 198983, 200761, 186511, 166697, 161090, 155814, 150208, 170289, 134541,
    132267, 138269, 121826, 103701, 116794, 69495,  73972,  96561,  71825,  71305,  68455,  57475,
    58911,  69183,  61544,  47362,  29628,  34947};

  G4double radial_error_2[102] = {
    7250,  2796,  4200,  4035,  7786,  9666,  11622, 19126, 11478, 13805, 12593, 13215, 16428,
    16835, 19401, 18540, 25199, 49984, 18332, 31378, 23938, 23782, 24629, 24555, 23622, 36273,
    31094, 24304, 23643, 23636, 31042, 56170, 34302, 34339, 33215, 26304, 49354, 22665, 29952,
    21531, 25864, 25428, 25025, 38030, 20431, 20274, 17357, 22920, 27045, 37569, 20274, 16006,
    37086, 22740, 19518, 19276, 13998, 19812, 21177, 19033, 14503, 16641, 21189, 16470, 22573,
    12532, 12043, 11123, 12376, 9735,  18756, 24040, 10334, 8444,  9300,  9080,  8998,  8635,
    7994,  11119, 12734, 9615,  8173,  6071,  9824,  6995,  5621,  5192,  5352,  8780,  4350,
    5095,  3481,  3458,  3110,  2947,  2877,  3199,  3585,  2193,  1587,  1689};

  for (G4int i = 0; i < 102; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_2[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[1]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_2[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_2[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[1] << " radius: " << radial_distance_2[i]
                   << " Data: " << radial_fluence_2[i] << G4endl;
  }

  G4double radial_distance_3[102] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -190.3, -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,  -64.1,  -60.5,
    -58.6,  -54.6,  -54.6,  -50.3,  -50.3,  -45.6,  -40.4,  -27.0,  -27.0,  27.0,   40.4,   45.6,
    46.0,   50.3,   54.6,   58.6,   60.5,   60.5,   60.5,   67.5,   67.5,   69.1,   70.8,   70.8,
    70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   85.2,   90.3,   90.3,   90.3,   90.3,   92.8,
    94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,  112.5,  117.4,  117.4,
    117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,  127.5,  131.0,  133.5,  137.7,
    143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,  154.6,  154.6,  164.5,
    164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_3[102] = {
    19654,  23255,  28522,  29988,  66240,  81051,  82117,  105704, 97672,  110736, 104929, 112948,
    139589, 135880, 144366, 148753, 173919, 151538, 152228, 183371, 192210, 178006, 171565, 190845,
    198896, 198543, 221945, 207810, 208064, 214427, 234634, 219351, 246342, 246591, 224161, 214558,
    178277, 200567, 204048, 187469, 195491, 190137, 193141, 159567, 171050, 172297, 151338, 171693,
    164591, 173333, 155287, 139681, 149211, 139398, 147546, 127203, 120764, 135269, 135597, 144734,
    119849, 124716, 103624, 109857, 104741, 104519, 97749,  93341,  94634,  82193,  85547,  85349,
    88025,  70531,  73397,  78090,  78385,  69544,  64191,  60962,  59521,  57224,  64979,  50557,
    47420,  52529,  46303,  39465,  43679,  26193,  27994,  36002,  27536,  26758,  26318,  21902,
    21728,  25603,  23063,  17841,  11053,  12680};

  G4double radial_error_3[102] = {
    3775,  1093,  1666,  1721,  3049,  3790,  4661,  9652,  4424,  5439,  4826,  5058, 6534,
    7070,  7872,  8028,  10862, 25055, 7128,  14961, 11922, 10140, 10552, 10178, 9008, 16385,
    12552, 10763, 9314,  9596,  13543, 29333, 14200, 15097, 13512, 11421, 34321, 9076, 12301,
    8673,  12339, 10444, 10628, 19067, 8002,  7824,  6762,  9521,  11145, 20852, 8340, 6195,
    15392, 9484,  7751,  9142,  5280,  7885,  8838,  8759,  5903,  6937,  10581, 7022, 10277,
    4950,  4845,  4274,  4839,  3672,  11095, 14371, 3608,  3157,  3470,  3603,  3459, 3096,
    2978,  4475,  4677,  3867,  3196,  2232,  4926,  2765,  2133,  2013,  1946,  3588, 1728,
    1892,  1363,  1340,  1193,  1139,  1135,  1154,  1463,  819,   672,   641};

  for (G4int i = 0; i < 102; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_3[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[2]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_3[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_3[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[2] << " radius: " << radial_distance_3[i]
                   << " Data: " << radial_fluence_3[i] << G4endl;
  }

  G4double radial_distance_4[102] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,  -64.1,  -60.5,
    -58.6,  -54.6,  -54.6,  -50.3,  -50.3,  -45.6,  -40.4,  -27.0,  -27.0,  27.0,   40.4,   45.6,
    46.0,   50.3,   54.6,   58.6,   60.5,   60.5,   60.5,   67.5,   67.5,   69.1,   70.8,   70.8,
    70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   85.2,   90.3,   90.3,   90.3,   90.3,   92.8,
    94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,  112.5,  117.4,  117.4,
    117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,  127.5,  131.0,  133.5,  137.7,
    143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,  154.6,  154.6,  164.5,
    164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_4[102] = {
    10787,  12973,  16029,  16823,  37060,  45469,  46469,  59793,  54722,  61839,  59306,  63596,
    79278,  77315,  82491,  84399,  97383,  87004,  83784,  103674, 109837, 101372, 98073,  109123,
    115080, 114736, 129676, 119080, 120122, 121708, 136590, 121178, 141880, 144317, 130661, 123695,
    103661, 115618, 114399, 110128, 112016, 109275, 109903, 91304,  97038,  100656, 85314,  97644,
    92960,  93183,  87705,  79562,  84859,  78723,  83401,  72673,  68564,  76771,  75967,  81769,
    68056,  69807,  58161,  61592,  58092,  58780,  54528,  52244,  52313,  45462,  47536,  47586,
    49418,  39064,  40854,  43521,  43647,  38148,  35901,  33429,  32990,  31792,  36060,  27950,
    25418,  28993,  25652,  21769,  24058,  14371,  15479,  19639,  15312,  14715,  14684,  12111,
    11859,  13995,  12605,  9860,   6029,   6884};

  G4double radial_error_4[102] = {
    2087, 682,   1040, 1044, 1781,  2252, 2928, 4664, 2772, 3443, 2869, 2944, 4101, 4362,  5048,
    5004, 6723,  7565, 4240, 8171,  6506, 6318, 6950, 6420, 6008, 9192, 8399, 6141, 5597,  5960,
    8962, 16016, 9266, 8974, 8990,  6426, 9005, 5571, 8418, 5254, 6757, 6325, 6871, 10339, 4866,
    4916, 4015,  5502, 6999, 12350, 5109, 3735, 9087, 6080, 5014, 5372, 3178, 5176, 5568,  4651,
    3515, 4261,  5694, 4372, 5792,  2908, 2858, 2581, 2905, 2131, 5508, 6566, 2269, 1821,  2149,
    2176, 2088,  1714, 1832, 2764,  3382, 2436, 1891, 1310, 2519, 1676, 1269, 1220, 1119,  1146,
    1026, 1159,  838,  845,  697,   703,  672,  670,  915,  472,  409,  386};

  for (G4int i = 0; i < 102; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_4[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[3]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_4[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_4[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[3] << " radius: " << radial_distance_4[i]
                   << " Data: " << radial_fluence_4[i] << G4endl;
  }

  G4double radial_distance_5[102] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,  -64.1,  -60.5,
    -58.6,  -54.6,  -54.6,  -50.3,  -50.3,  -45.6,  -40.4,  -27.0,  -27.0,  27.0,   40.4,   45.6,
    46.0,   50.3,   54.6,   58.6,   60.5,   60.5,   60.5,   67.5,   67.5,   69.1,   70.8,   70.8,
    70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   85.2,   90.3,   90.3,   90.3,   90.3,   92.8,
    94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,  112.5,  117.4,  117.4,
    117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,  127.5,  131.0,  133.5,  137.7,
    143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,  154.6,  154.6,  164.5,
    164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_5[102] = {
    6356,  7770,  9671,  10087, 22305, 27386, 28211, 36102, 32910, 37046, 36039, 38516, 48296,
    47140, 50574, 51043, 58376, 53488, 49520, 63378, 67202, 61844, 59825, 66840, 71393, 70743,
    81271, 73381, 74414, 73987, 85294, 71502, 87630, 90560, 81680, 76659, 64615, 71544, 68582,
    69490, 68951, 67424, 67016, 56597, 59077, 63197, 51658, 59608, 56113, 54894, 53143, 48723,
    51845, 47582, 50552, 44713, 41868, 46715, 45754, 49557, 41508, 41942, 34765, 37086, 34516,
    35491, 32622, 31437, 30984, 27074, 28172, 28223, 29812, 23289, 24423, 25922, 25999, 22529,
    21580, 19663, 19634, 18993, 21502, 16645, 14701, 17184, 15277, 12887, 14273, 8467,  9204,
    11513, 9145,  8705,  8803,  7190,  6976,  8241,  7400,  5866,  3533,  4036};

  G4double radial_error_5[102] = {
    961,  395,  648,  658,  1000, 1314, 1707, 2540, 1560, 1964, 1666, 1628, 2394, 2562, 3085,
    3071, 4215, 4327, 2336, 4274, 3622, 3630, 4102, 3822, 3463, 5261, 5004, 3678, 3313, 3483,
    5308, 7595, 5685, 5954, 5539, 3624, 5776, 3376, 5240, 3165, 3730, 3927, 4148, 6311, 2863,
    2896, 2209, 3596, 4570, 6303, 3046, 2063, 5663, 3704, 2925, 2580, 1840, 3037, 3523, 2689,
    2079, 2583, 3364, 2552, 3363, 1705, 1659, 1499, 1709, 1159, 2939, 3654, 1318, 1011, 1257,
    1190, 1173, 1026, 1046, 1822, 1954, 1550, 1090, 707,  1414, 973,  719,  753,  620,  671,
    676,  712,  491,  484,  403,  409,  379,  384,  506,  268,  242,  219};

  for (G4int i = 0; i < 102; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_5[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[4]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_5[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_5[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[4] << " radius: " << radial_distance_5[i]
                   << " Data: " << radial_fluence_5[i] << G4endl;
  }

  G4double radial_distance_6[99] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3,
    -97.5,  -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,
    -64.1,  -60.5,  -58.6,  -54.6,  -54.6,  -50.3,  -50.3,  -45.6,  -40.4,  -27.0,  40.4,
    45.6,   46.0,   50.3,   54.6,   58.6,   60.5,   60.5,   60.5,   67.5,   67.5,   69.1,
    70.8,   70.8,   70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   90.3,   90.3,   90.3,
    90.3,   92.8,   94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,
    112.5,  117.4,  117.4,  117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,
    127.5,  131.0,  133.5,  137.7,  143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,
    153.9,  154.6,  154.6,  154.6,  164.5,  164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_6[99] = {
    1257,  1637,  2095,  2100,  4818,  5894,  6242,  7688,  7037,  7789,  8061,  8519,  10821,
    10558, 11620, 10888, 12217, 12334, 9923,  14786, 15212, 13871, 13237, 15221, 17121, 16296,
    20293, 17330, 17794, 16357, 21007, 13977, 20184, 18474, 15769, 17059, 14299, 17890, 16155,
    15908, 15019, 14015, 13224, 15993, 11361, 13493, 11989, 11636, 11672, 11225, 11806, 10268,
    10546, 9553,  10405, 9890,  10890, 9383,  8936,  7069,  7980,  6998,  7709,  6816,  6772,
    6210,  5635,  5572,  5549,  6473,  4862,  5114,  5263,  5297,  4539,  4624,  3854,  4039,
    3982,  4465,  3469,  2756,  3480,  3162,  2592,  2932,  1676,  1897,  2242,  1907,  1761,
    1865,  1465,  1390,  1646,  1446,  1215,  687,   800};

  G4double radial_error_6[99] = {
    232,  101,  181,  180,  261,  351,  564,  951,  439,  545,  442,  438,  676,  712,  901,
    845,  1184, 1316, 631,  1543, 1217, 1084, 1239, 1191, 1055, 2031, 1740, 1135, 1034, 978,
    1554, 3914, 2036, 1227, 1721, 1049, 1467, 1071, 1277, 1284, 1125, 2598, 845,  944,  608,
    1001, 1244, 2673, 874,  600,  2153, 1143, 849,  510,  903,  1090, 838,  613,  733,  1136,
    774,  1132, 474,  451,  411,  471,  298,  1058, 1456, 332,  256,  327,  324,  303,  252,
    270,  472,  582,  453,  298,  181,  428,  266,  185,  202,  160,  185,  197,  190,  138,
    132,  108,  117,  96,   97,   149,  73,   63,   57};

  for (G4int i = 0; i < 99; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_6[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[5]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_6[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_6[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[5] << " radius: " << radial_distance_6[i]
                   << " Data: " << radial_fluence_6[i] << G4endl;
  }

  G4double radial_distance_7[92] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -70.8,  -67.5,  -67.5,  -64.1,  -64.1,  -60.5,
    -54.6,  -50.3,  -50.3,  45.6,   46.0,   50.3,   58.6,   60.5,   60.5,   60.5,   67.5,   69.1,
    70.8,   70.8,   70.8,   76.9,   76.9,   79.7,   81.1,   83.9,   85.2,   90.3,   90.3,   90.3,
    90.3,   92.8,   94.0,   95.2,   102.0,  102.0,  105.3,  105.3,  105.3,  111.5,  111.5,  112.5,
    117.4,  117.4,  117.4,  120.2,  124.8,  124.8,  124.8,  127.5,  127.5,  127.5,  127.5,  131.0,
    133.5,  137.7,  143.3,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,  154.6,
    154.6,  164.5,  164.5,  164.5,  165.2,  171.2,  187.5,  187.5};

  G4double radial_fluence_7[92] = {
    262,  372,  493,  464,  1136, 1377, 1509, 1745, 1628, 1760, 1988, 2072, 2660, 2587, 2955, 2468,
    2729, 3136, 2108, 3912, 3773, 3403, 3158, 3810, 4100, 4597, 4773, 5030, 4336, 4563, 5307, 4229,
    4197, 3683, 4033, 4644, 2733, 3368, 2742, 2802, 2795, 2879, 2976, 2387, 2653, 2803, 2422, 2535,
    2330, 2609, 2341, 2058, 1500, 1864, 1501, 1822, 1526, 1589, 1312, 1267, 1144, 1123, 1532, 1097,
    1153, 1124, 1140, 978,  1077, 793,  889,  902,  996,  783,  542,  750,  703,  553,  649,  351,
    420,  460,  427,  381,  427,  318,  296,  351,  297,  271,  141,  170};

  G4double radial_error_7[92] = {
    40,  31,  57,  57,  71,  92,  149, 266, 114, 149, 126, 115, 196, 219, 277, 268, 305, 388, 685,
    532, 381, 307, 346, 307, 658, 330, 347, 390, 540, 354, 344, 393, 360, 328, 862, 300, 183, 315,
    386, 784, 265, 172, 596, 292, 255, 310, 137, 264, 338, 239, 179, 203, 322, 213, 378, 125, 134,
    107, 142, 75,  309, 377, 82,  66,  90,  80,  76,  67,  82,  136, 203, 124, 83,  45,  139, 74,
    49,  49,  41,  55,  59,  50,  36,  33,  28,  37,  28,  24,  38,  18,  17,  15};

  for (G4int i = 0; i < 92; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_7[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[6]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_7[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_7[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[6] << " radius: " << radial_distance_7[i]
                   << " Data: " << radial_fluence_7[i] << G4endl;
  }

  G4double radial_distance_8[80] = {
    -187.5, -164.5, -154.6, -150.2, -132.7, -124.8, -121.2, -111.5, -111.5, -107.4, -105.3, -97.5,
    -95.2,  -94.0,  -90.3,  -90.3,  -76.9,  -70.8,  -67.5,  -67.5,  -64.1,  -54.6,  -50.3,  -50.3,
    -45.6,  45.6,   46.0,   50.3,   58.6,   60.5,   67.5,   69.1,   70.8,   70.8,   70.8,   76.9,
    76.9,   79.7,   81.1,   85.2,   90.3,   90.3,   90.3,   92.8,   94.0,   95.2,   102.0,  105.3,
    105.3,  111.5,  111.5,  112.5,  117.4,  117.4,  117.4,  120.2,  124.8,  124.8,  124.8,  127.5,
    127.5,  131.0,  133.5,  137.7,  143.3,  144.1,  144.1,  145.6,  150.2,  150.2,  153.9,  154.6,
    154.6,  154.6,  164.5,  164.5,  164.5,  165.2,  171.2,  187.5};

  G4double radial_fluence_8[80] = {
    122,  182,  247,  223,  571,  687,  769,  853,  808,  860,  1027, 1061, 1369, 1329, 1554, 1205,
    1323, 1647, 2115, 1952, 1749, 2137, 2484, 2594, 1997, 2761, 2392, 2476, 3072, 2266, 2299, 2653,
    1389, 1751, 1348, 1436, 1416, 1524, 1557, 1342, 1517, 1297, 1170, 1322, 1217, 1018, 930,  916,
    741,  795,  615,  619,  525,  511,  771,  537,  563,  530,  540,  536,  427,  442,  483,  383,
    356,  341,  261,  314,  163,  203,  212,  208,  181,  211,  152,  139,  166,  137,  132,  80};

  G4double radial_error_8[80] = {
    24,  19,  31,  33,  46,  63,  89,  163, 72,  90,  83,  68,  111, 128, 164, 161,
    192, 250, 311, 203, 182, 405, 194, 214, 192, 225, 336, 181, 198, 232, 461, 189,
    110, 185, 198, 475, 149, 101, 367, 156, 194, 172, 200, 148, 124, 127, 124, 85,
    88,  73,  82,  46,  162, 214, 51,  41,  58,  50,  46,  51,  102, 74,  53,  27,
    43,  31,  30,  27,  32,  36,  34,  24,  20,  18,  21,  17,  17,  25,  12,  10};

  for (G4int i = 0; i < 80; ++i)
  {
    analysisManager->FillNtupleDColumn(9, 0, radial_distance_8[i]);
    analysisManager->FillNtupleDColumn(9, 1, radial_mean_energy[7]);
    analysisManager->FillNtupleDColumn(9, 2, radial_fluence_8[i]);
    analysisManager->FillNtupleDColumn(9, 3, radial_error_8[i]);
    analysisManager->AddNtupleRow(9);

    fLowRadialFile << " Energy: " << radial_mean_energy[7] << " radius: " << radial_distance_8[i]
                   << " Data: " << radial_fluence_8[i] << G4endl;
  }

  G4double he3_radial_distance_2[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_2[10] = {665235, 593298, 602514, 387449, 420921,
                                       331789, 275497, 218988, 207045, 115093};

  G4double he3_radial_error_2[10] = {9979, 9474, 9483, 7658, 7966, 7105, 6362, 5755, 5607, 4100};

  std::ofstream fDataRadialFile;
  G4String name2("RadialData.dat");
  fDataRadialFile.open(name2);

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_2[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[1]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_2[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_2[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[1]
                    << " radius: " << he3_radial_distance_2[i]
                    << " Data: " << he3_radial_fluence_2[i] << G4endl;
  }

  G4double he3_radial_distance_3[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_3[10] = {284113, 245043, 266956, 168293, 180439,
                                       130339, 119252, 93661,  88648,  42240};

  G4double he3_radial_error_3[10] = {4660, 4334, 4501, 3576, 3725, 3166, 3013, 2677, 2604, 1747};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_3[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[2]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_3[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_3[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[2]
                    << " radius: " << he3_radial_distance_3[i]
                    << " Data: " << he3_radial_fluence_3[i] << G4endl;
  }

  G4double he3_radial_distance_4[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_4[10] = {167763, 143894, 157657, 96559, 107418,
                                       75270,  67229,  54325,  49735, 24108};

  G4double he3_radial_error_4[10] = {3135, 2910, 3064, 2392, 2504, 2105, 2005, 1795, 1721, 1196};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_4[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[3]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_4[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_4[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[3]
                    << " radius: " << he3_radial_distance_4[i]
                    << " Data: " << he3_radial_fluence_4[i] << G4endl;
  }

  G4double he3_radial_distance_5[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_5[10] = {105971, 90614, 100489, 59547, 66724,
                                       47353,  41776, 32850,  30290, 14574};

  G4double he3_radial_error_5[10] = {2071, 1921, 2019, 1556, 1647, 1383, 1301, 1151, 1107, 755};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_5[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[4]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_5[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_5[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[4]
                    << " radius: " << he3_radial_distance_5[i]
                    << " Data: " << he3_radial_fluence_5[i] << G4endl;
  }

  G4double he3_radial_distance_6[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_6[10] = {26393, 21271, 24639, 14014, 15459,
                                       10163, 9156,  7627,  6354,  2970};

  G4double he3_radial_error_6[10] = {656, 588, 634, 478, 502, 407, 387, 353, 322, 220};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_6[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[5]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_6[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_6[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[5]
                    << " radius: " << he3_radial_distance_6[i]
                    << " Data: " << he3_radial_fluence_6[i] << G4endl;
  }

  G4double he3_radial_distance_7[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_7[10] = {6932, 5332, 6363, 3500, 3989, 2635, 2244, 1833, 1480, 755};

  G4double he3_radial_error_7[10] = {247, 217, 240, 180, 188, 154, 140, 125, 113, 79};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_7[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[6]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_7[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_7[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[6]
                    << " radius: " << he3_radial_distance_7[i]
                    << " Data: " << he3_radial_fluence_7[i] << G4endl;
  }

  G4double he3_radial_distance_8[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_8[10] = {3746, 2736, 3440, 1835, 2070, 1359, 1092, 938, 705, 316};

  G4double he3_radial_error_8[10] = {144, 124, 139, 102, 108, 87, 78, 72, 61, 43};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_8[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[7]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_8[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_8[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[7]
                    << " radius: " << he3_radial_distance_8[i]
                    << " Data: " << he3_radial_fluence_8[i] << G4endl;
  }

  G4double he3_radial_distance_9[10] = {16.8, 40.4,  45.6,  69.1,  81.1,
                                        98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_9[10] = {496, 399, 603, 313, 316, 200, 148, 128, 82, 47};

  G4double he3_radial_error_9[10] = {32, 28, 36, 25, 26, 21, 16, 16, 9, 9};

  for (G4int i = 0; i < 10; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_9[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[8]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_9[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_9[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[8]
                    << " radius: " << he3_radial_distance_9[i]
                    << " Data: " << he3_radial_fluence_9[i] << G4endl;
  }

  G4double he3_radial_distance_10[6] = {81.1, 98.6, 105.3, 113.5, 124.8, 153.9};

  G4double he3_radial_fluence_10[6] = {81, 61, 27, 30, 16, 9};

  G4double he3_radial_error_10[6] = {10, 8, 4, 6, 5, 3};

  for (G4int i = 0; i < 6; ++i)
  {
    analysisManager->FillNtupleDColumn(10, 0, he3_radial_distance_10[i]);
    analysisManager->FillNtupleDColumn(10, 1, radial_mean_energy[9]);
    analysisManager->FillNtupleDColumn(10, 2, he3_radial_fluence_10[i]);
    analysisManager->FillNtupleDColumn(10, 3, he3_radial_error_10[i]);
    analysisManager->AddNtupleRow(10);

    fDataRadialFile << " Energy: " << radial_mean_energy[9]
                    << " radius: " << he3_radial_distance_10[i]
                    << " Data: " << he3_radial_fluence_10[i] << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::createNeutronFluxHisto(G4int events, const Test15Run* tarcRun)
{
  G4cout << " fill flux histograms for events: " << events << G4endl;
  G4double absolute_totalflux = (tarcRun->total_flux * (1000000000.0 / (G4double)events))
                                / 26130.008;  // dividing by surface area of the shell: 4*pi*r^2
  G4double radius = 32.0 * mm;
  G4double cylinder_length = 150.0 * mm;
  G4double sphere_volume = (4.0 / 3.0) * 3.14159265 * std::pow(radius, 3.0);
  G4double cylinder_volume = 3.14159265 * std::pow(radius, 2.0) * cylinder_length;
  G4double shell_volume = (4.0 / 3.0) * 3.14159265 * std::pow(shell_outer_radius, 3.0)
                          - (4.0 / 3.0) * 3.14159265 * std::pow(shell_inner_radius, 3.0);

  G4double tarc_integral = 0.0;
  G4double tarc_integral_E = 0.0;
  G4double tarc_lithium = 0.0;
  G4double tarc_lithium_E = 0.0;
  G4double tarc_helium = 0.0;
  G4double tarc_helium_E = 0.0;

  std::ofstream fHighEnergyFile;
  G4String name4("HighEnergyOutput.dat");
  fHighEnergyFile.open(name4);
  std::ofstream fLowEnergyFile;
  G4String name5("LowEnergyOutput.dat");
  fLowEnergyFile.open(name5);
  std::ofstream fLithiumEnergyFile;
  G4String name6("LithiumEnergyOutput.dat");
  fLithiumEnergyFile.open(name6);

  fHighEnergyFile << std::setw(12) << "Energy" << std::setw(12) << "Geant4" << std::setw(12)
                  << "error" << std::setw(12) << "Data" << std::setw(12) << "stat" << std::setw(12)
                  << "syst" << std::setw(12) << "bin_min" << std::setw(12) << "bin_max" << G4endl;
  fLowEnergyFile << std::setw(12) << "Energy" << std::setw(12) << "Geant4" << std::setw(12)
                 << "error" << std::setw(12) << "Data" << std::setw(12) << "stat" << std::setw(12)
                 << "syst" << std::setw(12) << "bin_min" << std::setw(12) << "bin_max" << G4endl;
  fLithiumEnergyFile << std::setw(12) << "Energy" << std::setw(12) << "Geant4" << std::setw(12)
                     << "error" << std::setw(12) << "Data" << std::setw(12) << "stat"
                     << std::setw(12) << "syst" << std::setw(12) << "bin_min" << std::setw(12)
                     << "bin_max" << G4endl;

  for (G4int i = 0; i < 32; ++i)
  {
    G4double mean_energy = (tarcRun->flux_energy[i] + tarcRun->flux_energy[i + 1]) / 2.;
    // 26130.008 is the surface area in cm2 of a sphere at 45.6cm from the lead centre
    G4double absolute_flux = (tarcRun->flux[i] * (1000000000.0 / (G4double)events)) / 26130.008;
    // xbug    G4double absolute_flux_perp =
    // (cos_flux[i]*(1000000000.0/(G4double)events))/26130.008;
    G4double temp_bin_width = tarcRun->flux_energy[i + 1] - tarcRun->flux_energy[i];
    G4double absolute_flux_perp =
      mean_energy
      * (((tarcRun->cos_flux[i] * (1000000000.0 / (G4double)events)) / 26130.008) / temp_bin_width);

    G4double absolute_fluence =
      mean_energy
      * (100.0 * (1000000000.0 / (G4double)events) * ((tarcRun->fluence_step[i]) / sphere_volume)
         / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_fluence_front =
      mean_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->fluence_front_step[i]) / sphere_volume)
         / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2
    // xbug    G4double absolute_fluence_cyl =
    // 100.0*(1000000000.0/(G4double)events)*((fluence_step_cyl[i])/cylinder_volume); // factor of
    // 100 to account for /cm2 from /mm2
    G4double absolute_fluence_cyl =
      mean_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->fluence_step_cyl[i]) / cylinder_volume)
         / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_fluence_shell =
      mean_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->fluence_step_shell[i]) / shell_volume)
         / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_fluence_shell_error = 0.0;
    if (tarcRun->fluence_step_shell[i] > 0.0)
    {
      absolute_fluence_shell_error =
        (std::sqrt(tarcRun->fluence_step_shell[i]) / tarcRun->fluence_step_shell[i])
        * absolute_fluence_shell;
    }

    tarc_helium += tarcRun->flux_data[i];
    tarc_helium_E += tarcRun->flux_data[i] * mean_energy;

    G4cout << " flux[i]: " << tarcRun->flux[i] << " events: " << events << G4endl;
    G4cout << " cos_flux[i]: " << tarcRun->cos_flux[i] << " events: " << events << G4endl;
    G4cout << " fluence_step[i]: " << tarcRun->fluence_step[i]
           << " fluence_step_cyl[i]: " << tarcRun->fluence_step_cyl[i] << G4endl;
    G4double absolute_eflux = (tarcRun->eflux[i] * (1000000000.0 / (G4double)events)) / 26130.008;
    G4double absolute_error = 0.0;
    if (tarcRun->flux[i] > 0.0)
      absolute_error = (std::sqrt(tarcRun->flux[i]) / tarcRun->flux[i]) * absolute_flux;

    analysisManager->FillNtupleDColumn(3, 0, mean_energy);
    analysisManager->FillNtupleDColumn(3, 1, tarcRun->flux_data[i]);
    analysisManager->FillNtupleDColumn(3, 2, tarcRun->flux_stat_error[i]);
    analysisManager->FillNtupleDColumn(3, 3, tarcRun->flux_syst_error[i]);
    analysisManager->FillNtupleDColumn(3, 4, absolute_flux);
    analysisManager->FillNtupleDColumn(3, 5, absolute_flux_perp);
    analysisManager->FillNtupleDColumn(3, 6, absolute_fluence);
    analysisManager->FillNtupleDColumn(3, 7, absolute_error);
    analysisManager->FillNtupleDColumn(3, 8, tarcRun->flux[i]);
    analysisManager->FillNtupleDColumn(3, 9, tarcRun->eflux_data[i]);
    analysisManager->FillNtupleDColumn(3, 10, absolute_eflux);
    analysisManager->FillNtupleDColumn(3, 11, tarcRun->fluence_step[i]);
    analysisManager->FillNtupleDColumn(3, 12, absolute_fluence_cyl);
    analysisManager->FillNtupleDColumn(3, 13, absolute_fluence_front);
    analysisManager->FillNtupleDColumn(3, 14, absolute_fluence_shell);
    analysisManager->FillNtupleDColumn(3, 15, absolute_fluence_shell_error);
    analysisManager->AddNtupleRow(3);

    G4cout << " Got here i:" << i << " mean_energy: " << mean_energy
           << " absolute_flux: " << absolute_flux << G4endl;
    G4cout << G4endl << G4endl << G4endl;
    G4cout << " Got here i:" << i << " mean_energy: " << mean_energy
           << " absolute_fluence: " << absolute_fluence << G4endl;
    G4cout << " Got here i:" << i << " mean_energy: " << mean_energy
           << " absolute_fluence cylinder: " << absolute_fluence_cyl << G4endl;
    G4cout << " and flux: " << tarcRun->flux[i] << " events: " << events
           << " and flux data: " << tarcRun->flux_data[i]
           << " and cos flux:" << tarcRun->cos_flux[i] << G4endl;
    G4cout << G4endl << G4endl << G4endl;
    local_energy_integral[2] += absolute_flux * mean_energy;
    local_energy_integral[3] += tarcRun->flux_data[i] * mean_energy;

    fHighEnergyFile << std::setw(12) << mean_energy << std::setw(12) << absolute_fluence_shell
                    << std::setw(12) << absolute_fluence_shell_error << std::setw(12)
                    << tarcRun->flux_data[i] << std::setw(12) << tarcRun->flux_stat_error[i]
                    << std::setw(12) << tarcRun->flux_syst_error[i] << std::setw(12)
                    << tarcRun->flux_energy[i] << std::setw(12) << tarcRun->flux_energy[i + 1]
                    << G4endl;
  }

  for (G4int i = 0; i < 100; ++i)
  {
    G4double mean_low_energy = std::sqrt((tarcRun->low_energy[i + 1]) * (tarcRun->low_energy[i]));
    G4double absolute_low_flux =
      (tarcRun->low_flux[i] * (1000000000.0 / (G4double)events)) / 26130.008;
    // xbug    G4double absolute_low_flux_perp =
    // (cos_low_flux[i]*(1000000000.0/(G4double)events))/26130.008;
    G4double temp_bin_width_low = tarcRun->low_energy[i + 1] - tarcRun->low_energy[i];
    G4double absolute_low_flux_perp =
      mean_low_energy
      * (((tarcRun->cos_low_flux[i] * (1000000000.0 / (G4double)events)) / 26130.008)
         / temp_bin_width_low);
    G4double absolute_low_fluence =
      100.0 * (1000000000.0 / (G4double)events)
      * ((tarcRun->low_fluence_step[i])
         / sphere_volume);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_low_fluence_front =
      100.0 * (1000000000.0 / (G4double)events)
      * ((tarcRun->low_fluence_front_step[i])
         / sphere_volume);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_low_fluence_cyl =
      mean_low_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->low_fluence_step_cyl[i]) / cylinder_volume)
         / temp_bin_width_low);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_low_fluence_shell =
      mean_low_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->low_fluence_step_shell[i]) / shell_volume)
         / temp_bin_width_low);  // factor of 100 to account for /cm2 from /mm2

    G4double absolute_low_fluence_shell_error = 0.0;
    if (tarcRun->low_fluence_step_shell[i] > 0)
    {
      absolute_low_fluence_shell_error =
        (std::sqrt(tarcRun->low_fluence_step_shell[i]) / tarcRun->low_fluence_step_shell[i])
        * absolute_low_fluence_shell;
    }

    G4double absolute_low_error = 0;
    if (tarcRun->low_flux[i] > 0.0)
    {
      absolute_low_error =
        (std::sqrt(tarcRun->low_flux[i]) / tarcRun->low_flux[i]) * absolute_low_flux;
    }
    //

    tarc_integral += tarcRun->low_flux_data[i];
    tarc_integral_E += tarcRun->low_flux_data[i] * mean_low_energy;

    analysisManager->FillNtupleDColumn(4, 0, mean_low_energy);
    analysisManager->FillNtupleDColumn(4, 1, tarcRun->low_flux_data[i]);
    analysisManager->FillNtupleDColumn(4, 2, tarcRun->low_stat[i]);
    analysisManager->FillNtupleDColumn(4, 3, tarcRun->low_syst[i]);
    analysisManager->FillNtupleDColumn(4, 4, absolute_low_flux);
    analysisManager->FillNtupleDColumn(4, 5, absolute_low_flux_perp);
    analysisManager->FillNtupleDColumn(4, 6, absolute_low_fluence);
    analysisManager->FillNtupleDColumn(4, 7, absolute_low_error);
    analysisManager->FillNtupleDColumn(4, 8, tarcRun->low_flux[i]);
    analysisManager->FillNtupleDColumn(4, 9, tarcRun->low_fluence_step[i]);
    analysisManager->FillNtupleDColumn(4, 10, absolute_low_fluence_cyl);
    analysisManager->FillNtupleDColumn(4, 11, absolute_low_fluence_front);
    analysisManager->FillNtupleDColumn(4, 12, absolute_low_fluence_shell);
    analysisManager->FillNtupleDColumn(4, 13, absolute_low_fluence_shell_error);
    analysisManager->AddNtupleRow(4);

    G4cout << " Got here i:" << i << " mean_energy: " << mean_low_energy
           << "lowE1: " << tarcRun->low_energy[i] << " lowE2: " << tarcRun->low_energy[i + 1]
           << " absolute_flux: " << absolute_low_flux << G4endl;
    G4cout << " and flux: " << tarcRun->low_flux[i] << " events: " << events
           << " and flux data: " << tarcRun->low_flux_data[i] << G4endl;
    G4cout << " and fluence: " << tarcRun->low_fluence[i] << " events: " << events
           << " and flux data: " << tarcRun->low_flux_data[i] << G4endl;

    fLowEnergyFile << std::setw(12) << mean_low_energy << std::setw(12)
                   << absolute_low_fluence_shell << std::setw(12)
                   << absolute_low_fluence_shell_error << std::setw(12) << tarcRun->low_flux_data[i]
                   << std::setw(12) << tarcRun->low_stat[i] << std::setw(12) << tarcRun->low_syst[i]
                   << std::setw(12) << tarcRun->low_energy[i] << std::setw(12)
                   << tarcRun->low_energy[i + 1] << G4endl;
  }

  for (G4int i = 0; i < 100; ++i)
  {
    G4double mean_lithium_energy =
      std::sqrt((tarcRun->lithium_energy[i + 1]) * (tarcRun->lithium_energy[i]));
    G4double absolute_lithium_flux =
      (tarcRun->lithium_flux[i] * (1000000000.0 / (G4double)events)) / 26130.008;
    // xbug    G4double absolute_lithium_flux_perp =
    // (cos_lithium_flux[i]*(1000000000.0/(G4double)events))/26130.008;
    G4double temp_bin_width_lithium = tarcRun->low_energy[i + 1] - tarcRun->low_energy[i];
    G4double absolute_lithium_flux_perp =
      mean_lithium_energy
      * (((tarcRun->cos_lithium_flux[i] * (1000000000.0 / (G4double)events)) / 26130.008)
         / temp_bin_width_lithium);
    //     G4double absolute_lithium_fluence =
    //     (lithium_fluence[i]*(1000000000.0/(G4double)events))*((lithium_fluence_step[i]/cm)/sphere_volume/cm3);
    G4double absolute_lithium_fluence =
      100.0 * (1000000000.0 / (G4double)events)
      * ((tarcRun->lithium_fluence_step[i])
         / sphere_volume);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_lithium_fluence_front =
      100.0 * (1000000000.0 / (G4double)events)
      * ((tarcRun->lithium_fluence_front_step[i])
         / sphere_volume);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_lithium_fluence_cyl =
      mean_lithium_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->lithium_fluence_step_cyl[i]) / cylinder_volume)
         / temp_bin_width_lithium);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_lithium_fluence_shell =
      mean_lithium_energy
      * (100.0 * (1000000000.0 / (G4double)events)
         * ((tarcRun->lithium_fluence_step_shell[i]) / shell_volume)
         / temp_bin_width_lithium);  // factor of 100 to account for /cm2 from /mm2
    G4double absolute_lithium_fluence_shell_error = 0.0;
    if (tarcRun->lithium_fluence_step_shell[i] > 0)
    {
      absolute_lithium_fluence_shell_error =
        (std::sqrt(tarcRun->lithium_fluence_step_shell[i]) / tarcRun->lithium_fluence_step_shell[i])
        * absolute_lithium_fluence_shell;
    }

    G4cout << G4endl << G4endl << G4endl
           << " lithium_fluence_step_shell = " << tarcRun->lithium_fluence_step_shell[i]
           << " for i: " << i << " mean energy: " << mean_lithium_energy << G4endl << G4endl
           << G4endl;

    G4double absolute_lithium_Zflux =
      (tarcRun->lithium_Zflux[i] * (1000000000.0 / (G4double)events)) / 859.539;
    G4double absolute_lithium_flux_5cm =
      (tarcRun->lithium_flux_5cm[i] * (1000000000.0 / (G4double)events)) / 12.57;
    G4double absolute_lithium_error = 0.0;
    if (tarcRun->lithium_flux[i] > 0.0)
    {
      absolute_lithium_error =
        (std::sqrt(tarcRun->lithium_flux[i]) / tarcRun->lithium_flux[i]) * absolute_lithium_flux;
    }
    G4double absolute_lithium_error_5cm = 0.0;
    if (tarcRun->lithium_flux_5cm[i] > 0.0)
    {
      absolute_lithium_error_5cm =
        (std::sqrt(tarcRun->lithium_flux_5cm[i]) / tarcRun->lithium_flux_5cm[i])
        * absolute_lithium_flux_5cm;
    }
    //

    tarc_lithium += tarcRun->low_flux_data[i];
    tarc_lithium_E += tarcRun->low_flux_data[i] * mean_lithium_energy;

    analysisManager->FillNtupleDColumn(5, 0, mean_lithium_energy);
    analysisManager->FillNtupleDColumn(5, 1, tarcRun->lithium_flux_data[i]);
    analysisManager->FillNtupleDColumn(5, 2, tarcRun->lithium_stat[i]);
    analysisManager->FillNtupleDColumn(5, 3, tarcRun->lithium_syst[i]);
    analysisManager->FillNtupleDColumn(5, 4, absolute_lithium_flux);
    analysisManager->FillNtupleDColumn(5, 5, absolute_lithium_flux_perp);
    analysisManager->FillNtupleDColumn(5, 6, absolute_lithium_fluence);
    analysisManager->FillNtupleDColumn(5, 7, absolute_lithium_Zflux);
    analysisManager->FillNtupleDColumn(5, 8, absolute_lithium_error);
    analysisManager->FillNtupleDColumn(5, 9, absolute_lithium_flux_5cm);
    analysisManager->FillNtupleDColumn(5, 10, absolute_lithium_error_5cm);
    analysisManager->FillNtupleDColumn(5, 11, tarcRun->lithium_flux[i]);
    analysisManager->FillNtupleDColumn(5, 12, tarcRun->lithium_fluence_step[i]);
    analysisManager->FillNtupleDColumn(5, 13, absolute_lithium_fluence_cyl);
    analysisManager->FillNtupleDColumn(5, 14, absolute_lithium_fluence_front);
    analysisManager->FillNtupleDColumn(5, 15, absolute_lithium_fluence_shell);
    analysisManager->FillNtupleDColumn(5, 16, absolute_lithium_fluence_shell_error);
    analysisManager->AddNtupleRow(5);

    G4cout << " LITHIUM Got here i:" << i << " lithium mean_energy: " << mean_lithium_energy
           << "lowE1: " << tarcRun->lithium_energy[i]
           << " lowE2: " << tarcRun->lithium_energy[i + 1]
           << " absolute_flux: " << absolute_lithium_flux << G4endl;
    G4cout << " and flux: " << tarcRun->lithium_flux[i] << " events: " << events
           << " and flux data: " << tarcRun->lithium_flux_data[i] << G4endl;
    G4cout << " and fluence: " << tarcRun->lithium_fluence[i] << " events: " << events
           << " and flux data: " << tarcRun->lithium_flux_data[i] << G4endl;

    fLithiumEnergyFile << std::setw(12) << mean_lithium_energy << std::setw(12)
                       << absolute_lithium_fluence_shell << std::setw(12)
                       << absolute_lithium_fluence_shell_error << std::setw(12)
                       << tarcRun->lithium_flux_data[i] << std::setw(12) << tarcRun->lithium_stat[i]
                       << std::setw(12) << tarcRun->lithium_syst[i] << std::setw(12)
                       << tarcRun->lithium_energy[i] << std::setw(12)
                       << tarcRun->lithium_energy[i + 1] << G4endl;
  }

  G4cout << G4endl << G4endl << G4endl << " total_flux: " << tarcRun->total_flux << G4endl
         << " absolute total_flux: " << absolute_totalflux << G4endl << G4endl << G4endl << G4endl;

  G4double integral_temp_radius = 45.6 * cm;
  G4double integral_surface_area = 4.0 * CLHEP::pi * integral_temp_radius * integral_temp_radius;
  G4double absolute_integral_scintillation =
    (tarcRun->integral_scintillation * (1.e9 / (G4double)events)) / (integral_surface_area / cm2);
  G4double absolute_integral_scintillation_E =
    (tarcRun->integral_scintillation_E * (1.e9 / (G4double)events)) / (integral_surface_area / cm2);

  G4double absolute_integral_lithium =
    (tarcRun->integral_lithium * (1000000000.0 / (G4double)events)) / 26130.008;
  G4double absolute_integral_lithium_E =
    (tarcRun->integral_lithium_E * (1000000000.0 / (G4double)events)) / 26130.008;
  G4double absolute_integral_helium =
    (tarcRun->integral_helium * (1000000000.0 / (G4double)events)) / 26130.008;
  G4double absolute_integral_helium_E =
    (tarcRun->integral_helium_E * (1000000000.0 / (G4double)events)) / 26130.008;

  if (absolute_integral_scintillation != 0)
    G4cout << G4endl << G4endl << G4endl << G4endl << G4endl
           << " Low Flux, Geant4: " << absolute_integral_scintillation << " data: " << tarc_integral
           << " ratio: " << tarc_integral / absolute_integral_scintillation
           << " surface area: " << G4BestUnit(integral_surface_area, "Surface") << G4endl
           << " Lithium Flux, Geant4: " << absolute_integral_lithium << " data: " << tarc_lithium
           << G4endl << " Helium Flux, Geant4: " << absolute_integral_helium
           << " data: " << tarc_helium << G4endl << G4endl << G4endl << G4endl
           << " ENERGY Low Flux, Geant4: " << absolute_integral_scintillation_E
           << " data: " << tarc_integral_E << G4endl
           << " ENERGY Lithium Flux, Geant4: " << absolute_integral_lithium_E
           << " data: " << tarc_lithium_E << G4endl
           << " ENERGY Helium Flux, Geant4: " << absolute_integral_helium_E
           << " data: " << tarc_helium_E << G4endl << G4endl << G4endl << G4endl;

  G4double surface_area;
  G4double temp_radius;

  for (G4int i = 0; i < 10; ++i)
  {
    for (G4int j = 0; j < 10; ++j)
    {
      temp_radius = tarcRun->radii[i];
      surface_area = 4.0 * CLHEP::pi * temp_radius * temp_radius;
      G4cout << " temp_radius: " << G4BestUnit(temp_radius, "Length")
             << " surface_area: " << G4BestUnit(surface_area, "Surface") << G4endl;
      G4double temp_denominator = tarcRun->fractional_bin_width * tarcRun->radii_energies[j] / eV;
      G4cout << " flux_radius[i][j]: " << tarcRun->flux_radius[i][j] << " events: " << events
             << " bin width in eV: " << temp_denominator << G4endl;

      analysisManager->FillNtupleDColumn(13, 0, tarcRun->radii[i]);
      analysisManager->FillNtupleDColumn(13, 1, tarcRun->radii_energies[j] / eV);
      analysisManager->FillNtupleDColumn(13, 2, copy_radial_fluence_1[i]);
      analysisManager->FillNtupleDColumn(13, 3, copy_radial_error_1[i]);
      analysisManager->AddNtupleRow(13);
    }
  }

  for (G4int i = 0; i < tarcRun->n_max; ++i)
  {
    temp_radius = 45.6 * cm;
    surface_area = 4.0 * CLHEP::pi * temp_radius * temp_radius;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::createRadialFluxHisto(G4int events, const Test15Run* aRun)
{
  //--------------------------------------------------
  // nRadialShellFluence = tpf->create("4999", "Radial Shell Fluence","float radius, energy,
  // fluence, true_e, true_f");

  const Test15Run* tarcRun = static_cast<const Test15Run*>(aRun);

  // loop over radii:
  std::ofstream fGeant4RadialFile;
  G4String name3("Geant4RadialData.dat");
  fGeant4RadialFile.open(name3);

  for (G4int i = 0; i < 26; ++i)
  {
    G4double shell_volume = (4.0 / 3.0) * 3.14159265 * std::pow(outer_radius[i], 3.0)
                            - (4.0 / 3.0) * 3.14159265 * std::pow(inner_radius[i], 3.0);
    G4double radius = 0.5 * (outer_radius[i] + inner_radius[i]);
    // loop over energies
    for (G4int j = 0; j < 10; ++j)
    {
      G4double temp_bin_width =
        tarcRun->lithium_radial_energy_upper[j] - tarcRun->lithium_radial_energy_lower[j];
      if (tarcRun->radial_fluence_step[i][j] != 0)
      {
        G4double absolute_radial_fluence =
          tarcRun->lithium_radial_mean[j]
          * (100.0 * (1000000000.0 / (G4double)events)
             * ((tarcRun->radial_fluence_step[i][j]) / shell_volume)
             / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2
        G4double absolute_radial_fluence_true =
          tarcRun->lithium_radial_true_mean[j]
          * (100.0 * (1000000000.0 / (G4double)events)
             * ((tarcRun->radial_fluence_step[i][j]) / shell_volume)
             / temp_bin_width);  // factor of 100 to account for /cm2 from /mm2

        G4cout << " radial_fluence_step = " << tarcRun->radial_fluence_step[i][j] << " for i: " << i
               << " and j: " << j << " mean energy: " << tarcRun->lithium_radial_mean[j]
               << " outer radius: " << outer_radius[i] << G4endl;

        analysisManager->FillNtupleDColumn(8, 0, radius);
        analysisManager->FillNtupleDColumn(8, 1, tarcRun->lithium_radial_mean[j]);
        analysisManager->FillNtupleDColumn(8, 2, absolute_radial_fluence);
        analysisManager->FillNtupleDColumn(
          8, 3, tarcRun->lithium_radial_true_mean[j]);  // true_e lithium_radial_true_mean[j]
        analysisManager->FillNtupleDColumn(8, 4, absolute_radial_fluence_true);  // true_f
        analysisManager->AddNtupleRow(8);

        fGeant4RadialFile << " Energy: " << tarcRun->lithium_radial_mean[j]
                          << " radius: " << radius / 10. << " Geant4: "
                          << absolute_radial_fluence / tarcRun->lithium_radial_mean[j] << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15RunAction::secondarySummary(G4int, const Test15Run*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
