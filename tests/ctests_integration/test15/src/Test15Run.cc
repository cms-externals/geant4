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
/// \file biasing/Test15/src/Test15Run.cc
/// \brief Implementation of the Test15Run class
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  (Description)
//    Test15Run Class is for accumulating scored quantities which is 
//  scored using G4MutiFunctionalDetector and G4VPrimitiveScorer.
//  Accumulation is done using G4THitsMap object.
//
//    The constructor Test15Run(const std::vector<G4String> mfdName)
//  needs a vector filled with MultiFunctionalDetector names which
//  was assigned at instantiation of MultiFunctionalDetector(MFD).
//  Then Test15Run constructor automatically scans primitive scorers
//  in the MFD, and obtains collectionIDs of all collections associated
//  to those primitive scorers. Futhermore, the G4THitsMap objects 
//  for accumulating during a RUN are automatically created too.
//  (*) Collection Name is same as primitive scorer name.
// 
//    The resultant information is kept inside Test15Run objects as
//  data members.
//  std::vector<G4String> fCollName;            // Collection Name,
//  std::vector<G4int> fCollID;                 // Collection ID,
//  std::vector<G4THitsMap<G4double>*> fRunMap; // HitsMap for RUN.
//
//  The resualtant HitsMap objects are obtain using access method,
//  GetHitsMap(..).
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Test15Run.hh"
#include "G4SDManager.hh"

#include "G4RunManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4SystemOfUnits.hh"

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Test15Run::Test15Run(): G4Run(),oldTrackID(-1)
{

  fractional_bin_width = 0.2;
  
  total_flux = 0;
  G4double flux_energy_in[] = {59500.,109000.,158500.,208000.,257500.,307000.,356500.,406000.,455500.,505000.,554500.,604000.,653500.,703000.,752500.,802000.,901000.,1000000.,1162308.,1350960.,1570232.,1825092.,2500000.,3500000.,5000000.,7000000.,10000000.,20000000.,50000000.,100000000.,500000000.,1000000000.,2000000000.};

  G4double flux_data_in[] = {4769990.,5513695.,7339415.,6917651.,10017312.,6167753.,9621671.,8704919.,8429564.,10091790.,7095639.,5036533.,5692027.,12173780.,2858499.,2481939.,2492927.,2233198.,2491258.,959095.,1814799.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  G4double flux_stat_in[] = {585185.,766895.,1079366.,1154588.,1627418.,1283921.,1816459.,1812074.,1884908.,2221605.,1910559.,1670429.,1887832.,3046648.,1457554.,1043202.,1142402.,938652.,1047119.,684930.,1063375.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  G4double flux_syst_in[] = {858598.,992465.,1321095.,1176001.,1702943.,1048518.,1635684.,1392787.,1348730.,1614686.,1135302.,805845.,910724.,1947805.,457360.,397110.,398868.,357312.,398601.,153455.,290368.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  G4double low_energy_in[] = 
    {0.38,
     0.44,
     0.52,
     0.61,
     0.72,
     0.84,
     0.99,
     1.16,
     1.37,
     1.6,
     1.88,
     2.21,
     2.6,
     3.06,
     3.59,
     4.22,
     4.96,
     5.82,
     6.84,
     8.04,
     9.44,
     11.09,
     13.03,
     15.31,
     17.99,
     21.14,
     24.83,
     29.17,
     34.28,
     40.27,
     47.32,
     55.59,
     65.31,
     76.74,
     90.16,
     105.93,
     124.45,
     146.22,
     171.79,
     201.84,
     237.14,
     278.61,
     327.34,
     384.59,
     451.86,
     530.89,
     623.74,
     732.83,
     861,
     1011.58,
     1188.5,
     1396.37,
     1640.59,
     1927.53,
     2264.64,
     2660.73,
     3126.08,
     3672.82,
     4315.19,
     5069.9,
     5956.62,
     6998.41,
     8222.42,
     9660.5,
     11350.09,
     13335.2,
     15667.5,
     18407.72,
     21627.19,
     25409.74,
     29853.85,
     35075.23,
     41209.82,
     48417.34,
     56885.43,
     66834.59,
     78523.82,
     92257.17};

  G4double low_flux_in[] = 
    {275278,
     301782,
     340272,
     356588,
     394892,
     429498,
     458291,
     484007,
     506543,
     544825,
     566643,
     616788,
     642024,
     676087,
     719414,
     742831,
     774353,
     827903,
     845510,
     885449,
     915560,
     941733,
     980162,
     1012818,
     1055418,
     1102804,
     1141888,
     1156113,
     1185135,
     1235172,
     1289648,
     1307669,
     1341033,
     1370895,
     1398922,
     1431318,
     1441628,
     1513853,
     1547521,
     1563220,
     1636756,
     1628540,
     1705225,
     1751752,
     1797444,
     1858746,
     1913164,
     1944251,
     1953811,
     2064894,
     2070736,
     2206711,
     2244527,
     2341376,
     2431927,
     2498864,
     2645949,
     2678028,
     2758465,
     2878502,
     2882548,
     2890662,
     2994883,
     3266327,
     3300371,
     3237443,
     3488670,
     3721156,
     3758134,
     3793469,
     4196991,
     4477528,
     4762069,
     4984865,
     5477323,
     6116790,
     6491265,
     7726631};

  G4double low_stat_in[] = 
    {3448,
     3756,
     4157,
     4421,
     4838,
     5256,
     5646,
     6029,
     6422,
     6938,
     7324,
     7971,
     8464,
     9038,
     9673,
     10263,
     10859,
     11755,
     12289,
     13134,
     13871,
     14636,
     15567,
     16446,
     17466,
     18600,
     19708,
     20572,
     21690,
     23087,
     24600,
     25735,
     27092,
     28509,
     29863,
     31571,
     32989,
     35157,
     36989,
     38597,
     41363,
     42734,
     45742,
     48197,
     50849,
     53944,
     56952,
     59634,
     62321,
     66872,
     69862,
     75060,
     79162,
     84198,
     89569,
     94963,
     101803,
     106813,
     113268,
     120584,
     126414,
     132323,
     141359,
     154144,
     162315,
     169011,
     183639,
     198814,
     209985,
     221589,
     244748,
     266618,
     290044,
     312674,
     347075,
     389450,
     424105,
     490180};

  G4double low_syst_in[] = 
    {37226,
     40810,
     46016,
     48223,
     53404,
     58085,
     61981,
     65460,
     68511,
     73691,
     76646,
     83433,
     86852,
     91466,
     97336,
     100514,
     104792,
     112054,
     114456,
     119885,
     123990,
     127567,
     132814,
     137288,
     143123,
     149624,
     155017,
     157055,
     161127,
     168088,
     175696,
     178382,
     183211,
     187623,
     191857,
     196778,
     198759,
     209411,
     214898,
     218059,
     229518,
     229762,
     242287,
     250941,
     259929,
     271730,
     283196,
     291933,
     298177,
     320998,
     328686,
     358571,
     374394,
     402087,
     431292,
     459100,
     505230,
     533181,
     574485,
     629077,
     663089,
     701981,
     769922,
     891239,
     958096,
     1002119,
     1153770,
     1317240,
     1426233,
     1545651,
     1838329,
     2110672,
     2418254,
     2729327,
     3235861,
     3901621,
     4472942,
     5754467};

  G4double lithium_energy_in[] =
    {0.021,	 
     0.024,	 
     0.029,	 
     0.033,	 
     0.039,	 
     0.046,	 
     0.054,	 
     0.064,	 
     0.075,	 
     0.088,	 
     0.104,	 
     0.122,	 
     0.143,	 
     0.168,	 
     0.197,	 
     0.232,	 
     0.272,	 
     0.320,	 
     0.376,	 
     0.442,	 
     0.519,	 
     0.610,	 
     0.716,	 
     0.841,	 
     0.989,	 
     1.161,	 
     1.365,	 
     1.603,	 
     1.884,	 
     2.213,	 
     2.600,	 
     3.055,	 
     3.589,	 
     4.217,	 
     4.955,	 
     5.821,	 
     6.839,	 
     8.035,	 
     9.441,	 
     11.092,	 
     13.032,	 
     15.311,	 
     17.989,	 
     21.135,	 
     24.831,	 
     29.174,	 
     34.277,	 
     40.272,	 
     47.315,	 
     55.590,	 
     65.313,	 
     76.736,	 
     90.157,	 
     105.925,	 
     124.451,	 
     146.217,	 
     171.791,	 
     201.836,	 
     237.137,	 
     278.612,	 
     327.340,	 
     384.591,	 
     451.855,	 
     530.884,	 
     623.734,	 
     732.823,	 
     860.992,	 
     1011.577, 
     1188.499,
     1396.365, 
     1640.586, 
     1927.521, 
     2264.639, 
     2660.720, 
     3126.072, 
     3672.816, 
     4315.182, 
     5069.899, 
     5956.609, 
     6998.403, 
     8222.412, 
     9660.496, 
     11350.080,
     13335.170,
     15667.480,
     18407.660,
     21627.140,
     25409.660,
     29853.750,
     35075.100,
     41209.670,
     48417.120,
     56885.140,
     66834.250,
     78523.430};

  G4double lithium_flux_in[] =
    {15908,
     17993,
     24624,
     38977,
     38936,
     56744,
     51858,
     78687,
     80212,
     92435,
     112138,
     138309,
     147323,
     178602,
     208492,
     231684,
     242909,
     306021,
     310552,
     349617,
     361640,
     372367,
     390333,
     453242,
     502635,
     493583,
     494861,
     482488,
     514136,
     543474,
     627154,
     656274,
     726710,
     724241,
     806797,
     763901,
     737945,
     867885,
     875563,
     829345,
     921048,
     965185,
     1014804,
     889206,
     973931,
     1114639,
     1058564,
     1232387,
     1529422,
     1345060,
     1210322,
     1453971,
     1418608,
     1362524,
     1227946,
     1511154,
     1429373,
     1627807,
     1582482,
     1628250,
     1723348,
     1848327,
     2178810,
     2137029,
     2204020,
     1845805,
     1936230,
     1998776,
     2311320,
     2387485,
     2268294,
     2483071,
     2455236,
     2450272,
     2766197,
     2755937,
     3095936,
     3118377,
     3160060,
     3445465,
     3091801,
     3706670,
     3717733,
     4654344,
     3894362,
     4508537,
     2939137,
     4439832,
     4631860,
     5931119,
     9105124,
     14277100,
     45261009,
     45816255,
     5828108};
  
  G4double lithium_stat_in[] =
    {2296,
     2545,
     3102,
     4069,
     4242,
     5340,
     5209,
     6631,
     6933,
     7721,
     8835,
     10199,
     10941,
     12514,
     14033,
     15343,
     16284,
     13826,
     14432,
     15860,
     16667,
     17347,
     18013,
     18952,
     18568,
     15919,
     12888,
     10717,
     10963,
     13377,
     19953,
     25914,
     28390,
     29453,
     30985,
     25506,
     24175,
     32372,
     29274,
     26276,
     30171,
     33120,
     32924,
     28541,
     33712,
     41872,
     45386,
     60071,
     70447,
     59887,
     57434,
     67287,
     66888,
     63152,
     66509,
     81408,
     78724,
     86350,
     87150,
     91315,
     99985,
     113129,
     132774,
     130999,
     130072,
     121329,
     130514,
     140505,
     158478,
     165673,
     166773,
     180262,
     186485,
     192890,
     214243,
     221442,
     244535,
     255526,
     387514,
     412492,
     400721,
     455571,
     482881,
     552433,
     517210,
     568926,
     469298,
     589659,
     606093,
     697767,
     885261,
     1135819,
     2044653,
     2084523,
     749593};

  G4double lithium_syst_in[] =
    {2204, 
     2492, 
     3411, 
     5399, 
     5393, 
     7860, 
     7183, 
     10900, 
     11111, 
     12804, 
     15533, 
     19159, 
     20407, 
     24740, 
     28881, 
     32093, 
     33648, 
     29984, 
     30420, 
     34292, 
     35494, 
     36485, 
     38490, 
     44399, 
     49262, 
     48453, 
     48492, 
     47270, 
     50449, 
     53391, 
     61800, 
     64432, 
     71219, 
     70982, 
     79082, 
     74898, 
     72464, 
     85110, 
     85930, 
     81397, 
     90669, 
     94758, 
     99699, 
     87559, 
     95789, 
     109678, 
     104235, 
     121456, 
     150878, 
     132850, 
     119781, 
     144068, 
     140843, 
     135945, 
     122624, 
     151257, 
     143870, 
     166021, 
     161064, 
     167587, 
     177132, 
     191310, 
     228913, 
     227725, 
     235396, 
     199994, 
     212977, 
     224305, 
     264033, 
     280557, 
     272401, 
     307337, 
     312542, 
     323612, 
     378190, 
     394280, 
     462686, 
     491293, 
     733096, 
     843086, 
     800244, 
     1017497, 
     1085017, 
     1447471, 
     1293239, 
     1601695, 
     1118906, 
     1813919, 
     2033568, 
     2801549, 
     4631794, 
     7828772, 
     26773503, 
     29256399, 
     4019784};

  G4double radii_temp[] = {16.8*cm,40.4*cm,45.6*cm,69.1*cm,81.1*cm,98.6*cm,105.3*cm,113.5*cm,124.8*cm,153.9*cm};
  //  G4double radii_energies_temp[] = {0.1*eV,1.5*eV,5.0*eV,10.0*eV,18.0*eV,100.0*eV,480.0*eV,1000.0*eV,10000.0*eV,50000.0*eV};
  G4double radii_energies_temp[] = {0.1,1.5,5.0,10.0,18.0,100.0,480.0,1000.0,10000.0,50000.0};

  for(G4int i=0; i<10; ++i) {
    radii[i] = radii_temp[i];
    radii_energies[i] = radii_energies_temp[i];
    for(G4int j=0; j<10; ++j) flux_radius[i][j] = 0;
  }

  for(G4int i=0; i<1000; ++i) fluence_spectrum[i] = 0.0;
  n_max = 0;

  //     radii_energies
  //     flux_radius

  energy_integral[0] = 0.0;
  energy_integral[1] = 0.0;
  energy_integral[2] = 0.0;
  energy_integral[3] = 0.0;
  //  flux_energy = flux_energy_in;

  neutflux[0] = 0;
  neutflux[1] = 0;
  neutflux[2] = 0;
  neutflux[3] = 0;

  enflux[0] = 0.0;
  enflux[1] = 0.0;
  enflux[2] = 0.0;
  enflux[3] = 0.0;


  //  flux_data = flux_data_in;
  G4double mean_energy; 
  eflux_integral = 0.0;
  G4int j = 0;
  for (G4int i=0; i<32; ++i) {
    flux[i] = 0.0;
    cos_flux[i] = 0.0;
    fluence[i] = 0.0;
    fluence_step[i] = 0.0;
    fluence_front_step[i] = 0.0;
    fluence_cyl[i] = 0.0;
    fluence_step_cyl[i] = 0.0;
    fluence_step_shell[i] = 0.0;
    eflux[i] = 0.0;
    flux_energy[i] = flux_energy_in[i];
    mean_energy = (flux_energy[i] + flux_energy_in[i+1])/2.;
    //BUG!    flux_data[i] = flux_data_in[i]/mean_energy;
    flux_data[i] = flux_data_in[i];
    eflux_data[i] = flux_data_in[i]*mean_energy;
    eflux_integral += flux_data_in[i]*mean_energy/1.e6;
    fine_energy[j] = flux_energy_in[i];
    fine_eflux[j] = 0.0;
    j++;
    fine_energy[j] = mean_energy;
    fine_eflux[j] = 0.0;
    j++;
    flux_stat_error[i] = flux_stat_in[i];
    flux_syst_error[i] = flux_syst_in[i];
  }
  flux_energy[32] = flux_energy_in[32];
  fine_energy[j] = flux_energy_in[32];

  for(G4int k=0; k<65; ++k) G4cout << " Flux energy, I: " << k << " and energy: " << fine_energy[k] << G4endl;

  G4double bin_width = (std::log(1.e5)-std::log(0.01))/100.0;
  G4double energy_start = 0.01;

  low_energy[0] = energy_start;
  //  lithium_energy[0] = energy_start*eV;
  lithium_energy[0] = energy_start;

  G4int k_idx = 0;
  G4int m_idx = 0;

  lithium_integral_data = 0.0;
  lithium_Eintegral_data = 0.0;
  integral_flux_5cm = 0.0;
  integral_flux_10cm = 0.0;
  integral_flux_46cm = 0.0;
  integral_Zflux_46cm = 0.0;
  integral_Eflux_46cm = 0.0;
  integral_Eflux_46cm_restricted = 0.0;
  integral_flux_70cm = 0.0;
  integral_flux_100cm = 0.0;
  integral_flux_120cm = 0.0;

  //  G4double radial_mean_energy[] = {0.1*eV, 1.5*eV, 5*eV, 10*eV, 18*eV, 100*eV, 480*eV, 1000*eV, 10*keV, 50*keV};
  G4double radial_mean_energy[] = {0.1, 1.5, 5.0, 10.0, 18.0, 100.0, 480.0, 1000.0, 10000.0, 50000.0};
  G4int radial_index = 0;

  //  for (G4int i=0; i<78; ++i) {    
  for (G4int i=0; i<100; ++i) {    
    low_flux[i] = 0.0;
    cos_low_flux[i] = 0.0;
    low_fluence[i] = 0.0;
    low_fluence_step[i] = 0.0;
    low_fluence_front_step[i] = 0.0;
    low_fluence_cyl[i] = 0.0;
    low_fluence_step_cyl[i] = 0.0;
    low_fluence_step_shell[i] = 0.0;
    lithium_flux[i] = 0.0;
    cos_lithium_flux[i] = 0.0;
    lithium_fluence[i] = 0.0;
    lithium_fluence_step[i] = 0.0;
    lithium_fluence_front_step[i] = 0.0;
    lithium_fluence_cyl[i] = 0.0;
    lithium_fluence_step_cyl[i] = 0.0;
    lithium_fluence_step_shell[i] = 0.0;
    lithium_Zflux[i] = 0.0;
    lithium_flux_5cm[i] = 0.0;
    //    low_energy[i+1] = std::exp(2*log(low_energy_in[i])-log(low_energy[i]));
    low_energy[i+1] = std::exp(bin_width+std::log(low_energy[i]));
    //    lithium_energy[i+1] = std::exp(bin_width+std::log(lithium_energy[i]))*eV;
    lithium_energy[i+1] = std::exp(bin_width+std::log(lithium_energy[i]));
    mean_energy = std::exp(0.5*(std::log(low_energy[i+1])+std::log(low_energy[i])));
    //    G4double lithium_mean_energy = std::exp(0.5*(std::log(lithium_energy[i+1])+std::log(lithium_energy[i])))*eV;
    G4double lithium_mean_energy = std::exp(0.5*(std::log(lithium_energy[i+1])+std::log(lithium_energy[i])));

    G4cout << " lithium_energy: " << lithium_energy[i] << " radial_index: " << radial_index << " mean_energy: " << radial_mean_energy[radial_index] << G4endl;
    if(lithium_energy[i] < radial_mean_energy[radial_index] && lithium_energy[i+1] > radial_mean_energy[radial_index]) {
      lithium_radial_energy_lower[radial_index] = lithium_energy[i];
      lithium_radial_energy_upper[radial_index] = lithium_energy[i+1];
      lithium_radial_mean[radial_index] = radial_mean_energy[radial_index];
      lithium_radial_true_mean[radial_index] = lithium_mean_energy;
      for(G4int k=0; k<26; ++k) radial_fluence_step[k][radial_index] = 0.0;
      radial_index++;
      if(radial_index > 9) radial_index=9;
    }

    G4double bug_ratio = mean_energy/(low_energy[i+1]-low_energy[i]);
    G4cout << " bug ratio: " << bug_ratio << " mean: " << mean_energy << "  lowE1: " << low_energy[i] << " lowE2: " << low_energy[i+1] << " and I: " << i << G4endl;
    G4cout << " lithium mean: " << lithium_mean_energy << "  lithium lowE1: " << lithium_energy[i] << " lithium lowE2: " << lithium_energy[i+1] << " and I: " << i << G4endl;
    low_flux_data[i] = 0.0;
    low_stat[i] = 0.0;
    low_syst[i] = 0.0;
    lithium_flux_data[i] = 0.0;
    lithium_stat[i] = 0.0;
    lithium_syst[i] = 0.0;
    if(std::abs(mean_energy/low_energy_in[k_idx]-1.0) < 0.05) {
      low_flux_data[i] = low_flux_in[k_idx];
      low_stat[i] = low_stat_in[k_idx];
      low_syst[i] = low_syst_in[k_idx];
      k_idx++;
    }
    if(m_idx < 95 && (std::abs(lithium_mean_energy/lithium_energy_in[m_idx]-1.0) < 0.05)) {
      lithium_flux_data[i] = lithium_flux_in[m_idx];
      lithium_stat[i] = lithium_stat_in[m_idx];
      lithium_syst[i] = lithium_syst_in[m_idx];
      lithium_integral_data += lithium_flux_in[m_idx];
      lithium_Eintegral_data += lithium_flux_in[m_idx]*lithium_mean_energy/1.e6;
//       std::ofstream hitsfile("energyintegral.out", std::ios::app);
//       hitsfile << " Lithium Data Energy Flux: " << G4BestUnit((lithium_flux_in[m]*lithium_mean_energy/1.e6),"Energy") << G4endl;
      m_idx++;
    }
    //    low_flux_data[i] = low_flux_in[i];
    // low_stat[i] = low_stat_in[i];
    // low_syst[i] = low_syst_in[i];
  }

  gamma_flux = 0;
  exiting_flux = 0;
  exitinggrichine_flux = 0;
  exiting_check_flux = 0;
  neutron_fluence = 0;
  neutron_fluence_cyl = 0;
  neutron_fluence_46cm = 0;
  integral_flux_46cm = 0;
  integral_Eflux_46cm = 0.0;
  neutron_flux = 0;
  neutron_check = 0;
  electron_flux = 0;
  piminus_flux = 0;
  piplus_flux = 0;
  pizero_flux = 0;
  positron_flux = 0;
  proton_flux = 0;
  muon_flux = 0;
  other_flux = 0;

  exiting_energy = 0.0;

  integral_scintillation = 0;
  integral_scintillation_E = 0.0;
  integral_lithium = 0;
  integral_lithium_E = 0.0;
  integral_helium = 0;
  integral_helium_E = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Test15Run::~Test15Run()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in an event,
//  is accumulated during a Run.
void Test15Run::RecordEvent(const G4Event* aEvent)
{

  numberOfEvent++;  // This is an original line.

  //=============================
  // HitsCollection of This Event
  //============================
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  if (!HCE) return;

  //=======================================================
  // Sum up HitsMap of this Event  into HitsMap of this RUN
  //=======================================================
  G4int Ncol = fCollID.size();
  for ( G4int i = 0; i < Ncol ; i++ ){  // Loop over HitsCollection
    G4THitsMap<G4double>* EvtMap=0;
    if ( fCollID[i] >= 0 ){           // Collection is attached to HCE
      EvtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID[i]));
    }else{
      G4cout <<" Error EvtMap Not Found "<< i << G4endl;
    }
    if ( EvtMap )  {
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *fRunMap[i] += *EvtMap;
      //======================================================
    }
  }

  // G4Run::RecordEvent(aEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Access method for HitsMap of the RUN
//
//-----
// Access HitsMap.
//  By  MultiFunctionalDetector name and Collection Name.
G4THitsMap<G4double>* Test15Run::GetHitsMap(const G4String& detName,
                                         const G4String& colName){
    G4String fullName = detName+"/"+colName;
    return GetHitsMap(fullName);
}

//-----
// Access HitsMap.
//  By full description of collection name, that is
//    <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4THitsMap<G4double>* Test15Run::GetHitsMap(const G4String& fullName){
    G4int Ncol = fCollName.size();
    for ( G4int i = 0; i < Ncol; i++){
        if ( fCollName[i] == fullName ){
            return fRunMap[i];
        }
    }
    return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// - Dump All HitsMap of this RUN. (for debuging and monitoring of quantity).
//   This method calls G4THisMap::PrintAll() for individual HitsMap.
void Test15Run::DumpAllScorer(){

  // - Number of HitsMap in this RUN.
  G4int n = GetNumberOfHitsMap();
  // - GetHitsMap and dump values.
  for ( G4int i = 0; i < n ; i++ ){
    G4THitsMap<G4double>* RunMap =GetHitsMap(i);
    if ( RunMap ) {
      G4cout << " PrimitiveScorer RUN " 
             << RunMap->GetSDname() <<","<< RunMap->GetName() << G4endl;
      G4cout << " Number of entries " << RunMap->entries() << G4endl;
      std::map<G4int,G4double*>::iterator itr = RunMap->GetMap()->begin();
      for(; itr != RunMap->GetMap()->end(); itr++) {
        G4cout << "  copy no.: " << itr->first
               << "  Run Value : " << *(itr->second) 
               << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Test15Run::Merge(const G4Run* aRun)
{
  const Test15Run * localRun = static_cast<const Test15Run *>(aRun);
  //=======================================================
  // Merge HitsMap of working threads
  //=======================================================
  G4int nCol = localRun->fCollID.size();
  for ( G4int i = 0; i < nCol ; i++ ){  // Loop over HitsCollection
    if ( localRun->fCollID[i] >= 0 ){
      *fRunMap[i] += *localRun->fRunMap[i];
    }
  }
  exiting_flux += localRun->exiting_flux;
  exitinggrichine_flux += localRun->exitinggrichine_flux;
  exiting_energy += localRun->exiting_energy;
  exiting_check_flux += localRun->exiting_check_flux;
  gamma_flux += localRun->gamma_flux;
  neutron_flux += localRun->neutron_flux;
  electron_flux += localRun->electron_flux;
  piminus_flux += localRun->piminus_flux;
  piplus_flux += localRun->piplus_flux;
  pizero_flux += localRun->pizero_flux;
  positron_flux += localRun->positron_flux;
  proton_flux += localRun->proton_flux;
  muon_flux += localRun->muon_flux;
  other_flux += localRun->other_flux;
  neutron_check += localRun->neutron_check;
  neutron_fluence += localRun->neutron_fluence;
  integral_flux_46cm += localRun->integral_flux_46cm;
  integral_Eflux_46cm += localRun->integral_Eflux_46cm;

  total_flux += localRun->total_flux;

  for(G4int i=0; i<32; ++i) {
    flux[i] += localRun->flux[i];  
    cos_flux[i] += localRun->cos_flux[i];  
    fluence_step[i] += localRun->fluence_step[i];  
    fluence_front_step[i] += localRun->fluence_front_step[i];  
    fluence_step_cyl[i] += localRun->fluence_step_cyl[i];  
    fluence_step_shell[i] += localRun->fluence_step_shell[i];  
    eflux[i] += localRun->eflux[i];  
    // G4cout << " Test15Run:: fluence_step_shell: " << fluence_step_shell[i] << G4endl;
    // getchar();
  }
  for(G4int i=0; i<60; ++i) {
    fine_eflux[i] += localRun->fine_eflux[i];  
  }
  for(G4int i=0; i<100; ++i) {
    low_flux[i] += localRun->low_flux[i];  
    cos_low_flux[i] += localRun->cos_low_flux[i];  
    low_fluence_step[i] += localRun->low_fluence_step[i];  
    low_fluence_front_step[i] += localRun->low_fluence_front_step[i];  
    low_fluence_step_cyl[i] += localRun->low_fluence_step_cyl[i];  
    low_fluence_step_shell[i] += localRun->low_fluence_step_shell[i];  

    lithium_flux[i] += localRun->lithium_flux[i];  
    cos_lithium_flux[i] += localRun->cos_lithium_flux[i];  
    lithium_fluence_step[i] += localRun->lithium_fluence_step[i];  
    lithium_fluence_front_step[i] += localRun->lithium_fluence_front_step[i];  
    lithium_fluence_step_cyl[i] += localRun->lithium_fluence_step_cyl[i];  
    lithium_fluence_step_shell[i] += localRun->lithium_fluence_step_shell[i];  
    lithium_Zflux[i] += localRun->lithium_Zflux[i];  
    lithium_flux_5cm[i] += localRun->lithium_flux_5cm[i];  

  }

  for(G4int i=0; i<10; ++i) {
    for(G4int k=0; k<26; ++k) radial_fluence_step[k][i] += localRun->radial_fluence_step[k][i];
  }


  G4Run::Merge(aRun);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//x2018 void Test15Run::analyseNeutronRadialFluence(G4double rad_energy, G4double rad_time, G4double rad_steplength, G4int radius_index)
void Test15Run::analyseNeutronRadialFluence(G4double rad_energy, G4double, G4double rad_steplength, G4int radius_index)
{
  if(radius_index < 0 || radius_index > 26) G4cout << " WARNING radius index is wrong!!!!!! " << radius_index << G4endl;

  G4double temp_energy = rad_energy/eV;

  if(temp_energy < lithium_radial_energy_upper[9] && temp_energy > lithium_radial_energy_lower[0]) {
    for (G4int i=0; i<10; ++i) {    
      if(temp_energy > lithium_radial_energy_lower[i] && temp_energy < lithium_radial_energy_upper[i]) {
	radial_fluence_step[radius_index][i] += rad_steplength/mm;      
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//x2018 void Test15Run::analyseNeutronFlux(G4double n_energy, G4double n_time, G4double startEnergy, G4int TrackID, G4int ParentID, G4double zMomentum, G4double startTime, G4double radius, G4double zPos, G4double parentEnergy, G4String parentParticle, G4double cos_angle, G4int number_generations, G4String Particle, G4bool reduced_tally)
void Test15Run::analyseNeutronFlux(G4double n_energy, G4double , G4double, G4int TrackID, G4int, G4double, G4double, G4double radius, G4double zPos, G4double, G4String , G4double cos_angle, G4int , G4String Particle, G4bool )
{

  if(Particle == "neutron") {
    if(TrackID == oldTrackID && std::abs(radius - 456.0*mm)<0.1) {
      duplicate_neutron++;
    }  else {
      duplicate_neutron = 0;
    }
    oldTrackID = TrackID;
  }
    
  G4double temp_energy = n_energy/eV;

  if(Particle == "neutron") {

    for (G4int i=0; i<10; ++i) {
      if( std::abs(radius-radii[i])<0.1) {
	for (G4int j=0; j<10; ++j) {
	  if ( std::abs(n_energy-radii_energies[j])<fractional_bin_width*radii_energies[j]) flux_radius[i][j] += 1.0/std::abs(cos_angle);
	}
      }
    }

    if( std::abs(radius-45.6*cm)<0.1) {
      if(n_energy > 0.345*eV && n_energy < 1.e5*eV) {
	integral_scintillation += 1.0/std::abs(cos_angle);
	integral_scintillation_E += n_energy/eV*1.0/std::abs(cos_angle);
      }
      
      G4int n = (G4int) ((2.0 + std::log10(n_energy/eV))/0.09);
      //BUG BUG BUG BUG BUG!!!!
      if(n<0) n = 0;
      if(n < 100) {
	fluence_spectrum[n] += 1.0/std::abs(cos_angle);
	if(n > n_max) n_max = n;
      }
      
      if(temp_energy > 0.0194 && temp_energy < 1.e5) {
	integral_lithium += 1.0/std::abs(cos_angle);
	integral_lithium_E += temp_energy*1.0/std::abs(cos_angle);
      }
      
      if(temp_energy > 59500 && temp_energy < 1825092) {
	integral_helium += 1.0/std::abs(cos_angle);
	integral_helium_E += temp_energy*1.0/std::abs(cos_angle);
      }
    }
  }

  if(Particle == "neutron") {

    if(std::abs(radius-50.0*mm) < 0.1) integral_flux_5cm++;
    if(std::abs(radius-100.0*mm) < 0.1) integral_flux_10cm++;
    if(std::abs(radius-456.0*mm) < 0.1) {
      integral_flux_46cm++;
      integral_Eflux_46cm += n_energy;
      if(n_energy > 0.1*eV && n_energy < 10.0*keV) neutron_fluence_46cm++;
      if(std::abs(zPos-75.0) < 15.0) integral_Zflux_46cm++;
    }
    if(std::abs(radius-700.0*mm) < 0.1) integral_flux_70cm++;
    if(std::abs(radius-1000.0*mm) < 0.1) integral_flux_100cm++;
    if(std::abs(radius-1200.0*mm) < 0.1) integral_flux_120cm++;
    
    if(std::abs(radius-456.0*mm) < 0.1) {
      total_flux++;
      
      if(n_energy/eV < flux_energy[0]) energy_integral[0]+=n_energy/eV;
      
      if(temp_energy<lithium_energy[100]) {
	for(G4int i = 0; i<100; ++i) {
	  if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) {
	    lithium_flux[i]++;
	    cos_lithium_flux[i] += 1.0/std::abs(cos_angle);
	    integral_Eflux_46cm_restricted += n_energy;
	    if(std::abs(zPos-75.0) < 15.0) lithium_Zflux[i]++;
	  }
	}
      }
      
      if(temp_energy<low_energy[100]) {
	for(G4int i = 0; i<100; ++i) {
	  if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1]) {
	    low_flux[i]++;
	    cos_low_flux[i] += 1.0/std::abs(cos_angle);
	  }
	}
      }
      else if(temp_energy>flux_energy[0]) {
	for(G4int i = 0; i<32; ++i) {
	  if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1]) {
	    energy_integral[1]+=n_energy/eV;
	    flux[i]++;
	    if(cos_angle != 0.0) cos_flux[i] += 1.0/std::abs(cos_angle);
	    eflux[i] += n_energy/eV;
	  }
	}
      } 
      for(G4int j = 0; j<60; ++j) {

	//BUG? 23/11/18 if(n_energy > fine_energy[j]/1000000.0 && energy < fine_energy[j+1]/1000000.0) fine_eflux[j] += n_energy/eV;
	if(n_energy > fine_energy[j]/1000000.0 && n_energy < fine_energy[j+1]/1000000.0) fine_eflux[j] += n_energy/eV;

      } 
    }
        
    if(std::abs(radius-50.0*mm) < 0.1) {
      
      if(temp_energy<lithium_energy[100]) {
	for(G4int i = 0; i<100; ++i) {
	  if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) lithium_flux_5cm[i]++;
	}
      }
      
    }
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//x2018 void Test15Run::analyseNeutronShellFluence(G4double energy, G4double time, G4double startEnergy, G4int TrackID, G4int ParentID, G4double zMomentum, G4double startTime, G4double radius, G4double zPos, G4double parentEnergy, G4String parentParticle, G4double steplength, G4bool enter_sph, G4bool enter_cyl, G4bool exit_sph, G4bool exit_cyl, G4String Volume, G4bool enter_sph_front, G4bool exit_sph_front, G4int preParentReplica, G4int postParentReplica, G4int preReplica, G4int postReplica)
void Test15Run::analyseNeutronShellFluence(G4double shell_energy, G4double , G4double , G4int , G4int , G4double , G4double , G4double , G4double, G4double , G4String , G4double steplength, G4bool , G4bool , G4bool , G4bool , G4String , G4bool , G4bool , G4int , G4int , G4int , G4int )
{

  G4double temp_energy = shell_energy/eV;

  if(temp_energy<lithium_energy[100]) {
    for(G4int i = 0; i<100; ++i) {
      if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1])
	lithium_fluence_step_shell[i] += steplength/mm;
    }
  }
  
  if(temp_energy<low_energy[100]) {
    for(G4int i = 0; i<100; ++i) {
      if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1])
	low_fluence_step_shell[i] += steplength/mm;
    }
  }
  if(temp_energy>flux_energy[0]) {
    for(G4int i = 0; i<32; ++i) {
      if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1])
	fluence_step_shell[i] += steplength/mm;
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

