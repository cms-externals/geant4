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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on density effect
//
// Authors:    A.Bagulya, A.Ivanchenko 28.10.2009
//
//----------------------------------------------------------------------------

#include "G4DensityEffectData.hh"

#include "G4SystemOfUnits.hh"

#include <iomanip>

G4DensityEffectData::G4DensityEffectData() { Initialize(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DensityEffectData::Initialize()
{
  // R.M. Sternheimer et al. Density Effect for the Ionization Loss of Charged
  // Particles in Various Substances. Atom. Data Nucl. Data Tabl. 30 (1984) 261-271.
  //
  // Data[10]:  Eplasma rho -C X_0 X_1 a m delta_0 DELTA_{max} I
  // Eplasma - Plasma energy (in eV)
  // rho     - Sternheimer adjustment factor for the atomic excitation energies
  // -C
  // X_0, X_1, m, a - parameters in fitting formulas
  // delta_0 - Density-effect value delta(X_0)
  // DELTA_{max} - Upper bound for the error inherent in the fitting procedure
  // I - Mean ionisation potential in keV used in this parameterisation

  for (G4int i = 0; i < NDENSELEM; ++i) {
    indexZ[i] = -1;
    state[i] = kStateSolid;
  }

  // G4_lH2  index=0
  G4double M1[NDENSARRAY] = {
    7.031, 1.546, 3.2632, 0.4759, 1.9215, 0.13483, 5.6249, 0., 0.021, 21.8};
  AddMaterial(M1, "G4_lH2");

  // G4_H  index=1
  G4double M0[NDENSARRAY] = {
    0.263, 1.412, 9.5835, 1.8639, 3.2718, 0.14092, 5.7273, 0.0, 0.024, 19.2};
  AddMaterial(M0, "G4_H");
  indexZ[1] = 1;
  state[1] = kStateGas;

  // G4_He  index=2
  G4double M2[NDENSARRAY] = {0.263, 1.7, 11.1393, 2.2017, 3.6122, 0.13443, 5.8347, 0, 0.024, 41.8};
  AddMaterial(M2, "G4_He");
  indexZ[2] = 2;
  state[2] = kStateGas;

  // G4_Li  index=3
  G4double M3[NDENSARRAY] = {
    13.844, 1.535, 3.1221, 0.1304, 1.6397, 0.95136, 2.4993, 0.14, 0.062, 40.0};
  AddMaterial(M3, "G4_Li");
  indexZ[3] = 3;

  // G4_Be  index=4
  G4double M4[NDENSARRAY] = {
    26.098, 1.908, 2.7847, 0.0392, 1.6922, 0.80392, 2.4339, 0.14, 0.029, 63.7};
  AddMaterial(M4, "G4_Be");
  indexZ[4] = 4;

  // G4_B  index=5
  G4double M5[NDENSARRAY] = {
    30.17, 2.32, 2.8477, 0.0305, 1.9688, 0.56224, 2.4512, 0.14, 0.024, 76.0};
  AddMaterial(M5, "G4_B");
  indexZ[5] = 5;

  // G4_C  index=6
  G4double M6[NDENSARRAY] = {
    28.803, 2.376, 2.9925, -0.0351, 2.486, 0.2024, 3.0036, 0.1, 0.038, 78.0};
  AddMaterial(M6, "G4_C");
  indexZ[6] = 6;

  // G4_N  index=7
  G4double M7[NDENSARRAY] = {
    0.695, 1.984, 10.5400, 1.7378, 4.1323, 0.15349, 3.2125, 0.0, 0.086, 82.};
  AddMaterial(M7, "G4_N");
  indexZ[7] = 7;
  state[7] = kStateGas;

  // G4_O  index=8
  G4double M8[NDENSARRAY] = {
    0.744, 2.314, 10.7004, 1.7541, 4.3213, 0.11778, 3.2913, 0.0, 0.101, 95.};
  AddMaterial(M8, "G4_O");
  indexZ[8] = 8;
  state[8] = kStateGas;

  // G4_F  index=9
  G4double M9[NDENSARRAY] = {
    0.788, 2.450, 10.9653, 1.8433, 4.4096, 0.11083, 3.2962, 0.0, 0.121, 115.};
  AddMaterial(M9, "G4_F");
  indexZ[9] = 9;
  state[9] = kStateGas;

  // G4_Ne  index=10
  G4double M10[NDENSARRAY] = {
    0.587, 2.577, 11.9041, 2.0735, 4.6421, 0.08064, 3.5771, 0.0, 0.110, 137.};
  AddMaterial(M10, "G4_Ne");
  indexZ[10] = 10;
  state[10] = kStateGas;

  // G4_Na  index=11
  G4double M11[NDENSARRAY] = {
    19.641, 2.648, 5.0526, 0.2880, 3.1962, 0.07772, 3.6452, 0.08, 0.098, 149.};
  AddMaterial(M11, "G4_Na");
  indexZ[11] = 11;

  // G4_Mg  index=12
  G4double M12[NDENSARRAY] = {
    26.708, 2.331, 4.5297, 0.1499, 3.0668, 0.08163, 3.6166, 0.08, 0.073, 156.};
  AddMaterial(M12, "G4_Mg");
  indexZ[12] = 12;

  // G4_Al  index=13
  G4double M13[NDENSARRAY] = {
    32.86, 2.18, 4.2395, 0.1708, 3.0127, 0.08024, 3.6345, 0.12, 0.061, 166.};
  AddMaterial(M13, "G4_Al");
  indexZ[13] = 13;

  // G4_Si  index=14
  G4double M14[NDENSARRAY] = {
    31.055, 2.103, 4.4351, 0.2014, 2.8715, 0.14921, 3.2546, 0.14, 0.059, 173.};
  AddMaterial(M14, "G4_Si");
  indexZ[14] = 14;

  // G4_P  index=15
  G4double M15[NDENSARRAY] = {
    29.743, 2.056, 4.5214, 0.1696, 2.7815, 0.2361, 2.9158, 0.14, 0.057, 173.};
  AddMaterial(M15, "G4_P");
  indexZ[15] = 15;

  // G4_S  index=16
  G4double M16[NDENSARRAY] = {
    28.789, 2.131, 4.6659, 0.158, 2.7159, 0.33992, 2.6456, 0.14, 0.059, 180.};
  AddMaterial(M16, "G4_S");
  indexZ[16] = 16;

  // G4_Cl  index=17
  G4double M17[NDENSARRAY] = {
    1.092, 1.734, 11.1421, 1.5555, 4.2994, 0.19849, 2.9702, 0.0, 0.041, 174.};
  AddMaterial(M17, "G4_Cl");
  indexZ[17] = 17;
  state[17] = kStateGas;

  // G4_Ar  index=18
  G4double M18[NDENSARRAY] = {
    0.789, 1.753, 11.9480, 1.7635, 4.4855, 0.19714, 2.9618, 0.0, 0.037, 188.};
  AddMaterial(M18, "G4_Ar");
  indexZ[18] = 18;
  state[18] = kStateGas;

  // G4_K  index=19
  G4double M19[NDENSARRAY] = {
    18.65, 1.830, 5.6423, 0.3851, 3.1724, 0.19827, 2.9233, 0.10, 0.035, 190.};
  AddMaterial(M19, "G4_K");
  indexZ[19] = 19;

  // G4_Ca  index=20
  G4double M20[NDENSARRAY] = {
    25.342, 1.666, 5.0396, 0.3228, 3.1191, 0.15643, 3.0745, 0.14, 0.031, 191.};
  AddMaterial(M20, "G4_Ca");
  indexZ[20] = 20;

  // G4_Sc  index=21
  G4double M21[NDENSARRAY] = {
    34.050, 1.826, 4.6949, 0.1640, 3.0593, 0.15754, 3.0517, 0.10, 0.027, 216.};
  AddMaterial(M21, "G4_Sc");
  indexZ[21] = 21;

  // G4_Ti  index=22
  G4double M22[NDENSARRAY] = {
    41.619, 1.969, 4.4450, 0.0957, 3.0386, 0.15662, 3.0302, 0.12, 0.025, 233.};
  AddMaterial(M22, "G4_Ti");
  indexZ[22] = 22;

  // G4_V  index=23
  G4double M23[NDENSARRAY] = {
    47.861, 2.070, 4.2659, 0.0691, 3.0322, 0.15436, 3.0163, 0.14, 0.024, 245.};
  AddMaterial(M23, "G4_V");
  indexZ[23] = 23;

  // G4_Cr  index=24
  G4double M24[NDENSARRAY] = {
    52.458, 2.181, 4.1781, 0.0340, 3.0451, 0.15419, 2.9896, 0.14, 0.023, 257.};
  AddMaterial(M24, "G4_Cr");
  indexZ[24] = 24;

  // G4_Mn  index=25
  G4double M25[NDENSARRAY] = {
    53.022, 2.347, 4.2702, 0.0447, 3.1074, 0.14973, 2.9796, 0.14, 0.021, 272.};
  AddMaterial(M25, "G4_Mn");
  indexZ[25] = 25;

  // G4_Fe  index=26
  G4double M26[NDENSARRAY] = {
    55.172, 2.504, 4.2911, -0.0012, 3.1531, 0.1468, 2.9632, 0.12, 0.021, 286.};
  AddMaterial(M26, "G4_Fe");
  indexZ[26] = 26;

  // G4_Co  index=27
  G4double M27[NDENSARRAY] = {
    58.188, 2.626, 4.2601, -0.0187, 3.1790, 0.14474, 2.9502, 0.12, 0.019, 297.};
  AddMaterial(M27, "G4_Co");
  indexZ[27] = 27;

  // G4_Ni  index=28
  G4double M28[NDENSARRAY] = {
    59.385, 2.889, 4.3115, -0.0566, 3.1851, 0.16496, 2.843, 0.10, 0.020, 311.};
  AddMaterial(M28, "G4_Ni");
  indexZ[28] = 28;

  // G4_Cu  index=29
  G4double M29[NDENSARRAY] = {
    58.270, 2.956, 4.4190, -0.0254, 3.2792, 0.14339, 2.9044, 0.08, 0.019, 322.};
  AddMaterial(M29, "G4_Cu");
  indexZ[29] = 29;

  // G4_Zn  index=30
  G4double M30[NDENSARRAY] = {
    52.132, 3.142, 4.6906, 0.0049, 3.3668, 0.14714, 2.8652, 0.08, 0.019, 330.};
  AddMaterial(M30, "G4_Zn");
  indexZ[30] = 30;

  // G4_Ga  index=31
  G4double M31[NDENSARRAY] = {
    46.688, 2.747, 4.9353, 0.2267, 3.5434, 0.09440, 3.1314, 0.14, 0.019, 334.};
  AddMaterial(M31, "G4_Ga");
  indexZ[31] = 31;

  // G4_Ge  index=32
  G4double M32[NDENSARRAY] = {
    44.141, 2.461, 5.1411, 0.3376, 3.6096, 0.07188, 3.3306, 0.14, 0.025, 350.};
  AddMaterial(M32, "G4_Ge");
  indexZ[32] = 32;

  // G4_As  index=33
  G4double M33[NDENSARRAY] = {
    45.779, 2.219, 5.0510, 0.1767, 3.5702, 0.06633, 3.4176, 0.00, 0.030, 347.};
  AddMaterial(M33, "G4_As");
  indexZ[33] = 33;

  // G4_Se  index=34
  G4double M34[NDENSARRAY] = {
    40.112, 2.104, 5.3210, 0.2258, 3.6264, 0.06568, 3.4317, 0.10, 0.024, 348.};
  AddMaterial(M34, "G4_Se");
  indexZ[34] = 34;

  // G4_Br  index=35
  G4double M35[NDENSARRAY] = {
    1.604, 1.845, 11.7307, 1.5262, 4.9899, 0.06335, 3.467, 0, 0.022, 343.};
  AddMaterial(M35, "G4_Br");
  indexZ[35] = 35;
  state[35] = kStateGas;

  // G4_Kr  index=36
  G4double M36[NDENSARRAY] = {
    1.114, 1.77, 12.5115, 1.7158, 5.0748, 0.07446, 3.4051, 0, 0.025, 352.};
  AddMaterial(M36, "G4_Kr");
  indexZ[36] = 36;
  state[36] = kStateGas;

  // G4_Rb  index=37
  G4double M37[NDENSARRAY] = {
    23.467, 1.823, 6.4776, 0.5737, 3.7995, 0.07261, 3.4177, 0.14, 0.026, 363.};
  AddMaterial(M37, "G4_Rb");
  indexZ[37] = 37;

  // G4_Sr  index=38
  G4double M38[NDENSARRAY] = {
    30.244, 1.707, 5.9867, 0.4585, 3.6778, 0.07165, 3.4435, 0.14, 0.026, 366.};
  AddMaterial(M38, "G4_Sr");
  indexZ[38] = 38;

  // G4_Y  index=39
  G4double M39[NDENSARRAY] = {
    40.346, 1.649, 5.4801, 0.3608, 3.5542, 0.07138, 3.4585, 0.14, 0.027, 379.};
  AddMaterial(M39, "G4_Y");
  indexZ[39] = 39;

  // G4_Zr  index=40
  G4double M40[NDENSARRAY] = {
    48.671, 1.638, 5.1774, 0.2957, 3.489, 0.07177, 3.4533, 0.14, 0.028, 393.};
  AddMaterial(M40, "G4_Zr");
  indexZ[40] = 40;

  // G4_Nb  index=41
  G4double M41[NDENSARRAY] = {
    56.039, 1.734, 5.0141, 0.1785, 3.2201, 0.13883, 3.093, 0.14, 0.036, 417.};
  AddMaterial(M41, "G4_Nb");
  indexZ[41] = 41;

  // G4_Mo  index=42
  G4double M42[NDENSARRAY] = {
    60.951, 1.658, 4.8793, 0.2267, 3.2784, 0.10525, 3.2549, 0.14, 0.03, 424.};
  AddMaterial(M42, "G4_Mo");
  indexZ[42] = 42;

  // G4_Tc  index=43
  G4double M43[NDENSARRAY] = {
    64.760, 1.727, 4.7769, 0.0949, 3.1253, 0.16572, 2.9738, 0.14, 0.040, 428.};
  AddMaterial(M43, "G4_Tc");
  indexZ[43] = 43;

  // G4_Ru  index=44
  G4double M44[NDENSARRAY] = {
    66.978, 1.780, 4.7694, 0.0599, 3.0834, 0.19342, 2.8707, 0.14, 0.046, 441.};
  AddMaterial(M44, "G4_Ru");
  indexZ[44] = 44;

  // G4_Rh  index=45
  G4double M45[NDENSARRAY] = {
    67.128, 1.804, 4.8008, 0.0576, 3.1069, 0.19205, 2.8633, 0.14, 0.046, 449.};
  AddMaterial(M45, "G4_Rh");
  indexZ[45] = 45;

  // G4_Pd  index=46
  G4double M46[NDENSARRAY] = {
    65.683, 1.911, 4.9358, 0.0563, 3.0555, 0.24178, 2.7239, 0.14, 0.047, 470.};
  AddMaterial(M46, "G4_Pd");
  indexZ[46] = 46;

  // G4_Ag  index=47
  G4double M47[NDENSARRAY] = {
    61.635, 1.933, 5.0630, 0.0657, 3.1074, 0.24585, 2.6899, 0.14, 0.052, 470.};
  AddMaterial(M47, "G4_Ag");
  indexZ[47] = 47;

  // G4_Cd  index=48
  G4double M48[NDENSARRAY] = {
    55.381, 1.895, 5.2727, 0.1281, 3.1667, 0.24609, 2.6772, 0.14, 0.051, 469.};
  AddMaterial(M48, "G4_Cd");
  indexZ[48] = 48;

  // G4_In  index=49
  G4double M49[NDENSARRAY] = {
    50.896, 1.851, 5.5211, 0.2406, 3.2032, 0.23879, 2.7144, 0.14, 0.044, 488.};
  AddMaterial(M49, "G4_In");
  indexZ[49] = 49;

  // G4_Sn  index=50
  G4double M50[NDENSARRAY] = {
    50.567, 1.732, 5.5340, 0.2879, 3.2959, 0.18689, 2.8576, 0.14, 0.037, 488.};
  AddMaterial(M50, "G4_Sn");
  indexZ[50] = 50;

  // G4_Sb  index=51
  G4double M51[NDENSARRAY] = {
    48.242, 1.645, 5.6241, 0.3189, 3.3489, 0.16652, 2.9319, 0.14, 0.034, 487.};
  AddMaterial(M51, "G4_Sb");
  indexZ[51] = 51;

  // G4_Te  index=52
  G4double M52[NDENSARRAY] = {
    45.952, 1.577, 5.7131, 0.3296, 3.4418, 0.13815, 3.0354, 0.14, 0.033, 485.};
  AddMaterial(M52, "G4_Te");
  indexZ[52] = 52;

  // G4_I  index=53
  G4double M53[NDENSARRAY] = {
    41.348, 1.498, 5.9488, 0.0549, 3.2596, 0.23766, 2.7276, 0.0, 0.045, 491.};
  AddMaterial(M53, "G4_I");
  indexZ[53] = 53;

  // G4_Xe  index=54
  G4double M54[NDENSARRAY] = {
    1.369, 1.435, 12.7281, 1.563, 4.7371, 0.23314, 2.7414, 0, 0.043, 482.};
  AddMaterial(M54, "G4_Xe");
  indexZ[54] = 54;
  state[54] = kStateGas;

  // G4_Cs  index=55
  G4double M55[NDENSARRAY] = {
    25.37, 1.462, 6.9135, 0.5473, 3.5914, 0.18233, 2.8866, 0.14, 0.035, 488.};
  AddMaterial(M55, "G4_Cs");
  indexZ[55] = 55;

  // G4_Ba  index=56
  G4double M56[NDENSARRAY] = {
    34.425, 1.410, 6.3153, 0.4190, 3.4547, 0.18268, 2.8906, 0.14, 0.035, 491.};
  AddMaterial(M56, "G4_Ba");
  indexZ[56] = 56;

  // G4_La  index=57
  G4double M57[NDENSARRAY] = {
    45.792, 1.392, 5.7850, 0.3161, 3.3293, 0.18591, 2.8828, 0.14, 0.036, 501.};
  AddMaterial(M57, "G4_La");
  indexZ[57] = 57;

  // G4_Ce  index=58
  G4double M58[NDENSARRAY] = {
    47.834, 1.461, 5.7837, 0.2713, 3.3432, 0.18885, 2.8592, 0.14, 0.040, 523.};
  AddMaterial(M58, "G4_Ce");
  indexZ[58] = 58;

  // G4_Pr  index=59
  G4double M59[NDENSARRAY] = {
    48.301, 1.520, 5.8096, 0.2333, 3.2773, 0.23265, 2.7331, 0.14, 0.041, 535.};
  AddMaterial(M59, "G4_Pr");
  indexZ[59] = 59;

  // G4_Ne  index=60
  G4double M60[NDENSARRAY] = {
    48.819, 1.588, 5.8290, 0.1984, 3.3063, 0.23530, 2.7050, 0.14, 0.044, 546.};
  AddMaterial(M60, "G4_Ne");
  indexZ[60] = 60;

  // G4_Pr  index=61
  G4double M61[NDENSARRAY] = {
    50.236, 1.672, 5.8224, 0.1627, 3.3199, 0.24280, 2.6674, 0.14, 0.048, 560.};
  AddMaterial(M61, "G4_Pr");
  indexZ[61] = 61;

  // G4_Sa  index=62
  G4double M62[NDENSARRAY] = {
    50.540, 1.749, 5.8597, 0.1520, 3.3460, 0.24698, 2.6403, 0.14, 0.053, 574.};
  AddMaterial(M62, "G4_Sa");
  indexZ[62] = 62;

  // G4_Eu  index=63
  G4double M63[NDENSARRAY] = {
    42.484, 1.838, 6.2278, 0.1888, 3.4633, 0.24448, 2.6245, 0.14, 0.06, 580.};
  AddMaterial(M63, "G4_Eu");
  indexZ[63] = 63;

  // G4_Gd  index=64
  G4double M64[NDENSARRAY] = {
    51.672, 1.882, 5.8738, 0.1058, 3.3932, 0.25109, 2.5977, 0.14, 0.061, 591.};
  AddMaterial(M64, "G4_Gd");
  indexZ[64] = 64;

  // G4_Tb  index=65
  G4double M65[NDENSARRAY] = {
    52.865, 1.993, 5.9045, 0.0947, 3.4224, 0.24453, 2.6056, 0.14, 0.063, 614.};
  AddMaterial(M65, "G4_Tb");
  indexZ[65] = 65;

  // G4_Dy  index=66
  G4double M66[NDENSARRAY] = {
    53.698, 2.081, 5.9183, 0.0822, 3.4474, 0.24665, 2.5849, 0.14, 0.061, 628.};
  AddMaterial(M66, "G4_Dy");
  indexZ[66] = 66;

  // G4_Ho  index=67
  G4double M67[NDENSARRAY] = {
    54.467, 2.197, 5.9587, 0.0761, 3.4782, 0.24638, 2.5726, 0.14, 0.062, 650.};
  AddMaterial(M67, "G4_Ho");
  indexZ[67] = 67;

  // G4_Er  index=68
  G4double M68[NDENSARRAY] = {
    55.322, 2.26, 5.9521, 0.0648, 3.4922, 0.24823, 2.5573, 0.14, 0.061, 658.};
  AddMaterial(M68, "G4_Er");
  indexZ[68] = 68;

  // G4_Tm  index=69
  G4double M69[NDENSARRAY] = {
    56.225, 2.333, 5.9677, 0.0812, 3.5085, 0.24889, 2.5469, 0.14, 0.062, 674.};
  AddMaterial(M69, "G4_Tm");
  indexZ[69] = 69;

  // G4_Yb  index=70
  G4double M70[NDENSARRAY] = {
    47.546, 2.505, 6.3325, 0.1199, 3.6246, 0.25295, 2.5141, 0.14, 0.071, 684.};
  AddMaterial(M70, "G4_Yb");
  indexZ[70] = 70;

  // G4_Lu  index=71
  G4double M71[NDENSARRAY] = {
    57.581, 2.348, 5.9785, 0.1560, 3.5218, 0.24033, 2.5643, 0.14, 0.054, 694.};
  AddMaterial(M71, "G4_Lu");
  indexZ[71] = 71;

  // G4_Hf  index=72
  G4double M72[NDENSARRAY] = {
    66.770, 2.174, 5.7139, 0.1965, 3.4337, 0.22918, 2.6155, 0.14, 0.035, 705.};
  AddMaterial(M72, "G4_Hf");
  indexZ[72] = 72;

  // G4_Ta  index=73
  G4double M73[NDENSARRAY] = {
    74.692, 2.07, 5.5262, 0.2117, 3.4805, 0.17798, 2.7623, 0.14, 0.03, 718.};
  AddMaterial(M73, "G4_Ta");
  indexZ[73] = 73;

  // G4_W  index=74
  G4double M74[NDENSARRAY] = {
    80.315, 1.997, 5.4059, 0.2167, 3.496, 0.15509, 2.8447, 0.14, 0.027, 727.};
  AddMaterial(M74, "G4_W");
  indexZ[74] = 74;

  // G4_Re  index=75
  G4double M75[NDENSARRAY] = {
    83.846, 1.976, 5.3445, 0.0559, 3.4845, 0.15184, 2.8627, 0.08, 0.026, 736.};
  AddMaterial(M75, "G4_Re");
  indexZ[75] = 75;

  // G4_Os  index=76
  G4double M76[NDENSARRAY] = {
    86.537, 1.947, 5.3083, 0.0891, 3.5414, 0.12751, 2.9608, 0.10, 0.023, 746.};
  AddMaterial(M76, "G4_Os");
  indexZ[76] = 76;

  // G4_Ir  index=77
  G4double M77[NDENSARRAY] = {
    86.357, 1.927, 5.3418, 0.0819, 3.5480, 0.12690, 2.9658, 0.10, 0.023, 757.};
  AddMaterial(M77, "G4_Ir");
  indexZ[77] = 77;

  // G4_Pt  index=78
  G4double M78[NDENSARRAY] = {
    84.389, 1.965, 5.4732, 0.1484, 3.6212, 0.11128, 3.0417, 0.12, 0.021, 790.};
  AddMaterial(M78, "G4_Pt");
  indexZ[78] = 78;

  // G4_Au  index=79
  G4double M79[NDENSARRAY] = {
    80.215, 1.926, 5.5747, 0.2021, 3.6979, 0.09756, 3.1101, 0.14, 0.020, 790.};
  AddMaterial(M79, "G4_Au");
  indexZ[79] = 79;

  // G4_Hg  index=80
  G4double M80[NDENSARRAY] = {
    66.977, 1.904, 5.9605, 0.2756, 3.7275, 0.11014, 3.0519, 0.14, 0.021, 800.};
  AddMaterial(M80, "G4_Hg");
  indexZ[80] = 80;

  // G4_Tl  index=81
  G4double M81[NDENSARRAY] = {
    62.104, 1.814, 6.1365, 0.3491, 3.8044, 0.09455, 3.1450, 0.14, 0.019, 810.};
  AddMaterial(M81, "G4_Tl");
  indexZ[81] = 81;

  // G4_Pb  index=82
  G4double M82[NDENSARRAY] = {
    61.072, 1.755, 6.2018, 0.3776, 3.8073, 0.09359, 3.1608, 0.14, 0.019, 823.};
  AddMaterial(M82, "G4_Pb");
  indexZ[82] = 82;

  // G4_Bi  index=83
  G4double M83[NDENSARRAY] = {
    56.696, 1.684, 6.3505, 0.4152, 3.8248, 0.0941, 3.1671, 0.14, 0.02, 823.};
  AddMaterial(M83, "G4_Bi");
  indexZ[83] = 83;

  // G4_Po  index=84
  G4double M84[NDENSARRAY] = {
    55.773, 1.637, 6.4003, 0.4267, 3.8293, 0.09282, 3.183, 0.14, 0.02, 830.};
  AddMaterial(M84, "G4_Po");
  indexZ[84] = 84;

  // G4_Rn  index=85
  G4double M85[NDENSARRAY] = {
    1.708, 1.458, 13.2839, 1.5368, 4.9889, 0.20798, 2.7409, 0, 0.057, 794.};
  AddMaterial(M85, "G4_Rn");
  indexZ[86] = 85;
  state[86] = kStateGas;

  // G4_Ra  index=86
  G4double M86[NDENSARRAY] = {
    40.205, 1.403, 7.0452, 0.5991, 3.9428, 0.08804, 3.2454, 0.14, 0.022, 826.};
  AddMaterial(M86, "G4_Ra");
  indexZ[88] = 86;

  // G4_Ac  index=87
  G4double M87[NDENSARRAY] = {
    57.254, 1.380, 6.3742, 0.4559, 3.7966, 0.08567, 3.2683, 0.14, 0.023, 841.};
  AddMaterial(M87, "G4_Ac");
  indexZ[89] = 87;

  // G4_Th  index=88
  G4double M88[NDENSARRAY] = {
    61.438, 1.363, 6.2473, 0.4202, 3.7681, 0.08655, 3.2610, 0.14, 0.025, 847.};
  AddMaterial(M88, "G4_Th");
  indexZ[90] = 88;

  // G4_Pa  index=89
  G4double M89[NDENSARRAY] = {
    70.901, 1.42, 6.0327, 0.3144, 3.5079, .14770, 2.9845, 0.14, 0.036, 878.};
  AddMaterial(M89, "G4_Pa");
  indexZ[91] = 89;

  // G4_U  index=90
  G4double M90[NDENSARRAY] = {
    77.986, 1.447, 5.8694, 0.2260, 3.3721, .19677, 2.8171, 0.14, 0.043, 890.};
  AddMaterial(M90, "G4_U");
  indexZ[92] = 90;

  // G4_Np  index=91
  G4double M91[NDENSARRAY] = {
    81.221, 1.468, 5.8149, 0.1869, 3.369, 0.19741, 2.8082, 0.14, 0.043, 902.};
  AddMaterial(M91, "G4_Np");
  indexZ[93] = 91;

  // G4_Pu  index=92
  G4double M92[NDENSARRAY] = {
    80.486, 1.519, 5.8748, 0.1557, 3.3981, 0.20419, 2.7679, 0.14, 0.057, 921.};
  AddMaterial(M92, "G4_Pu");
  indexZ[94] = 92;

  // G4_Am  index=93
  G4double M93[NDENSARRAY] = {
    66.607, 1.552, 6.2813, 0.2274, 3.5021, 0.20308, 2.7615, 0.14, 0.056, 934.};
  AddMaterial(M93, "G4_Am");
  indexZ[95] = 93;

  // G4_Cm  index=94
  G4double M94[NDENSARRAY] = {
    66.022, 1.559, 6.3097, 0.2484, 3.516, .20257, 2.7579, 0.14, 0.056, 939.};
  AddMaterial(M94, "G4_Cm");
  indexZ[96] = 94;

  // G4_Bk  index=95
  G4double M95[NDENSARRAY] = {
    67.557, 1.574, 6.2912, 0.2378, 3.5186, .20192, 2.7560, 0.14, 0.062, 952.};
  AddMaterial(M95, "G4_Bk");
  indexZ[97] = 95;

  // G4_A-150_TISSUE  index=96
  G4double M96[NDENSARRAY] = {
    22.667, 1.950, 3.1100, 0.1329, 2.6234, 0.10783, 3.4442, 0, 0.048, 65.1};
  AddMaterial(M96, "G4_A-150_TISSUE");

  // G4_ACETONE  index=97
  G4double M97[NDENSARRAY] = {
    19.010, 1.976, 3.4341, 0.2197, 2.6928, 0.11100, 3.4047, 0, 0.069, 64.2};
  AddMaterial(M97, "G4_ACETONE");

  // G4_ACETYLENE  index=98
  G4double M98[NDENSARRAY] = {
    0.700, 1.784, 9.8419, 1.6017, 4.0074, 0.12167, 3.4277, 0, 0.080, 58.2};
  AddMaterial(M98, "G4_ACETYLENE");

  // G4_ADENINE  index=99
  G4double M99[NDENSARRAY] = {
    24.098, 1.892, 3.1724, 0.1295, 2.4219, 0.20908, 3.0271, 0, 0.052, 71.4};
  AddMaterial(M99, "G4_ADENINE");

  // G4_ADIPOSE_TISSUE_ICRP  index=100
  G4double M100[NDENSARRAY] = {
    20.655, 1.987, 3.2367, 0.1827, 2.6530, 0.10278, 3.4817, 0, 0.060, 63.2};
  AddMaterial(M100, "G4_ADIPOSE_TISSUE_ICRP");

  // G4_AIR  index=101
  G4double M101[NDENSARRAY] = {
    0.707, 2.054, 10.5961, 1.7418, 4.2759, 0.10914, 3.3994, 0, 0.090, 85.7};
  AddMaterial(M101, "G4_AIR");

  // G4_ALANINE  index=102
  G4double M102[NDENSARRAY] = {
    25.204, 2.074, 3.0965, 0.1354, 2.6336, 0.11484, 3.3526, 0, 0.056, 71.9};
  AddMaterial(M102, "G4_ALANINE");

  // G4_ALUMINIM_OXIDE  index=103
  G4double M103[NDENSARRAY] = {
    40.206, 2.394, 3.5682, 0.0402, 2.8665, 0.08500, 3.5458, 0, 0.031, 145.2};
  AddMaterial(M103, "G4_ALUMINIM_OXIDE");

  // G4_AMBER  index=104
  G4double M104[NDENSARRAY] = {
    22.450, 1.946, 3.0701, 0.1335, 2.5610, 0.11934, 3.4098, 0, 0.053, 63.2};
  AddMaterial(M104, "G4_AMBER");

  // G4_AMMONIA  index=105
  G4double M105[NDENSARRAY] = {
    0.635, 1.814, 9.8763, 1.6822, 4.1158, 0.08315, 3.6464, 0, 0.102, 53.7};
  AddMaterial(M105, "G4_AMMONIA");

  // G4_ANILINE  index=106
  G4double M106[NDENSARRAY] = {
    21.361, 1.938, 3.2622, 0.1618, 2.5805, 0.13134, 3.3434, 0, 0.052, 66.2};
  AddMaterial(M106, "G4_ANILINE");

  // G4_ANTHRACENE  index=107
  G4double M107[NDENSARRAY] = {
    23.704, 1.954, 3.1514, 0.1146, 2.5213, 0.14677, 3.2831, 0, 0.042, 69.5};
  AddMaterial(M107, "G4_ANTHRACENE");

  // G4_B-100_BONE  index=108
  G4double M108[NDENSARRAY] = {
    25.199, 2.013, 3.4528, 0.1252, 3.042, 0.05268, 3.7365, 0, 0.043, 85.9};
  AddMaterial(M108, "G4_B-100_BONE");

  // G4_BAKELITE  index=109
  G4double M109[NDENSARRAY] = {
    23.408, 2.046, 3.2582, 0.1471, 2.6055, 0.12713, 3.347, 0, 0.052, 72.4};
  AddMaterial(M109, "G4_BAKELITE");

  // G4_BARIUM_FLUORIDE  index=110
  G4double M110[NDENSARRAY] = {
    41.398, 1.727, 5.4122, -0.0098, 3.3871, 0.15991, 2.8867, 0, 0.034, 375.9};
  AddMaterial(M110, "G4_BARIUM_FLUORIDE");

  // G4_BARIUM_SULFATE  index=111
  G4double M111[NDENSARRAY] = {
    40.805, 1.893, 4.8923, -0.0128, 3.4069, 0.11747, 3.0427, 0, 0.03, 285.7};
  AddMaterial(M111, "G4_BARIUM_SULFATE");

  // G4_BENZENE  index=112
  G4double M112[NDENSARRAY] = {
    19.806, 1.873, 3.3269, 0.171, 2.5091, 0.16519, 3.2174, 0, 0.052, 63.4};
  AddMaterial(M112, "G4_BENZENE");

  // G4_BERYLLIUM_OXIDE  index=113
  G4double M113[NDENSARRAY] = {
    34.629, 2.296, 2.9801, 0.0241, 2.5846, 0.10755, 3.4927, 0, 0.031, 93.2};
  AddMaterial(M113, "G4_BERYLLIUM_OXIDE");

  // G4_BGO  index=114
  G4double M114[NDENSARRAY] = {
    49.904, 2.121, 5.7409, 0.0456, 3.7816, 0.09569, 3.0781, 0, 0.023, 534.1};
  AddMaterial(M114, "G4_BGO");

  // G4_BLOOD_ICRP  index=115
  G4double M115[NDENSARRAY] = {
    22.001, 2.184, 3.4581, 0.2239, 2.8017, 0.08492, 3.5406, 0, 0.088, 75.2};
  AddMaterial(M115, "G4_BLOOD_ICRP");

  // G4_BONE_COMPACT_ICRU  index=116
  G4double M116[NDENSARRAY] = {
    28.536, 2.091, 3.3390, 0.0944, 3.0201, 0.05822, 3.6419, 0, 0.042, 91.9};
  AddMaterial(M116, "G4_BONE_COMPACT_ICRU");

  // G4_BONE_CORTICAL_ICRP  index=117
  G4double M117[NDENSARRAY] = {
    28.298, 2.118, 3.6488, 0.1161, 3.0919, 0.06198, 3.5919, 0, 0.04, 106.4};
  AddMaterial(M117, "G4_BONE_CORTICAL_ICRP");

  // G4_BORON_CARBIDE  index=118
  G4double M118[NDENSARRAY] = {
    31.38, 2.14, 2.9859, 0.0093, 2.1006, 0.37087, 2.8076, 0, 0.022, 84.7};
  AddMaterial(M118, "G4_BORON_CARBIDE");

  // G4_BORON_OXIDE  index=119
  G4double M119[NDENSARRAY] = {
    27.107, 2.446, 3.6027, 0.1843, 2.7379, 0.11548, 3.3832, 0, 0.053, 99.6};
  AddMaterial(M119, "G4_BORON_OXIDE");

  // G4_BRAIN_ICRP  index=120
  G4double M120[NDENSARRAY] = {
    21.772, 2.162, 3.4279, 0.2206, 2.8021, 0.08255, 3.5585, 0, 0.086, 73.3};
  AddMaterial(M120, "G4_BRAIN_ICRP");

  // G4_BUTANE  index=121
  G4double M121[NDENSARRAY] = {1.101, 1.727, 8.5633, 1.3788, 3.7524, 0.10852, 3.4884, 0, 0.1, 48.3};
  AddMaterial(M121, "G4_BUTANE");

  // G4_N-BUTYL_ALCOHOL  index=122
  G4double M122[NDENSARRAY] = {
    19.52, 1.942, 3.2425, 0.1937, 2.6439, 0.10081, 3.5139, 0, 0.065, 59.9};
  AddMaterial(M122, "G4_N-BUTYL_ALCOHOL");

  // G4_C-552  index=123
  G4double M123[NDENSARRAY] = {
    27.023, 2.128, 3.3338, 0.151, 2.7083, 0.10492, 3.4344, 0, 0.053, 86.8};
  AddMaterial(M123, "G4_C-552");

  // G4_CADMIUM_TELLURIDE  index=124
  G4double M124[NDENSARRAY] = {
    46.314, 1.935, 5.9096, 0.0438, 3.2836, 0.2484, 2.6665, 0, 0.057, 539.3};
  AddMaterial(M124, "G4_CADMIUM_TELLURIDE");

  // G4_CADMIUM_TUNGSTATE  index=125
  G4double M125[NDENSARRAY] = {
    52.954, 2.289, 5.3594, 0.0123, 3.5941, 0.12861, 2.915, 0, 0.027, 468.3};
  AddMaterial(M125, "G4_CADMIUM_TUNGSTATE");

  // G4_CALCIUM_CARBONATE  index=126
  G4double M126[NDENSARRAY] = {
    34.08, 2.141, 3.7738, 0.0492, 3.0549, 0.08301, 3.412, 0, 0.037, 136.4};
  AddMaterial(M126, "G4_CALCIUM_CARBONATE");

  // G4_CALCIUM_FLUORIDE  index=127
  G4double M127[NDENSARRAY] = {
    35.849, 2.127, 4.0653, 0.0676, 3.1683, 0.06942, 3.5263, 0, 0.044, 166.0};
  AddMaterial(M127, "G4_CALCIUM_FLUORIDE");

  // G4_CALCIUM_OXIDE  index=128
  G4double M128[NDENSARRAY] = {
    36.988, 1.973, 4.1209, -0.0172, 3.0171, 0.12128, 3.1936, 0, 0.024, 176.1};
  AddMaterial(M128, "G4_CALCIUM_OXIDE");

  // G4_CALCIUM_SULFATE  index=129
  G4double M129[NDENSARRAY] = {
    35.038, 2.179, 3.9388, 0.0587, 3.1229, 0.07708, 3.4495, 0, 0.021, 152.3};
  AddMaterial(M129, "G4_CALCIUM_SULFATE");

  // G4_CALCIUM_TUNGSTATE  index=130
  G4double M130[NDENSARRAY] = {
    46.934, 2.262, 5.2603, 0.0323, 3.8932, 0.0621, 3.2649, 0, 0.021, 395.0};
  AddMaterial(M130, "G4_CALCIUM_TUNGSTATE");

  // G4_CARBON_DIOXIDE  index=131
  G4double M131[NDENSARRAY] = {
    0.874, 2.118, 10.1537, 1.6294, 4.1825, 0.11768, 3.3227, 0, 0.091, 85.0};
  AddMaterial(M131, "G4_CARBON_DIOXIDE");

  // G4_CARBON_TETRACHLORIDE  index=132
  G4double M132[NDENSARRAY] = {
    25.234, 1.742, 4.7712, 0.1773, 2.9165, 0.19018, 3.0116, 0, 0.041, 166.3};
  AddMaterial(M132, "G4_CARBON_TETRACHLORIDE");

  // G4_CELLULOSE_CELLOPHANE  index=133
  G4double M133[NDENSARRAY] = {25.008, 2.17, 3.2647, 0.158, 2.6778, 0.11151, 3.381, 0, 0.06, 77.6};
  AddMaterial(M133, "G4_CELLULOSE_CELLOPHANE");

  // G4_CELLULOSE_BUTYRATE  index=134
  G4double M134[NDENSARRAY] = {
    23.041, 2.128, 3.3497, 0.1794, 2.6809, 0.11444, 3.3738, 0, 0.056, 74.6};
  AddMaterial(M134, "G4_CELLULOSE_BUTYRATE");

  // G4_CELLULOSE_NITRATE  index=135
  G4double M135[NDENSARRAY] = {
    25.224, 2.252, 3.4762, 0.1897, 2.7253, 0.11813, 3.3237, 0, 0.063, 87.0};
  AddMaterial(M135, "G4_CELLULOSE_NITRATE");

  // G4_CERIC_SULFATE  index=136
  G4double M136[NDENSARRAY] = {
    21.743, 2.205, 3.5212, 0.2363, 2.8769, 0.07666, 3.5607, 0, 0.095, 76.7};
  AddMaterial(M136, "G4_CERIC_SULFATE");

  // G4_CESIUM_FLUORIDE  index=137
  G4double M137[NDENSARRAY] = {
    37.942, 1.714, 5.9046, 0.0084, 3.3374, 0.22052, 2.728, 0, 0.044, 440.7};
  AddMaterial(M137, "G4_CESIUM_FLUORIDE");

  // G4_CESIUM_IODIDE  index=138
  G4double M138[NDENSARRAY] = {
    39.455, 1.672, 6.2807, 0.0395, 3.3353, 0.25381, 2.6657, 0, 0.067, 553.1};
  AddMaterial(M138, "G4_CESIUM_IODIDE");

  // G4_CHLOROBENZENE  index=139
  G4double M139[NDENSARRAY] = {
    21.752, 1.889, 3.8201, 0.1714, 2.9272, 0.09856, 3.3797, 0, 0.031, 89.1};
  AddMaterial(M139, "G4_CHLOROBENZENE");

  // G4_CHLOROFORM  index=140
  G4double M140[NDENSARRAY] = {
    24.462, 1.734, 4.7055, 0.1786, 2.9581, 0.16959, 3.0627, 0, 0.038, 156.};
  AddMaterial(M140, "G4_CHLOROFORM");

  // G4_CONCRETE  index=141
  G4double M141[NDENSARRAY] = {
    30.986, 2.322, 3.9464, 0.1301, 3.0466, 0.07515, 3.5467, 0, 0.024, 135.2};
  AddMaterial(M141, "G4_CONCRETE");

  // G4_CYCLOHEXANE  index=142
  G4double M142[NDENSARRAY] = {
    19.207, 1.861, 3.1544, 0.1728, 2.5549, 0.12035, 3.4278, 0, 0.057, 56.4};
  AddMaterial(M142, "G4_CYCLOHEXANE");

  // G4_1,2-DICHLOROBENZENE  index=143
  G4double M143[NDENSARRAY] = {
    23.354, 1.862, 4.0348, 0.1587, 2.8276, 0.1601, 3.0836, 0, 0.029, 106.5};
  AddMaterial(M143, "G4_1,2-DICHLOROBENZENE");

  // G4_DICHLORODIETHYL_ETHER  index=144
  G4double M144[NDENSARRAY] = {
    22.894, 1.903, 4.0135, 0.1773, 3.1586, 0.06799, 3.525, 0, 0.026, 103.3};
  AddMaterial(M144, "G4_DICHLORODIETHYL_ETHER");

  // G4_1,2-DICHLOROETHANE  index=145
  G4double M145[NDENSARRAY] = {
    22.764, 1.618, 4.1849, 0.1375, 2.9529, 0.13383, 3.1675, 0, 0.03, 111.9};
  AddMaterial(M145, "G4_1,2-DICHLOROETHANE");

  // G4_DIETHYL_ETHER  index=146
  G4double M146[NDENSARRAY] = {
    18.326, 1.951, 3.3721, 0.2231, 2.6745, 0.1055, 3.4586, 0, 0.07, 60.0};
  AddMaterial(M146, "G4_DIETHYL_ETHER");

  // G4_N,N-DIMETHYL_FORMAMIDE  index=147
  G4double M147[NDENSARRAY] = {
    20.763, 2.005, 3.3311, 0.1977, 2.6686, 0.1147, 3.371, 0, 0.065, 66.6};
  AddMaterial(M147, "G4_N,N-DIMETHYL_FORMAMIDE");

  // G4_DIMETHYL_SULFOXIDE  index=148
  G4double M148[NDENSARRAY] = {
    22.173, 2.075, 3.9844, 0.2021, 3.1263, 0.06619, 3.5708, 0, 0.03, 98.6};
  AddMaterial(M148, "G4_DIMETHYL_SULFOXIDE");

  // G4_ETHANE  index=149
  G4double M149[NDENSARRAY] = {
    0.789, 1.69, 9.1043, 1.5107, 3.8743, 0.09627, 3.6095, 0, 0.097, 45.4};
  AddMaterial(M149, "G4_ETHANE");

  // G4_ETHYL_ALCOHOL  index=150
  G4double M150[NDENSARRAY] = {
    19.232, 2.013, 3.3699, 0.2218, 2.7052, 0.09878, 3.4834, 0, 0.071, 62.9};
  AddMaterial(M150, "G4_ETHYL_ALCOHOL");

  // G4_ETHYL_CELLULOSE  index=151
  G4double M151[NDENSARRAY] = {
    22.594, 2.065, 3.2415, 0.1683, 2.6527, 0.11077, 3.4098, 0, 0.057, 69.3};
  AddMaterial(M151, "G4_ETHYL_CELLULOSE");

  // G4_ETHYLENE  index=152
  G4double M152[NDENSARRAY] = {
    0.746, 1.733, 9.438, 1.5528, 3.9327, 0.10636, 3.5387, 0, 0.085, 50.7};
  AddMaterial(M152, "G4_ETHYLENE");

  // G4_EYE_LENS_ICRP  index=153
  G4double M153[NDENSARRAY] = {22.388, 2.154, 3.372, 0.207, 2.7446, 0.0969, 3.455, 0, 0.077, 73.3};
  AddMaterial(M153, "G4_EYE_LENS_ICRP");

  // G4_FERRIC_OXIDE  index=154
  G4double M154[NDENSARRAY] = {
    45.331, 2.747, 4.2245, -0.0074, 3.2573, 0.10478, 3.1313, 0, 0.026, 227.3};
  AddMaterial(M154, "G4_FERRIC_OXIDE");

  // G4_FERROBORIDE  index=155
  G4double M155[NDENSARRAY] = {
    52.546, 2.726, 4.2057, -0.0988, 3.1749, 0.12911, 3.024, 0, 0.022, 261.0};
  AddMaterial(M155, "G4_FERROBORIDE");

  // G4_FERROUS_OXIDE  index=156
  G4double M156[NDENSARRAY] = {
    47.327, 2.769, 4.3175, -0.0279, 3.2002, 0.12959, 3.0168, 0, 0.022, 248.6};
  AddMaterial(M156, "G4_FERROUS_OXIDE");

  // G4_FERROUS_SULFATE  index=157
  G4double M157[NDENSARRAY] = {
    21.69, 2.208, 3.5183, 0.2378, 2.8254, 0.08759, 3.4923, 0, 0.096, 76.4};
  AddMaterial(M157, "G4_FERROUS_SULFATE");

  // G4_FREON-12  index=158
  G4double M158[NDENSARRAY] = {
    21.121, 1.974, 4.8251, 0.3035, 3.2659, 0.07978, 3.4626, 0, 0.025, 143.0};
  AddMaterial(M158, "G4_FREON-12");

  // G4_FREON-12B2  index=159
  G4double M159[NDENSARRAY] = {
    25.877, 2.195, 5.7976, 0.3406, 3.7956, 0.05144, 3.5565, 0, 0.021, 284.9};
  AddMaterial(M159, "G4_FREON-12B2");

  // G4_FREON-13  index=160
  G4double M160[NDENSARRAY] = {
    19.432, 2.116, 4.7483, 0.3659, 3.2337, 0.07238, 3.5551, 0, 0.05, 126.6};
  AddMaterial(M160, "G4_FREON-13");

  // G4_FREON-13B1  index=161
  G4double M161[NDENSARRAY] = {
    23.849, 2.233, 5.3555, 0.3522, 3.7554, 0.03925, 3.7194, 0, 0.036, 210.5};
  AddMaterial(M161, "G4_FREON-13B1");

  // G4_FREON-13I1  index=162
  G4double M162[NDENSARRAY] = {
    25.615, 1.924, 5.8774, 0.2847, 3.728, 0.09112, 3.1658, 0, 0.025, 293.5};
  AddMaterial(M162, "G4_FREON-13I1");

  // G4_GADOLINIUM_OXYSULFIDE  index=163
  G4double M163[NDENSARRAY] = {
    51.099, 2.179, 5.5347, -0.1774, 3.4045, 0.22161, 2.63, 0, 0.056, 493.3};
  AddMaterial(M163, "G4_GADOLINIUM_OXYSULFIDE");

  // G4_GALLIUM_ARSENIDE  index=164
  G4double M164[NDENSARRAY] = {
    44.17, 2.652, 5.3299, 0.1764, 3.642, 0.07152, 3.3356, 0, 0.027, 384.9};
  AddMaterial(M164, "G4_GALLIUM_ARSENIDE");

  // G4_GEL_PHOTO_EMULSION  index=165
  G4double M165[NDENSARRAY] = {
    24.058, 2.156, 3.2687, 0.1709, 2.7058, 0.10102, 3.4418, 0, 0.06, 74.8};
  AddMaterial(M165, "G4_GEL_PHOTO_EMULSION");

  // G4_Pyrex_Glass  index=166
  G4double M166[NDENSARRAY] = {
    30.339, 2.369, 3.9708, 0.1479, 2.9933, 0.0827, 3.5224, 0, 0.022, 134.0};
  AddMaterial(M166, "G4_Pyrex_Glass");

  // G4_GLASS_LEAD  index=167
  G4double M167[NDENSARRAY] = {
    46.631, 2.085, 5.8476, 0.0614, 3.8146, 0.09544, 3.074, 0, 0.025, 526.4};
  AddMaterial(M167, "G4_GLASS_LEAD");

  // G4_GLASS_PLATE  index=168
  G4double M168[NDENSARRAY] = {
    31.481, 2.329, 4.0602, 0.1237, 3.0649, 0.07678, 3.5381, 0, 0.025, 145.4};
  AddMaterial(M168, "G4_GLASS_PLATE");

  // G4_GLUCOSE  index=169
  G4double M169[NDENSARRAY] = {
    26.153, 2.174, 3.1649, 0.1411, 2.67, 0.10783, 3.3946, 0, 0.061, 77.2};
  AddMaterial(M169, "G4_GLUCOSE");

  // G4_GLUTAMINE  index=170
  G4double M170[NDENSARRAY] = {
    25.437, 2.077, 3.1167, 0.1347, 2.6301, 0.11931, 3.3254, 0, 0.055, 73.3};
  AddMaterial(M170, "G4_GLUTAMINE");

  // G4_GLYCEROL  index=171
  G4double M171[NDENSARRAY] = {
    23.846, 2.12, 3.2267, 0.1653, 2.6862, 0.10168, 3.4481, 0, 0.067, 72.6};
  AddMaterial(M171, "G4_GLYCEROL");

  // G4_GUANINE  index=172
  G4double M172[NDENSARRAY] = {
    26.022, 1.97, 3.1171, 0.1163, 2.4296, 0.2053, 3.0186, 0, 0.069, 75.0};
  AddMaterial(M172, "G4_GUANINE");

  // G4_GYPSUM  index=173
  G4double M173[NDENSARRAY] = {
    31.379, 2.187, 3.8382, 0.0995, 3.1206, 0.06949, 3.5134, 0, 0.038, 129.7};
  AddMaterial(M173, "G4_GYPSUM");

  // G4_N-HEPTANE  index=174
  G4double M174[NDENSARRAY] = {
    18.128, 1.848, 3.1978, 0.1928, 2.5706, 0.11255, 3.4885, 0, 0.059, 54.4};
  AddMaterial(M174, "G4_N-HEPTANE");

  // G4_N-HEXANE  index=175
  G4double M175[NDENSARRAY] = {
    17.836, 1.843, 3.2156, 0.1984, 2.5757, 0.11085, 3.5027, 0, 0.061, 54.0};
  AddMaterial(M175, "G4_N-HEXANE");

  // G4_KAPTON  index=176
  G4double M176[NDENSARRAY] = {
    24.586, 2.109, 3.3497, 0.1509, 2.5631, 0.15972, 3.1921, 0, 0.05, 79.6};
  AddMaterial(M176, "G4_KAPTON");

  // G4_LANTHANUM_OXYBROMIDE  index=177
  G4double M177[NDENSARRAY] = {
    47.125, 1.831, 5.4666, -0.035, 3.3288, 0.1783, 2.8457, 0, 0.04, 439.7};
  AddMaterial(M177, "G4_LANTHANUM_OXYBROMIDE");

  // G4_LANTHANUM_OXYSULFIDE  index=178
  G4double M178[NDENSARRAY] = {
    45.394, 1.681, 5.6151, -0.0934, 3.2741, 0.22579, 2.7075, 0, 0.065, 456.2};
  AddMaterial(M178, "G4_LANTHANUM_OXYSULFIDE");

  // G4_LEAD_OXIDE  index=179
  G4double M179[NDENSARRAY] = {
    56.488, 2.012, 6.2162, 0.0356, 3.5456, 0.19645, 2.7299, 0, 0.039, 766.7};
  AddMaterial(M179, "G4_LEAD_OXIDE");

  // G4_LITHIUM_AMIDE  index=180
  G4double M180[NDENSARRAY] = {22.609, 1.74, 2.7961, 0.0198, 2.5152, 0.0874, 3.7534, 0, 0.05, 55.5};
  AddMaterial(M180, "G4_LITHIUM_AMIDE");

  // G4_LITHIUM_CARBONATE index=181
  G4double M181[NDENSARRAY] = {
    29.217, 2.246, 3.2029, 0.0551, 2.6598, 0.09936, 3.5417, 0, 0.062, 87.9};
  AddMaterial(M181, "G4_LITHIUM_CARBONATE");

  // G4_LITHIUM_FLUORIDE  index=182
  G4double M182[NDENSARRAY] = {
    31.815, 2.197, 3.1667, 0.0171, 2.7049, 0.07593, 3.7478, 0, 0.084, 94.0};
  AddMaterial(M182, "G4_LITHIUM_FLUORIDE");

  // G4_LITHIUM_HYDRIDE  index=183
  G4double M183[NDENSARRAY] = {
    18.51, 1.482, 2.358, -0.0988, 1.4515, 0.90567, 2.5849, 0, 0.035, 36.5};
  AddMaterial(M183, "G4_LITHIUM_HYDRIDE");

  // G4_LITHIUM_IODIDE  index=184
  G4double M184[NDENSARRAY] = {
    34.841, 1.706, 6.2671, 0.0892, 3.3702, 0.23274, 2.7146, 0, 0.043, 485.1};
  AddMaterial(M184, "G4_LITHIUM_IODIDE");

  // G4_LITHIUM_OXIDE  index=185
  G4double M185[NDENSARRAY] = {
    27.984, 2.039, 2.934, -0.0511, 2.5874, 0.08035, 3.7878, 0, 0.043, 73.6};
  AddMaterial(M185, "G4_LITHIUM_OXIDE");

  // G4_LITHIUM_TETRABORATE  index=186
  G4double M186[NDENSARRAY] = {
    31.343, 2.36, 3.2093, 0.0737, 2.6502, 0.11075, 3.4389, 0, 0.048, 94.6};
  AddMaterial(M186, "G4_LITHIUM_TETRABORATE");

  // G4_LUNG_ICRP  index=187
  G4double M187[NDENSARRAY] = {
    21.891, 2.184, 3.4708, 0.2261, 2.8001, 0.08588, 3.5353, 0, 0.089, 75.3};
  AddMaterial(M187, "G4_LUNG_ICRP");

  // G4_M3_WAX  index=188
  G4double M188[NDENSARRAY] = {22, 1.975, 3.254, 0.1523, 2.7529, 0.07864, 3.6412, 0, 0.044, 67.9};
  AddMaterial(M188, "G4_M3_WAX");

  // G4_MAGNESIUM_CARBONATE  index=189
  G4double M189[NDENSARRAY] = {
    34.979, 2.388, 3.4319, 0.086, 2.7997, 0.09219, 3.5003, 0, 0.045, 118.0};
  AddMaterial(M189, "G4_MAGNESIUM_CARBONATE");

  // G4_MAGNESIUM_FLUORIDE  index=190
  G4double M190[NDENSARRAY] = {
    34.634, 2.33, 3.7105, 0.1369, 2.863, 0.07934, 3.6485, 0, 0.085, 134.3};
  AddMaterial(M190, "G4_MAGNESIUM_FLUORIDE");

  // G4_MAGNESIUM_OXIDE  index=191
  G4double M191[NDENSARRAY] = {
    38.407, 2.412, 3.6404, 0.0575, 2.858, 0.08313, 3.5968, 0, 0.055, 143.8};
  AddMaterial(M191, "G4_MAGNESIUM_OXIDE");

  // G4_MAGNESIUM_TETRABORATE  index=192
  G4double M192[NDENSARRAY] = {
    32.089, 2.43, 3.4328, 0.1147, 2.7635, 0.09703, 3.4893, 0, 0.044, 108.3};
  AddMaterial(M192, "G4_MAGNESIUM_TETRABORATE");

  // G4_MERCURIC_IODIDE  index=193
  G4double M193[NDENSARRAY] = {
    46.494, 1.892, 6.3787, 0.104, 3.4728, 0.21513, 2.7264, 0, 0.047, 684.5};
  AddMaterial(M193, "G4_MERCURIC_IODIDE");

  // G4_METHANE  index=194
  G4double M194[NDENSARRAY] = {
    0.588, 1.662, 9.5243, 1.6263, 3.9716, 0.09253, 3.6257, 0, 0.112, 41.7};
  AddMaterial(M194, "G4_METHANE");

  // G4_METHANOL  index=195
  G4double M195[NDENSARRAY] = {19.214, 2.125, 3.516, 0.2529, 2.7639, 0.0897, 3.5477, 0, 0.08, 67.6};
  AddMaterial(M195, "G4_METHANOL");

  // G4_MIX_D_WAX  index=196
  G4double M196[NDENSARRAY] = {
    21.547, 1.905, 3.078, 0.1371, 2.7145, 0.0749, 3.6823, 0, 0.047, 60.9};
  AddMaterial(M196, "G4_MIX_D_WAX");

  // G4_MS20_TISSUE  index=197
  G4double M197[NDENSARRAY] = {
    21.153, 2.07, 3.5341, 0.1997, 2.8033, 0.08294, 3.6061, 0, 0.053, 75.1};
  AddMaterial(M197, "G4_MS20_TISSUE");

  // G4_MUSCLE_SCELETAL_ICRP  index=198
  G4double M198[NDENSARRAY] = {
    21.781, 2.185, 3.4809, 0.2282, 2.7999, 0.08636, 3.533, 0, 0.089, 75.3};
  AddMaterial(M198, "G4_MUSCLE_SCELETAL_ICRP");

  // G4_MUSCLE_STRIATED_ICRU  index=199
  G4double M199[NDENSARRAY] = {
    21.795, 2.174, 3.4636, 0.2249, 2.8032, 0.08507, 3.5383, 0, 0.086, 74.7};
  AddMaterial(M199, "G4_MUSCLE_STRIATED_ICRU");

  // G4_MUSCLE_WITH_SUCROSE  index=200
  G4double M200[NDENSARRAY] = {22.48, 2.169, 3.391, 0.2098, 2.755, 0.09481, 3.4699, 0, 0.08, 74.3};
  AddMaterial(M200, "G4_MUSCLE_WITH_SUCROSE");

  // G4_MUSCLE_WITHOUT_SUCROSE  index=201
  G4double M201[NDENSARRAY] = {
    22.109, 2.173, 3.4216, 0.2187, 2.768, 0.09143, 3.4982, 0, 0.056, 74.2};
  AddMaterial(M201, "G4_MUSCLE_WITHOUT_SUCROSE");

  // G4_NAPHTALENE  index=202
  G4double M202[NDENSARRAY] = {
    22.459, 1.956, 3.2274, 0.1374, 2.5429, 0.14766, 3.2654, 0, 0.051, 68.4};
  AddMaterial(M202, "G4_NAPHTALENE");

  // G4_NITROBENZENE  index=203
  G4double M203[NDENSARRAY] = {
    22.747, 2.065, 3.4073, 0.1777, 2.663, 0.12727, 3.3091, 0, 0.051, 75.8};
  AddMaterial(M203, "G4_NITROBENZENE");

  // G4_NITROUS_OXIDE  index=204
  G4double M204[NDENSARRAY] = {
    0.872, 2.059, 10.1575, 1.6477, 4.1565, 0.11992, 3.3318, 0, 0.086, 84.9};
  AddMaterial(M204, "G4_NITROUS_OXIDE");

  // G4_NYLON-8062  index=205
  G4double M205[NDENSARRAY] = {
    22.221, 1.967, 3.125, 0.1503, 2.6004, 0.11513, 3.4044, 0, 0.054, 64.3};
  AddMaterial(M205, "G4_NYLON-8062");

  // G4_NYLON-6/6  index=206
  G4double M206[NDENSARRAY] = {
    22.774, 1.931, 3.0634, 0.1336, 2.5834, 0.11818, 3.3826, 0, 0.051, 63.9};
  AddMaterial(M206, "G4_NYLON-6-6");

  // G4_NYLON-6/10  index=207
  G4double M207[NDENSARRAY] = {
    22.866, 1.942, 3.0333, 0.1304, 2.5681, 0.11852, 3.3912, 0, 0.05, 63.2};
  AddMaterial(M207, "G4_NYLON-6-10");

  // G4_NYLON-11_RILSAN  index=208
  G4double M208[NDENSARRAY] = {
    25.661, 1.902, 2.7514, 0.0678, 2.4281, 0.14868, 3.2576, 0, 0.044, 61.6};
  AddMaterial(M208, "G4_NYLON-11_RILSAN");

  // G4_OCTANE  index=209
  G4double M209[NDENSARRAY] = {
    18.36, 1.851, 3.1834, 0.1882, 2.5664, 0.11387, 3.4776, 0, 0.057, 54.7};
  AddMaterial(M209, "G4_OCTANE");

  // G4_PARAFFIN  index=210
  G4double M210[NDENSARRAY] = {
    21.031, 1.844, 2.9551, 0.1289, 2.5084, 0.12087, 3.4288, 0, 0.052, 55.9};
  AddMaterial(M210, "G4_PARAFFIN");

  // G4_N-PENTANE  index=211
  G4double M211[NDENSARRAY] = {
    17.398, 1.842, 3.2504, 0.2086, 2.5855, 0.10809, 3.5265, 0, 0.064, 53.6};
  AddMaterial(M211, "G4_N-PENTANE");

  // G4_PHOTO_EMULSION  index=212
  G4double M212[NDENSARRAY] = {
    37.946, 2.264, 5.3319, 0.1009, 3.4866, 0.12399, 3.0094, 0, 0.028, 331.0};
  AddMaterial(M212, "G4_PHOTO_EMULSION");

  // G4_PLASTIC_SC_VINYLTOLUENE  index=213
  G4double M213[NDENSARRAY] = {
    21.54, 1.929, 3.1997, 0.1464, 2.4855, 0.16101, 3.2393, 0, 0.05, 64.7};
  AddMaterial(M213, "G4_PLASTIC_SC_VINYLTOLUENE");

  // G4_PLUTONIUM_DIOXIDE  index=214
  G4double M214[NDENSARRAY] = {
    62.143, 1.846, 5.9719, -0.2311, 3.5554, 0.20594, 2.6522, 0, 0.111, 746.5};
  AddMaterial(M214, "G4_PLUTONIUM_DIOXIDE");

  // G4_POLYACRYLONITRILE  index=215
  G4double M215[NDENSARRAY] = {
    22.642, 1.955, 3.2459, 0.1504, 2.5159, 0.16275, 3.1975, 0, 0.05, 69.6};
  AddMaterial(M215, "G4_POLYACRYLONITRILE");

  // G4_POLYCARBONATE  index=216
  G4double M216[NDENSARRAY] = {
    22.915, 2.06, 3.3201, 0.1606, 2.6225, 0.1286, 3.3288, 0, 0.049, 73.1};
  AddMaterial(M216, "G4_POLYCARBONATE");

  // G4_POLYCHLOROSTYRENE  index=217
  G4double M217[NDENSARRAY] = {
    23.81, 1.902, 3.4659, 0.1238, 2.9241, 0.0753, 3.5441, 0, 0.029, 81.7};
  AddMaterial(M217, "G4_POLYCHLOROSTYRENE");

  // G4_POLYETHYLENE  index=218
  G4double M218[NDENSARRAY] = {
    21.099, 1.882, 3.0016, 0.137, 2.5177, 0.12108, 3.4292, 0, 0.051, 57.4};
  AddMaterial(M218, "G4_POLYETHYLENE");

  // G4_MYLAR  index=219
  G4double M219[NDENSARRAY] = {
    24.595, 2.144, 3.3262, 0.1562, 2.6507, 0.12679, 3.3076, 0, 0.052, 78.7};
  AddMaterial(M219, "G4_MYLAR");

  // G4_LUCITE  index=220
  G4double M220[NDENSARRAY] = {
    23.086, 2.173, 3.3297, 0.1824, 2.6681, 0.11433, 3.3836, 0, 0.056, 74.0};
  AddMaterial(M220, "G4_LUCITE");

  // G4_POLYOXOMETHYLENE  index=221
  G4double M221[NDENSARRAY] = {
    25.11, 2.175, 3.2514, 0.1584, 2.6838, 0.10808, 3.4002, 0, 0.063, 77.4};
  AddMaterial(M221, "G4_POLYOXOMETHYLENE");

  // G4_POLYPROPYLENE  index=222
  G4double M222[NDENSARRAY] = {
    20.457, 1.884, 3.1252, 0.1534, 2.4822, 0.15045, 3.2855, 0, 0.055, 59.2};
  AddMaterial(M222, "G4_POLYPROPYLENE");

  // G4_POLYSTYRENE  index=223
  G4double M223[NDENSARRAY] = {
    21.754, 2.027, 3.2999, 0.1647, 2.5031, 0.16454, 3.2224, 0, 0.051, 68.7};
  AddMaterial(M223, "G4_POLYSTYRENE");

  // G4_TEFLON  index=224
  G4double M224[NDENSARRAY] = {
    29.609, 2.142, 3.4161, 0.1648, 2.7404, 0.10606, 3.4046, 0, 0.073, 99.1};
  AddMaterial(M224, "G4_TEFLON");

  // G4_POLYTRIFLUOROCHLOROETHYLENE  index=225
  G4double M225[NDENSARRAY] = {
    28.955, 2.094, 3.8551, 0.1714, 3.0265, 0.07727, 3.5085, 0, 0.035, 120.7};
  AddMaterial(M225, "G4_POLYTRIFLUOROCHLOROETHYLENE");

  // G4_POLYVINYL_ACETATE  index=226
  G4double M226[NDENSARRAY] = {
    22.978, 2.116, 3.3309, 0.1769, 2.6747, 0.11442, 3.3762, 0, 0.055, 73.7};
  AddMaterial(M226, "G4_POLYVINYL_ACETATE");

  // G4_PLOYVINYL_ALCOHOL  index=227
  G4double M227[NDENSARRAY] = {
    24.251, 2.071, 3.1115, 0.1401, 2.6315, 0.11178, 3.3893, 0, 0.056, 69.7};
  AddMaterial(M227, "G4_PLOYVINYL_ALCOHOL");

  // G4_POLYVINYL_BUTYRAL  index=228
  G4double M228[NDENSARRAY] = {
    22.521, 2.021, 3.1865, 0.1555, 2.6186, 0.11544, 3.3983, 0, 0.054, 67.2};
  AddMaterial(M228, "G4_POLYVINYL_BUTYRAL");

  // G4_POLYVINYL_CHLORIDE  index=229
  G4double M229[NDENSARRAY] = {
    23.51, 1.84, 4.0532, 0.1559, 2.9415, 0.12438, 3.2104, 0, 0.027, 108.2};
  AddMaterial(M229, "G4_POLYVINYL_CHLORIDE");

  // G4_POLYVINYLIDENE_CHLORIDE  index=230
  G4double M230[NDENSARRAY] = {
    26.437, 1.814, 4.2506, 0.1314, 2.9009, 0.15466, 3.102, 0, 0.034, 134.3};
  AddMaterial(M230, "G4_POLYVINYLIDENE_CHLORIDE");

  // G4_POLYVINYLIDENE_FLUORIDE  index=231
  G4double M231[NDENSARRAY] = {27.024, 2.16, 3.3793, 0.1717, 2.7375, 0.10316, 3.42, 0, 0.067, 88.8};
  AddMaterial(M231, "G4_POLYVINYLIDENE_FLUORIDE");

  // G4_POLYVINYL_PYRROLIDONE  index=232
  G4double M232[NDENSARRAY] = {
    23.671, 1.989, 3.1017, 0.1324, 2.5867, 0.12504, 3.3326, 0, 0.031, 67.7};
  AddMaterial(M232, "G4_POLYVINYL_PYRROLIDONE");

  // G4_POTASSIUM_IODIDE  index=233
  G4double M233[NDENSARRAY] = {
    33.575, 1.784, 6.1088, 0.1044, 3.3442, 0.22053, 2.7558, 0, 0.042, 431.9};
  AddMaterial(M233, "G4_POTASSIUM_IODIDE");

  // G4_POTASSIUM_OXIDE  index=234
  G4double M234[NDENSARRAY] = {
    30.672, 2.065, 4.6463, 0.048, 3.011, 0.16789, 3.0121, 0, 0.027, 189.9};
  AddMaterial(M234, "G4_POTASSIUM_OXIDE");

  // G4_PROPANE  index=235
  G4double M235[NDENSARRAY] = {
    0.959, 1.708, 8.7878, 1.4326, 3.7998, 0.09916, 3.592, 0, 0.093, 47.1};
  AddMaterial(M235, "G4_PROPANE");

  // G4_lPROPANE  index=236
  G4double M236[NDENSARRAY] = {
    14.509, 1.844, 3.5529, 0.2861, 2.6568, 0.10329, 3.562, 0, 0.068, 52.0};
  AddMaterial(M236, "G4_lPROPANE");

  // G4_N-PROPYL_ALCOHOL index=237
  G4double M237[NDENSARRAY] = {
    19.429, 1.972, 3.2915, 0.2046, 2.6681, 0.09644, 3.5415, 0, 0.07, 61.1};
  AddMaterial(M237, "N-PROPYL_ALCOHOL");

  // G4_PYRIDINE  index=238
  G4double M238[NDENSARRAY] = {
    20.807, 1.895, 3.3148, 0.167, 2.5245, 0.16399, 3.1977, 0, 0.051, 66.2};
  AddMaterial(M238, "G4_PYRIDINE");

  // G4_RUBBER_BUTYL  index=239
  G4double M239[NDENSARRAY] = {
    20.873, 1.852, 2.9915, 0.1347, 2.5154, 0.12108, 3.4296, 0, 0.051, 56.5};
  AddMaterial(M239, "G4_RUBBER_BUTYL");

  // G4_RUBBER_NATURAL  index=240
  G4double M240[NDENSARRAY] = {
    20.644, 1.889, 3.1272, 0.1512, 2.4815, 0.15058, 3.2879, 0, 0.053, 59.8};
  AddMaterial(M240, "G4_RUBBER_NATURAL");

  // G4_RUBBER_NEOPRENE  index=241
  G4double M241[NDENSARRAY] = {
    23.036, 1.874, 3.7911, 0.1501, 2.9461, 0.09763, 3.3632, 0, 0.026, 93.0};
  AddMaterial(M241, "G4_RUBBER_NEOPRENE");

  // G4_SILICON_DIOXIDE  index=242
  G4double M242[NDENSARRAY] = {
    31.014, 2.335, 4.0029, 0.1385, 3.0025, 0.08408, 3.5064, 0, 0.018, 139.2};
  AddMaterial(M242, "G4_SILICON_DIOXIDE");

  // G4_SILVER_BROMIDE  index=243
  G4double M243[NDENSARRAY] = {
    48.448, 2.271, 5.6139, 0.0352, 3.2109, 0.24582, 2.682, 0, 0.043, 486.6};
  AddMaterial(M243, "G4_SILVER_BROMIDE");

  // G4_SILVER_CHLORIDE  index=244
  G4double M244[NDENSARRAY] = {
    45.405, 2.096, 5.3437, -0.0139, 3.2022, 0.22968, 2.7041, 0, 0.062, 398.4};
  AddMaterial(M244, "G4_SILVER_CHLORIDE");

  // G4_SILVER_HALIDES  index=245
  G4double M245[NDENSARRAY] = {
    48.433, 2.27, 5.6166, 0.0353, 3.2117, 0.24593, 2.6814, 0, 0.043, 487.1};
  AddMaterial(M245, "G4_SILVER_HALIDES");

  // G4_SILVER_IODIDE  index=246
  G4double M246[NDENSARRAY] = {
    46.105, 1.945, 5.9342, 0.0148, 3.2908, 0.25059, 2.6572, 0, 0.071, 543.5};
  AddMaterial(M246, "G4_SILVER_IODIDE");

  // G4_SKIN_ICRP  index=247
  G4double M247[NDENSARRAY] = {22.4, 2.14, 3.3546, 0.2019, 2.7526, 0.09459, 3.4643, 0, 0.076, 72.7};
  AddMaterial(M247, "G4_SKIN_ICRP");

  // G4_SODIUM_CARBONATE  index=248
  G4double M248[NDENSARRAY] = {
    32.117, 2.557, 3.7178, 0.1287, 2.8591, 0.08715, 3.5638, 0, 0.074, 125.0};
  AddMaterial(M248, "G4_SODIUM_CARBONATE");

  // G4_SODIUM_IODIDE  index=249
  G4double M249[NDENSARRAY] = {
    36.057, 1.857, 6.0572, 0.1203, 3.592, 0.12516, 3.0398, 0, 0.031, 452.0};
  AddMaterial(M249, "G4_SODIUM_IODIDE");

  // G4_SODIUM_MONOXIDE  index=250
  G4double M250[NDENSARRAY] = {
    30.205, 2.689, 4.1892, 0.1652, 2.9793, 0.07501, 3.6943, 0, 0.097, 148.8};
  AddMaterial(M250, "G4_SODIUM_MONOXIDE");

  // G4_SODIUM_NITRATE  index=251
  G4double M251[NDENSARRAY] = {
    30.459, 2.456, 3.6502, 0.1534, 2.8221, 0.09391, 3.5097, 0, 0.081, 114.6};
  AddMaterial(M251, "G4_SODIUM_NITRATE");

  // G4_STILBENE  index=252
  G4double M252[NDENSARRAY] = {
    20.719, 1.963, 3.368, 0.1734, 2.5142, 0.16659, 3.2168, 0, 0.052, 67.7};
  AddMaterial(M252, "G4_STILBENE");

  // G4_SUCROSE  index=253
  G4double M253[NDENSARRAY] = {
    26.416, 2.167, 3.1526, 0.1341, 2.6558, 0.11301, 3.363, 0, 0.057, 77.5};
  AddMaterial(M253, "G4_SUCROSE");

  // G4_TERPHENYL  index=254
  G4double M254[NDENSARRAY] = {
    23.116, 1.976, 3.2639, 0.1322, 2.5429, 0.14964, 3.2685, 0, 0.043, 71.7};
  AddMaterial(M254, "G4_TERPHENYL");

  // G4_TESTES_ICRP  index=255
  G4double M255[NDENSARRAY] = {
    21.815, 2.185, 3.4698, 0.2274, 2.7988, 0.08533, 3.5428, 0, 0.091, 75.0};
  AddMaterial(M255, "G4_TESTES_ICRP");

  // G4_TETRACHLOROETHYLENE  index=256
  G4double M256[NDENSARRAY] = {
    25.513, 1.79, 4.6619, 0.1713, 2.9083, 0.18595, 3.0156, 0, 0.038, 159.2};
  AddMaterial(M256, "G4_TETRACHLOROETHYLENE");

  // G4_THALIUM_CHLORIDE  index=257
  G4double M257[NDENSARRAY] = {
    48.749, 1.997, 6.3009, 0.0705, 3.5716, 0.18599, 2.769, 0, 0.04, 690.3};
  AddMaterial(M257, "G4_THALIUM_CHLORIDE");

  // G4_TISSUE_SOFT_ICRP  index=258
  G4double M258[NDENSARRAY] = {
    21.394, 2.144, 3.4354, 0.2211, 2.7799, 0.08926, 3.511, 0, 0.077, 72.3};
  AddMaterial(M258, "G4_TISSUE_SOFT_ICRP");

  // G4_TISSUE_SOFT_ICRU-4  index=259
  G4double M259[NDENSARRAY] = {
    21.366, 2.192, 3.5087, 0.2377, 2.7908, 0.09629, 3.4371, 0, 0.092, 74.9};
  AddMaterial(M259, "G4_TISSUE_SOFT_ICRU-4");

  // G4_TISSUE-METHANE  index=260
  G4double M260[NDENSARRAY] = {0.697, 1.89, 9.95, 1.6442, 4.1399, 0.09946, 3.4708, 0, 0.098, 61.2};
  AddMaterial(M260, "G4_TISSUE-METHANE");

  // G4_TISSUE-PROPANE  index=261
  G4double M261[NDENSARRAY] = {
    0.913, 1.856, 9.3529, 1.5139, 3.9916, 0.09802, 3.5159, 0, 0.092, 59.5};
  AddMaterial(M261, "G4_TISSUE-PROPANE");

  // G4_TITANIUM_DIOXIDE  index=262
  G4double M262[NDENSARRAY] = {
    41.022, 2.307, 3.9522, -0.0119, 3.1647, 0.08569, 3.3267, 0, 0.027, 179.5};
  AddMaterial(M262, "G4_TITANIUM_DIOXIDE");

  // G4_TOLUENE  index=263
  G4double M263[NDENSARRAY] = {
    19.764, 1.88, 3.3026, 0.1722, 2.5728, 0.13284, 3.3558, 0, 0.052, 62.5};
  AddMaterial(M263, "G4_TOLUENE");

  // G4_TRICHLOROETHYLENE  index=264
  G4double M264[NDENSARRAY] = {
    24.301, 1.789, 4.6148, 0.1803, 2.914, 0.18272, 3.0137, 0, 0.036, 148.1};
  AddMaterial(M264, "G4_TRICHLOROETHYLENE");

  // G4_TRIETHYL_PHOSPHATE  index=265
  G4double M265[NDENSARRAY] = {
    21.863, 2.1, 3.6242, 0.2054, 2.9428, 0.06922, 3.6302, 0, 0.049, 81.2};
  AddMaterial(M265, "G4_TRIETHYL_PHOSPHATE");

  // G4_TUNGSTEN_HEXAFLUORIDE  index=266
  G4double M266[NDENSARRAY] = {
    29.265, 2.325, 5.9881, 0.302, 4.2602, 0.03658, 3.5134, 0, 0.055, 354.4};
  AddMaterial(M266, "G4_TUNGSTEN_HEXAFLUORIDE");

  // G4_URANIUM_DICARBIDE  index=267
  G4double M267[NDENSARRAY] = {
    60.969, 1.703, 6.0247, -0.2191, 3.5208, 0.2112, 2.6577, 0, 0.12, 752.0};
  AddMaterial(M267, "G4_URANIUM_DICARBIDE");

  // G4_URANIUM_MONOCARBIDE  index=268
  G4double M268[NDENSARRAY] = {
    66.602, 1.68, 6.121, -0.2524, 3.4941, 0.22972, 2.6169, 0, 0.132, 862.0};
  AddMaterial(M268, "G4_URANIUM_MONOCARBIDE");

  // G4_URANIUM_OXIDE  index=269
  G4double M269[NDENSARRAY] = {
    60.332, 1.76, 5.9605, -0.1938, 3.5292, 0.20463, 2.6711, 0, 0.098, 720.6};
  AddMaterial(M269, "G4_URANIUM_OXIDE");

  // G4_UREA  index=270
  G4double M270[NDENSARRAY] = {
    24.194, 2.022, 3.2032, 0.1603, 2.6525, 0.11609, 3.3461, 0, 0.06, 72.8};
  AddMaterial(M270, "G4_UREA");

  // G4_VALINE  index=271
  G4double M271[NDENSARRAY] = {
    23.622, 2.024, 3.1059, 0.1441, 2.6227, 0.11386, 3.3774, 0, 0.056, 67.7};
  AddMaterial(M271, "G4_VALINE");

  // G4_VITON  index=272
  G4double M272[NDENSARRAY] = {
    26.948, 2.227, 3.5943, 0.2106, 2.7874, 0.09965, 3.4556, 0, 0.07, 98.6};
  AddMaterial(M272, "G4_VITON");

  // G4_WATER  index=273
  G4double M273[NDENSARRAY] = {
    21.469, 2.203, 3.5017, 0.24, 2.8004, 0.09116, 3.4773, 0, 0.097, 75.0};
  AddMaterial(M273, "G4_WATER");

  // G4_WATER_VAPOR  index=274
  G4double M274[NDENSARRAY] = {
    0.59, 2.175, 10.5962, 1.7952, 4.3437, 0.08101, 3.5901, 0, 0.121, 71.6};
  AddMaterial(M274, "G4_WATER_VAPOR");

  // G4_XYLENE  index=275
  G4double M275[NDENSARRAY] = {
    19.866, 1.882, 3.2698, 0.1695, 2.5675, 0.13216, 3.3564, 0, 0.051, 61.8};
  AddMaterial(M275, "G4_XYLENE");

  // G4_GRAPHITE  index=276
  G4double M276[NDENSARRAY] = {
    30.652, 2.29, 2.868, -0.0178, 2.3415, 0.26142, 2.8697, 0.12, 0.038, 78.0};
  AddMaterial(M276, "G4_GRAPHITE");

  // G4_GRAPHITE_POROUS  index=277
  G4double M277[NDENSARRAY] = {
    26.555, 2.49, 3.155, 0.048, 2.5387, 0.20762, 2.9532, 0.14, 0.038, 78.0};
  AddMaterial(M277, "G4_GRAPHITE_POROUS");
}

G4int G4DensityEffectData::GetElementIndex(G4int Z, G4State) const
{
  return (Z >= 0 && Z < NDENSELEM) ? indexZ[Z] : -1;
}

G4int G4DensityEffectData::GetIndex(const G4String& matName) const
{
  G4int idx = -1;

  for (G4int i = 0; i < NDENSDATA; ++i) {
    if (names[i] == matName) {
      idx = i;
      break;
    }
  }
  return idx;
}

void G4DensityEffectData::AddMaterial(G4double* val, const G4String& matName)
{
  for (G4int i = 0; i < NDENSARRAY; ++i) {
    data[index][i] = val[i];
  }
  data[index][0] *= CLHEP::eV;
  data[index][9] *= CLHEP::eV;
  names.push_back(matName);
  ++index;
}

void G4DensityEffectData::PrintData(const G4String& matName) const
{
  if (matName.empty() || "all" == matName) {
    DumpData();
    return;
  }
  G4int idx = GetIndex(matName);
  if (idx >= 0) {
    G4cout << "G4DensityEffectData for <" << matName << "> index= " << idx << G4endl;
    G4cout << "I(eV)= " << data[idx][9] / CLHEP::eV << "Eplasma(eV)= " << data[idx][0] / CLHEP::eV
           << " rho= " << data[idx][1] << " -C= " << data[idx][2] << " x0= " << data[idx][3]
           << " x1= " << data[idx][4] << " a= " << data[idx][5] << " m= " << data[idx][6]
           << " d0= " << data[idx][7] << " err= " << data[idx][8] << G4endl;
  }
  else {
    G4cout << "G4DensityEffectData does not have <" << matName << ">" << G4endl;
  }
}

void G4DensityEffectData::DumpData() const
{
  G4cout << "======================================================================" << G4endl;
  G4cout << "     Material        Eplasma(eV)  rho  -C   x0   x1   a   m  d0  err" << G4endl;
  G4cout << "======================================================================" << G4endl;
  for (G4int i = 0; i < NDENSDATA; ++i) {
    G4cout << std::setw(3) << i << ". " << std::setw(25) << names[i] << std::setw(8)
           << data[i][0] / eV;
    for (G4int j = 1; j < NDENSARRAY; ++j) {
      G4cout << std::setw(8) << data[i][j];
    }
    G4cout << G4endl;
  }
  G4cout << "======================================================================" << G4endl;
}
