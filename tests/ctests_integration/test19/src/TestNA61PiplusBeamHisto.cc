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

#include "TestNA61PiplusBeamHisto.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

TestNA61PiplusBeamHisto::TestNA61PiplusBeamHisto( G4String htitle )
   : TstHistoSet(htitle)
{

  G4String title;
  G4String outcome; 
    
  // comparison vs recent data from Pub-2019
  //

  // The momentum binning per theta-bin needs to be properly exported
  // from https://www.hepdata.net/record/ins1754136
  // For now, just use something generic, starting from 0. and up to 40.
/*  
  // secondary proton 
  //
  double tmp_p[] = {1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 
                    5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 
		    9.2, 9.6, 10, 10.8, 11.6, 12.4, 13.2, 14, 14.8, 15.6, 
		    16.4, 17.2, 18, 18.8, 19.6, 20.4, 21.2, 22};

  // 0-20
  fNProtonBinsP2019[0] = 30;
  fProtonBinsP2019[0] = new double[fNProtonBinsP2019[0]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[0]; ++i )
  {
     fProtonBinsP2019[0][i] = tmp_p[i];
  }
  // 20-40
  fNProtonBinsP2019[2] = 37;
  fProtonBinsP2019[2] = new double[fNProtonBinsP2019[2]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[2]; ++i )
  {
     fProtonBinsP2019[2][i] = tmp_p[i];
  }
  // 40-60
  fNProtonBinsP2019[3] = 35;
  fProtonBinsP2019[3] = new double[fNProtonBinsP2019[3]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[3]; ++i )
  {
     fProtonBinsP2019[3][i] = tmp_p[i];
  }
  // 60-100
  fNProtonBinsP2019[4] = 32;
  fProtonBinsP2019[4] = new double[fNProtonBinsP2019[4]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[4]; ++i )
  {
     fProtonBinsP2019[4][i] = tmp_p[i];
  }
  // 100-140
  fNProtonBinsP2019[5] = 26;
  fProtonBinsP2019[5] = new double[fNProtonBinsP2019[5]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[5]; ++i )
  {
     fProtonBinsP2019[5][i] = tmp_p[i];
  }
  // 140-180
  fNProtonBinsP2019[6] = 22;
  fProtonBinsP2019[6] = new double[fNProtonBinsP2019[6]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[6]; ++i )
  {
     fProtonBinsP2019[6][i] = tmp_p[i];
  }
  // 180-240
  fNProtonBinsP2019[7] = 16;
  fProtonBinsP2019[7] = new double[fNProtonBinsP2019[7]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[7]; ++i )
  {
     fProtonBinsP2019[7][i] = tmp_p[i];
  }
  // 240-300
  fNProtonBinsP2019[8] = 10;
  fProtonBinsP2019[8] = new double[fNProtonBinsP2019[8]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[8]; ++i )
  {
     fProtonBinsP2019[8][i] = tmp_p[i];
  }
  // 300-360
  fNProtonBinsP2019[9] = 6;
  fProtonBinsP2019[9] = new double[fNProtonBinsP2019[9]+1]; 
  for ( int i=0; i<=fNProtonBinsP2019[9]; ++i )
  {
     fProtonBinsP2019[9][i] = tmp_p[i];
  }
*/
  outcome = " -> X + proton ";

  title = htitle + outcome + " (0<theta<20 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_0_20", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[0], fProtonBinsP2019[0] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_20_40", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[2], fProtonBinsP2019[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_40_60", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[3], fProtonBinsP2019[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_60_100", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[4], fProtonBinsP2019[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_100_140", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[5], fProtonBinsP2019[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_140_180", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[6], fProtonBinsP2019[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_180_240", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[7], fProtonBinsP2019[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_240_300", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[8], fProtonBinsP2019[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_300_360", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[9], fProtonBinsP2019[9] ) );
  title = htitle + outcome + " (360<theta<420 [mrad]) ";
  fHistoProtonP2019.push_back( new TH1F( "proton_360_420", title.c_str(), 100, 0., 40. ) ); //fNProtonBinsP2019[10], fProtonBinsP2019[10] ) );

  // secondary charged pions
  // NOTE: sometimes there's -1 or -2 offset and/or comment about "pi+ beam only"
  //       those are left over from when we "hacked" from analysis of the 31GeV/c p+C,
  //       and certain things were kind of common for both pi+ or p beam
  //       going forward, it needs to be checked more carefully and cleaned up
  
  // secondary pi+
  //

  double tmp_pi_0[] = {      
      0.4,
      1.5,
      3.0,
      4.0,
      5.0,
      6.0,
      7.0,
      8.0,
      9.0,
      10.0,
      11.0,
      12.0,
      13.0,
      14.0,
      15.0,
      17.0,
      19.0,
      21.0,
      23.0,
      25.0,
      27.0,
      29.0,
      31.0, 
      33.0 };

  double tmp_pi_1[] = {       
      0.4,
      1.5,
      2.2,
      2.6,
      3.0,
      3.4,
      3.8,
      4.2,
      4.6,
      5.0,
      5.4,
      5.8,
      6.2,
      6.6,
      7.0,
      7.8,
      8.6,
      9.4,
      10.2,
      11.0,
      11.8,
      12.6,
      13.4,
      14.2,
      15.0,
      17.0,
      19.0,
      21.0,
      23.0,
      25.0,
      27.0,
      29.0,
      31.0,
      33.0,
      35.0,
      38.0,
      41.0,
      44.0,
      48.0 };

   double tmp_pi_2[] = {       
      0.4,
      0.8,
      1.5,
      2.2,
      2.6,
      3.0,
      3.4,
      3.8,
      4.2,
      4.6,
      5.0,
      5.4,
      5.8,
      6.2,
      6.6,
      7.0,
      7.8,
      8.6,
      9.4,
      10.2,
      11.0,
      11.8,
      12.6,
      13.4,
      14.2,
      15.0,
      17.0,
      19.0,
      21.0,
      23.0,
      25.0,
      29.0,
      35.0 };


  
  // 0-10
  fNPiplusBinsP2019[0] = 23;  
  fPiplusBinsP2019[0] = new double[fNPiplusBinsP2019[0]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[0]; ++i )
  {
     fPiplusBinsP2019[0][i] = tmp_pi_0[i]; 
  }
  // 10-20
  fNPiplusBinsP2019[1] = 38;  
  fPiplusBinsP2019[1] = new double[fNPiplusBinsP2019[1]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[1]; ++i )  
  {
     fPiplusBinsP2019[1][i] = tmp_pi_1[i]; 
  }
  // 20-40
  fNPiplusBinsP2019[2] = 32; 
  fPiplusBinsP2019[2] = new double[fNPiplusBinsP2019[2]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[2]; ++i )   
  {
     fPiplusBinsP2019[2][i] = tmp_pi_2[i];
  }
  // 40-60
  fNPiplusBinsP2019[3] = 28; 
  fPiplusBinsP2019[3] = new double[fNPiplusBinsP2019[3]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[3]-1; ++i ) 
  {
     fPiplusBinsP2019[3][i] = tmp_pi_2[i]; 
  }
  fPiplusBinsP2019[3][fNPiplusBinsP2019[3]] = 25.0; 
  // 60-100
  fNPiplusBinsP2019[4] = 22; 
  fPiplusBinsP2019[4] = new double[fNPiplusBinsP2019[4]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[4]-2; ++i )   
  {
     fPiplusBinsP2019[4][i] = tmp_pi_2[i]; 
  }
  fPiplusBinsP2019[4][fNPiplusBinsP2019[4]-1] = 13.0; 
  fPiplusBinsP2019[4][fNPiplusBinsP2019[4]] = 16.0; 
  
  double tmp_pi_3[] = {       
      0.4,
      0.8, 
      1.5,
      2.2,
      2.6,
      3.0,
      3.4,
      3.8,
      4.2,
      4.6,
      5.0,
      5.4,
      5.8,
      6.2,
      8.0,
      12.0 };

    // 100-140
  fNPiplusBinsP2019[5] = 15; 
  fPiplusBinsP2019[5] = new double[fNPiplusBinsP2019[5]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[5]; ++i )
  {
     fPiplusBinsP2019[5][i] = tmp_pi_3[i]; 
  }
  
  double tmp_pi_4[] = { 
      1.5,
      2.2,
      2.6,
      3.0,
      3.4,
      3.8,
      4.2,
      4.6,
      5.0,
      6.0,
      9.0 };
      
  // 140-180
  fNPiplusBinsP2019[6] = 10; 
  fPiplusBinsP2019[6] = new double[fNPiplusBinsP2019[6]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[6]; ++i )
  {
     fPiplusBinsP2019[6][i] = tmp_pi_4[i]; 
  }

  double tmp_pi_5[] ={ 
      0.4,
      0.8,
      1.5,
      2.2,
      2.6,
      3.0,
      3.4,
      3.8,
      4.2,
      5.0,
      7.0 };
      
  // 180-240
  fNPiplusBinsP2019[7] = 10; 
  fPiplusBinsP2019[7] = new double[fNPiplusBinsP2019[7]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[7]; ++i )
  {
     fPiplusBinsP2019[7][i] = tmp_pi_5[i]; 
  }
  
  double tmp_pi_6[] = {       
      0.4,
      0.8,
      1.5,
      2.2,
      2.6,
      3.0,
      3.8,
      5.0 };
  
  // 240-300
  fNPiplusBinsP2019[8] = 7; 
  fPiplusBinsP2019[8] = new double[fNPiplusBinsP2019[8]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[8]; ++i )
  {
     fPiplusBinsP2019[8][i] = tmp_pi_6[i]; 
  }
  
  double tmp_pi_7[] = { 
      0.4,
      0.8,
      1.5,
      2.2,
      2.6,
      3.8 };
  
  // 300-360
  fNPiplusBinsP2019[9] = 5; 
  fPiplusBinsP2019[9] = new double[fNPiplusBinsP2019[9]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[9]; ++i )
  {
     fPiplusBinsP2019[9][i] = tmp_pi_7[i]; 
  }
  
  double tmp_pi_8[] = { 
      0.4,
      0.8,
      1.5,
      2.2,
      2.6,
      3.4 };
  
  // 360-420
  fNPiplusBinsP2019[10] = 5; 
  fPiplusBinsP2019[10] = new double[fNPiplusBinsP2019[10]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2019[10]; ++i )
  {
     fPiplusBinsP2019[10][i] = tmp_pi_8[i]; 
  }

  outcome = " -> X + pi+ ";

  title = htitle + outcome + " (0<theta<10 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_0_10", title.c_str(), fNPiplusBinsP2019[0], fPiplusBinsP2019[0] ) );
  title = htitle + outcome + " (10<theta<20 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_10_20", title.c_str(), fNPiplusBinsP2019[1], fPiplusBinsP2019[1] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_20_40", title.c_str(), fNPiplusBinsP2019[2], fPiplusBinsP2019[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_40_60", title.c_str(), fNPiplusBinsP2019[3], fPiplusBinsP2019[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_60_100", title.c_str(), fNPiplusBinsP2019[4], fPiplusBinsP2019[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_100_140", title.c_str(), fNPiplusBinsP2019[5], fPiplusBinsP2019[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_140_180", title.c_str(), fNPiplusBinsP2019[6], fPiplusBinsP2019[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_180_240", title.c_str(), fNPiplusBinsP2019[7], fPiplusBinsP2019[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_240_300", title.c_str(), fNPiplusBinsP2019[8], fPiplusBinsP2019[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_300_360", title.c_str(), fNPiplusBinsP2019[9], fPiplusBinsP2019[9] ) );
  title = htitle + outcome + " (360<theta<420 [mrad]) ";
  fHistoPiPlusP2019.push_back( new TH1F( "piplus_360_420", title.c_str(), fNPiplusBinsP2019[10], fPiplusBinsP2019[10] ) );
    
  // secondary pi-
  //
  
  // 0-10
  fNPiminusBinsP2019[0] = 23; 
  fPiminusBinsP2019[0] = new double[fNPiminusBinsP2019[0]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[0]; ++i )
  {
     fPiminusBinsP2019[0][i] = tmp_pi_0[i]; // tmp_pim_0[i];
  }
  // 10-20
  fNPiminusBinsP2019[1] = 38; 
  fPiminusBinsP2019[1] = new double[fNPiminusBinsP2019[1]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[1]; ++i )
  {
     fPiminusBinsP2019[1][i] = tmp_pi_1[i]; 
  }
  // 20-40
  fNPiminusBinsP2019[2] = 32;  
  fPiminusBinsP2019[2] = new double[fNPiminusBinsP2019[2]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[2]; ++i )
  {
     fPiminusBinsP2019[2][i] = tmp_pi_2[i]; 
  }
  // 40-60
  fNPiminusBinsP2019[3] = 28; 
  fPiminusBinsP2019[3] = new double[fNPiminusBinsP2019[3]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[3]-1; ++i ) 
  {
     fPiminusBinsP2019[3][i] = tmp_pi_2[i]; 
  }
  fPiminusBinsP2019[3][fNPiminusBinsP2019[3]] = 25.0; 
  // 60-100
  fNPiminusBinsP2019[4] = 22;  
  fPiminusBinsP2019[4] = new double[fNPiminusBinsP2019[4]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[4]-2; ++i ) // -2 for pi+ beam only
  {
     fPiminusBinsP2019[4][i] = tmp_pi_2[i]; // tmp_pim_1[i];
  }
  fPiminusBinsP2019[4][fNPiminusBinsP2019[4]-1] = 13.0; // pi+ beam only
  fPiminusBinsP2019[4][fNPiminusBinsP2019[4]] = 16.0;  // pi+ beam only
  
  // 100-140
  fNPiminusBinsP2019[5] = 15; 
  fPiminusBinsP2019[5] = new double[fNPiminusBinsP2019[5]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[5]; ++i )
  {
     fPiminusBinsP2019[5][i] = tmp_pi_3[i]; // tmp_pim_1[i];
  }
  // 140-180
  fNPiminusBinsP2019[6] = 10; 
  fPiminusBinsP2019[6] = new double[fNPiminusBinsP2019[6]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[6]; ++i )
  {
     fPiminusBinsP2019[6][i] = tmp_pi_4[i]; 
  }
  // 180-240
  fNPiminusBinsP2019[7] = 10;
  fPiminusBinsP2019[7] = new double[fNPiminusBinsP2019[7]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[7]; ++i )
  {
     fPiminusBinsP2019[7][i] = tmp_pi_5[i]; 
  }
  // 240-300
  fNPiminusBinsP2019[8] = 7; 
  fPiminusBinsP2019[8] = new double[fNPiminusBinsP2019[8]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[8]; ++i )
  {
     fPiminusBinsP2019[8][i] = tmp_pi_6[i];
  }
  // 300-360
  fNPiminusBinsP2019[9] = 5; 
  fPiminusBinsP2019[9] = new double[fNPiminusBinsP2019[9]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[9]; ++i )
  {
     fPiminusBinsP2019[9][i] = tmp_pi_7[i]; 
  }
  // 360-420
  fNPiminusBinsP2019[10] = 5; 
  fPiminusBinsP2019[10] = new double[fNPiminusBinsP2019[10]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2019[10]; ++i )
  {
     fPiminusBinsP2019[10][i] = tmp_pi_8[i]; 
  }
  
  outcome = " -> X + pi- ";

  title = htitle + outcome + " (0<theta<10 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_0_10", title.c_str(), fNPiminusBinsP2019[0], fPiminusBinsP2019[0] ) );
  title = htitle + outcome + " (10<theta<20 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_10_20", title.c_str(), fNPiminusBinsP2019[1], fPiminusBinsP2019[1] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_20_40", title.c_str(), fNPiminusBinsP2019[2], fPiminusBinsP2019[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_40_60", title.c_str(), fNPiminusBinsP2019[3], fPiminusBinsP2019[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_60_100", title.c_str(), fNPiminusBinsP2019[4], fPiminusBinsP2019[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_100_140", title.c_str(), fNPiminusBinsP2019[5], fPiminusBinsP2019[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_140_180", title.c_str(), fNPiminusBinsP2019[6], fPiminusBinsP2019[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_180_240", title.c_str(), fNPiminusBinsP2019[7], fPiminusBinsP2019[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_240_300", title.c_str(), fNPiminusBinsP2019[8], fPiminusBinsP2019[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_300_360", title.c_str(), fNPiminusBinsP2019[9], fPiminusBinsP2019[9] ) );
  title = htitle + outcome + " (360<theta<420 [mrad]) ";
  fHistoPiMinusP2019.push_back( new TH1F( "piminus_360_420", title.c_str(), fNPiminusBinsP2019[10], fPiminusBinsP2019[10] ) );

  // secondary K+/-
  //
  
//  double tmp_k[] = {0.8, 1.6, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.6, 8.8, 
//                     10, 11.6, 13.2, 14.8, 16.4, 18, 19.6};
 
  double tmp_k_1[] = { 4.0,  6.0,  9.0, 12.0, 15.0, 18.0, 
                      21.0, 25.0, 29.0, 33.0, 37.0, 41.0 };
 
  double tmp_k_2[] ={ 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 18.0, 23.0, 35.0 };
  
  double tmp_k_3[] = { 4.0, 6.0, 8.0, 11.0, 15.0, 25.0 };
  
  double tmp_k_4[] = { 4.0, 5.0, 6.5, 8.0, 11.0, 16.0 };
  
  double tmp_kplus_5[] = { 4.0, 6.5, 12.0 };

  // K+, 0-20
  fNKplusBinsP2019[0] = 9;
  fKplusBinsP2019[0] = new double[fNKplusBinsP2019[0]+1];
  for ( int i=0; i<=fNKplusBinsP2019[0]; ++i )
  {
     fKplusBinsP2019[0][i] = tmp_k_1[i];
  } 
  // K+, 20-40
  fNKplusBinsP2019[1] = 8;
  fKplusBinsP2019[1] = new double[fNKplusBinsP2019[1]+1];
  for ( int i=0; i<=fNKplusBinsP2019[1]; ++i )
  {
     fKplusBinsP2019[1][i] = tmp_k_2[i];
  } 
  // K+, 40-60
  fNKplusBinsP2019[2] = 5;
  fKplusBinsP2019[2] = new double[fNKplusBinsP2019[2]+1];
  for ( int i=0; i<=fNKplusBinsP2019[2]; ++i )
  {
     fKplusBinsP2019[2][i] = tmp_k_3[i];
  } 
  // K+, 60-100
  fNKplusBinsP2019[3] = 5;
  fKplusBinsP2019[3] = new double[fNKplusBinsP2019[3]+1];
  for ( int i=0; i<=fNKplusBinsP2019[3]; ++i )
  {
     fKplusBinsP2019[3][i] = tmp_k_4[i];
  } 
  // K+, 100-140
  fNKplusBinsP2019[4] = 2;
  fKplusBinsP2019[4] = new double[fNKplusBinsP2019[4]+1];
  for ( int i=0; i<=fNKplusBinsP2019[4]; ++i )
  {
     fKplusBinsP2019[4][i] = tmp_kplus_5[i];
  } 
  
  outcome = " -> X + K+ ";

  title = htitle + outcome + " (0<theta<20 [mrad]) ";
  fHistoKPlusP2019.push_back( new TH1F( "kplus_0_20", title.c_str(), fNKplusBinsP2019[0], fKplusBinsP2019[0] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoKPlusP2019.push_back( new TH1F( "kplus_20_40", title.c_str(), fNKplusBinsP2019[0], fKplusBinsP2019[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoKPlusP2019.push_back( new TH1F( "kplus_40_60", title.c_str(), fNKplusBinsP2019[2], fKplusBinsP2019[2] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoKPlusP2019.push_back( new TH1F( "kplus_60_100", title.c_str(), fNKplusBinsP2019[3], fKplusBinsP2019[3] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoKPlusP2019.push_back( new TH1F( "kplus_100_140", title.c_str(), fNKplusBinsP2019[4], fKplusBinsP2019[4] ) );


  double tmp_kminus_5[] = { 4.0, 6.0, 12.0 };

  // K-, 0-20 
  fNKminusBinsP2019[0] = 11;
  fKminusBinsP2019[0] = new double[fNKminusBinsP2019[0]+1];
  for ( int i=0; i<=fNKminusBinsP2019[0]; ++i )
  {
     fKminusBinsP2019[0][i] = tmp_k_1[i];
  }
  // K-, 20-40
  fNKminusBinsP2019[1] = 8;
  fKminusBinsP2019[1] = new double[fNKminusBinsP2019[1]+1];
  for ( int i=0; i<=fNKminusBinsP2019[1]; ++i )
  {
     fKminusBinsP2019[1][i] = tmp_k_2[i];
  }
  // K-, 40-60
  fNKminusBinsP2019[2] = 5;
  fKminusBinsP2019[2] = new double[fNKminusBinsP2019[2]+1];
  for ( int i=0; i<=fNKminusBinsP2019[2]; ++i )
  {
     fKminusBinsP2019[2][i] = tmp_k_3[i];
  }
  // K-, 60-100
  fNKminusBinsP2019[3] = 5;
  fKminusBinsP2019[3] = new double[fNKminusBinsP2019[3]+1];
  for ( int i=0; i<=fNKminusBinsP2019[3]; ++i )
  {
     fKminusBinsP2019[3][i] = tmp_k_4[i];
  }
  // K-, 100-140
  fNKminusBinsP2019[4] = 2;
  fKminusBinsP2019[4] = new double[fNKminusBinsP2019[4]+1];
  for ( int i=0; i<=fNKminusBinsP2019[4]; ++i )
  {
     fKminusBinsP2019[4][i] = tmp_kminus_5[i];
  }
  
  outcome = " -> X + K- ";

  title = htitle + outcome + " (0<theta<20 [mrad]) ";
  fHistoKMinusP2019.push_back( new TH1F( "kminus_0_20", title.c_str(), fNKminusBinsP2019[0], fKminusBinsP2019[0] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoKMinusP2019.push_back( new TH1F( "kminus_20_40", title.c_str(), fNKminusBinsP2019[1], fKminusBinsP2019[1] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoKMinusP2019.push_back( new TH1F( "kminus_40_60", title.c_str(), fNKminusBinsP2019[2], fKminusBinsP2019[2] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoKMinusP2019.push_back( new TH1F( "kminus_60_100", title.c_str(), fNKminusBinsP2019[3], fKminusBinsP2019[3] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoKMinusP2019.push_back( new TH1F( "kminus_100_140", title.c_str(), fNKminusBinsP2019[4], fKminusBinsP2019[4] ) );

/*
  // secondary K0s
  //
  
  double tmp_k0s_0[] = {0.4, 3.4, 6.4, 9.4, 12.4};
  double tmp_k0s_1[] = {0.4, 2.4, 4.4, 6.4, 8.4, 10.4, 12.4};
  double tmp_k0s_2[] = {0.4, 2.4, 4.4, 6.4, 9.4};
  
  // 0-40
  fNK0sBinsP2019[0] = 4;
  fK0sBinsP2019[0] = new double[fNK0sBinsP2019[0]+1];
  for ( int i=0; i<=fNK0sBinsP2019[0]; ++i )
  {
     fK0sBinsP2019[0][i] = tmp_k0s_0[i];
  }
  // 40-60
  fNK0sBinsP2019[1] = 4;
  fK0sBinsP2019[1] = new double[fNK0sBinsP2019[1]+1];
  for ( int i=0; i<=fNK0sBinsP2019[1]; ++i )
  {
     fK0sBinsP2019[1][i] = tmp_k0s_0[i];
  }
  // 60-100
  fNK0sBinsP2019[2] = 6;
  fK0sBinsP2019[2] = new double[fNK0sBinsP2019[2]+1];
  for ( int i=0; i<=fNK0sBinsP2019[2]; ++i )
  {
     fK0sBinsP2019[2][i] = tmp_k0s_1[i];
  }
  // 100-140
  fNK0sBinsP2019[3] = 5;
  fK0sBinsP2019[3] = new double[fNK0sBinsP2019[3]+1];
  for ( int i=0; i<=fNK0sBinsP2019[3]; ++i )
  {
     fK0sBinsP2019[3][i] = tmp_k0s_1[i];
  }
  // 140-180
  fNK0sBinsP2019[4] = 4;
  fK0sBinsP2019[4] = new double[fNK0sBinsP2019[4]+1];
  for ( int i=0; i<=fNK0sBinsP2019[4]; ++i )
  {
     fK0sBinsP2019[4][i] = tmp_k0s_2[i];
  }
  // 180-240
  fNK0sBinsP2019[5] = 3;
  fK0sBinsP2019[5] = new double[fNK0sBinsP2019[5]+1];
  for ( int i=0; i<=fNK0sBinsP2019[5]; ++i )
  {
     fK0sBinsP2019[5][i] = tmp_k0s_2[i];
  }
  // 240-300
  fNK0sBinsP2019[6] = 2;
  fK0sBinsP2019[6] = new double[fNK0sBinsP2019[6]+1];
  for ( int i=0; i<=fNK0sBinsP2019[6]; ++i )
  {
     fK0sBinsP2019[6][i] = tmp_k0s_2[i];
  }
  
  outcome = " -> X + K0s ";

  title = htitle + outcome + " (0<theta<40 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_0_40", title.c_str(), fNK0sBinsP2019[0], fK0sBinsP2019[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_40_60", title.c_str(), fNK0sBinsP2019[1], fK0sBinsP2019[1] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_60_100", title.c_str(), fNK0sBinsP2019[2], fK0sBinsP2019[2] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_100_140", title.c_str(), fNK0sBinsP2019[3], fK0sBinsP2019[3] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_140_180", title.c_str(), fNK0sBinsP2019[4], fK0sBinsP2019[4] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_180_240", title.c_str(), fNK0sBinsP2019[5], fK0sBinsP2019[5] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoK0sP2019.push_back( new TH1F( "k0s_240_300", title.c_str(), fNK0sBinsP2019[6], fK0sBinsP2019[6] ) );

*/

/*
  // secondary Lambda
  //
  
  double tmp_l_0[] = {0.4, 4.4, 7.4, 10.4, 13.4, 17.4, 21.4};
  double tmp_l_1[] = {0.4, 2.9, 4.4, 6.4, 8.4, 10.4, 12.4, 18};
  double tmp_l_2[] = {0.4, 2.9, 4.4, 6.4, 8, 11.4};
  double tmp_l_3[] = {0.4, 2.4, 4.4, 6.4, 11};
  double tmp_l_4[] = {0.4, 2.4, 4.4, 10.4};
  
  // 0-40
  fNLambdaBinsP2019[0] = 6;
  fLambdaBinsP2019[0] = new double[fNLambdaBinsP2019[0]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[0]; ++i )
  {
     fLambdaBinsP2019[0][i] = tmp_l_0[i];
  }
  // 40-60
  fNLambdaBinsP2019[1] = 6;
  fLambdaBinsP2019[1] = new double[fNLambdaBinsP2019[1]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[1]; ++i )
  {
     fLambdaBinsP2019[1][i] = tmp_l_0[i];
  }
  // 60-100
  fNLambdaBinsP2019[2] = 7;
  fLambdaBinsP2019[2] = new double[fNLambdaBinsP2019[2]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[2]; ++i )
  {
     fLambdaBinsP2019[2][i] = tmp_l_1[i];
  }
  // 100-140
  fNLambdaBinsP2019[3] = 6;
  fLambdaBinsP2019[3] = new double[fNLambdaBinsP2019[3]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[3]; ++i )
  {
     fLambdaBinsP2019[3][i] = tmp_l_1[i];
  }
  // 140-180
  fNLambdaBinsP2019[4] = 5;
  fLambdaBinsP2019[4] = new double[fNLambdaBinsP2019[4]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[4]; ++i )
  {
     fLambdaBinsP2019[4][i] = tmp_l_2[i];
  }
  // 180-240
  fNLambdaBinsP2019[5] = 4;
  fLambdaBinsP2019[5] = new double[fNLambdaBinsP2019[5]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[5]; ++i )
  {
     fLambdaBinsP2019[5][i] = tmp_l_3[i];
  }
  // 240-300
  fNLambdaBinsP2019[6] = 3;
  fLambdaBinsP2019[6] = new double[fNLambdaBinsP2019[6]+1];
  for ( int i=0; i<=fNLambdaBinsP2019[6]; ++i )
  {
     fLambdaBinsP2019[6][i] = tmp_l_4[i];
  }

  outcome = " -> X + Lambda ";

  title = htitle + outcome + " (0<theta<40 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_0_40", title.c_str(), fNLambdaBinsP2019[0], fLambdaBinsP2019[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_40_60", title.c_str(), fNLambdaBinsP2019[1], fLambdaBinsP2019[1] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_60_100", title.c_str(), fNLambdaBinsP2019[2], fLambdaBinsP2019[2] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_100_140", title.c_str(), fNLambdaBinsP2019[3], fLambdaBinsP2019[3] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_140_180", title.c_str(), fNLambdaBinsP2019[4], fLambdaBinsP2019[4] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_180_240", title.c_str(), fNLambdaBinsP2019[5], fLambdaBinsP2019[5] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoLambdaP2019.push_back( new TH1F( "lambda_240_300", title.c_str(), fNLambdaBinsP2019[6], fLambdaBinsP2019[6] ) );
*/

}

TestNA61PiplusBeamHisto::~TestNA61PiplusBeamHisto()
{

   for (size_t i=0; i<fHistoProtonP2019.size(); ++i )
   {
      delete fHistoProtonP2019[i];
   }

   for (size_t i=0; i<fHistoPiPlusP2019.size(); ++i )
   {
      delete fHistoPiPlusP2019[i];
   }
   for (size_t i=0; i<fHistoPiMinusP2019.size(); ++i )
   {
      delete fHistoPiMinusP2019[i];
   }
   for (size_t i=0; i<fHistoKPlusP2019.size(); ++i )
   {
      delete fHistoKPlusP2019[i];
   }
   for (size_t i=0; i<fHistoKMinusP2019.size(); ++i )
   {
      delete fHistoKMinusP2019[i];
   }
   for (size_t i=0; i<fHistoK0sP2019.size(); ++i )
   {
      delete fHistoK0sP2019[i];
   }
   for (size_t i=0; i<fHistoLambdaP2019.size(); ++i )
   {
      delete fHistoLambdaP2019[i];
   }

}

void TestNA61PiplusBeamHisto::FillEvt( G4VParticleChange* aChange, 
                                       const G4LorentzVector&, 
				       const G4LorentzVector& ) 
{

   G4int NSec = aChange->GetNumberOfSecondaries();
     
   const G4DynamicParticle* sec = 0;
      
   int NHProP2019Size = fHistoProtonP2019.size();
   int NHPipP2019Size = fHistoPiPlusP2019.size();
   int NHPimP2019Size = fHistoPiMinusP2019.size();
   int NHKpP2019Size = fHistoKPlusP2019.size();
   int NHKmP2019Size = fHistoKMinusP2019.size();
   int NHK0sP2019Size = fHistoK0sP2019.size();
   int NHLamP2019Size = fHistoLambdaP2019.size();
     
   for (G4int i=0; i<NSec; i++) 
   {

      sec = aChange->GetSecondary(i)->GetDynamicParticle();
			
      const G4String& pname = sec->GetDefinition()->GetParticleName();
	
      // G4ThreeVector mom = sec->GetMomentum();
	
      double theta = sec->GetMomentum().theta();
      theta /= mrad;
	
      double pmom = sec->GetTotalMomentum() / GeV ;
	
      int id = -1;
      double dth1 = 1.;

      if ( pname == "pi+" || pname == "pi-" ) 
      {
	if ( theta >= 0. && theta < 10. )
	{
	   id = 0;
	   dth1 = 10. * mrad;
	}
	else if ( theta >=10. && theta < 20. )
	{
	   id = 1;
//	   id1 = 1;
	   dth1 = 10. * mrad;
	}
	else if ( theta >= 20. && theta < 40. )
	{
	   id = 2; // 1;
//	   id1 = 2;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 40. && theta < 60. )
	{
	   id = 3; // 2;
//	   id1 = 3;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 60. && theta < 100. )
	{
	   id = 4; // 3;
//	   id1 = 4;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 100. && theta < 140. )
	{
	   id = 5; // 4;
//	   id1 = 5;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 140. && theta < 180. )
	{
	   id = 6; // 5;
//	   id1 = 6;
	   dth1 = 40. * mrad;
	}
	else if ( theta > 180. && theta < 240. )
	{
	   id = 7; // 6;
//	   id1 = 7;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 240. && theta < 300. )
	{
	   id = 8; // 7;
//	   id1 = 8;
	   dth1 = 60. * mrad;
	} 
	else if ( theta >= 300. && theta < 360. )
	{
	   id = 9; // 8;
//	   id1 = 9;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 360. && theta < 420. )
	{
	   id = 10; // 9;
//	   id1 = 10;
	   dth1 = 60. * mrad;
	}
      }
	
      if ( pname == "proton" || pname == "kaon+" || pname == "kaon-")
      {
	if ( theta >= 0. && theta < 20. )
	{
	   id = 0;
//	   id1 = 0;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 20. && theta < 40. )
	{
	   id = 1;
//	   id1 = 2;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 40. && theta < 60. )
	{
	   id = 2;
//	   id1 = 3;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 60. && theta < 100. )
	{
	   id = 3;
//	   id1 = 4;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 100. && theta < 140. )
	{
	   id = 4;
//	   id1 = 5;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 140. && theta < 180. )
	{
	   id = 5;
//	   id1 = 6;
	   dth1 = 40. * mrad;
	}
	else if ( theta > 180. && theta < 240. )
	{
	   id = 6;
//	   id1 = 7;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 240. && theta < 300. )
	{
	   id = 7;
//	   id1 = 8;
	   dth1 = 60. * mrad;
	} 
	else if ( theta >= 300. && theta < 360. )
	{
	   id = 8;
//	   id1 = 9;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 360. && theta < 420. )
	{
	   id = 9;
//	   id1 = 10;
	   dth1 = 60. * mrad;
	}
      }
      
      if ( pname == "proton" )
      {	   
	   if ( id >= 0 && id < NHProP2019Size ) fHistoProtonP2019[id]->Fill( pmom, 1./dth1 );
      }
      if ( pname == "pi-" )
      {
	   if ( id >= 0 && id < NHPimP2019Size ) fHistoPiMinusP2019[id]->Fill( pmom, 1./dth1 );
      }
      if ( pname == "pi+" )
      {
	   if ( id >= 0 && id < NHPipP2019Size ) fHistoPiPlusP2019[id]->Fill( pmom, 1./dth1 );
      }
      if ( pname == "kaon+" )
      {
	   if ( id >=0 && id < NHKpP2019Size ) fHistoKPlusP2019[id]->Fill( pmom, 1./dth1 );
      }
      if ( pname == "kaon-" )
      {
	   if ( id >=0 && id < NHKmP2019Size ) fHistoKMinusP2019[id]->Fill( pmom, 1./dth1 );
      }
      
      if ( pname == "proton" || pname == "pi-" || pname == "kaon-" ) continue;
      
      if ( pname == "pi+" || pname == "kaon+" ) continue;

      id = -1;
      dth1 = 1.;
      
      if ( pname == "kaon0S" || pname == "lambda" )
      {
         if ( theta >= 0. && theta < 40. )
	 {
	    id = 0;
	    dth1 = 40. * mrad;
	 }
	 else if ( theta >= 40. && theta < 60. )
	 {
	    id = 1;
	    dth1 = 20. * mrad;
	 }
	 else if ( theta >= 60. && theta < 100. )
	 {
	    id = 2;
	    dth1 = 40. * mrad;
	 }
	 else if ( theta >= 100. && theta < 140. )
	 {
	    id = 3;
	    dth1 = 40. * mrad;
	 }
	 else if ( theta >= 140. && theta < 180. )
	 {
	    id = 4;
	    dth1 = 40. * mrad;
	 }
	 else if ( theta >= 180. && theta < 240. )
	 {
	    id = 5;
	    dth1 = 60. * mrad;
	 }
	 else if ( theta >= 240. && theta < 300. )
	 {
	    id = 7;
	    dth1 = 60. * mrad;
	 }
	 
	 if ( pname == "kaon0S" )
	 {
	    if ( id >= 0 && id < NHK0sP2019Size ) fHistoK0sP2019[id]->Fill(pmom,1./dth1);
	 }
	 if ( pname == "lambda" )
	 {
	    if ( id >= 0 && id < NHLamP2019Size ) fHistoLambdaP2019[id]->Fill(pmom, 1./dth1);
	 }
	 
      }

   }
      
   return;
   
}

void TestNA61PiplusBeamHisto::Write( G4int stat, G4double ) 
{

   double scale = 1. / ((double)stat) ;
   
   for ( size_t i=0; i<fHistoProtonP2019.size(); ++i )
   {
      fHistoProtonP2019[i]->Scale(scale,"width");
      fHistoProtonP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoPiPlusP2019.size(); ++i )
   {
      fHistoPiPlusP2019[i]->Scale(scale,"width");
      fHistoPiPlusP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoPiMinusP2019.size(); ++i )
   {
      fHistoPiMinusP2019[i]->Scale(scale,"width");
      fHistoPiMinusP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoKPlusP2019.size(); ++i )
   {
      fHistoKPlusP2019[i]->Scale(scale,"width");
      fHistoKPlusP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoKMinusP2019.size(); ++i )
   {
      fHistoKMinusP2019[i]->Scale(scale,"width");
      fHistoKMinusP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoK0sP2019.size(); ++i )
   {
      fHistoK0sP2019[i]->Scale(scale,"width");
      fHistoK0sP2019[i]->Write();
   }
   for ( size_t i=0; i<fHistoLambdaP2019.size(); ++i )
   {
      fHistoLambdaP2019[i]->Scale(scale,"width");
      fHistoLambdaP2019[i]->Write();
   }
  
   return;

}
