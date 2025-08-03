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



#include "Tst23NA61Histo.hh"
#include "Tst23ParticleChange.hh"

#include "G4TrackVector.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"

#include "G4SystemOfUnits.hh"

// #include "G4Proton.hh" // tmp needed for Ecm/beta/gamma cross-checks


Tst23NA61Histo::Tst23NA61Histo( G4String htitle )
   : TstHistoSet( htitle )
{

  G4String title;
  G4String outcome;

  fHistoNSec = new TH1D( "NSec", htitle.c_str(), 100, 0., 100. );
  
  // secondary proton case
  //
    
  outcome = " -> X + proton ";
       
  title = htitle + outcome + " (0<theta<20 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_0_20",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_20_40",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_40_60",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_60_100",title.c_str(),250,0.,25.) ); 
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_100_140",title.c_str(),150,0.,15.) ); 
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_140_180",title.c_str(),150,0.,15.) ); 
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
  fHistoSecProton.push_back( new TH1F("protonMult_180_240",title.c_str(),150,0.,15.) ); 
//
// no exp.data for this bin
//
//  title = htitle + outcome + " (240<theta<300 (mrad)) ";
//  fHistoSecProton.push_back( new TH1F("protonMult_240_300",title.c_str(),150,0.,15.) ); 
      
  // secondary pi+ case
  //
  fNSecPiplusThetaBins = new int[10];
  fNSecPiplusThetaBins[0] = 25; // 24+1 - due to the "hole" at p=1.0-1.2
  fNSecPiplusThetaBins[1] = 43; // 42+1 - due to the "hole" at p=1.0-1.2
  fNSecPiplusThetaBins[2] = 50; // 49+1 - due to the "hole" at p=1.0-1.2
  fNSecPiplusThetaBins[3] = 44;
  fNSecPiplusThetaBins[4] = 38;
  fNSecPiplusThetaBins[5] = 31;
  fNSecPiplusThetaBins[6] = 27;
  fNSecPiplusThetaBins[7] = 13;
  fNSecPiplusThetaBins[8] =  8;
  fNSecPiplusThetaBins[9] =  8;
  
  fNSecPiminusThetaBins = new int[10];
  fNSecPiminusThetaBins[0] = 29;
  fNSecPiminusThetaBins[1] = 46;
  fNSecPiminusThetaBins[2] = 49;
  fNSecPiminusThetaBins[3] = 46;
  fNSecPiminusThetaBins[4] = 40;
  fNSecPiminusThetaBins[5] = 35;
  fNSecPiminusThetaBins[6] = 32;
  fNSecPiminusThetaBins[7] = 25;
  fNSecPiminusThetaBins[8] = 20;
  fNSecPiminusThetaBins[9] = 18;
      
  fSecPiplusBins = new double*[10];
  fSecPiminusBins = new double*[10];
  
  for ( int ii=0; ii<10; ++ii )
  {
     fSecPiplusBins[ii] = new double[fNSecPiplusThetaBins[ii]+1]; // N+1 to handle bin low/high edges
     fSecPiminusBins[ii] = new double[fNSecPiminusThetaBins[ii]+1];
     
     // I know that the first 9 points it's all the same (from 0.2 to 1.0)
     // it's common for pi+ & pi-
     //
     for ( int jj=0; jj<9; ++jj )
     {
        fSecPiplusBins[ii][jj]  = 0.1*jj + 0.2;
	fSecPiminusBins[ii][jj] = 0.1*jj + 0.2;
     }
  }
  
   //
   // NOTE: there's a "hole" in the p-spectra at 1.0-1.2 
   // for theta 0-20, 20-40, and 40-60
   //
   for ( int ii=0; ii<10; ++ii )
   {
      fSecPiminusBins[ii][9]  = 1.2;
      if ( ii > 7 ) continue;
      fSecPiplusBins[ii][9] = 1.2;
   }
   
   
   // pi+ specific
   //
   // binning is 0.4 from 1.2 and up to 10.0 for 0-20 & 20-40
   //
   for ( int ii=0; ii<2; ++ii )
   {
      int limit = std::min( 22, fNSecPiplusThetaBins[ii]-10+1 ); // 22 = (10.-1.2)/0.4
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiplusBins[ii][jj+10] = 0.4*jj + 1.6; // 1.6 = 1.2 + 0.4
      }
   }
   
   // binning is 0.8 from 10 and up to 19.6 (where applicable)
   //
   for ( int ii=0; ii<2; ++ii )
   {
      int limit = std::min( 12, fNSecPiplusThetaBins[ii]-32+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiplusBins[ii][jj+32] = 0.8*jj + 10.8;
      }
   }
   
   //
   // binning is 0.2 from 1.2 and up to 5.2 for all theta bins from 40 up
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 20, fNSecPiplusThetaBins[ii]-10+1 ); // 20 = (5.2 - 1.2)/0.2
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiplusBins[ii][jj+10] = 0.2*jj + 1.4; // 1.4 = 1.2 + 0.2
      }
   }
   //
   // binning is 0.4 from 5.2 and up to 10. for all theta bins from 40 up
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 12, fNSecPiplusThetaBins[ii]-30+1 ); // 12 = (10.-5.2)/0.4; 30=20+10, see above
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiplusBins[ii][jj+30] = 0.4*jj + 5.6;
      }
   }
   
   // from p=10. and up to 19.6 the binning is 0.8 (where applicable), for both pi+ and pi-,
   // but the p-bin numbering is shifted, so it needs to be done separately
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 12, fNSecPiplusThetaBins[ii]-42+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiplusBins[ii][jj+42] = 0.8*jj + 10.8;
      }
   }


   // now pi- specific
   //
   // no "holes" in the spectra
   //
   // binning is 0.2 from 1.2 and up to at least 2.8 for all theta bins
   //
   for ( int ii=0; ii<10; ++ii )
   {
      int limit = std::min( 8, fNSecPiminusThetaBins[ii]-10+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+10] = 0.2*jj + 1.4;
      }
   }
   
   // binning is 0.4 from 2.8 and up to 10., where applicable, for theta = 0-20, 20-40
   //
   for ( int ii=0; ii<2; ++ii )
   {
      int limit = std::min( 18, fNSecPiminusThetaBins[ii]-18+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+18] = 0.4*jj + 3.2; // 3.2=2.8+0.4
      }
   }
   
   // binning is 0.8 from 10 and up to 18.8 
   //
   for ( int ii=0; ii<2; ++ii )
   {
      int limit = std::min( 11, fNSecPiminusThetaBins[ii]-36+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+36] = 0.8*jj + 10.8;
      }
   }
      
   // binning is also 0.2 from 2.8 up to 5.2,  from theta=40 up (where applicable)
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 12, fNSecPiminusThetaBins[ii]-18+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+18] = 0.2*jj + 3.0; // 3.0=2.8+0.2
      }
   }
   
   // binning is 0.4 from p=5.2 and up to 10.0, for theta=40 and up
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 12, fNSecPiminusThetaBins[ii]-30+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+30] = 0.4*jj + 5.6;
      }
   }
         
   // from p=10. onwards the binning is 0.8 (where applicable), for both pi+ and pi-,
   // but the p-bin numbering is different, that's why it nedds to be done separately
   //
   for ( int ii=2; ii<10; ++ii )
   {
      int limit = std::min( 11, fNSecPiminusThetaBins[ii]-42+1 );
      for ( int jj=0; jj<limit; ++jj )
      {
         fSecPiminusBins[ii][jj+42] = 0.8*jj + 10.8;
      }
   }
   
  
/*
// cross-check printouts   
  
  std::cout << " PI+ BUSINESS " << std::endl;
  
  for ( int ii=0; ii<10; ++ii )
  {
     
     std::cout << " ii = " << ii << std::endl; 
     
     for ( int jj=0; jj<fNSecPiplusThetaBins[ii]+1; ++jj )
     {
        std::cout << fSecPiplusBins[ii][jj] << " " ;
     }
     std::cout << std::endl;
  }
  
  std::cout << " PI- NUSINESS " << std::endl;
  
  for ( int ii=0; ii<10; ++ii )
  {
     
     std::cout << " ii = " << ii << std::endl; 
     
     for ( int jj=0; jj<fNSecPiminusThetaBins[ii]+1; ++jj )
     {
        std::cout << fSecPiminusBins[ii][jj] << " " ;
     }
     std::cout << std::endl;
  }
*/  
  
    
  outcome = " -> X + pi- ";

  title = htitle + outcome + " (0<theta<20 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_0_20", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_0_20", title.c_str(), fNSecPiminusThetaBins[0], fSecPiminusBins[0] ) );
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_20_40", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_20_40", title.c_str(), fNSecPiminusThetaBins[1], fSecPiminusBins[1] ) );
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_40_60", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_40_60", title.c_str(), fNSecPiminusThetaBins[2], fSecPiminusBins[2] ) );
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_60_100", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_60_100", title.c_str(), fNSecPiminusThetaBins[3], fSecPiminusBins[3] ) );
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_100_140", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_100_140", title.c_str(), fNSecPiminusThetaBins[4], fSecPiminusBins[4] ) );
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_140_180", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_140_180", title.c_str(), fNSecPiminusThetaBins[5], fSecPiminusBins[5] ) );
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_180_240", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_180_240", title.c_str(), fNSecPiminusThetaBins[6], fSecPiminusBins[6] ) );
  title = htitle + outcome + " (240<theta<300 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_240_300", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_240_300", title.c_str(), fNSecPiminusThetaBins[7], fSecPiminusBins[7] ) );
  title = htitle + outcome + " (300<theta<360 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_300_360", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_300_360", title.c_str(), fNSecPiminusThetaBins[8], fSecPiminusBins[8] ) );
  title = htitle + outcome + " (360<theta<420 (mrad)) ";
//  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_360_420", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiMinus.push_back( new TH1F( "piminusMult_360_420", title.c_str(), fNSecPiminusThetaBins[9], fSecPiminusBins[9] ) );

  outcome = " -> X + pi+ ";

  title = htitle + outcome + " (0<theta<20 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_0_20", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_0_20", title.c_str(), fNSecPiplusThetaBins[0], fSecPiplusBins[0] ) );
  title = htitle + outcome + " (20<theta<40 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_20_40", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_20_40", title.c_str(), fNSecPiplusThetaBins[1], fSecPiplusBins[1] ) );
  title = htitle + outcome + " (40<theta<60 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_40_60", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_40_60", title.c_str(), fNSecPiplusThetaBins[2], fSecPiplusBins[2] ) );
  title = htitle + outcome + " (60<theta<100 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_60_100", title.c_str(), 250, 0., 25. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_60_100", title.c_str(), fNSecPiplusThetaBins[3], fSecPiplusBins[3] ) );
  title = htitle + outcome + " (100<theta<140 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_100_140", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_100_140", title.c_str(), fNSecPiplusThetaBins[4], fSecPiplusBins[4] ) );
  title = htitle + outcome + " (140<theta<180 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_140_180", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_140_180", title.c_str(), fNSecPiplusThetaBins[5], fSecPiplusBins[5] ) );
  title = htitle + outcome + " (180<theta<240 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_180_240", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_180_240", title.c_str(), fNSecPiplusThetaBins[6], fSecPiplusBins[6] ) );
  title = htitle + outcome + " (240<theta<300 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_240_300", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_240_300", title.c_str(), fNSecPiplusThetaBins[7], fSecPiplusBins[7] ) );
  title = htitle + outcome + " (300<theta<360 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_300_360", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_300_360", title.c_str(), fNSecPiplusThetaBins[8], fSecPiplusBins[8] ) );
  title = htitle + outcome + " (360<theta<420 (mrad)) ";
//  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_360_420", title.c_str(), 150, 0., 15. ) );
  fHistoSecPiPlus.push_back( new TH1F( "piplusMult_360_420", title.c_str(), fNSecPiplusThetaBins[9], fSecPiplusBins[9] ) );
  
  outcome = " -> X + K+ ";
  
  title = htitle + outcome + " (20<theta<140 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplusMult_20_140",     title.c_str(), 100, 0., 10. ) );
  title = htitle + outcome + " (140<theta<240 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplusMult_140_240",    title.c_str(), 100, 0., 10. ) );
  title = htitle + ", K+/pi+ (20<theta<140 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplus2piplus_20_140",  title.c_str(), 100, 0., 10. ) );
  title = htitle + ", K+/pi+ (140<theta<240 (mrad)) ";
  fHistoSecKPlus.push_back( new TH1F( "kplus2piplus_140_240",  title.c_str(), 100, 0., 10. ) );
  
  fHistoSecPiPlus2.push_back( new TH1F( "piplus_20_140",  " ", 100, 0., 10. ) );
  fHistoSecPiPlus2.push_back( new TH1F( "piplus_140_240", " ", 100, 0., 10. ) );

  // now, for comparisons vs most recent data from Pub-2015
  //
  
  // secondary proton 
  //
  double tmp_p[] = {1.2, 1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 
                    5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 
		    9.2, 9.6, 10, 10.8, 11.6, 12.4, 13.2, 14, 14.8, 15.6, 
		    16.4, 17.2, 18, 18.8, 19.6, 20.4, 21.2, 22};

  // 0-10
  fNProtonBinsP2015[0] = 30;
  fProtonBinsP2015[0] = new double[fNProtonBinsP2015[0]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[0]; ++i )
  {
     fProtonBinsP2015[0][i] = tmp_p[i];
  }
  // 10-20
  fNProtonBinsP2015[1] = 37;
  fProtonBinsP2015[1] = new double[fNProtonBinsP2015[1]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[1]; ++i )
  {
     fProtonBinsP2015[1][i] = tmp_p[i];
  }
  // 20-40
  fNProtonBinsP2015[2] = 37;
  fProtonBinsP2015[2] = new double[fNProtonBinsP2015[2]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[2]; ++i )
  {
     fProtonBinsP2015[2][i] = tmp_p[i];
  }
  // 40-60
  fNProtonBinsP2015[3] = 35;
  fProtonBinsP2015[3] = new double[fNProtonBinsP2015[3]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[3]; ++i )
  {
     fProtonBinsP2015[3][i] = tmp_p[i];
  }
  // 60-100
  fNProtonBinsP2015[4] = 32;
  fProtonBinsP2015[4] = new double[fNProtonBinsP2015[4]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[4]; ++i )
  {
     fProtonBinsP2015[4][i] = tmp_p[i];
  }
  // 100-140
  fNProtonBinsP2015[5] = 26;
  fProtonBinsP2015[5] = new double[fNProtonBinsP2015[5]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[5]; ++i )
  {
     fProtonBinsP2015[5][i] = tmp_p[i];
  }
  // 140-180
  fNProtonBinsP2015[6] = 22;
  fProtonBinsP2015[6] = new double[fNProtonBinsP2015[6]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[6]; ++i )
  {
     fProtonBinsP2015[6][i] = tmp_p[i];
  }
  // 180-240
  fNProtonBinsP2015[7] = 16;
  fProtonBinsP2015[7] = new double[fNProtonBinsP2015[7]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[7]; ++i )
  {
     fProtonBinsP2015[7][i] = tmp_p[i];
  }
  // 240-300
  fNProtonBinsP2015[8] = 10;
  fProtonBinsP2015[8] = new double[fNProtonBinsP2015[8]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[8]; ++i )
  {
     fProtonBinsP2015[8][i] = tmp_p[i];
  }
  // 300-360
  fNProtonBinsP2015[9] = 6;
  fProtonBinsP2015[9] = new double[fNProtonBinsP2015[9]+1]; 
  for ( int i=0; i<=fNProtonBinsP2015[9]; ++i )
  {
     fProtonBinsP2015[9][i] = tmp_p[i];
  }

  outcome = " -> X + proton ";


  title = htitle + outcome + " (0<theta<10 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_0_10", title.c_str(), fNProtonBinsP2015[0], fProtonBinsP2015[0] ) );
  title = htitle + outcome + " (10<theta<20 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_10_20", title.c_str(), fNProtonBinsP2015[1], fProtonBinsP2015[1] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_20_40", title.c_str(), fNProtonBinsP2015[2], fProtonBinsP2015[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_40_60", title.c_str(), fNProtonBinsP2015[3], fProtonBinsP2015[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_60_100", title.c_str(), fNProtonBinsP2015[4], fProtonBinsP2015[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_100_140", title.c_str(), fNProtonBinsP2015[5], fProtonBinsP2015[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_140_180", title.c_str(), fNProtonBinsP2015[6], fProtonBinsP2015[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_180_240", title.c_str(), fNProtonBinsP2015[7], fProtonBinsP2015[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_240_300", title.c_str(), fNProtonBinsP2015[8], fProtonBinsP2015[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoProtonP2015.push_back( new TH1F( "proton_300_360", title.c_str(), fNProtonBinsP2015[9], fProtonBinsP2015[9] ) );

  // secondary pi+
  //
  double tmp_pip_0[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 
                       1.6, 2, 2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 
		       5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 
		       9.6, 10, 10.8, 11.6, 12.4, 13.2, 14, 14.8, 15.6, 16.4, 
		       17.2, 18, 18.8, 19.6, 20.4, 21.2, 22};  
  double tmp_pip_1[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 
                       1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 
		       3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.2, 
		       5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 
		       9.6, 10, 10.8, 11.6, 12.4, 13.2, 14, 14.8, 15.6, 16.4, 
		       17.2, 18, 18.8, 19.6, 20.4};

  // 0-10
  fNPiplusBinsP2015[0] = 39; 
  fPiplusBinsP2015[0] = new double[fNPiplusBinsP2015[0]+1]; 
  for ( int i=1; i<=fNPiplusBinsP2015[0]+1; ++i )
  {
     fPiplusBinsP2015[0][i] = tmp_pip_0[i];
  }
  // 10-20
  fNPiplusBinsP2015[1] = 46; 
  fPiplusBinsP2015[1] = new double[fNPiplusBinsP2015[1]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[1]; ++i )
  {
     fPiplusBinsP2015[1][i] = tmp_pip_0[i];
  }
  // 20-40
  fNPiplusBinsP2015[2] = 46; 
  fPiplusBinsP2015[2] = new double[fNPiplusBinsP2015[2]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[2]; ++i )
  {
     fPiplusBinsP2015[2][i] = tmp_pip_0[i];
  }
  // 40-60
  fNPiplusBinsP2015[3] = 54; 
  fPiplusBinsP2015[3] = new double[fNPiplusBinsP2015[3]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[3]; ++i )
  {
     fPiplusBinsP2015[3][i] = tmp_pip_1[i];
  }
  // 60-100
  fNPiplusBinsP2015[4] = 51; 
  fPiplusBinsP2015[4] = new double[fNPiplusBinsP2015[4]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[4]; ++i )
  {
     fPiplusBinsP2015[4][i] = tmp_pip_1[i];
  }
  // 100-140
  fNPiplusBinsP2015[5] = 46; 
  fPiplusBinsP2015[5] = new double[fNPiplusBinsP2015[5]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[5]; ++i )
  {
     fPiplusBinsP2015[5][i] = tmp_pip_1[i];
  }
  // 140-180
  fNPiplusBinsP2015[6] = 42; 
  fPiplusBinsP2015[6] = new double[fNPiplusBinsP2015[6]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[6]; ++i )
  {
     fPiplusBinsP2015[6][i] = tmp_pip_1[i];
  }
  // 180-240
  fNPiplusBinsP2015[7] = 36; 
  fPiplusBinsP2015[7] = new double[fNPiplusBinsP2015[7]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[7]; ++i )
  {
     fPiplusBinsP2015[7][i] = tmp_pip_1[i];
  }
  // 240-300
  fNPiplusBinsP2015[8] = 28; 
  fPiplusBinsP2015[8] = new double[fNPiplusBinsP2015[8]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[8]; ++i )
  {
     fPiplusBinsP2015[8][i] = tmp_pip_1[i];
  }
  // 300-360
  fNPiplusBinsP2015[9] = 18; 
  fPiplusBinsP2015[9] = new double[fNPiplusBinsP2015[9]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[9]; ++i )
  {
     fPiplusBinsP2015[9][i] = tmp_pip_1[i];
  }
  // 360-420
  fNPiplusBinsP2015[10] = 8; 
  fPiplusBinsP2015[10] = new double[fNPiplusBinsP2015[10]+1]; 
  for ( int i=0; i<=fNPiplusBinsP2015[10]; ++i )
  {
     fPiplusBinsP2015[10][i] = tmp_pip_1[i];
  }

  outcome = " -> X + pi+ ";

  title = htitle + outcome + " (0<theta<10 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_0_10", title.c_str(), fNPiplusBinsP2015[0], fPiplusBinsP2015[0] ) );
  title = htitle + outcome + " (10<theta<20 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_10_20", title.c_str(), fNPiplusBinsP2015[1], fPiplusBinsP2015[1] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_20_40", title.c_str(), fNPiplusBinsP2015[2], fPiplusBinsP2015[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_40_60", title.c_str(), fNPiplusBinsP2015[3], fPiplusBinsP2015[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_60_100", title.c_str(), fNPiplusBinsP2015[4], fPiplusBinsP2015[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_100_140", title.c_str(), fNPiplusBinsP2015[5], fPiplusBinsP2015[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_140_180", title.c_str(), fNPiplusBinsP2015[6], fPiplusBinsP2015[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_180_240", title.c_str(), fNPiplusBinsP2015[7], fPiplusBinsP2015[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_240_300", title.c_str(), fNPiplusBinsP2015[8], fPiplusBinsP2015[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_300_360", title.c_str(), fNPiplusBinsP2015[9], fPiplusBinsP2015[9] ) );
  title = htitle + outcome + " (360<theta<420 [mrad]) ";
  fHistoPiPlusP2015.push_back( new TH1F( "piplus_360_420", title.c_str(), fNPiplusBinsP2015[10], fPiplusBinsP2015[10] ) );

  // secondary pi-
  //
  double tmp_pim_0[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 
                        1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3.2, 3.6, 
			4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 
			8, 8.4, 8.8, 9.2, 9.6, 10, 10.8, 11.6, 12.4, 13.2, 
			14, 14.8, 15.6, 16.4, 17.2, 18, 18.8, 19.6, 20.4, 21.2};
  double tmp_pim_1[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 
                        1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 
			3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.2, 
			5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 
			9.6, 10, 10.8, 11.6, 12.4, 13.2, 14, 14.8, 15.6, 16.4, 
			17.2, 18, 18.8, 19.6, 20.4};
  
  // 0-10
  fNPiminusBinsP2015[0] = 42; 
  fPiminusBinsP2015[0] = new double[fNPiminusBinsP2015[0]+1]; 
  for ( int i=1; i<=fNPiminusBinsP2015[0]; ++i )
  {
     fPiminusBinsP2015[0][i] = tmp_pim_0[i];
  }
  // 10-20
  fNPiminusBinsP2015[1] = 49; 
  fPiminusBinsP2015[1] = new double[fNPiminusBinsP2015[1]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[1]; ++i )
  {
     fPiminusBinsP2015[1][i] = tmp_pim_0[i];
  }
  // 20-40
  fNPiminusBinsP2015[2] = 49; 
  fPiminusBinsP2015[2] = new double[fNPiminusBinsP2015[2]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[2]; ++i )
  {
     fPiminusBinsP2015[2][i] = tmp_pim_0[i];
  }
  // 40-60
  fNPiminusBinsP2015[3] = 54; 
  fPiminusBinsP2015[3] = new double[fNPiminusBinsP2015[3]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[3]; ++i )
  {
     fPiminusBinsP2015[3][i] = tmp_pim_1[i];
  }
  // 60-100
  fNPiminusBinsP2015[4] = 51; 
  fPiminusBinsP2015[4] = new double[fNPiminusBinsP2015[4]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[4]; ++i )
  {
     fPiminusBinsP2015[4][i] = tmp_pim_1[i];
  }
  // 100-140
  fNPiminusBinsP2015[5] = 46; 
  fPiminusBinsP2015[5] = new double[fNPiminusBinsP2015[5]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[5]; ++i )
  {
     fPiminusBinsP2015[5][i] = tmp_pim_1[i];
  }
  // 140-180
  fNPiminusBinsP2015[6] = 43; 
  fPiminusBinsP2015[6] = new double[fNPiminusBinsP2015[6]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[6]; ++i )
  {
     fPiminusBinsP2015[6][i] = tmp_pim_1[i];
  }
  // 180-240
  fNPiminusBinsP2015[7] = 38; 
  fPiminusBinsP2015[7] = new double[fNPiminusBinsP2015[7]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[7]; ++i )
  {
     fPiminusBinsP2015[7][i] = tmp_pim_1[i];
  }
  // 240-300
  fNPiminusBinsP2015[8] = 33; 
  fPiminusBinsP2015[8] = new double[fNPiminusBinsP2015[8]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[8]; ++i )
  {
     fPiminusBinsP2015[8][i] = tmp_pim_1[i];
  }
  // 300-360
  fNPiminusBinsP2015[9] = 23; 
  fPiminusBinsP2015[9] = new double[fNPiminusBinsP2015[9]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[9]; ++i )
  {
     fPiminusBinsP2015[9][i] = tmp_pim_1[i];
  }
  // 360-420
  fNPiminusBinsP2015[10] = 19; 
  fPiminusBinsP2015[10] = new double[fNPiminusBinsP2015[10]+1]; 
  for ( int i=0; i<=fNPiminusBinsP2015[10]; ++i )
  {
     fPiminusBinsP2015[10][i] = tmp_pim_1[i];
  }
  
  outcome = " -> X + pi- ";

  title = htitle + outcome + " (0<theta<10 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_0_10", title.c_str(), fNPiminusBinsP2015[0], fPiminusBinsP2015[0] ) );
  title = htitle + outcome + " (10<theta<20 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_10_20", title.c_str(), fNPiminusBinsP2015[1], fPiminusBinsP2015[1] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_20_40", title.c_str(), fNPiminusBinsP2015[2], fPiminusBinsP2015[2] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_40_60", title.c_str(), fNPiminusBinsP2015[3], fPiminusBinsP2015[3] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_60_100", title.c_str(), fNPiminusBinsP2015[4], fPiminusBinsP2015[4] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_100_140", title.c_str(), fNPiminusBinsP2015[5], fPiminusBinsP2015[5] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_140_180", title.c_str(), fNPiminusBinsP2015[6], fPiminusBinsP2015[6] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_180_240", title.c_str(), fNPiminusBinsP2015[7], fPiminusBinsP2015[7] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_240_300", title.c_str(), fNPiminusBinsP2015[8], fPiminusBinsP2015[8] ) );
  title = htitle + outcome + " (300<theta<360 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_300_360", title.c_str(), fNPiminusBinsP2015[9], fPiminusBinsP2015[9] ) );
  title = htitle + outcome + " (360<theta<420 [mrad]) ";
  fHistoPiMinusP2015.push_back( new TH1F( "piminus_360_420", title.c_str(), fNPiminusBinsP2015[10], fPiminusBinsP2015[10] ) );

  // secondary K+/-
  //
  
  double tmp_k[] = {0.8, 1.6, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.6, 8.8, 
                     10, 11.6, 13.2, 14.8, 16.4, 18, 19.6};
 
  // K+, 0-20
  fNKplusBinsP2015[0] = 14;
  fKplusBinsP2015[0] = new double[fNKplusBinsP2015[0]+1];
  for ( int i=0; i<=fNKplusBinsP2015[0]; ++i )
  {
     fKplusBinsP2015[0][i] = tmp_k[i];
  } 
  // K+, 20-40
  fNKplusBinsP2015[1] = 14;
  fKplusBinsP2015[1] = new double[fNKplusBinsP2015[1]+1];
  for ( int i=0; i<=fNKplusBinsP2015[1]; ++i )
  {
     fKplusBinsP2015[1][i] = tmp_k[i];
  } 
  // K+, 40-60
  fNKplusBinsP2015[2] = 16;
  fKplusBinsP2015[2] = new double[fNKplusBinsP2015[2]+1];
  for ( int i=0; i<=fNKplusBinsP2015[2]; ++i )
  {
     fKplusBinsP2015[2][i] = tmp_k[i];
  } 
  // K+, 60-100
  fNKplusBinsP2015[3] = 14;
  fKplusBinsP2015[3] = new double[fNKplusBinsP2015[3]+1];
  for ( int i=0; i<=fNKplusBinsP2015[3]; ++i )
  {
     fKplusBinsP2015[3][i] = tmp_k[i];
  } 
  // K+, 100-140
  fNKplusBinsP2015[4] = 12;
  fKplusBinsP2015[4] = new double[fNKplusBinsP2015[4]+1];
  for ( int i=0; i<=fNKplusBinsP2015[4]; ++i )
  {
     fKplusBinsP2015[4][i] = tmp_k[i];
  } 
  // K+, 140-180
  fNKplusBinsP2015[5] = 10;
  fKplusBinsP2015[5] = new double[fNKplusBinsP2015[5]+1];
  for ( int i=0; i<=fNKplusBinsP2015[5]; ++i )
  {
     fKplusBinsP2015[5][i] = tmp_k[i];
  } 
  // K+, 180-240
  fNKplusBinsP2015[6] = 8;
  fKplusBinsP2015[6] = new double[fNKplusBinsP2015[6]+1];
  for ( int i=0; i<=fNKplusBinsP2015[6]; ++i )
  {
     fKplusBinsP2015[6][i] = tmp_k[i];
  } 
  // K+, 240-300
  fNKplusBinsP2015[7] = 6;
  fKplusBinsP2015[7] = new double[fNKplusBinsP2015[7]+1];
  for ( int i=0; i<=fNKplusBinsP2015[7]; ++i )
  {
     fKplusBinsP2015[7][i] = tmp_k[i];
  }

  
  outcome = " -> X + K+ ";

  title = htitle + outcome + " (0<theta<20 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_0_20", title.c_str(), fNKplusBinsP2015[0], fKplusBinsP2015[0] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_20_40", title.c_str(), fNKplusBinsP2015[0], fKplusBinsP2015[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_40_60", title.c_str(), fNKplusBinsP2015[2], fKplusBinsP2015[2] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_60_100", title.c_str(), fNKplusBinsP2015[3], fKplusBinsP2015[3] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_100_140", title.c_str(), fNKplusBinsP2015[4], fKplusBinsP2015[4] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_140_180", title.c_str(), fNKplusBinsP2015[5], fKplusBinsP2015[5] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_180_240", title.c_str(), fNKplusBinsP2015[6], fKplusBinsP2015[6] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoKPlusP2015.push_back( new TH1F( "kplus_240_300", title.c_str(), fNKplusBinsP2015[7], fKplusBinsP2015[7] ) );

  // K-, 0-20 
  fNKminusBinsP2015[0] = 14;
  fKminusBinsP2015[0] = new double[fNKminusBinsP2015[0]+1];
  for ( int i=0; i<=fNKminusBinsP2015[0]; ++i )
  {
     fKminusBinsP2015[0][i] = tmp_k[i];
  }
  // K-, 20-40
  fNKminusBinsP2015[1] = 16;
  fKminusBinsP2015[1] = new double[fNKminusBinsP2015[1]+1];
  for ( int i=0; i<=fNKminusBinsP2015[1]; ++i )
  {
     fKminusBinsP2015[1][i] = tmp_k[i];
  }
  // K-, 40-60
  fNKminusBinsP2015[2] = 15;
  fKminusBinsP2015[2] = new double[fNKminusBinsP2015[2]+1];
  for ( int i=0; i<=fNKminusBinsP2015[2]; ++i )
  {
     fKminusBinsP2015[2][i] = tmp_k[i];
  }
  // K-, 60-100
  fNKminusBinsP2015[3] = 15;
  fKminusBinsP2015[3] = new double[fNKminusBinsP2015[3]+1];
  for ( int i=0; i<=fNKminusBinsP2015[3]; ++i )
  {
     fKminusBinsP2015[3][i] = tmp_k[i];
  }
  // K-, 100-140
  fNKminusBinsP2015[4] = 12;
  fKminusBinsP2015[4] = new double[fNKminusBinsP2015[4]+1];
  for ( int i=0; i<=fNKminusBinsP2015[4]; ++i )
  {
     fKminusBinsP2015[4][i] = tmp_k[i];
  }
  // K-, 140-180
  fNKminusBinsP2015[5] = 9;
  fKminusBinsP2015[5] = new double[fNKminusBinsP2015[5]+1];
  for ( int i=0; i<=fNKminusBinsP2015[5]; ++i )
  {
     fKminusBinsP2015[5][i] = tmp_k[i];
  }
  // K-, 180-240
  fNKminusBinsP2015[6] = 7;
  fKminusBinsP2015[6] = new double[fNKminusBinsP2015[6]+1];
  for ( int i=0; i<=fNKminusBinsP2015[6]; ++i )
  {
     fKminusBinsP2015[6][i] = tmp_k[i];
  }
  // K-, 240-300
  fNKminusBinsP2015[7] = 6;
  fKminusBinsP2015[7] = new double[fNKminusBinsP2015[7]+1];
  for ( int i=0; i<=fNKminusBinsP2015[7]; ++i )
  {
     fKminusBinsP2015[7][i] = tmp_k[i];
  }
  
  outcome = " -> X + K- ";

  title = htitle + outcome + " (0<theta<20 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_0_20", title.c_str(), fNKminusBinsP2015[0], fKminusBinsP2015[0] ) );
  title = htitle + outcome + " (20<theta<40 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_20_40", title.c_str(), fNKminusBinsP2015[1], fKminusBinsP2015[1] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_40_60", title.c_str(), fNKminusBinsP2015[2], fKminusBinsP2015[2] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_60_100", title.c_str(), fNKminusBinsP2015[3], fKminusBinsP2015[3] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_100_140", title.c_str(), fNKminusBinsP2015[4], fKminusBinsP2015[4] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_140_180", title.c_str(), fNKminusBinsP2015[5], fKminusBinsP2015[5] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_180_240", title.c_str(), fNKminusBinsP2015[6], fKminusBinsP2015[6] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoKMinusP2015.push_back( new TH1F( "kminus_240_300", title.c_str(), fNKminusBinsP2015[7], fKminusBinsP2015[7] ) );


  // secondary K0s
  //
  
  double tmp_k0s_0[] = {0.4, 3.4, 6.4, 9.4, 12.4};
  double tmp_k0s_1[] = {0.4, 2.4, 4.4, 6.4, 8.4, 10.4, 12.4};
  double tmp_k0s_2[] = {0.4, 2.4, 4.4, 6.4, 9.4};
  
  // 0-40
  fNK0sBinsP2015[0] = 4;
  fK0sBinsP2015[0] = new double[fNK0sBinsP2015[0]+1];
  for ( int i=0; i<=fNK0sBinsP2015[0]; ++i )
  {
     fK0sBinsP2015[0][i] = tmp_k0s_0[i];
  }
  // 40-60
  fNK0sBinsP2015[1] = 4;
  fK0sBinsP2015[1] = new double[fNK0sBinsP2015[1]+1];
  for ( int i=0; i<=fNK0sBinsP2015[1]; ++i )
  {
     fK0sBinsP2015[1][i] = tmp_k0s_0[i];
  }
  // 60-100
  fNK0sBinsP2015[2] = 6;
  fK0sBinsP2015[2] = new double[fNK0sBinsP2015[2]+1];
  for ( int i=0; i<=fNK0sBinsP2015[2]; ++i )
  {
     fK0sBinsP2015[2][i] = tmp_k0s_1[i];
  }
  // 100-140
  fNK0sBinsP2015[3] = 5;
  fK0sBinsP2015[3] = new double[fNK0sBinsP2015[3]+1];
  for ( int i=0; i<=fNK0sBinsP2015[3]; ++i )
  {
     fK0sBinsP2015[3][i] = tmp_k0s_1[i];
  }
  // 140-180
  fNK0sBinsP2015[4] = 4;
  fK0sBinsP2015[4] = new double[fNK0sBinsP2015[4]+1];
  for ( int i=0; i<=fNK0sBinsP2015[4]; ++i )
  {
     fK0sBinsP2015[4][i] = tmp_k0s_2[i];
  }
  // 180-240
  fNK0sBinsP2015[5] = 3;
  fK0sBinsP2015[5] = new double[fNK0sBinsP2015[5]+1];
  for ( int i=0; i<=fNK0sBinsP2015[5]; ++i )
  {
     fK0sBinsP2015[5][i] = tmp_k0s_2[i];
  }
  // 240-300
  fNK0sBinsP2015[6] = 2;
  fK0sBinsP2015[6] = new double[fNK0sBinsP2015[6]+1];
  for ( int i=0; i<=fNK0sBinsP2015[6]; ++i )
  {
     fK0sBinsP2015[6][i] = tmp_k0s_2[i];
  }
  
  outcome = " -> X + K0s ";

  title = htitle + outcome + " (0<theta<40 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_0_40", title.c_str(), fNK0sBinsP2015[0], fK0sBinsP2015[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_40_60", title.c_str(), fNK0sBinsP2015[1], fK0sBinsP2015[1] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_60_100", title.c_str(), fNK0sBinsP2015[2], fK0sBinsP2015[2] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_100_140", title.c_str(), fNK0sBinsP2015[3], fK0sBinsP2015[3] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_140_180", title.c_str(), fNK0sBinsP2015[4], fK0sBinsP2015[4] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_180_240", title.c_str(), fNK0sBinsP2015[5], fK0sBinsP2015[5] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoK0sP2015.push_back( new TH1F( "k0s_240_300", title.c_str(), fNK0sBinsP2015[6], fK0sBinsP2015[6] ) );

  // secondary Lambda
  //
  
  double tmp_l_0[] = {0.4, 4.4, 7.4, 10.4, 13.4, 17.4, 21.4};
  double tmp_l_1[] = {0.4, 2.9, 4.4, 6.4, 8.4, 10.4, 12.4, 18};
  double tmp_l_2[] = {0.4, 2.9, 4.4, 6.4, 8, 11.4};
  double tmp_l_3[] = {0.4, 2.4, 4.4, 6.4, 11};
  double tmp_l_4[] = {0.4, 2.4, 4.4, 10.4};
  
  // 0-40
  fNLambdaBinsP2015[0] = 6;
  fLambdaBinsP2015[0] = new double[fNLambdaBinsP2015[0]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[0]; ++i )
  {
     fLambdaBinsP2015[0][i] = tmp_l_0[i];
  }
  // 40-60
  fNLambdaBinsP2015[1] = 6;
  fLambdaBinsP2015[1] = new double[fNLambdaBinsP2015[1]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[1]; ++i )
  {
     fLambdaBinsP2015[1][i] = tmp_l_0[i];
  }
  // 60-100
  fNLambdaBinsP2015[2] = 7;
  fLambdaBinsP2015[2] = new double[fNLambdaBinsP2015[2]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[2]; ++i )
  {
     fLambdaBinsP2015[2][i] = tmp_l_1[i];
  }
  // 100-140
  fNLambdaBinsP2015[3] = 6;
  fLambdaBinsP2015[3] = new double[fNLambdaBinsP2015[3]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[3]; ++i )
  {
     fLambdaBinsP2015[3][i] = tmp_l_1[i];
  }
  // 140-180
  fNLambdaBinsP2015[4] = 5;
  fLambdaBinsP2015[4] = new double[fNLambdaBinsP2015[4]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[4]; ++i )
  {
     fLambdaBinsP2015[4][i] = tmp_l_2[i];
  }
  // 180-240
  fNLambdaBinsP2015[5] = 4;
  fLambdaBinsP2015[5] = new double[fNLambdaBinsP2015[5]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[5]; ++i )
  {
     fLambdaBinsP2015[5][i] = tmp_l_3[i];
  }
  // 240-300
  fNLambdaBinsP2015[6] = 3;
  fLambdaBinsP2015[6] = new double[fNLambdaBinsP2015[6]+1];
  for ( int i=0; i<=fNLambdaBinsP2015[6]; ++i )
  {
     fLambdaBinsP2015[6][i] = tmp_l_4[i];
  }

  outcome = " -> X + Lambda ";

  title = htitle + outcome + " (0<theta<40 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_0_40", title.c_str(), fNLambdaBinsP2015[0], fLambdaBinsP2015[0] ) );
  title = htitle + outcome + " (40<theta<60 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_40_60", title.c_str(), fNLambdaBinsP2015[1], fLambdaBinsP2015[1] ) );
  title = htitle + outcome + " (60<theta<100 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_60_100", title.c_str(), fNLambdaBinsP2015[2], fLambdaBinsP2015[2] ) );
  title = htitle + outcome + " (100<theta<140 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_100_140", title.c_str(), fNLambdaBinsP2015[3], fLambdaBinsP2015[3] ) );
  title = htitle + outcome + " (140<theta<180 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_140_180", title.c_str(), fNLambdaBinsP2015[4], fLambdaBinsP2015[4] ) );
  title = htitle + outcome + " (180<theta<240 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_180_240", title.c_str(), fNLambdaBinsP2015[5], fLambdaBinsP2015[5] ) );
  title = htitle + outcome + " (240<theta<300 [mrad]) ";
  fHistoLambdaP2015.push_back( new TH1F( "lambda_240_300", title.c_str(), fNLambdaBinsP2015[6], fLambdaBinsP2015[6] ) );

}

Tst23NA61Histo::~Tst23NA61Histo()
{

   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for ( size_t i=0; i<fHistoSecKPlus.size(); i++ )
   {
      delete fHistoSecKPlus[i];
   }
   for ( size_t i=0; i<fHistoSecPiPlus2.size(); i++ )
   {
      delete fHistoSecPiPlus2[i];
   }
   for (size_t i=0; i<fHistoProtonP2015.size(); ++i )
   {
      delete fHistoProtonP2015[i];
   }

   for (size_t i=0; i<fHistoPiPlusP2015.size(); ++i )
   {
      delete fHistoPiPlusP2015[i];
   }
   for (size_t i=0; i<fHistoPiMinusP2015.size(); ++i )
   {
      delete fHistoPiMinusP2015[i];
   }
   for (size_t i=0; i<fHistoKPlusP2015.size(); ++i )
   {
      delete fHistoKPlusP2015[i];
   }
   for (size_t i=0; i<fHistoKMinusP2015.size(); ++i )
   {
      delete fHistoKMinusP2015[i];
   }
   for (size_t i=0; i<fHistoK0sP2015.size(); ++i )
   {
      delete fHistoK0sP2015[i];
   }
   for (size_t i=0; i<fHistoLambdaP2015.size(); ++i )
   {
      delete fHistoLambdaP2015[i];
   }

}


void Tst23NA61Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

   G4int NSec = aChange->GetNumberOfSecondaries();
 
   if ( NSec <= 0 ) return;
  
   bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
   
   if ( isFirst && fDoResDecay ) AccountForResDecay( aChange );
      
   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );
     
   const G4Track*           trk = 0;  
   const G4DynamicParticle* sec = 0;
       
   
   for (G4int i=0; i<NSec; i++) 
   {
        trk = aChange->GetSecondary(i);
        sec = trk->GetDynamicParticle();			
      const G4String& pname = sec->GetDefinition()->GetParticleName();
	
      // G4ThreeVector mom = sec->GetMomentum();
	
      double theta = sec->GetMomentum().theta();
      theta /= mrad;
	
      double pmom = sec->GetTotalMomentum() / GeV ;
	
      int id = -1;
      int id1 = -1;
      double dth1 = 1.;
	
      int NHProSize = fHistoSecProton.size();
      int NHPimSize = fHistoSecPiMinus.size();
      int NHPipSize = fHistoSecPiPlus.size();
      
      int NHProP2015Size = fHistoProtonP2015.size();
      int NHPipP2015Size = fHistoPiPlusP2015.size();
      int NHPimP2015Size = fHistoPiMinusP2015.size();
      int NHKpP2015Size = fHistoKPlusP2015.size();
      int NHKmP2015Size = fHistoKMinusP2015.size();
      int NHK0sP2015Size = fHistoK0sP2015.size();
      int NHLamP2015Size = fHistoLambdaP2015.size();

      if ( pname == "proton" || 
           pname == "pi+" || pname == "pi-" || 
	   pname == "kaon+" || pname == "kaon-" )
      {
	if ( theta >= 0. && theta < 10. )
	{
	   id = 0;
	   id1 = 0;
	   dth1 = 10. * mrad;
	}
	else if ( theta >=0. && theta < 20. )
	{
	   id = 0;
	   id1 = 1;
	   dth1 = 10. * mrad;
	}
	else if ( theta >= 20. && theta < 40. )
	{
	   id = 1;
	   id1 = 2;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 40. && theta < 60. )
	{
	   id = 2;
	   id1 = 3;
	   dth1 = 20. * mrad;
	}
	else if ( theta >= 60. && theta < 100. )
	{
	   id = 3;
	   id1 = 4;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 100. && theta < 140. )
	{
	   id = 4;
	   id1 = 5;
	   dth1 = 40. * mrad;
	}
	else if ( theta >= 140. && theta < 180. )
	{
	   id = 5;
	   id1 = 6;
	   dth1 = 40. * mrad;
	}
	else if ( theta > 180. && theta < 240. )
	{
	   id = 6;
	   id1 = 7;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 240. && theta < 300. )
	{
	   id = 7;
	   id1 = 8;
	   dth1 = 60. * mrad;
	} 
	else if ( theta >= 300. && theta < 360. )
	{
	   id = 8;
	   id1 = 9;
	   dth1 = 60. * mrad;
	}
	else if ( theta >= 360. && theta < 420. )
	{
	   id = 9;
	   id1 = 10;
	   dth1 = 60. * mrad;
	}
	
	if ( pname == "proton" )
	{	   
	   if ( id >=0 && id < NHProSize ) fHistoSecProton[id]->Fill( pmom );
	   if ( id1 >= 0 && id1 < NHProP2015Size ) fHistoProtonP2015[id1]->Fill( pmom, 1./dth1 );
	}
	if ( pname == "pi-" )
	{
	   if ( id >=0 && id < NHPimSize ) fHistoSecPiMinus[id]->Fill( pmom );
	   if ( id1 >= 0 && id1 < NHPimP2015Size ) fHistoPiMinusP2015[id1]->Fill( pmom, 1./dth1 );
	}
	if ( pname == "pi+" )
	{
	   if ( id >= 0 && id < NHPipSize ) fHistoSecPiPlus[id]->Fill( pmom );
	   if ( id1 >= 0 && id1 < NHPipP2015Size ) fHistoPiPlusP2015[id1]->Fill( pmom, 1./dth1 );
	}
	if ( pname == "kaon+" )
	{
	   if ( id >=0 && id < NHKpP2015Size ) fHistoKPlusP2015[id]->Fill( pmom, 1./dth1 );
	}
	if ( pname == "kaon-" )
	{
	   if ( id >=0 && id < NHKmP2015Size ) fHistoKMinusP2015[id]->Fill( pmom, 1./dth1 );
	}
      }
      
      if ( pname == "proton" || pname == "pi-" || pname == "kaon-" ) continue;

      id = -1;

      if ( pname == "kaon+" || pname == "pi+" )
      {
         if ( theta >=20. && theta < 140. )
	 {
	    id = 0;
	 }
	 else if ( theta >= 140. && theta < 240. )
	 {
	    id = 1;
	 }
	 if ( pname == "kaon+" )
	 {
	    if ( id >= 0 )
	    {
	       fHistoSecKPlus[id]->Fill( pmom ) ;
	       fHistoSecKPlus[id+2]->Fill( pmom ) ;
	    }
	 }
	 else if ( pname == "pi+" )
	 {
	    if ( id >= 0 ) fHistoSecPiPlus2[id]->Fill( pmom ) ;
	 }	 
      }

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
	    if ( id >= 0 && id < NHK0sP2015Size ) fHistoK0sP2015[id]->Fill(pmom,1./dth1);
	 }
	 if ( pname == "lambda" )
	 {
	    if ( id >= 0 && id < NHLamP2015Size ) fHistoLambdaP2015[id]->Fill(pmom, 1./dth1);
	 }
	 
      }

   }
      
   return;
   
}

void Tst23NA61Histo::Write( G4int, G4double xsec )
{

   double norm = 1.;
   
   // NOTE: bear in mind that fHistoNSec->Integral() is is the number of entries WITHIN histo xmin-xmax
   //       while fHistoNSec->GetEntries() is the total number of entries, incl. under/overflow (i.e. the actual NEvents)
   //       however, AFTER the NORMALIZATION of fHistoSec, GetEntries() and Integral() give the same result
   //
   norm = (double)(fHistoNSec->GetEntries());
   
   fHistoNSec->Scale( 1./(double)norm, "width" );
   fHistoNSec->Write();
   
   for ( size_t i=0; i<fHistoSecProton.size(); ++i )
   {
      double xbin = fHistoSecProton[i]->GetBinWidth(1);
      double scale = 1. / ((double)norm * xbin) ;
      fHistoSecProton[i]->Scale(scale);
      fHistoSecProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      double scale = 1. / ((double)norm) ;
      fHistoSecPiMinus[i]->Scale(scale,"width");
      fHistoSecPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      double scale = 1. / ((double)norm) ;
      fHistoSecPiPlus[i]->Scale(scale,"width");
      fHistoSecPiPlus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecKPlus.size()-2; i++ )
   {
      double xbin = fHistoSecKPlus[i]->GetBinWidth(1);
      double scale = 1. / ((double)norm * xbin) ;
      fHistoSecKPlus[i]->Scale(scale);
      fHistoSecKPlus[i]->Write();
      fHistoSecKPlus[i+2]->Divide( fHistoSecPiPlus2[i] );
      fHistoSecKPlus[i+2]->Write();
   }

   double scale1 = xsec / ((double)norm);

   for ( size_t i=0; i<fHistoProtonP2015.size(); ++i )
   {
      fHistoProtonP2015[i]->Scale(scale1,"width");
      fHistoProtonP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoPiPlusP2015.size(); ++i )
   {
      fHistoPiPlusP2015[i]->Scale(scale1,"width");
      fHistoPiPlusP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoPiMinusP2015.size(); ++i )
   {
      fHistoPiMinusP2015[i]->Scale(scale1,"width");
      fHistoPiMinusP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoKPlusP2015.size(); ++i )
   {
      fHistoKPlusP2015[i]->Scale(scale1,"width");
      fHistoKPlusP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoKMinusP2015.size(); ++i )
   {
      fHistoKMinusP2015[i]->Scale(scale1,"width");
      fHistoKMinusP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoK0sP2015.size(); ++i )
   {
      fHistoK0sP2015[i]->Scale(scale1,"width");
      fHistoK0sP2015[i]->Write();
   }
   for ( size_t i=0; i<fHistoLambdaP2015.size(); ++i )
   {
      fHistoLambdaP2015[i]->Scale(scale1,"width");
      fHistoLambdaP2015[i]->Write();
   }
   
   return;

}

