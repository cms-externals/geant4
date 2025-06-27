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

#include "TestNA49Histo.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"


TestNA49Histo::TestNA49Histo( G4String htitle )
   : TstHistoSet( htitle )   
{

  std::string title;
  std::string outcome; 
  std::string hname; 
    
  fHistoNSec = new TH1D( "NSec", htitle.c_str(), 50, 0., 50. );
  
  outcome = " -> X + proton ";       
  title = htitle + outcome;

  // in principle, the exp.data are binned by 0.05 from -0.8 to 0.95, 
  // except in xF=(-0.1,0.1) where the binning is finer (0.025) 
  // so for now we do 0.025 everywhere, for simplicity
  //
  
  fNPBinsXF = 41; // 40+1 for bin edges  
  fPBinsXF = new double[fNPBinsXF]; 
  
  std::vector<double> bincenter;
  double center = 0.;
  
  for ( int ii=0; ii<15; ++ii )
  {
     center = -0.8 + 0.05*ii;
     bincenter.push_back( center );
  } 
  for (int ii=0; ii<8; ++ii )
  {
     center = -0.1 + 0.025*(ii+1);
     bincenter.push_back( center );
  }
  for ( int ii=0; ii<17; ++ii ) // 17 = fNPBinsXF-1-23, 23 - from previous 2 loops
  {
     center = 0.1 + 0.05*(ii+1);
     bincenter.push_back( center );
  }
    
  fPBinsXF[0] = -0.825;
  for ( int ii=1; ii<fNPBinsXF-1; ++ii )
  {
     fPBinsXF[ii] = (bincenter[ii-1]+bincenter[ii])/2.;
  }
  fPBinsXF[fNPBinsXF-1] = fPBinsXF[fNPBinsXF-2] + 0.05;
  
  bincenter.clear();
     
  fNPBinsPT = 20; // 19+1 for edges
  for ( int ii=0; ii<10; ++ii )
  {
     center = 0.05 + ii*0.05;
     bincenter.push_back( center );
  }
  for ( int ii=0; ii<4; ++ii )
  {
     center = 0.6 + 0.1*ii;
     bincenter.push_back( center );
  }
  for ( int ii=0; ii<5; ++ii )
  {
     center = 1.1 + ii*0.2;
     bincenter.push_back( center ); 
  }
    
  fPBinsPT = new double[fNPBinsPT];
  fPBinsPT[0] = 0.025;
  for ( int ii=1; ii<fNPBinsPT-1; ++ii )
  {
     fPBinsPT[ii] = (bincenter[ii-1]+bincenter[ii])/2.;
  }
  fPBinsPT[fNPBinsPT-1] = fPBinsPT[fNPBinsPT-2] + 0.2;
    
  bincenter.clear();

  std::cout << " booking histo for proton integrated spectra " << std::endl;
  
  
   fHistoSecProton.push_back( new TH1D( "proton_dXSecdxF", title.c_str(), fNPBinsXF-1, fPBinsXF ) );	
   fHistoSecProton.push_back( new TH1D( "proton_dNdxF",    title.c_str(), fNPBinsXF-1, fPBinsXF ) );	
   fHistoSecProton.push_back( new TProfile( "proton_pT",   title.c_str(), fNPBinsXF-1, fPBinsXF, 0., 10. ) );	
   fHistoSecProton.push_back( new TProfile( "proton_pT2",  title.c_str(), fNPBinsXF-1, fPBinsXF, 0., 10. ) );	

   std::cout << " book proton histo for ddiff pt spectra " << std::endl;

  for ( int nb=0; nb<fNPBinsXF-1; ++nb ) // FIXME !!! fNPBinsXF is the actual # of bins
                                          // while all other are N+1 !!!
  {
     double xF = (fPBinsXF[nb]+fPBinsXF[nb+1])/2.;
     std::ostringstream osCount1;
     osCount1 << xF;
     std::ostringstream osCount2;
     osCount2 << nb;
     G4String subtitle = ", xF = " + osCount1.str();
     hname = "pTpro" + osCount2.str();
     fHistoPTProton.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), fNPBinsPT-1, fPBinsPT ) );     
  }
       
  outcome = " -> X + antiproton ";
  title = htitle + outcome;
  
//  double pbarbins[] = { -0.3, -0.2, -0.15, -0.1, -0.075, -0.05, -0.025, 
//                         0.0,
//			 0.025, 0.05, 0.1, 0.15, 0.2, 0.3 }
  // we specify bin left edge, to make the bin center match the exp.data
  //
  double pbarbins[] = { -0.55, -0.45, -0.35, -0.25, -0.175, -0.125, -0.0875, -0.0625, -0.0375, 
                         -0.0125,
			 0.0125, 0.0375, 0.075, 0.125, 0.175, 0.25, 0.35, 0.45, 0.55 };
  int npbarbins = sizeof(pbarbins) / sizeof(double) - 1;
  
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dXSecdxF", title.c_str(), 80, -1., 1. ) );	
  fHistoSecAntiProton.push_back( new TH1D( "antiproton_dNdxF",    title.c_str(), npbarbins, pbarbins ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT",   title.c_str(), npbarbins, pbarbins, 0., 10. ) );	
  fHistoSecAntiProton.push_back( new TProfile( "antiproton_pT2",  title.c_str(), npbarbins, pbarbins, 0., 10. ) );	
  
   fNPbarBinsXF = 13; 
   fPbarBinsXF = new double[fNPbarBinsXF+1]; 
   fPbarBinsXF[0] = -0.225;
   fPbarBinsXF[1] = -0.175;
   fPbarBinsXF[2] = -0.125;
   fPbarBinsXF[3] = -0.0875;
   fPbarBinsXF[4] = -0.0625;
   fPbarBinsXF[5] = -0.0375;
   fPbarBinsXF[6] = -0.0125;
   fPbarBinsXF[7] =  0.0125;
   fPbarBinsXF[8] =  0.0375;
   fPbarBinsXF[9] =  0.075;
   fPbarBinsXF[10]=  0.125;
   fPbarBinsXF[11]=  0.175;
   fPbarBinsXF[12]=  0.25;
   fPbarBinsXF[13]=  0.35;

   std::vector<std::string> PbarXFLabel;
   PbarXFLabel.push_back("-0.200");
   PbarXFLabel.push_back("-0.150");
   PbarXFLabel.push_back("-0.100");
   PbarXFLabel.push_back("-0.075");
   PbarXFLabel.push_back("-0.050");
   PbarXFLabel.push_back("-0.025");
   PbarXFLabel.push_back("0.000");
   PbarXFLabel.push_back("0.025");
   PbarXFLabel.push_back("0.050");
   PbarXFLabel.push_back("0.100");
   PbarXFLabel.push_back("0.150");
   PbarXFLabel.push_back("0.200");
   PbarXFLabel.push_back("0.300");
   
   fNPbarBinsPT = new int[fNPbarBinsXF];
   fNPbarBinsPT[0] = 9; // -0.2
   fNPbarBinsPT[1] = 10; // -0.15
   fNPbarBinsPT[2] = 10; // -0.1
   fNPbarBinsPT[3] = 8; // -0.075
   fNPbarBinsPT[4] = 10; // -0.05
   fNPbarBinsPT[5] = 8; // -0.025
   fNPbarBinsPT[6] = 11; // 0.
   fNPbarBinsPT[7] = 8; // 0.025
   fNPbarBinsPT[8] = 10; // 0.05
   fNPbarBinsPT[9] = 11; // 0.1
   fNPbarBinsPT[10] = 9; // 0.15
   fNPbarBinsPT[11] = 11; // 0.2
   fNPbarBinsPT[12] = 6; // 0.3
   
   fPbarBinsPT = new double*[fNPbarBinsXF];
   for ( int i=0; i<fNPbarBinsXF; ++i )
   {
      fPbarBinsPT[i] = new double[fNPbarBinsPT[i]+1];
      if ( i == 12 )
      {
         fPbarBinsPT[i][0] = 0.0;
	 fPbarBinsPT[i][1] = 0.2;
	 fPbarBinsPT[i][2] = 0.4;
	 fPbarBinsPT[i][3] = 0.6;
	 fPbarBinsPT[i][4] = 0.8;
	 fPbarBinsPT[i][5] = 1.0;
	 fPbarBinsPT[i][6] = 1.2;
	 continue;
      }
      else
      {
         fPbarBinsPT[i][0] = 0.05;
         fPbarBinsPT[i][1] = 0.15;
         fPbarBinsPT[i][2] = 0.25;
         fPbarBinsPT[i][3] = 0.35;
         fPbarBinsPT[i][4] = 0.45;
	 fPbarBinsPT[i][5] = 0.55;
	 fPbarBinsPT[i][6] = 0.65;
	 fPbarBinsPT[i][7] = 0.75;
	 fPbarBinsPT[i][8] = 1.0;
      }
      if ( i == 0 ) // override
      {
	 fPbarBinsPT[i][6] = 0.8;
	 fPbarBinsPT[i][7] = 1.0;
	 fPbarBinsPT[i][8] = 1.2;
	 fPbarBinsPT[i][9] = 1.4; 
      }
      else
      {
         for ( int j=8; j<fNPbarBinsPT[i]; ++j )
	 {
	    fPbarBinsPT[i][j+1] = fPbarBinsPT[i][j] + 0.2;
	 }     
      }
   }
   
   for ( int nb=0; nb<fNPbarBinsXF; ++nb ) 
   {
     std::ostringstream osBin;
     osBin << nb;
     std::string subtitle = ", xF = " + PbarXFLabel[nb];
     hname = "pTpbar" + osBin.str();
     fHistoPTAntiProton.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), fNPbarBinsPT[nb], fPbarBinsPT[nb] ) );     
   }

   // pions

  fNPiBinsXF = 30;
  fPiBinsXF = new double[fNPiBinsXF];  
  fPiBinsXF[0] = -0.550;
  fPiBinsXF[1] = -0.450;
  fPiBinsXF[2] = -0.350; 
  fPiBinsXF[3] = -0.275;
  fPiBinsXF[4] = -0.225;
  fPiBinsXF[5] = -0.175;
  fPiBinsXF[6] = -0.1375;
  fPiBinsXF[7] = -0.1125;
  fPiBinsXF[8] = -0.0875; 
  fPiBinsXF[9] = -0.0625;
  fPiBinsXF[10] = -0.045;
  fPiBinsXF[11] = -0.035;
  fPiBinsXF[12] = -0.025;
  fPiBinsXF[13] = -0.015; 
  fPiBinsXF[14] = -0.005;
  fPiBinsXF[15] =  0.005;  
  fPiBinsXF[16] =  0.015;
  fPiBinsXF[17] =  0.025;
  fPiBinsXF[18] =  0.035;  
  fPiBinsXF[19] =  0.045;
  fPiBinsXF[20] =  0.0625;  
  fPiBinsXF[21] =  0.0875;
  fPiBinsXF[22] =  0.1125;
  fPiBinsXF[23] =  0.1375;
  fPiBinsXF[24] =  0.175; 
  fPiBinsXF[25] =  0.225;
  fPiBinsXF[26] =  0.275;
  fPiBinsXF[27] =  0.350;
  fPiBinsXF[28] =  0.450;
  fPiBinsXF[29] =  0.55;  
  
  fNPiBinsPT = 17;
  fPiBinsPT = new double[fNPiBinsPT];
  fPiBinsPT[0] = 0.025;
  fPiBinsPT[1] = 0.075;
  fPiBinsPT[2] = 0.125;
  fPiBinsPT[3] = 0.175;
  fPiBinsPT[4] = 0.225;
  fPiBinsPT[5] = 0.275;
  fPiBinsPT[6] = 0.35;
  fPiBinsPT[7] = 0.45;
  fPiBinsPT[8] = 0.55;
  fPiBinsPT[9] = 0.65;
  fPiBinsPT[10] = 0.75;
  fPiBinsPT[11] = 0.85;
  fPiBinsPT[12] = 0.95; 
  fPiBinsPT[13] = 1.1;
  fPiBinsPT[14] = 1.3;
  fPiBinsPT[15] = 1.5;
  fPiBinsPT[16] = 1.7;

  outcome = " -> X + pi- ";  
  title = htitle + outcome;
//  double pionbins[] = { -0.500, -0.400, -0.300, -0.250, -0.200, -0.150, -0.125, -0.100, -0.075, 
//                     -0.050, -0.040, -0.030, -0.020, -0.010,  
//		      0.0,
//		      0.010,  0.020,  0.030,  0.040,  0.050,
//		      0.075,  0.100,  0.125,  0.150,  0.200,  0.250,  0.300,  0.400,  0.500 }; 
		       
  // the idea is to make the bin center correspond to the number in the NA49 data file(s)
  //
/*
  double pibins[] = { -0.550, -0.450, -0.350, -0.275, -0.225, -0.175, -0.1375, -0.1125, -0.0875, 
                     -0.0625, -0.045, -0.035, -0.025, -0.015,  
		      -0.005,
		      0.005,  0.015,  0.025,  0.035,  0.045,
		      0.0625,  0.0875,  0.1125,  0.1375,  0.175, 0.225,  0.275,  0.350,  0.450, 0.55 };  
  int npibins = sizeof(pibins) / sizeof(double) - 1;
*/

//  double pibins_pt = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8 }
/*
  double pibins_pt[] = { 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 
                         0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 
		         1.1, 1.3, 1.5, 1.7 };
  int npibins_pt = sizeof(pibins_pt)/sizeof(double) - 1; 
*/
  std::cout << " book pi- histo for integrated spectra " << std::endl; 
  
  fHistoSecPiMinus.push_back( new TH1D( "pimimus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiMinus.push_back( new TH1D( "piminus_dNdxF",    title.c_str(), fNPiBinsXF-1, fPiBinsXF ) );
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT",   title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	
//  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), npibins, pibins, 0., 10. ) );	
  fHistoSecPiMinus.push_back( new TProfile( "piminus_pT2",  title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	

   std::cout << " book pi- histo for ddiff pt spectra " << std::endl;

  for ( int nb=0; nb<fNPiBinsXF-1; ++nb )
  {
     double xF = (fPiBinsXF[nb]+fPiBinsXF[nb+1])/2.;
     std::ostringstream osCount1;
     osCount1 << xF;
     std::ostringstream osCount2;
     osCount2 << nb;
     G4String subtitle = ", xF = " + osCount1.str();
     hname = "pTpim" + osCount2.str();
//     fHistoPTPiMinus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), npibins_pt, pibins_pt) );     
     fHistoPTPiMinus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), fNPiBinsPT-1, fPiBinsPT ) );     
  }

  outcome = " -> X + pi+ ";
  title = htitle + outcome;

  std::cout << " book pi+ histo for integrated spectra " << std::endl; 

  fHistoSecPiPlus.push_back( new TH1D( "piplus_dXSecdxF", title.c_str(), 40, -1., 1. ) );	
  fHistoSecPiPlus.push_back( new TH1D( "piplus_dNdxF",    title.c_str(), fNPiBinsXF-1, fPiBinsXF ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT",   title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	
  fHistoSecPiPlus.push_back( new TProfile( "piplus_pT2",  title.c_str(), fNPiBinsXF-1, fPiBinsXF, 0., 10. ) );	

  std::cout << " book pi+ histo for ddiff spectra " << std::endl; 

  for ( int nb=0; nb<fNPiBinsXF-1; ++nb )
  {
     double xF = (fPiBinsXF[nb]+fPiBinsXF[nb+1])/2.;
     std::ostringstream osCount1;
     osCount1 << xF;
     std::ostringstream osCount2;
     osCount2 << nb;
     G4String subtitle = ", xF = " + osCount1.str();
     hname = "pTpip" + osCount2.str();
//     fHistoPTPiPlus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), npibins_pt, pibins_pt ) );     
     fHistoPTPiPlus.push_back( new TH1D( hname.c_str(), (title+subtitle).c_str(), fNPiBinsPT-1, fPiBinsPT ) );     
  }

  outcome = " -> X + neutron ";
  title = htitle + outcome;

  // this binning is reasonable because the data go by 0.1 in (0.1-0.6), then just 2 bins of 0.15
  //
  fHistoSecNeutron.push_back( new TH1D( "neutron_dNdxF",    title.c_str(), 10, 0.05, 1.05 ) );	
    
  fInteraction = new G4VParticleChange();

}

TestNA49Histo::~TestNA49Histo()
{

   if (fHistoNSec) delete fHistoNSec;
   
   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecAntiProton.size(); i++ )
   {
      delete fHistoSecAntiProton[i];
   }
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for ( size_t i=0; i<fHistoSecNeutron.size(); i++ )
   {
      delete fHistoSecNeutron[i];
   }

   for ( size_t i=0; i<fHistoPTProton.size(); ++i )
   {
      delete fHistoPTProton[i];
   }
   
   delete [] fPBinsXF;
   delete [] fPBinsPT;
   
   for ( size_t i=0; i<fHistoPTAntiProton.size(); ++i )
   {
      delete fHistoPTAntiProton[i];
   }
   
   for ( int i=0; i<fNPbarBinsXF; ++i )
   {
      delete [] fPbarBinsPT[i];
   }
   
   delete [] fPbarBinsPT;
   delete [] fPbarBinsXF;
   
   for ( size_t i=0; i<fHistoPTPiMinus.size(); ++i )
   {
      delete fHistoPTPiMinus[i];
   }
   for ( size_t i=0; i<fHistoPTPiPlus.size(); ++i )
   {
      delete fHistoPTPiPlus[i];
   }
   
   delete [] fPiBinsXF;
   delete [] fPiBinsPT;
   
   if ( fInteraction )
   {
      int nsc = fInteraction->GetNumberOfSecondaries();
      if ( nsc > 0  )
      {
         for ( int i=0; i<nsc; i++) 
         {   
            delete fInteraction->GetSecondary(i);
         } 
      }
      fInteraction->Clear();
      delete fInteraction;    
   }
      
}

void TestNA49Histo::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& labp ) 
{

   G4int NSec = fInteraction->GetNumberOfSecondaries();

   if ( NSec > 0 )
   {
      for ( int i=0; i<NSec; i++ )
      {
         delete fInteraction->GetSecondary(i);
      }
      fInteraction->Clear();
   }
   
   NSec = aChange->GetNumberOfSecondaries();
   
   if ( NSec <= 0 ) return;
   
   for ( int i=0; i<NSec; i++ )
   {
      G4Track* trk = new G4Track( *(aChange->GetSecondary(i)) );
      fInteraction->AddSecondary( trk );
   }

   if ( fDoResDecay ) AccountForResDecay( fInteraction );
   
   double SQRT_S = 17.2;
   
   G4ThreeVector boostLabp = labp.boostVector();
      
   NSec = fInteraction->GetNumberOfSecondaries();
   
   fHistoNSec->Fill( (double)NSec );
   
   const G4Track*           trk = 0; 
   const G4DynamicParticle* sec = 0;
       
   for (G4int i=0; i<NSec; i++) 
   {
        trk = fInteraction->GetSecondary(i);
	sec = trk->GetDynamicParticle();			
	const G4String& pname = sec->GetDefinition()->GetParticleName();
	
	G4ThreeVector mom = sec->GetMomentum() / GeV ;
	double mass = sec->GetDefinition()->GetPDGMass() / GeV;
	double ekin = sec->GetKineticEnergy() / GeV ;
	
	G4LorentzVector boostSec( mom, ekin+mass );
	boostSec.boost(-boostLabp);
	
	double xF  = 2 * (boostSec.z()) / SQRT_S;
	
/*
	if ( xF < -1. || xF > 1. )
	{
	   std::cout << " xF = " << xF << std::endl;
	   std::cout << " mom: " << mom.x() << " " << mom.y() << " " << mom.z() << " ekin+mass = " << ekin+mass << std::endl;
	} 
*/	
	double pT  = mom.perp() ;
	double pT2 = pT * pT ;
        
	if ( pname == "neutron" )
	{
	   fHistoSecNeutron[0]->Fill( xF );
	}
	else if ( pname == "pi-" )
	{
	   fHistoSecPiMinus[1]->Fill( xF );
	   fHistoSecPiMinus[2]->Fill( xF, pT );
	   fHistoSecPiMinus[3]->Fill( xF, pT2 );
	   int nb = -1;
	   for ( int ib=0; ib<fNPiBinsXF-1; ++ib )
	   {
	         if ( xF >= fPiBinsXF[ib] && xF < fPiBinsXF[ib+1] )
		 {
		    nb = ib;
		    break;
		 }
	   }
	   if ( nb == -1 ) continue;
	   // calculate weight as Epart/dP3 (see NA49 papers on the SPS hadrons site at CERN) 
	   double wei = CalculateBinWeight( labp, pname, ekin, mass, pT, nb, SQRT_S );
	   fHistoPTPiMinus[nb]->Fill( pT, wei ); 
	}
	else if ( pname == "pi+" )
	{
	   fHistoSecPiPlus[1]->Fill( xF );
	   fHistoSecPiPlus[2]->Fill( xF, pT );
	   fHistoSecPiPlus[3]->Fill( xF, pT2 );
	   int nb = -1;
	   for ( int ib=0; ib<fNPiBinsXF-1; ++ib )
	   {
	         if ( xF > fPiBinsXF[ib] && xF <= fPiBinsXF[ib+1] )
		 {
		    nb = ib;
		    break;
		 }
	   }
	   if ( nb == -1 ) continue;
	   // calculate weight as Epart/dP3 (see NA49 papers on the SPS hadrons site at CERN) 
	   double wei = CalculateBinWeight( labp, pname, ekin, mass, pT, nb, SQRT_S );
	   fHistoPTPiPlus[nb]->Fill( pT, wei ); 
	}
	else if ( pname == "proton" )
	{
	   fHistoSecProton[1]->Fill( xF );
	   fHistoSecProton[2]->Fill( xF, pT );
	   fHistoSecProton[3]->Fill( xF, pT2 );
	   int nb = -1;
	   for ( int ib=0; ib<fNPBinsXF; ++ib ) // # number of proton bins are N and the size is N+1, while for pions # is N+1
	                                        // need to make it uniform !!!
	   {
	         if ( xF > fPBinsXF[ib] && xF <= fPBinsXF[ib+1] )
		 {
		    nb = ib;
		    break;
		 }
	   }
	   if ( nb == -1 ) continue;
	   // calculate weight as Epart/dP3 (see NA49 papers on the SPS hadrons site at CERN) 
	   double wei = CalculateBinWeight( labp, pname, ekin, mass, pT, nb, SQRT_S );
	   fHistoPTProton[nb]->Fill( pT, wei ); 
	}
	else if ( pname == "anti_proton" )
	{
	   fHistoSecAntiProton[1]->Fill( xF );
	   fHistoSecAntiProton[2]->Fill( xF, pT );
	   fHistoSecAntiProton[3]->Fill( xF, pT2 );
	   int nb = -1;
	   for ( int ib=0; ib<fNPbarBinsXF; ++ib )
	   {
	      if ( xF > fPbarBinsXF[ib] && xF <= fPbarBinsXF[ib+1] )
	      {
	         nb = ib;
		 break;
	      }
	   }
	   if ( nb == -1 ) continue;
	   double wei = CalculateBinWeight( labp, pname, ekin, mass, pT, nb, SQRT_S );
	   fHistoPTAntiProton[nb]->Fill( pT, wei );
	}	
   }
      
   return;
   
}

void TestNA49Histo::Write( G4int stat, G4double )
{

   fHistoNSec->Scale( 1./((double)stat),"width" ) ;
   fHistoNSec->Write();
   
   for ( size_t i=0; i<fHistoSecProton.size(); ++i )   
   {
      if ( i <2 ) fHistoSecProton[i]->Scale(1./((double)stat),"width");
      fHistoSecProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoPTProton.size(); ++i )
   {
      fHistoPTProton[i]->Scale( 1./((double)stat) ); // Note: NO scaling with "width" because it had to be E/dP3
                                                     //       which is already taken into account as a weight 
      fHistoPTProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecAntiProton.size(); ++i )   
   {
      if ( i < 2 ) fHistoSecAntiProton[i]->Scale(1./((double)stat),"width");
      fHistoSecAntiProton[i]->Write();
   }
   for ( size_t i=0; i<fHistoPTAntiProton.size(); ++i )
   {
      fHistoPTAntiProton[i]->Scale( 1./((double)stat) ); // Note: NO scaling with "width" because it had to be E/dP3
                                                         //       which is already taken into account as a weight
      fHistoPTAntiProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      if ( i < 2 ) fHistoSecPiMinus[i]->Scale(1./((double)stat),"width");
      fHistoSecPiMinus[i]->Write();
   }
   for ( size_t i=0; i<fHistoPTPiMinus.size(); ++i )
   {
      fHistoPTPiMinus[i]->Scale( 1./((double)stat) ); // Note: NO scaling with "width" because it had to be E/dP3
                                                      //       which is already taken into account as a weight 
      fHistoPTPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      if ( i < 2 ) fHistoSecPiPlus[i]->Scale(1./((double)stat),"width");
      fHistoSecPiPlus[i]->Write();
   }
   for ( size_t i=0; i<fHistoPTPiPlus.size(); ++i )
   {
      fHistoPTPiPlus[i]->Scale( 1./((double)stat) ); // NO scaling with "width" - see earlier comment for fHistoPTPiMinus
      fHistoPTPiPlus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecNeutron.size(); ++i )
   {
     if ( i < 2 ) fHistoSecNeutron[i]->Scale(1./((double)stat),"width");
     fHistoSecNeutron[i]->Write();
   }

   return;

}

double TestNA49Histo::CalculateBinWeight( const G4LorentzVector& labp, const G4String& pname, 
                                          double Ekin, double mass, double pT, int xFbin, double sqrt_s )
{

   double wei = 1.;
   int pTbin = -1;
   
// FIZME !!! Need to redo this method and make it more generic !!!

   if ( pname == "pi+" || pname == "pi-" )
   {
      for (  int ib=0; ib<fNPiBinsPT; ++ib )
      {
         if ( pT >= fPiBinsPT[ib] && pT < fPiBinsPT[ib+1] )
         {
            pTbin = ib;
	    break;
         }
      }
   }
   else if ( pname == "proton" )
   {
      for (  int ib=0; ib<fNPBinsPT; ++ib )
      {
         if ( pT >= fPBinsPT[ib] && pT < fPBinsPT[ib+1] )
         {
            pTbin = ib;
	    break;
         }
      }      
   }
   else if ( pname == "anti_proton" )
   {
      for ( int ib=0; ib<fNPbarBinsPT[xFbin]; ++ib )
      {
         if ( pT >= fPbarBinsPT[xFbin][ib] && pT < fPbarBinsPT[xFbin][ib+1] )
	 {
	    pTbin = ib;
	    break;
	 }
      }
   }
   
   if ( pTbin == -1 ) return wei;

// NOTE: This fragment of code draws inspiration in a similar application
//       originally implemented by Mike Kordosky (W&M/MINERvA/NuMI).
//       Credits go to Mike !!!
//
   double dpT2  = 0.; // fPiBinsPT[pTbin+1]*fPiBinsPT[pTbin+1] - fPiBinsPT[pTbin]*fPiBinsPT[pTbin];
   double pLmin = 0.; // fPiBinsXF[xFbin]*sqrt_s/2.;
   double pLmax = 0.; // fPiBinsXF[xFbin+1]*sqrt_s/2.;

   if ( pname == "pi+" || pname =="pi-" )
   {
      dpT2 = fPiBinsPT[pTbin+1]*fPiBinsPT[pTbin+1] - fPiBinsPT[pTbin]*fPiBinsPT[pTbin];
      pLmin = fPiBinsXF[xFbin]*sqrt_s/2.;
      pLmax = fPiBinsXF[xFbin+1]*sqrt_s/2.;
   }
   else if ( pname == "proton" )
   {
      dpT2 = fPBinsPT[pTbin+1]*fPBinsPT[pTbin+1] - fPBinsPT[pTbin]*fPBinsPT[pTbin];
      pLmin = fPBinsXF[xFbin]*sqrt_s/2.;
      pLmax = fPBinsXF[xFbin+1]*sqrt_s/2.;
   }
   else if ( pname == "anti_proton" )
   {
      dpT2 = fPbarBinsPT[xFbin][pTbin+1]*fPbarBinsPT[xFbin][pTbin+1] - fPbarBinsPT[xFbin][pTbin]*fPbarBinsPT[xFbin][pTbin];
      pLmin = fPbarBinsXF[xFbin]*sqrt_s/2.;
      pLmax = fPbarBinsXF[xFbin+1]*sqrt_s/2.;
   }

   double EPCM1 = sqrt( (pLmin*pLmin+pT*pT) + mass*mass );
   double EPCM2 = sqrt( (pLmax*pLmax+pT*pT) + mass*mass );
   double pZmin = labp.boostVector().gamma()*(labp.beta()*EPCM1 + pLmin );
   double pZmax = labp.boostVector().gamma()*(labp.beta()*EPCM2 + pLmax );
   double dP3   = CLHEP::pi * dpT2 *(pZmax-pZmin);

   wei = (Ekin+mass)/dP3;

   return wei;
   
}
