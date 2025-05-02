#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
// #include <list>

#include <math.h>
// #include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
// #include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "NA61.h"
#include "ReadNA61Data.C"
#include "Chi2Calc.C"

#include "REGRESSION_TEST.h"

// #include "TMatrixD.h"
// #include "TDecompSVD.h"

#ifndef G4VAL_CHI2CALCNA61_C
#define G4VAL_CHI2CALCNA61_C

double Chi2MomSpectrumNA61ThetaBin( std::string beam, std::string target, std::string secondary, 
                                    std::string sec_angle_min, std::string sec_angle_max,
			            std::string model,
				    int& NDF,
				    std::string version=".",
				    bool pub2015=false,
				    bool useStatEr=true, bool useSysEr=true )
{

   double chi2 = 0.;

//   std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 
      std::string location = "";
      if ( version == CurrentVersion || version == "." )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + version;
         location = regre_test_dir + "/" + TEST_NAME + "/" + version;
      }
      std::string histofile = location + "/na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 
      
   TFile* f = TFile::Open( histofile.c_str() );
      
   // std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
   std::string histoname = secondary;
   if ( pub2015 )
   {
         histoname += "_";
   }
   else
   {
         histoname += "Mult_"; 
   }
   histoname += ( sec_angle_min + "_" + sec_angle_max );
   
   TH1F* h = (TH1F*)f->Get( histoname.c_str() );
   int NX = h->GetNbinsX();

   /*
   for ( int k=0; k<=NX; k++ ) 
   {
      double xx1 = h->GetBinLowEdge(k);
      double xx2 = h->GetBinWidth(k);
      for (int kk=0; kk<NPointsNA61; kk++ )
      {
	 if ( xx1 < SecMom[kk] && xx1+xx2 > SecMom[kk] )
	 {
	    double yy1  = h->GetBinContent(k);
	    double eyy1 = h->GetBinError(k);
	    // Chi2 += ( yy1 - SecSigma[kk] )*( yy1 - SecSigma[kk] ) / SecSigma[kk];
	    if ( ( eyy1*eyy1 + SecESigma[kk]*SecESigma[kk] ) > 0 )
	    {
	       Chi2 += ( yy1 - SecSigma[kk] )*( yy1 - SecSigma[kk] ) / ( eyy1*eyy1 + SecESigma[kk]*SecESigma[kk] );
	       ++NDF;
	    }
	    break; 
	 }
      }
   }
   */
   TGraphErrors* gdata = getNA61AsGraph();
   chi2 += Chi2( gdata, h, NDF );
   
   // basically, I seem to be using shape errors only, so... subtract 1 degree of freedom ?
   // NDF--;   

   f->Close();
   
   return chi2; 

}

double Chi2NA61IntegralData2015( std::string beam, std::string target, std::string secondary,
                                 std::string model,
                                 int& NDF,
				 std::string version=".",
                                 bool useStatEr=true, bool useSysEr=true )
{

   double chi2 = 0.;
   
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      readMomSpectrum( "proton", "C", secondary, "0", "10", "../test23/na61-exp-data/Pub-2015/" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "0", "10",
				           model,
					   NDF, version, true );
      readMomSpectrum( "proton", "C", secondary, "10", "20", "../test23/na61-exp-data/Pub-2015/" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "10", "20",
				           model,
					   NDF, version, true );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      readMomSpectrum( "proton", "C", secondary, "0", "20", "../test23/na61-exp-data/Pub-2015/" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "0", "20",
				           model,
					   NDF, version, true );
   
   }
   else if ( secondary == "k0s" )
   {
      readMomSpectrum( "proton", "C", "K0s", "0", "40", "../test23/na61-exp-data/Pub-2015/" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "0", "40",
				           model,
					   NDF, version, true );
   
   }
   else if ( secondary == "lambda" )
   {
      readMomSpectrum( "proton", "C", "lambda", "0", "40", "../test23/na61-exp-data/Pub-2015/" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "0", "40",
				           model,
					   NDF, version, true );
   
   }
   
   readMomSpectrum( "proton", "C", secondary, "40", "60", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "40", "60",
				           model,
				           NDF, version, true );
   readMomSpectrum( "proton", "C", secondary, "60", "100", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "60", "100",
				           model,
					   NDF, version, true );
   readMomSpectrum( "proton", "C", secondary, "100", "140", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "100", "140",
				           model,
					   NDF, version, true );
   readMomSpectrum( "proton", "C", secondary, "140", "180", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "140", "180",
				           model,
					   NDF, version, true );
   readMomSpectrum( "proton", "C", secondary, "180", "240", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "180", "240",
				           model,
					   NDF, version, true );

   readMomSpectrum( "proton", "C", secondary, "240", "300", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "240", "300",
				           model,
					   NDF, version, true );

   if ( secondary == "lambda" || secondary == "k0s" || 
        secondary == "kplus"  || secondary == "kminus" ) return chi2;
   
   readMomSpectrum( "proton", "C", secondary, "300", "360", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "300", "360",
				           model,
					   NDF, version, true );
   
   if ( secondary == "proton" ) return chi2;

   readMomSpectrum( "proton", "C", secondary, "360", "420", "../test23/na61-exp-data/Pub-2015/" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "360", "420",
				           model,
					   NDF, version, true );
   
   return chi2;

}

double Chi2MomSpectrumNA61Integral( std::string beam, std::string target, std::string secondary,
                                std::string model,
                                int& NDF,
				std::string version=".",
                                bool useStatEr=true, bool useSysEr=true )
{

   double chi2 = 0.;

   readMomSpectrum( "proton", "C", secondary, "0", "20" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "0", "20",
				           model,
					   NDF, version );
   readMomSpectrum( "proton", "C", secondary, "20", "40" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "20", "40",
				           model,
				           NDF, version );
   readMomSpectrum( "proton", "C", secondary, "40", "60" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "40", "60",
				           model,
				           NDF, version );
   readMomSpectrum( "proton", "C", secondary, "60", "100" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "60", "100",
				           model,
					   NDF, version );
   readMomSpectrum( "proton", "C", secondary, "100", "140" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "100", "140",
				           model,
					   NDF, version );
   readMomSpectrum( "proton", "C", secondary, "140", "180" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "140", "180",
				           model,
					   NDF, version );
   readMomSpectrum( "proton", "C", secondary, "180", "240" );
   chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "180", "240",
				           model,
					   NDF, version );
   if ( secondary != "proton" )
   {
         readMomSpectrum( "proton", "C", secondary, "240", "300" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
                                           "240", "300",
				           model,
					   NDF, version );
	 // exclude this theta-bin for now as we only plot 8 theta-bins for pions 
	 //         readMomSpectrum( "proton", "C", secondary, "300", "360" );
	 //         Chi2 += Chi2MomSpectrumNA61ThetaBin( "proton", "C", secondary, 
	 //                                           "300", "360",
					     //				           model,
					     //					   NDF );
   }
   
   return chi2;

}

#endif
