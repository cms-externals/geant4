#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "NA49.h"

// #include "TMatrixD.h"
// #include "TDecompSVD.h"

#include "Chi2Calc.C"
#include "ReadNA49Data.C"

#include "REGRESSION_TEST.h"


#ifndef G4VAL_CHI2CALCNA49_C
#define G4VAL_CHI2CALCNA49_C

double Chi2DDiffXSecNA49( std::string beam, std::string target, 
                          std::string secondary, 
			  std::string model, 
                          int& NDF,
                          std::string version = ".",
			  bool useStatEr=true, bool useSysEr=true )
{
   
   double chi2 = 0;

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
      std::string histofile = location + "/na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 
      
//   std::string histofile = "";
//   histofile += version;
//   histofile += ("/na49-histo/" + beam + target + "158.0GeV" + model + ".root"); 
   TFile* f = TFile::Open( histofile.c_str() );

   for ( int icount=0; icount<NSubSetsNA49; ++icount )
   {
      std::ostringstream cnt;
      if ( secondary == "piplus" || secondary == "piminus" )
      {
         cnt << (icount+7);
      }
      else
         cnt << icount;

      std::string hname = "pT";
      if ( secondary == "piplus" )
      {
	    hname +="pip";
      }
      else if ( secondary == "piminus" )
      {
	    hname += "pim";
      }
      else if ( secondary == "proton" )
      {
            hname += "pro";
      }
      else if ( secondary == "antiproton" )
      {
         hname += "pbar";
      }

      hname += cnt.str();
      
      TH1F* h = (TH1F*)f->Get( hname.c_str() );
      h->Scale(226.); // or should it be xsec=251. as in Geant4 ???   

      TGraphErrors* gdata = get1DDiffXSecAsGraph( icount );
      chi2 += Chi2( gdata, h, NDF );  
         
   } 
   
   //  subtract 1 degree of freedom ?
   // NDF--;   
   
   f->Close();
   
   return chi2;
   
}

double Chi2IntegratedSpectrumNA49( std::string beam, std::string target, 
                                   std::string secondary, 
                                   std::string histo,
			           std::string model,  
                                   int& NDF,
                                   std::string version = ".",
			           bool useStatEr=true, bool useSysEr=true )
{

   double chi2 = 0.;

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
      
      std::string histofile = location + "/na49-histo/" + beam + target + "158.0GeV" + model + ".root";
      
//   std::string histofile = "";
//   histofile += version;
//   histofile += ("/na49-histo/" + beam + target + "158.0GeV" + model + ".root"); 
   TFile* f = TFile::Open( histofile.c_str() );

   std::string histoname = secondary + "_" + histo;
   TH1F* h = (TH1F*)f->Get( histoname.c_str() );

   TGraphErrors* gdata = getIntegratedSpectrumAsGraph( histo );
   chi2 += Chi2( gdata, h, NDF ); 

   f->Close();
   
   return chi2;

}

#endif
