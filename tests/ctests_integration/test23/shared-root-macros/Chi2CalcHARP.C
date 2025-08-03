#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "HARP.h"
#include "Chi2Calc.C"

#include "REGRESSION_TEST.h"

#ifndef G4VAL_CHI2CALCHARP_C
#define G4VAL_CHI2CALCHARP_C

double Chi2MomSpectrumHARP( std::string beam, std::string target, std::string energy,
                            std::string secondary, 
                            std::string region, // must be "FW" or "LA"
			    std::string model,
			    int& NDF,
			    std::string version = ".",
			    bool useStatEr=true, bool useSysEr=true )
{

// NOTE (JVY): At present, I'm NOT using (under/over)flow in the MC histo because
//             the X-range (P,GeV/c) is standatd for all of them, but the data 
//             may start at hight value of P.


   double Chi2 = 0.;
      
   std::string location = "";
   if ( version == CurrentVersion || version == "." )
   {
         location = ".";
   }
   else
   {
         location = regre_test_dir + "/" + TEST_NAME + "/" + version;
   }

   std::string histofile = location + "/harp-histo/" + beam + target + energy + "GeV" + model + ".root";
                     
   TFile* f = TFile::Open( histofile.c_str() );
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      
      char buf[5];
      sprintf( buf, "%i", i );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );
            
      TH1F* h = (TH1F*)f->Get( histoname.c_str() );
      
      int NX = h->GetNbinsX();      
      
      int first = -1;
      int last  = -1;
      // find the first bin that matches data & MC
      for ( int k=0; k<=NX; ++k )
      {
         double xx1 = h->GetBinLowEdge(k);
	 double xx2 = h->GetBinWidth(k);
	 for ( int kk=0; kk<NPointsHARP[i]; ++kk )
	 {
	    double xcenter = 0.5 * (XMinHARP[i][kk]+XMaxHARP[i][kk]);
	    if ( xx1 < xcenter && xx1+xx2 > xcenter )
	    {
	       if ( first == - 1 ) first=k;
	       break;
	    }
	 }
	 if ( first != - 1 ) break;
      }
      
      // from here on, we assume that the binning is the same in MC and Data
      last = first + NPointsHARP[i] ;
      if ( last > h->GetXaxis()->GetLast() )
      {
         last = h->GetXaxis()->GetLast();
      }
            
      int nmtx = last - first;
      
      for ( int k=0; k<nmtx; k++ ) 
      {
            if ( YHARP[i][k] <=0. || EYHARP[i][k] <= 0. ) continue; // in some datasets -1 marks missing/non-existing data
	    
	    double xx1 = h->GetBinLowEdge(k+first);
	    double xx2 = h->GetBinWidth(k+first);
	    double xcenter = 0.5 * (XMinHARP[i][k]+XMaxHARP[i][k]);
            double yy1  = h->GetBinContent(k+first);
	    double eyy1 = h->GetBinError(k+first);
	    if ( ( eyy1*eyy1 + EYHARP[i][k]*EYHARP[i][k] ) > 0 )
	    {
	       Chi2 += ( yy1 - YHARP[i][k] )*( yy1 - YHARP[i][k] ) / ( eyy1*eyy1 + EYHARP[i][k]*EYHARP[i][k] );
	       ++NDF;
	    }
      }
               
   }
   
   // should I do subtract 1 degree of freedom at the end ? 
   // NDF--;     
   
   f->Close();
   
   return Chi2;

}

#endif
