#include <iostream>

#include <math.h>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"

// #include "TMatrixD.h"
// #include "TDecompSVD.h"
      

#ifndef G4VAL_CHI2_CALC_C
#define G4VAL_CHI2_CALC_C

double Chi2( TH1F* hdata, TH1F* hsim, int& NDF )
{

   double Chi2 = 0.;
   
   int NXData = hdata->GetNbinsX();
   int NXSim = hsim->GetNbinsX();

   for ( int k=1; k<=NXSim; ++k ) 
   { 
         double xx1 = hsim->GetBinLowEdge(k);
	 double xx2 = hsim->GetBinWidth(k);
	 for (int kk=1; kk<=NXData; ++kk )
	 {
	   
	   double xd1 = hdata->GetBinLowEdge(kk);
	   double xd2 = hdata->GetBinWidth(kk);
	   
	   if ( xx1 < (xd1+xd2/2.) && xx1+xx2 > (xd1+xd2/2.) ) // in principle, this is not safe because (in general) 
                                                               // the bin width maybe different..
	    {
	       double yy1  = hsim->GetBinContent(k);
	       double eyy1 = hsim->GetBinError(k);
	       double yd1  = hdata->GetBinContent(kk); 
	       double eyd1 = hdata->GetBinError(kk);
	       if ( ( eyy1*eyy1 + eyd1*eyd1 ) > 0 )
	       {
	          Chi2 += ( yy1 - yd1 )*( yy1 - yd1 ) / ( eyy1*eyy1 + eyd1*eyd1 );
	          ++NDF;
	       } 
	       break;
	    }
	    
	 }
   }
               
   return Chi2;
   
}

double Chi2( TGraphErrors* gdata, TH1F* hsim, int& NDF )
{

   double Chi2=0.;

   int     NXData = gdata->GetN();

   // std::cout << "NXData = " << NXData << std::endl; 

   double* XData  = gdata->GetX();
   double* YData  = gdata->GetY();
   double* EYData = gdata->GetEY();
   
            
   int NXSim  = hsim->GetNbinsX();

   // std::cout << " NXSim = " << NXSim << std::endl;

      for ( int k=1; k<=NXSim; ++k ) 
      { 
         double xx1 = hsim->GetBinLowEdge(k);
	 double xx2 = hsim->GetBinWidth(k);
	 for (int kk=0; kk<NXData; ++kk )
	 {
	   if ( xx1 < XData[kk] && xx1+xx2 > XData[kk] ) // in principle, this is not safe because (in general) 
                                                         // the bin width maybe different..
	    {
	       double yy1  = hsim->GetBinContent(k);
	       double eyy1 = hsim->GetBinError(k);
	       if ( ( eyy1*eyy1 + EYData[kk]*EYData[kk] ) > 0 )
	       {
	          Chi2 += ( yy1 - YData[kk] )*( yy1 - YData[kk] ) / ( eyy1*eyy1 + EYData[kk]*EYData[kk] );
	          ++NDF;
	       } 
	       break;
	    }
	 }
      }      

      return Chi2;

}

#endif
