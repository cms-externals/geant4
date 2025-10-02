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

#include "../test23/shared-root-macros/NA49.h"
#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/HARP.h"

#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/ReadNA49Data.C"
#include "../test23/shared-root-macros/Chi2CalcNA49.C"
#include "../test23/shared-root-macros/ReadNA61Data.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"
#include "../test23/shared-root-macros/ReadHARPData.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"

const int NModels = 1;
std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };  
// std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[2]  = { "qgsp_bert-with-decays", "ftfp_bert-with-decays" };  
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

// static int isNA49UtilLoaded = 0;
// static int isNA61UtilLoaded = 0;
// static int isHARPUtilLoaded = 0;
// static int isChi2UtilLoaded = 0;

void calcChi2DDiffXSecNA49( std::string beam, std::string target, std::string secondary )
{

   readDDiffSpectra( beam, target, secondary );
      
   std::cout << " chi2 is calculated over double diff. pt spectra for " << secondary << ", as measured by NA49: " << std::endl;
   for ( int m=0; m<NModels; ++m )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, secondary, ModelName[m], NDF );
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;
   }   

   return;

}

void calcChi2IntegratedSpectrumNA49( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );  

   std::cout << " chi2 is calculated over " << histo << " spectrum for " << secondary << ", as measured by NA49:" << std::endl; 
   for ( int m=0; m<NModels; ++m )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, histo, ModelName[m], NDF );
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;
   } 

   return;

}

void calcChi2MomSpectrumNA61( std::string beam, std::string target, std::string secondary )
{
  
   std::cout << " chi2 is calculated over momentum spectra in different angular bin, as measured by NA61: " << std::endl;

   for ( int m=0; m<NModels; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 = Chi2MomSpectrumNA61Integral( beam, target, secondary, ModelName[m], NDF );

      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;

   }   

   return;

} 


void calcChi2MomSpectrumHARP( std::string beam, std::string target, std::string energy, std::string secondary )
{

  /*
   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      isHARPUtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
  */

   for ( int m=0; m<NModels; ++m )
   {
      
      double chi2 = 0.;
      double chi2_fw = 0.;
      double chi2_la = 0.; 
      int NDF = 0;
      int NDF_fw = 0;
      int NDF_la = 0;

      ReadHARPData( beam, target, energy, secondary, "FW" );   
      chi2_fw += Chi2MomSpectrumHARP( beam, target, energy, secondary, "FW", ModelName[m], NDF_fw );
      
      std::cout << "FW only:  chi2/NDF = " << chi2_fw << "/" << NDF_fw << " = " << chi2_fw/NDF_fw << " for model " << ModelName[m] << std::endl;
      
      if ( secondary == "proton" ) continue;
      
      ReadHARPData( beam, target, energy, secondary, "LA" );   
      chi2_la += Chi2MomSpectrumHARP( beam, target, energy, secondary, "LA", ModelName[m], NDF_la );

      std::cout << "LA only:  chi2/NDF = " << chi2_la << "/" << NDF_la << " = " << chi2_la/NDF_la << " for model " << ModelName[m] << std::endl;
      
      NDF = NDF_fw + NDF_la;
      chi2 = chi2_fw + chi2_la;

      std::cout << "   Total:  chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;

   }
   
   return;

}

