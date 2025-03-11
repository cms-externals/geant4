
/* remove DbReader, revert to using local ASCII files
R__LOAD_LIBRARY(libDbReader.so);
*/
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
#include "TGraph.h"

/* remove DbReader, revert to using local ASCII files
#include "../DbReader/DbReader.h"
*/
std::string TEST_NAME="test23";

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/G4PHYSLIST_HighEnergy.h"

#include "../test23/shared-root-macros/NA49.h"
#include "../test23/shared-root-macros/ReadNA49Data.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA49.C"
#include "../test23/shared-root-macros/DrawNA49Spectra.C"

#include "../test23/shared-root-macros/PlotNA49Regre.C"


void NA49Regre()
{

   for ( int m=0; m<NModels_HE; ++m )
   {
      plot_dNdxF_pT_Regre( "proton", "C", "piplus", ModelName_HE[m] ); 
      plot_dNdxF_pT_Regre( "proton", "C", "piminus", ModelName_HE[m] ); 
      plot_dNdxF_pT_Regre( "proton", "C", "proton", ModelName_HE[m] ); 
      plot_dNdxF_pT_Regre( "proton", "C", "antiproton", ModelName_HE[m] ); 
      plot_dNdxF_pT_Regre( "proton", "C", "neutron", ModelName_HE[m] );       
   }
   
   return;
   
}

