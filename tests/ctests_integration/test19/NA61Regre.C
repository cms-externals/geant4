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

std::string TEST_NAME="test19";

#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_HighEnergy.h"

#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/ReadNA61Data.C"
#include "../test23/shared-root-macros/DrawNA61Spectra.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"

#include "../test23/shared-root-macros/PlotNA61Regre.C"


void NA61Regre()
{

  for ( int m=0; m<NModels_HE; ++m )
  {
/*
     plotSecondarySumCombinedRegre( "piplus", ModelName_HE[m] ); 
     plotSecondarySumCombinedRegre( "piminus", ModelName_HE[m] ); 
     plotSecondarySumCombinedRegre( "proton", ModelName_HE[m] ); 
*/
     plotRegreData2015( "piplus", ModelName_HE[m] );
     plotRegreData2015( "piminus", ModelName_HE[m] );
     plotRegreData2015( "kplus", ModelName_HE[m] );
     plotRegreData2015( "kminus", ModelName_HE[m] );
     plotRegreData2015( "k0s", ModelName_HE[m] );
     plotRegreData2015( "lambda", ModelName_HE[m] );
     plotRegreData2015( "proton", ModelName_HE[m] );

     plotKPlus2PiPlusRatioRegre( "proton", "C", ModelName_HE[m] );
  }

}

