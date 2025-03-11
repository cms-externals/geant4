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

std::string TEST_NAME="test19";
#include "../test19/G4MODELS_HighEnergy.h"

#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/ReadNA61Data.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"
#include "../test23/shared-root-macros/DrawNA61Spectra.C"

#include "../test23/shared-root-macros/PlotNA61Models.C"

void NA61Models()
{

/*
   plotSecondarySumCombined2pages("piplus");
   plotSecondarySumCombined2pages("piminus");
   plotSecondarySumCombined2pages("proton");
*/   
   plotModels2Data2015( "piplus" );
   plotModels2Data2015( "piminus" );
   plotModels2Data2015( "kplus" );
   plotModels2Data2015( "kminus" );
   plotModels2Data2015( "k0s" );
   plotModels2Data2015( "lambda" );
   plotModels2Data2015( "proton" );
   
   
   plotKPlus2PiPlusRatio( "proton", "C" );
   
   return;
   
}

