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
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TStyle.h"
#include "TGraph.h"

// #include "../test23/shared-root-macros/Chi2Calc.C"

#include "../test47/PlotITEPAnalysis.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

// void ITEPModels( char beam[8]="piminus" )
void ITEPModels( std::string beampart="piminus" )
{

   // std::string beampart(beam);
   
   std::string en = "5.00";
   if ( beampart == "proton" )
   {
      en = "7.50";
   }
   
   plotModelsMC2DataSummary( beampart, "C", "1.40" );
   plotModelsMC2DataSummary( beampart, "C", en );
   plotModelsMC2DataSummary( beampart, "U", "1.40" );
   plotModelsMC2DataSummary( beampart, "U", en );

   return ;

}
