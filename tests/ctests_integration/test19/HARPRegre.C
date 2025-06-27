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
#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_IntermediateEnergy.h"

#include "../test23/shared-root-macros/HARP.h"
#include "../test23/shared-root-macros/ReadHARPData.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"
#include "../test23/shared-root-macros/DrawHARPSpectra.C"

#include "../test23/shared-root-macros/PlotHARPRegre.C"


//void HARPRegre( char beam[8]="piminus", char target[3]="Be", char energy[5]="3.0" )
void HARPRegre( std::string beampart="piminus", std::string tg="Be", std::string en="3.0" )
{

   // std::string beampart(beam);
   // std::string tg(target);
   // std::string en(energy);
   
   for ( int m=0; m<NModels_IE; ++m )
   {

      PlotHARPAnalysisRegre( beampart, tg, en, "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, en, "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, en, "piminus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, en, "piminus", "LA", ModelName_IE[m] ); 

      if ( beampart == "proton" && tg == "Ta" && en == "8.0" )
      {	 
	 PlotHARPRegreForMu2e( "piplus", ModelName_IE[m] );
	 PlotHARPRegreForMu2e( "piminus", ModelName_IE[m] );
      }

   }
   
   return;

}

/*
void HARPRegre( char beam[8]="piminus", char target[3]="Be" )
{

   std::string beampart(beam);
   std::string tg(target);
   std::string en(energy);


   for ( int m=0; m<NModels_IE; ++m )
   {

      PlotHARPAnalysisRegre( beampart, tg, "3.0", "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "3.0", "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "3.0", "piminus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "3.0", "piminus", "LA", ModelName_IE[m] ); 

      PlotHARPAnalysisRegre( beampart, tg, "5.0", "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "5.0", "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "5.0", "piminus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "5.0", "piminus", "LA", ModelName_IE[m] ); 

      PlotHARPAnalysisRegre( beampart, tg, "8.0", "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "8.0", "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "8.0", "piminus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "8.0", "piminus", "LA", ModelName_IE[m] ); 

      PlotHARPAnalysisRegre( beampart, tg, "12.0", "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "12.0", "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "12.0", "piminus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, "12.0", "piminus", "LA", ModelName_IE[m] ); 
            
//      if ( beampart == "proton" && tg == "Be" )
//      {
//         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piplus", "FW", ModelName_IE[m] ); 
//         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piplus", "LA", ModelName_IE[m] ); 
//         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piminus", "FW", ModelName_IE[m] ); 
//         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piminus", "LA", ModelName_IE[m] );          
//      }

      if ( tg == "Ta" )
      {
         PlotHARPForMu2eRegre( "piplus" );
	 PlotHARPForMu2eRegre( "piminus" );
      }

   }
   
   return;
   
}
*/
