
#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test23/G4PHYSLIST_IntermediateEnergy.h"
#include "../test23/shared-root-macros/PlotHARPRegre.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"

std::string TEST_NAME="test23";

void HARPRegreAllEnergies( char beam[8]="piminus", char target[3]="Be" )
{
      
   std::string beampart(beam);
   std::string tg(target);

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
/*
      if ( beampart == "proton" && tg == "Be" )
      {
         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piplus", "FW", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piplus", "LA", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piminus", "FW", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, "8.9", "piminus", "LA", ModelName_IE[m] );          
      }
*/
   }
   
   return;
   
}

