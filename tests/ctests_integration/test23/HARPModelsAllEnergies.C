
std::string TEST_NAME="test23";

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/G4PHYSLIST_IntermediateEnergy.h"
#include "../test23/shared-root-macros/PlotHARPModels.C"

// void HARPModelsAllEnergies( char beam[8]="piminus", char target[3]="Be" )
void HARPModelsAllEnergies( std::string beampart="piminus", std::string tg="Be" )
{

   // std::string beampart(beam);
   // std::string tg(target);
   
   PlotHARPAnalysis( beampart, tg, "3.0", "piplus",  "FW");
   PlotHARPAnalysis( beampart, tg, "3.0", "piplus",  "LA");
   PlotHARPAnalysis( beampart, tg, "3.0", "piminus", "FW");
   PlotHARPAnalysis( beampart, tg, "3.0", "piminus", "LA");

   PlotHARPAnalysis( beampart, tg, "5.0", "piplus",  "FW");
   PlotHARPAnalysis( beampart, tg, "5.0", "piplus",  "LA");
   PlotHARPAnalysis( beampart, tg, "5.0", "piminus", "FW");
   PlotHARPAnalysis( beampart, tg, "5.0", "piminus", "LA");

   PlotHARPAnalysis( beampart, tg, "8.0", "piplus",  "FW");
   PlotHARPAnalysis( beampart, tg, "8.0", "piplus",  "LA");
   PlotHARPAnalysis( beampart, tg, "8.0", "piminus", "FW");
   PlotHARPAnalysis( beampart, tg, "8.0", "piminus", "LA");
   
   PlotHARPAnalysis( beampart, tg, "12.0", "piplus",  "FW");
   PlotHARPAnalysis( beampart, tg, "12.0", "piplus",  "LA");
   PlotHARPAnalysis( beampart, tg, "12.0", "piminus", "FW");
   PlotHARPAnalysis( beampart, tg, "12.0", "piminus", "LA");

/*
   if ( beampart == "proton" && tg == "Be" )
   {
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","LA");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","LA");
   }
*/
   return;   

}
