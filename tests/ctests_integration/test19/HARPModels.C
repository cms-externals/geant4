std::string TEST_NAME="test19";
#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_IntermediateEnergy.h"
#include "../test23/shared-root-macros/PlotHARPModels.C"

//void HARPModels( char beam[8]="piminus", char target[3]="Be", double energy=3 )
void HARPModels( std::string beampart="piminus", std::string tg="Be", double energy=3 )
{

   //std::string beampart(beam);
   //std::string tg(target);
   
   std::ostringstream os;
   os << energy;
   std::string en = os.str(); 
   if ( en.find(".") == std::string::npos )
   {
      en += ".0";
   }
   
   PlotHARPAnalysis(beampart,tg,en,"piplus","FW");
   PlotHARPAnalysis(beampart,tg,en,"piplus","LA");
   PlotHARPAnalysis(beampart,tg,en,"piminus","FW");
   PlotHARPAnalysis(beampart,tg,en,"piminus","LA");

   if ( beampart == "proton" && ( tg == "Ta" || tg == "C" ) && en == "8.0" )
   {
      PlotHARPForMu2e( "piplus", tg );
      PlotHARPForMu2e( "piminus" , tg );
   }   

   return;   

}

/* old implementation before we had to split jobs by target and energies

void HARPModels( char beam[8]="piminus", char target[3]="Be" )
{

   std::string beampart(beam);
   std::string tg(target);
   
   PlotHARPAnalysis(beampart,tg,"3.0","piplus","FW");
   PlotHARPAnalysis(beampart,tg,"3.0","piplus","LA");
   PlotHARPAnalysis(beampart,tg,"3.0","piminus","FW");
   PlotHARPAnalysis(beampart,tg,"3.0","piminus","LA");

   PlotHARPAnalysis(beampart,tg,"5.0","piplus","FW");
   PlotHARPAnalysis(beampart,tg,"5.0","piplus","LA");
   PlotHARPAnalysis(beampart,tg,"5.0","piminus","FW");
   PlotHARPAnalysis(beampart,tg,"5.0","piminus","LA");

   PlotHARPAnalysis(beampart,tg,"8.0","piplus","FW");
   PlotHARPAnalysis(beampart,tg,"8.0","piplus","LA");
   PlotHARPAnalysis(beampart,tg,"8.0","piminus","FW");
   PlotHARPAnalysis(beampart,tg,"8.0","piminus","LA");

   PlotHARPAnalysis(beampart,tg,"12.0","piplus","FW");
   PlotHARPAnalysis(beampart,tg,"12.0","piplus","LA");
   PlotHARPAnalysis(beampart,tg,"12.0","piminus","FW");
   PlotHARPAnalysis(beampart,tg,"12.0","piminus","LA");
   
   if ( beampart == "proton"  && tg == "Be" )
   {
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","LA");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","LA");
   }

   return;   

}
*/
