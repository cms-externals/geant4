
std::string TEST_NAME="test23";

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/G4PHYSLIST_IntermediateEnergy.h"
#include "../test23/shared-root-macros/PlotHARPModels.C"

// void HARPModels( char beam[8]="piminus", char target[3]="Be", int energy=12000 )
void HARPModels( std::string beampart="piminus", std::string tg="Be", int energy=12000 )
{

   // std::string beampart(beam);
   // std::string tg(target);
   // std::cout << " energy = " << energy << std::endl;
      
   std::ostringstream os;
   os.precision(4);   
   double denergy = energy;
   os << (denergy/1000);
   std::string ene = os.str();
   if ( (energy%1000) == 0 )
   {
      ene += ".0";
   }
   
   if ( beampart == "proton" && tg == "Ta" && ene != "8.0" )
   {
      std::cout << "proton+Ta is available only at 8.0GeV/c" << std::endl;
      return;
   } 

   if ( ene == "8.9" && !( beampart=="proton" && tg=="Be" ) )
   {
      std::cout << "Only proton+Be is avalable at 8.9GeV/c" << std::endl;
      return;
   }
   
   PlotHARPAnalysis( beampart, tg, ene, "piplus",  "FW");
   PlotHARPAnalysis( beampart, tg, ene, "piplus",  "LA");
   PlotHARPAnalysis( beampart, tg, ene, "piminus", "FW");
   PlotHARPAnalysis( beampart, tg, ene, "piminus", "LA");
   
   if ( beampart == "proton" && tg == "Be" )
   {
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piplus","LA");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","FW");
      PlotHARPAnalysis(beampart,tg,"8.9","piminus","LA");
   }

   if ( beampart == "proton" && tg == "Ta" && ene == "8.0" )
   {	 
	 PlotHARPForMu2e( "piplus" );
	 PlotHARPForMu2e( "piminus" );
   }

   return;   

}
