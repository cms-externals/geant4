
std::string TEST_NAME="test23";

#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test23/G4PHYSLIST_IntermediateEnergy.h"

#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"

#include "../test23/shared-root-macros/PlotHARPRegre.C"

// void HARPRegre( char beam[8]="piminus", char target[3]="Be", int energy=12000 )
void HARPRegre( std::string beampart="piminus", std::string tg="Be", int energy=12000 )
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
   
   if ( ene == "8.9" && !( beampart=="proton" && tg=="Be" ) )
   {
      std::cout << "Only proton+Be is avalable at 8.9GeV/c" << std::endl;
      return;
   }
   
   if ( beampart == "proton" && tg == "Ta" && ene != "8.0" )
   {
      std::cout << "proton+Ta is available only at 8.0GeV/c" << std::endl;
      return;
   } 

   for ( int m=0; m<NModels_IE; ++m )
   {

      PlotHARPAnalysisRegre( beampart, tg, ene, "piplus", "FW", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, ene, "piplus", "LA", ModelName_IE[m] ); 
      PlotHARPAnalysisRegre( beampart, tg, ene, "piminus", "FW", ModelName_IE[m] );
      PlotHARPAnalysisRegre( beampart, tg, ene, "piminus", "LA", ModelName_IE[m] ); 

/*
      if ( beampart == "proton" && tg == "Be" && ene == "8.9" )
      {
         PlotHARPAnalysisRegre( beampart, tg, ene, "piplus", "FW", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, ene, "piplus", "LA", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, ene, "piminus", "FW", ModelName_IE[m] ); 
         PlotHARPAnalysisRegre( beampart, tg, ene, "piminus", "LA", ModelName_IE[m] );          
      }
*/
      if ( beampart == "proton" && tg == "Ta" && ene == "8.0" )
      {	 
	 PlotHARPRegreForMu2e( "piplus", ModelName_IE[m] );
	 PlotHARPRegreForMu2e( "piminus", ModelName_IE[m] );
      }

   }
   
   return;
   
}

