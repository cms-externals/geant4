
#include "../../test23/shared-upload-scripts/SQLCompose.C"

// NOTE(JVY): "using namespace std" is declared in SQLCompose.C

void test48_gen_histo_sql( string g4version, string status )
{

   // pi- capture part
   //
   int ntgt4pimin = 7;
   string target[7] = { "C", "N", "O", "Al", "Cu", "Ta", "Pb" };

   for ( int it=0; it<ntgt4pimin; ++it )
   {
      SQLCompose( 48, "piminus", target[it], "0.GeV", "", "neutron",
	          "BertiniPreCo", 
		  g4version,
		  "NvsT",
		  "",
		  "neutron kinetic energy [MeV]", "neutron yield",
		  "",
		  "neutron yield versus kinetic energy",
		  status, 0 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
   }
   
   // mu- capture part
   //
   int ntgt4muminSingerExp = 8; 
   string targetSingerExp[8] = { "Al", "Si", "Ca", "Fe", "Ag", "I", "Au", "Pb" };
   int ntgt4muminSundelin = 3;
   string targetSundelin[3] = { "Si", "S", "Ca" };

   for ( int it=0; it<ntgt4muminSingerExp; ++it )
   {
      SQLCompose( 48, "muminus", targetSingerExp[it], "0.GeV", "", "neutron",
	          "captureUpdate", 
		  g4version,
		  "NNeutrons",
		  "",
		  "Number of neutrons per mu- capture", "event rate",
		  "",
		  "Yield of neutrons per mu- capture",
		  status, 0 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)      
   }

   for ( int it=0; it<ntgt4muminSundelin; ++it )
   {
      SQLCompose( 48, "muminus", targetSundelin[it], "0.GeV", "", "neutron",
	          "captureUpdate", 
		  g4version,
		  "NeutronKineticEnergy",
		  "",
		  "Neutron kinetic energy [MeV]", "Number of neutrons per capture per MeV",
		  "",
		  "Neutron yield vs kinetic energy",
		  status, 0 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)         
   }

   return;

}
