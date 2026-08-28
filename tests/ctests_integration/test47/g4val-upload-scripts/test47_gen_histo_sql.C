
#include "../../test23/shared-upload-scripts/SQLCompose.C"

// NOTE(JVY): "using namespace std" is declared in SQLCompose.C

void test47_gen_histo_sql( string g4version, string status )
{

   int nbeams = 3;
   string beam[3] = { "proton", "piplus", "piminus" };

   int ntgtITEP = 2;
   string target[2] = { "C", "U" };
   
   int nsecITEP = 2;
   string secondary[2] = { "proton", "neutron" };

   int nangITEP = 4;
   string angle[4] = { "59.1", "89.0", "119.0", "159.6" };
      
   int nmodels = 4;
   string model[4]  = { "binary", "bertini", "inclxx", "ftfp" }; // NOTE: inclxx results are NOT available for 9.6.p04 !!!
//   int nmodels = 3;
//   string model[3]  = { "binary", "bertini", "ftfp" }; // NOTE: inclxx results are NOT available for 9.6.p04 !!!

   for ( int ib=0; ib<nbeams; ++ib )
   {

      for ( int it=0; it<ntgtITEP; ++it )
      {

         for ( int is=0; is<nsecITEP; ++is )
         {

            for ( int ia=0; ia<nangITEP; ++ia )
	    {
               // 1.4GeV/c beam for the first 3 models (except ftfp)
	       //
	       for ( int im=0; im<nmodels-1; ++im )
	       {
	          string histoname = "KE" + secondary[is] + "0" + beam[ib] + target[it] + model[im] + "1.40GeV";
	          int counter = 5 - angle[ia].length();
	          for ( int ic=0; ic<counter; ++ic )
	          {
	             histoname += " "; // pad up for the fact that initially sec.angle 
		                       // was supposed to be char[6] - no more, no less...
	          }
	          histoname += angle[ia];
		  string reaction_details = ", at " + angle[ia] +"deg";
	          SQLCompose( 47, beam[ib], target[it], "", "1.40GeV", secondary[is],
	                      model[im],
			      g4version,
			      histoname,
			      reaction_details, 
			      "kinetic energy [GeV]", "Lorentz invariant cross section [mb/GeV**2]",
			      "", // no exp-specific subdir for MC results (Root files)
			      "Lorentz invariant cross section vs kinetic energy",
			      status, 1 );
	       }

	       // either 5GeV/c (pion beam) or 7.5GeV/c(proton beam) 
	       // which is the case for all except the 1st (binary)
	       //
	       for ( int im=1; im<nmodels; ++im )
	       {
	          string mom = "5.00GeV";
		  if ( beam[ib] == "proton" )
		  {
		     mom = "7.50GeV";
		  }		  
	          string histoname = "KE" + secondary[is] + "0" + beam[ib] + target[it] + model[im] + mom;
	          int couter = 5 - angle[ia].length();
	          for ( int ic=0; ic<counter; ++ic )
	          {
	             histoname += " "; // pad up for the fact that initially sec.angle 
		                       // was supposed to be char[6] - no more, no less...
	          }
	          histoname += angle[ia];
		  string reaction_details = ", at " + angle[ia] +"deg";
	          SQLCompose( 47, beam[ib], target[it], "", mom, secondary[is],
	                      model[im],
			      g4version,
			      histoname,
			      reaction_details, 
			      "kinetic energy [GeV]", "Lorentz invariant cross section [mb/GeV**2]",
			      "", // no exp-specific subdir for MC results (Root files)
			      "Lorentz invariant cross section vs kinetic energy",
			      status, 1 );
	       }
	    }
         }
      }
   }

// NOTE(JVY): Still need to add a utility to upload p+Ta->pbar results !!!

   return;

}
