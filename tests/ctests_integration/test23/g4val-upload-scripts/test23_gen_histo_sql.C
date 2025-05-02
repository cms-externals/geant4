
#include "../../test23/shared-upload-scripts/SQLCompose.C"

// NOTE(JVY): "using namespace std" is declared in SQLCompose.C

void test23_gen_histo_sql( string g4version, string status )
{


   int nsecNA49 = 5;
   int nsecNA61 = 3;
   string secondaries[5] = { "piplus", "piminus", "proton", "antiproton", "neutron" };
   int nmodels = 3;
   string HEmodels[3]  = { "ftfp_bert", "qgsp_bert", "NuBeam" };
   int nbinsNA61 = 10;
   string thetaminNA61[10] = {   "0",  "20",  "40",  "60", "100",
                               "140", "180", "240", "300", "360" };
   string thetamaxNA61[10] = {  "20",  "40",  "60", "100", "140",
                               "180", "240", "300", "360", "420" };    
   // NA49 part
   //
   for ( int is=0; is<nsecNA49; ++is )
   {
      for ( int im=0; im<nmodels; ++im )
      {
         string histoname = secondaries[is] + "_" + "dNdxF";
	 SQLCompose( 23, "proton", "C", "", "158.0GeV", secondaries[is],
	         HEmodels[im], 
		 // "geant4-10-01-patch-02",
		 g4version,
		 histoname,
		 "",
		 "xF", "average multiplicity",
		 "NA49",
		 "average multiplicity vs xF",
		 status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
         if ( secondaries[is] == "neutron" ) continue;
	 histoname = secondaries[is] + "_" + "pT";
	 SQLCompose( 23, "proton", "C", "", "158.0GeV", secondaries[is],
	         HEmodels[im], 
		 // "geant4-10-01-patch-02",
		 g4version,
		 histoname,
		 "",
		 "xF", "average pT",
		 "NA49",
		 "average pT vs xF",
		 status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
      }
   }

   //
   // NA61 part
   //
   // NOTE(JVY): in the case of NA61-like validation, 
   //              the theta-bin is coded in the name of each histo
   //
 
   for ( int is=0; is<nsecNA61; ++is )
   {
      for ( int im=0; im<nmodels; ++im )
      {
         for ( int ith=0; ith<nbinsNA61; ++ith )
	 {
	    if ( secondaries[is] == "proton" && ith>=(nbinsNA61-3) ) continue;
	    string histoname = secondaries[is] + "Mult_" + thetaminNA61[ith] + "_" + thetamaxNA61[ith];
	    string reaction_details = " (theta_secondary=" + thetaminNA61[ith];
	    reaction_details += ( "-" + thetamaxNA61[ith] + "mrad)" ); 
	    SQLCompose( 23, "proton", "C", "", "31.0GeV", secondaries[is],
	                HEmodels[im], 
		        g4version,
		        histoname,
			reaction_details,
		        "momentum(GeV/c) of secondary", "production rate",
		        "NA61",
		        "production rate vs momentum of secondary",
		        status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
         }
      }
   }
   // K+/pi+ ratio
   //
   for ( int im=0; im<nmodels; ++im )
   {
      SQLCompose( 23, "proton", "C", "", "31.0GeV", "K+,pi+",
	          HEmodels[im], 
		  g4version,
		  "kplus2piplus_20_140",
		  "(theta_secondary=20-140mrad)",
		  "momentum(GeV/c) of secondary", "K+/pi+ ratio",
		  "NA61",
		  "K+/pi+ vs momentum of secondaries",
		  status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
      SQLCompose( 23, "proton", "C", "", "31.0GeV", "K+,pi+",
	          HEmodels[im], 
		  g4version,
		  "kplus2piplus_140_240",
		  "(theta_secondary=140-240mrad)",
		  "momentum(GeV/c) of secondary", "K+/pi+ ratio",
		  "NA61",
		  "K+/pi+ vs momentum of secondaries",
		  status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
   }   
   
   // NOTE-3(JVY): This is towards HARP-like validation
   //              Unlike NA61, the theta-bin is NOT hardcoded, 
   //              but is "indexed" in either FW (0-3) or LA (0-8) region  

   // NOTE-4(JVY): OK, at the conceptual level, everything appears to be operational but...
   //              we bump (again) into Root/Cint limitation - too mant open files... 
   //              so we need to make it into a compiled code !

   int nsecHARP = 2; // technically speaking, it should be 3 since there're HARP data on proton production...
   int ntgtHARP = 2;
   string targetHARP[2] = { "Be", "C" }; // obviously, there're more targets but this is what we do at present
                                         // in addition, p+Ta at 8GeV/c will be treated as a special case for now
   int nbeamsHARP = 3;
   string beamHARP[3] = { "proton", "piplus", "piminus" };
   
   int nmomHARP = 4; // technically speaking, there's an extra point: p+Be at 8.9GeV/c
   string momHARP[4] = { "3.0GeV", "5.0GeV", "8.0GeV", "12.0GeV" };
   
   int nmodelsIE = 2;
   string IEmodels[2]  = { "ftfp_bert", "NuBeam" }; // no MC results for inclxx for until 10.1.p02
   
   int nHARPFW = 4;
   string thetaminHARPFW[4] = { "0.05", "0.1", "0.15", "0.2" };
   string thetamaxHARPFW[4] = { "0.1", "0.15", "0.2", "0.25" }; 
   int nHARPLA = 9;
   string thetaminHARPLA[9] = { "0.35", "0.55", "0.75", "0.95", "1.15", "1.35", "1.55", "1.75", "1.95" };
   string thetamaxHARPLA[9] = { "0.55", "0.75", "0.95", "1.15", "1.35", "1.55", "1.75", "1.95", "2.15" };
   
   for ( int im=0; im<nmodelsIE; ++im ) // loop over models
   {
      for ( int ib=0; ib<nbeamsHARP; ++ib )
      {
         for ( int it=0; it<ntgtHARP; ++it )
	 {
	    for ( int imm=0; imm<nmomHARP; ++imm )
	    {
	       for ( int is=0; is<nsecHARP; ++is )
               {
                  // the FW part
	          for ( int ia=0; ia<4; ++ia )
	          {
                     ostringstream os;
                     os << ia;
	             string histoname = secondaries[is] + "_FW_" + os.str();
	             string reaction_details = " (theta_secondary=" + thetaminHARPFW[ia];
	             reaction_details += ( "-" + thetamaxHARPFW[ia] + "rad)" );
	             SQLCompose( 23, beamHARP[ib], targetHARP[it], "", momHARP[imm], secondaries[is],
	                         IEmodels[im], 
		                 g4version,
		                 histoname,
			         reaction_details,
		                 "momentum(GeV/c) of secondary", "production rate",
		                 "HARP",
		                 "production rate vs momentum of secondary",
		                 status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
	          }
	          // the LA part
		  for ( int ia=0; ia<nHARPLA; ++ia )
		  {
                     ostringstream os;
                     os << ia;
	             string histoname = secondaries[is] + "_LA_" + os.str();
	             string reaction_details = " (theta_secondary=" + thetaminHARPLA[ia];
	             reaction_details += ( "-" + thetamaxHARPLA[ia] + "rad)" );
	             SQLCompose( 23, beamHARP[ib], targetHARP[it], "", momHARP[imm], secondaries[is],
	                         IEmodels[im], 
		                 g4version,
		                 histoname,
			         reaction_details,
		                 "momentum(GeV/c) of secondary", "production rate",
		                 "HARP",
		                 "production rate vs momentum of secondary",
		                 status, 1 ); // last arg mask=1 --> beam momentum (mask=0 --> energy)
		  }
		  //
		  // FIXME !!! Need to add theta spectra in diff. mom bins !!!!! 
		  //
               } // end loop over secondaries
	    } // end loop over momenta
	 } // end loop over targets
      } // end loop over beam types
   } // end loop over models
   
   // NOTE(JVY): Here I also have to add 2 special cases:
   //            a) p+Ta at 8GeV/c 
   //            b) p+Be at 8.9GeV/c


   return;

}
