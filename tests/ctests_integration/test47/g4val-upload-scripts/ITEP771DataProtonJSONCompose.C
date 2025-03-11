/*
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
*/

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"

#include "../PlotITEPAnalysis.C"

#include "DumpITEP771.C"

void ITEP771DataJSONCompose( std::ofstream&, std::string, int, std::string );

void ITEP771DataJSONMaster( std::string tgt, int tgtlnk )
{

   std::string fjsonname = "ITEP771-spectra-proton-";
   fjsonname += tgt;
   fjsonname += ".json";
   
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;
   
   ITEP771DataJSONCompose( jsonfile, tgt, tgtlnk, "proton" );
      
   ITEP771DataJSONCompose( jsonfile, tgt, tgtlnk, "neutron" );
      
   jsonfile << "]}" << std::endl;
    
   return;

}

void ITEP771DataJSONCompose( std::ofstream& jsonfile, std::string tgt, int tgtlnk, std::string secondary )
{
   

   std::string sec = secondary;
   if ( secondary == "piplus" )
   {
      sec = "pi+";
   }
   else if ( secondary == "piminus" )
   {
      sec = "pi-";
   }
   
   int seclink = -1;
   if ( secondary == "proton" )
   {
      seclink = 2212;
   }
   else if ( secondary == "neutron" )
   {
      seclink = 2112;
   } 
   
   std::string angles[4] = { "59.1", "89.0", "119.0", "159.6" };

   //
   // 1GeV/c p+A
   //

   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "1.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=50 is 1 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":50,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 1GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 1.4GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "1.40", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=51 is 1.4 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":51,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 1.4GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 2GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "2.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=52 is 2 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":52,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 2GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }

   // 3GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "3.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=37 is 3 GeV/c proton (overlap with HARP !!!) 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":37,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 3GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 5GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "5.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=38 is 5 GeV/c proton (overlap with HARP !!!)
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":38,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 5GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 6GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "6.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=53 is 6 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":53,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 6GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 6.25GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "6.25", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=54 is 6.25 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":54,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 6.25GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 6.5GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "6.50", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=55 is 6.5 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":55,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 6.5GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 7GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "7.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=56 is 7 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":56,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 7GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 8.25GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "8.25", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=58 is 8.25 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":58,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 8.25GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 8.5GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "8.50", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=59 is 8.5 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":59,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 8.5GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // 9GeV/c p+A
   //
   for ( int ia=0; ia<4; ++ia )
   {
   
      bool ok = readKEData( "proton", tgt, "9.00", secondary, angles[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=60 is 9 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":60,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 9GeV/c, theta=" 
	          << angles[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
         jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	 jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << angles[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }
   
   // special case:
   // at 7.50GeV, more anglilar measurements for C, Cu, Pb, U targets
   //
   std::vector<std::string> ThetaBin;
   ThetaBin.push_back("10.1");
   ThetaBin.push_back("15.0");
   ThetaBin.push_back("19.8");
   ThetaBin.push_back("24.8");
   ThetaBin.push_back("29.5");
   ThetaBin.push_back("34.6");
   ThetaBin.push_back("39.6");
   ThetaBin.push_back("44.3");
   ThetaBin.push_back("49.3");
   ThetaBin.push_back("54.2");
   ThetaBin.push_back("59.1");
   ThetaBin.push_back("64.1");
   ThetaBin.push_back("69.1");
   ThetaBin.push_back("74.1");
   ThetaBin.push_back("79.1");
   ThetaBin.push_back("84.1");
   ThetaBin.push_back("89.0");
   ThetaBin.push_back("98.9");
   ThetaBin.push_back("108.9");
   ThetaBin.push_back("119.0");
   ThetaBin.push_back("129.1");
   ThetaBin.push_back("139.1");
   ThetaBin.push_back("149.3");
   ThetaBin.push_back("159.6");
   ThetaBin.push_back("161.4");
   ThetaBin.push_back("165.5");
   ThetaBin.push_back("169.5");
   ThetaBin.push_back("173.5");
   ThetaBin.push_back("177.0");
   
   int na = ThetaBin.size();
   for ( int ia=0; ia<na; ++ia )
   {

      bool ok = readKEData( "proton", tgt, "7.50", secondary, ThetaBin[ia] ); 
   
      if ( ok )
      {      
         jsonfile << "{" << std::endl;
      
         // NOTE: referencelnk=17 Yu.D.Bayukov in Sov.J.Nucl.Phys 1985 (in English)
         //
         jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":17,\"mcdetaillnk\":1," << std::endl;
   
         // NOTE: observablelnk=9 means Lorentz-inv E(d3sigma/dp3)
         //
         // NOTE: beamlnk=57 is 7.5 GeV/c proton 
         //
         // NOTE: reactionlnk=1 means "particle production"
         //
         jsonfile << "\"beamlnk\":57,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":9,\"secondarylnk\":" 
                  << seclink << ",\"reactionlnk\":1," << std::endl;

         jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
         jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
         jsonfile << "\"title\":\"Production of " << sec << " in proton-"<< tgt << " interactions at 7.5GeV/c, theta=" 
	          << ThetaBin[ia] << " [deg]\"," << std::endl;
         jsonfile << "\"npoints\":0," << std::endl;
         jsonfile << "\"nbins\":[" << NPtKE << "]," << std::endl;
         jsonfile << "\"axisTitle\":[\"kinetic energy of " << sec << " [GeV]\",\"E(d3sigma/dp3) [mb/GeV**2]\"]," << std::endl;
      
         DumpITEP771( jsonfile, secondary );

         jsonfile << "}," << std::endl; // end of datatable 
      
        jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl; 
	jsonfile << "\"parnames\":[\"theta [deg] of secondary particle\"],\"parvalues\":[\"" << ThetaBin[ia]
	          << "\"]" << std::endl;

         jsonfile << "}," << std::endl; 
      }
   }

         
   return;

}

