#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"

#include "../shared-root-macros/NA49.h"
#include "../shared-root-macros/ReadNA49Data.C"

void NA49DataAsGraphCompose( std::ofstream&, std::string );
void DumpXF( std::ofstream& );
void DumpPT( std::ofstream& );

void NA49Master()
{

   std::string fjsonname = "NA49-proton-H-integrated-spectra.json";
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 

   jsonfile << "{\"ResultList\":[" << std::endl;
   
   NA49DataAsGraphCompose( jsonfile, "piplus" );
  
   jsonfile << "," << std::endl;
   
   NA49DataAsGraphCompose( jsonfile, "piminus" );
   
   jsonfile << "," << std::endl;

   NA49DataAsGraphCompose( jsonfile, "kplus" );
   
   jsonfile << "," << std::endl;
   
   NA49DataAsGraphCompose( jsonfile, "kminus" );
   
   jsonfile << "," << std::endl;

   NA49DataAsGraphCompose( jsonfile, "neutron" ); // no comman afterwards

   NA49DataAsGraphCompose( jsonfile, "proton" );
   
   jsonfile << "," << std::endl;

   NA49DataAsGraphCompose( jsonfile, "antiproton" );
   
   jsonfile << "]}" << std::endl;

   return;

}


void NA49DataAsGraphCompose( std::ofstream& jsonfile, std::string secondary )
{

   // 158 GeV/c p+C -> pi+

   readIntegratedSpectra( "proton", "H", secondary, true, false ); // last arg=false means "no sys errors added !"
                                                                   // will do it explicitly here: sys.err=2% 
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
   else if ( secondary == "antiproton" || secondary == "anti_proton" )
   {
      seclink = -2212;
   }
   else if ( secondary == "piplus" )
   {
      seclink = 211;
   } 
   else if ( secondary == "piminus" )
   {
      seclink = -211;
   }
   else if ( secondary == "kplus" )
   {
      seclink = 321;
   }
   else if ( secondary == "kminus" )
   {
      seclink = -321;
   }
      
   if ( NPointsNA49 <= 0 ) 
   {
      std::cout << " Invalid record, NPointsNA49 = " << NPointsNA49 << std::endl; 
      return;
   }

   // dN/dxF
   //
   
   jsonfile << "{" << std::endl;
   
   jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":53,\"mcdetaillnk\":1," << std::endl;
   
   // NOTE-1: observablelnk=4 means "average mult."
   //
   jsonfile << "\"beamlnk\":7,\"targetlnk\":1,\"observablelnk\":4,\"secondarylnk\":" << seclink << ",\"reactionlnk\":1," << std::endl;
   
   jsonfile << "\"datatable\":{" << std::endl;
   jsonfile << "\"dtid\":1,\"datatypeslnk\":1000," << std::endl;
   jsonfile << "\"title\":\"Production of " << sec << " in proton-H interactions at 158GeV/c\",\"npoints\":" << NPointsNA49 << "," << std::endl;
   jsonfile << "\"nbins\":[],\"axisTitle\":[\"xF\",\"dN/dxF\"]," << std::endl;

   DumpXF( jsonfile );

   jsonfile << "}," << std::endl; // end of datatable 
   
   jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parnames\":[], \"parvalues\":[]" << std::endl;
   
   jsonfile << "}," << std::endl;

   if ( secondary == "neutron" ) return;

   // d<pT>/dxF
   //
   jsonfile << "{" << std::endl;
   
   jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":53,\"mcdetaillnk\":1," << std::endl;
   
   // NOTE-1: observablelnk=4 means "average mult."
   //
   jsonfile << "\"beamlnk\":7,\"targetlnk\":1,\"observablelnk\":5,\"secondarylnk\":" << seclink << ",\"reactionlnk\":1," << std::endl;
   
   jsonfile << "\"datatable\":{" << std::endl;
   jsonfile << "\"dtid\":1,\"datatypeslnk\":1000," << std::endl;
   jsonfile << "\"title\":\"Production of " << sec << " in proton-H interactions at 158GeV/c\",\"npoints\":" << NPointsNA49 << "," << std::endl;
   jsonfile << "\"nbins\":[],\"axisTitle\":[\"xF\",\"<pT>\"]," << std::endl;

   DumpPT( jsonfile );

   jsonfile << "}," << std::endl; // end of datatable 
   
   jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parnames\":[],\"parvalues\":[]" << std::endl;
   
   jsonfile << "}" << std::endl;

   return;

}

void DumpXF( std::ofstream& jsonfile )
{

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << xF[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xF[i];
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << dNdxF[i];
   }
   jsonfile << "]," << std::endl;

   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_dNdxF[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // stat err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_dNdxF[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.02*dNdxF[i])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.02*dNdxF[i])/2.;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"binMin\":[], \"binmax\":[]" << std::endl; // no comma after the last square bracket

   return;
   
}

void DumpPT( std::ofstream& jsonfile )
{

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << xF[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xF[i];
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << pT[i];
   }
   jsonfile << "]," << std::endl;

   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_pT[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // stat err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_pT[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.02*pT[i])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.02*pT[i])/2.;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"binMin\":[], \"binmax\":[]" << std::endl; // no comma after the last square bracket

   return;

}
