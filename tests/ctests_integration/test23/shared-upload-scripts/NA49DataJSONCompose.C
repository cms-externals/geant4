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

void NA49DataJSONCompose( std::ofstream&, std::string );
void DumpJSON_xF(   std::ofstream& );
void DumpJSON_pT(   std::ofstream& );

double xFmin_pi[] = { -0.055, -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035,
                       0.045,  0.0675, 0.0875, 0.1125, 0.1375, 0.175, 0.225, 0.275, 0.359, 0.450 } ;  
double xFmax_pi[] = { -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045,
                       0.0675, 0.0875, 0.1125, 0.1375, 0.175, 0.225, 0.275, 0.359, 0.450, 0.550 } ;  

double xFmin_pbar[] = { -0.225, -0.175, -0.125, -0.0875, -0.0625, -0.0375, -0.0125,
			 0.0125, 0.0375, 0.075, 0.125, 0.175, 0.25 };
double xFmax_pbar[] = { -0.175, -0.125, -0.0875, -0.0625, -0.0375, -0.0125,
			 0.0125, 0.0375, 0.075, 0.125, 0.175, 0.25, 0.35 };
double xFmin_neut[] = { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.85 };
double xFmax_neut[] = { 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.85, 0.95 };

double* xmin = 0;
double* xmax = 0;

void NA49DataJSONMaster()
{

   std::string fjsonname = "NA49-proton-C-integrated-spectra.json";
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;
   
   NA49DataJSONCompose( jsonfile, "piplus" );
   
   jsonfile << "," << std::endl;
   
   NA49DataJSONCompose( jsonfile, "piminus" );
   
   jsonfile << "," << std::endl;

   NA49DataJSONCompose( jsonfile, "proton" );
   
   jsonfile << "," << std::endl;

   NA49DataJSONCompose( jsonfile, "antiproton" );

   jsonfile << "," << std::endl;

   NA49DataJSONCompose( jsonfile, "neutron" );
   
   jsonfile << "]}" << std::endl;
    
   return;

}

void NA49DataJSONCompose( std::ofstream& jsonfile, std::string secondary )
{
   

   // 158 GeV/c p+C -> pi+

   readIntegratedSpectra( "proton", "C", secondary, true, false ); // last arg=false means "no sys errors added !"
                                                                   // will do it explicitly here: sys.err=2.5%
   
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
      
   if ( NPointsNA49 <= 0 ) 
   {
      std::cout << " Invalid record, NPointsNA49 = " << NPointsNA49 << std::endl; 
      return;
   }
   
   xmin = new double[NPointsNA49];
   xmax = new double[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; ++i )
   {
      if ( secondary == "piplus" || secondary == "piminus" )
      {
         xmin[i] = xFmin_pi[i];
	 xmax[i] = xFmax_pi[i];
      }
      else if ( secondary == "antiproton" || secondary == "anti_proton" )
      {
         xmin[i] = xFmin_pbar[i];
	 xmax[i] = xFmax_pbar[i];
      }
      else if ( secondary == "neutron" )
      {
         xmin[i] = xFmin_neut[i];
	 xmax[i] = xFmax_neut[i];
      }
   }
   
   // special case - too much to fill up by hands
   //
   if ( secondary == "proton" )
   {
      for ( int i=0; i<14; ++i )
      {
         xmin[i] = -0.825 + 0.05*i;
	 xmax[i] = xmin[i] + 0.05;
//	 std::cout << i << " " << xmin[i] << " " << xmax[i] << std::endl;
      }
      xmin[14] = -0.125;
      xmax[14] = -0.0875;
//      std::cout << 14 << " " << xmin[14] << " " << xmax[14] << std::endl;
       

      for ( int i=0; i<7; ++i )
      {
         xmin[i+15] = -0.0875 + 0.025*i;
	 xmax[i+15] = xmin[i+15] + 0.025;
//	 std::cout << (i+15) << " " << xmin[i+15] << " " << xmax[i+15] << std::endl;
      }
      xmin[22] = 0.0875;
      xmax[22] = 0.125;
//      std::cout << 14 << " " << xmin[14] << " " << xmax[14] << std::endl;

      for ( int i=0; i<17; ++i )
      {
         xmin[23+i] = 0.125 + 0.05*i;
	 xmax[23+i] = xmin[i+23] + 0.05;
//	 std::cout << (i+23) << " " << xmin[i+23] << " " << xmax[i+23] << std::endl;
      }
   }
   
   jsonfile << "{" << std::endl;
   
   jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":53,\"mcdetaillnk\":1," << std::endl;
   
   // NOTE-1: observablelnk=4 means "average mult."
   //
   jsonfile << "\"beamlnk\":7,\"targetlnk\":6,\"observablelnk\":4,\"secondarylnk\":" << seclink << ",\"reactionlnk\":1," << std::endl;
   
   jsonfile << "\"datatable\":{" << std::endl;
   jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
   jsonfile << "\"title\":\"Production of " << sec << " in proton-Carbon interactions at 158GeV/c\",\"npoints\":0," << std::endl;
   jsonfile << "\"nbins\":[" << NPointsNA49 << "],\"axisTitle\":[\"xF\",\"dN/dxF\"]," << std::endl;
   
   DumpJSON_xF( jsonfile );
      
   jsonfile << "}," << std::endl; // end of datatable 
   
   jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parnames\":[], \"parvalues\":[]" << std::endl;
   
   jsonfile << "}," << std::endl;

   if ( secondary == "neutron" ) return;

   jsonfile << "{" << std::endl;
   
   jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":53,\"mcdetaillnk\":1," << std::endl;
   
   // NOTE-1: observablelnk=5 means "avereage pT"
   //
   jsonfile << "\"beamlnk\":7,\"targetlnk\":6,\"observablelnk\":5,\"secondarylnk\":" << seclink << ",\"reactionlnk\":1," << std::endl;
   
   jsonfile << "\"datatable\":{" << std::endl;
   jsonfile << "\"dtid\":1,\"datatypeslnk\":1," << std::endl;
   jsonfile << "\"title\":\"Production of " << sec << " in proton-Carbon interactions at 158GeV/c\",\"npoints\":0," << std::endl;
   jsonfile << "\"nbins\":[" << NPointsNA49 << "],\"axisTitle\":[\"xF\",\"<pT>\"]," << std::endl;

   DumpJSON_pT( jsonfile );

   jsonfile << "}," << std::endl; // end of datatable 
   
   jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parnames\":[],\"parvalues\":[]" << std::endl;
   
   jsonfile << "}" << std::endl;

   return;

}

void DumpJSON_xF(   std::ofstream& jsonfile )
{

   // dN/dxF
   //
   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << dNdxF[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << dNdxF[i];
   }
   jsonfile << "]," << std::endl;
   
   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << err_dNdxF[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_dNdxF[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // stat err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << err_dNdxF[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_dNdxF[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "up" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << (0.025*dNdxF[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.025*dNdxF[i])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << (0.025*dNdxF[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.025*dNdxF[i])/2.;
   }
   jsonfile << "]," << std::endl;
   
   // xf bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << xmin[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xmin[i];
   }
   jsonfile << "]," << std::endl;

   // xf bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << xmax[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xmax[i];
   }
   jsonfile << "]" << std::endl; // NOTE: NO comma after this last square bracket 
                                 // since it's the last one in the list of data arrays

   return;

}

void DumpJSON_pT(   std::ofstream& jsonfile )
{

   // <pT> vs xF
   //
   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << pT[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << pT[i];
   }
   jsonfile << "]," << std::endl;
   
   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << err_pT[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_pT[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // stat err "minus" ("down") - half of the total value
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << err_pT[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << err_pT[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "up" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << (0.025*pT[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.025*pT[i])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << (0.025*pT[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << (0.025*pT[i])/2.;
   }
   jsonfile << "]," << std::endl;
   
   // xf bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << xmin[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xmin[i];
   }
   jsonfile << "]," << std::endl;

   // xf bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << xmax[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      jsonfile << ", " << xmax[i];
   }
   jsonfile << "]" << std::endl; // NOTE: NO comma after this last square bracket 
                                 // since it's the last one in the list of data arrays

   return;

}
