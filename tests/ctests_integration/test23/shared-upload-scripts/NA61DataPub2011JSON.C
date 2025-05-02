#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"

#include "../shared-root-macros/NA61.h"
#include "../shared-root-macros/ReadNA61Data.C"

const int RefPub2011ChgPi=49;
const int RefPub2012K2Pi=50;
// const int RefPub2015=69;
const int NThetaBins = 10;
std::string thetamin[10] = { "0", "20", "40", "60", "100", "140", "180", "240", "300", "360" };
std::string thetamax[10] = { "20", "40", "60", "100", "140", "180", "240", "300", "360", "420" };
const int NThetaBinsK2Pi = 2;
std::string K2Pithetamin[2] = { "20", "140" };
std::string K2Pithetamax[2] = { "140", "240" };

float sigmaKPlus[2][7] = { 1.94, 2.25, 2.39, 2.10, 1.94, 1.49, 1.17, 2.89, 1.32, 0.55, 0., 0., 0., 0. };
float statKPlus[2][7] = { 0.26, 0.21, 0.22, 0.20, 0.18, 0.17, 0.17, 0.41, 0.17, 0.12, 0., 0., 0., 0. };
float sysKPlus[2][7] = { 0.07, 0.08, 0.08, 0.08, 0.07, 0.05, 0.08, 0.11, 0.06, 0.06, 0., 0., 0., 0. };

void NA61DataPub2011Master()
{

   std::string fjsonname = "NA61-spectra-Pub2011-proton-C.json";
   
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;

   NAA61DataPub2011JSONCompose( jsonfile, "piplus" );
   NAA61DataPub2011JSONCompose( jsonfile, "piminus" );
   NAA61DataPub2012JSONCompose( jsonfile );

   jsonfile << "]}" << std::endl;
   
   return;

}

void NAA61DataPub2011JSONCompose( std::ofstream& jsonfile,  std::string secondary )
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
   else if ( secondary == "piplus" )
   {
      seclink = 211;
   } 
   else if ( secondary == "piminus" )
   {
      seclink = -211;
   }
   
   for ( int ith=0; ith<NThetaBins; ++ith )
   {
      readMomSpectrum( "proton", "C", secondary, thetamin[ith], thetamax[ith], "../na61-exp-data/" );
      if ( NPointsNA61 <= 0 )
      {
         std::cout << " Invalid record for " << secondary << 
	              " at " << thetamin[ith] << "<theta<" << thetamax[ith] << std::endl;
	 continue;
      }
       
      jsonfile << "{" << std::endl;
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":" << RefPub2011ChgPi <<",\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE: beamlnk=36 --> 31GeVc/ proton (CERN SPS)
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":36,\"targetlnk\":6,\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;

      jsonfile << "\"title\":\"Production of " << sec << " in proton-C interactions at 31GeV/c, " 
               << thetamin[ith] << "<theta<" << thetamax[ith] << " [mrad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsNA61 << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c)]\"]," << std::endl;
      
      DumpNA61( jsonfile );
      
      jsonfile << "}," << std::endl; // end of datatable 

      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle [mrad]\"], \"parvalues\":[\"" << thetamin[ith] 
               << "-" << thetamax[ith] << "\"]" << std::endl;

//      if ( ith == NThetaBins-1 )
//      {
//         jsonfile << "}" << std::endl; 
//      }
//      else
//      {
         jsonfile << "}," << std::endl; 
//      } 
   
   }

   return;

}

void NAA61DataPub2012JSONCompose( std::ofstream& jsonfile )
{

   for ( int ith=0; ith<NThetaBinsK2Pi; ++ith )
   {
      readKPlus2PiPlusRatio( "proton", "C", K2Pithetamin[ith], K2Pithetamax[ith] );
      if ( NPointsNA61 <= 0 )
      {
         std::cout << " Invalid record for K+/pi+ at " << K2Pithetamin[ith] << "<theta<" << K2Pithetamax[ith] << std::endl;
	 continue;
      }
      
      jsonfile << "{" << std::endl;
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":" << RefPub2012K2Pi <<",\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=20 --> particle production ratio (e.g. (sigma(K+)/sigma(pi+)) / dtheta dp)
      //
      // NOTE: beamlnk=36 --> 31GeVc/ proton (CERN SPS)
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":36,\"targetlnk\":6,\"observablelnk\":20,\"secondarylnk\":0,\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;

      jsonfile << "\"title\":\"K+/pi+ ratio in proton-C interactions at 31GeV/c, " 
               << K2Pithetamin[ith] << "<theta<" << K2Pithetamax[ith] << " [mrad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsNA61 << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum[GeV/c]\",\"K+/pi+ ratio\"]," << std::endl;
      
      DumpNA61K2PiRatio( jsonfile );

      jsonfile << "}," << std::endl; // end of datatable 

      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particles [mrad]\"], \"parvalues\":[\"" << K2Pithetamin[ith] 
               << "-" << K2Pithetamax[ith] << "\"]" << std::endl;

/*
      if ( ith == NThetaBinsK2Pi-1 )
      {
         jsonfile << "}" << std::endl; 
      }
      else
      {
*/
         jsonfile << "}," << std::endl; 
/*
      } 
*/
      // K+ production, do it explicitly
      
      jsonfile << "{" << std::endl;
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":" << RefPub2012K2Pi <<",\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 
      //
      // NOTE: beamlnk=36 --> 31GeVc/ proton (CERN SPS)
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":36,\"targetlnk\":6,\"observablelnk\":8,\"secondarylnk\":321,\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;

      jsonfile << "\"title\":\"Production of K+ in proton-C interactions at 31GeV/c, " 
               << K2Pithetamin[ith] << "<theta<" << K2Pithetamax[ith] << " [mrad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsNA61 << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of K+ [GeV/c]\",\"Cross section [mb/(GeV/c)]\"]," << std::endl;

      DumpNA61KPlus( jsonfile, ith ); 
      
      
      jsonfile << "}," << std::endl; // end of datatable 

      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particles [mrad]\"], \"parvalues\":[\"" << K2Pithetamin[ith] 
               << "-" << K2Pithetamax[ith] << "\"]" << std::endl;

      if ( ith == NThetaBinsK2Pi-1 )
      {
         jsonfile << "}" << std::endl; 
      }
      else
      {
         jsonfile << "}," << std::endl; 
      } 

   }
   
   return;

}

void DumpNA61( std::ofstream& jsonfile )
{

   // float XSex = 1.;
   float XSec = 229.3;
   
   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << SecSigma[0]*XSec;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << SecSigma[i]*XSec;
   }

   jsonfile << "]," << std::endl;

   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << (SecESigma[0]*XSec)/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESigma[i]*XSec)/2.;
   }
   jsonfile << "]," << std::endl;
   
   // stat err "minus" ("down")
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << (SecESigma[0]*XSec)/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESigma[i]*XSec)/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err up
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << (SecESys[0]*XSec)/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESys[0]*XSec)/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err down
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << (SecESys[0]*XSec)/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESys[0]*XSec)/2.;
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << PMin[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMin[i];
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << PMax[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMax[i];
   }
   jsonfile << "]" << std::endl; // no comma since it's the last array in the table

   return;

}

void DumpNA61K2PiRatio( std::ofstream& jsonfile )
{

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << K2PiRatio[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << K2PiRatio[i];
   }

   jsonfile << "]," << std::endl;

   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << K2PiERatio[0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << K2PiERatio[i]/2.;
   }
   jsonfile << "]," << std::endl;
   
   // stat err "minus" ("down")
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << K2PiERatio[0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << K2PiERatio[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err up
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", 0.";
   }
   jsonfile << "]," << std::endl;

   // sys err down
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", 0.";
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << PMin[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMin[i];
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << PMax[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMax[i];
   }
   jsonfile << "]" << std::endl;

   return;
   
} 

void DumpNA61KPlus( std::ofstream& jsonfile, int count )
{

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << sigmaKPlus[count][0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << sigmaKPlus[count][i];
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << statKPlus[count][0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << statKPlus[count][i]/2.;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << statKPlus[count][0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << statKPlus[count][i]/2.;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << sysKPlus[count][0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << sysKPlus[count][i]/2.;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << sysKPlus[count][0]/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << sysKPlus[count][i]/2.;
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << PMin[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMin[i];
   }
   jsonfile << "]," << std::endl;

   // momentum bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << PMax[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << PMax[i];
   }
   jsonfile << "]" << std::endl;

   return;

}
