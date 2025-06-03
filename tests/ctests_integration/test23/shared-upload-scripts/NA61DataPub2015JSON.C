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

const int RefPub2015=69;

// for pions and protons...
// ... although for protons it goes only up to 360mrad
//
const int NThetaBins = 11;
std::string thetamin[11] = { "0", "10", "20", "40", "60", "100", "140", "180", "240", "300", "360" };
std::string thetamax[11] = { "10", "20", "40", "60", "100", "140", "180", "240", "300", "360", "420" };

// for charged kaons...
// ... although for K- it goes only to 240mrad
//
const int NThetaBinsK = 8;
std::string Kthetamin[8] = { "0", "20", "40", "60", "100", "140", "180", "240" };
std::string Kthetamax[8] = { "20", "40", "60", "100", "140", "180", "240", "300" };
  
const int NThetaBinsK0s = 7;
std::string K0sthetamin[7] = { "0", "40", "60", "100", "140", "180", "240" };
std::string K0sthetamax[7] = { "40", "60", "100", "140", "180", "240", "300" };

const int NThetaBinsLambda = 8;
std::string Lthetamin[8] = { "0", "40", "60", "100", "140", "180", "240", "300" };
std::string Lthetamax[8] = { "40", "60", "100", "140", "180", "240", "300", "420" };

int NBins = 0;
std::string* thmin = 0;
std::string* thmax = 0;

void NA61DataPub2015JSONCompose( std::ofstream&, std::string );
void DumpNA61( std::ofstream& );

void NA61DataPub2015Master( std::string secondary )
{

   std::string fjsonname = "NA61-spectra-Pub2015-proton-C-" + secondary + ".json";
   
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;

   NA61DataPub2015JSONCompose( jsonfile, secondary );

   jsonfile << "]}" << std::endl;
   
   return;

}

void NA61DataPub2015JSONCompose( std::ofstream& jsonfile,  std::string secondary )
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
   else if ( secondary == "kplus" )
   {
      sec = "K+";
   }
   else if ( secondary == "kminus" )
   {
      sec = "K-";
   }
   
   int seclink = -1;
   if ( secondary == "proton" )
   {
      seclink = 2212;
      NBins = NThetaBins;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = thetamin[i];
	 thmax[i] = thetamax[i];
      }
   }
   else if ( secondary == "piplus" )
   {
      seclink = 211;
      NBins = NThetaBins;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = thetamin[i];
	 thmax[i] = thetamax[i];
      }      
   } 
   else if ( secondary == "piminus" )
   {
      seclink = -211;
      NBins = NThetaBins;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = thetamin[i];
	 thmax[i] = thetamax[i];
      }
   }
   else if ( secondary == "kplus" )
   {
      seclink = 321;
      NBins = NThetaBinsK;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = Kthetamin[i];
	 thmax[i] = Kthetamax[i];
      }
   }
   else if ( secondary == "kminus" )
   {
      seclink = -321;
      NBins = NThetaBinsK;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = Kthetamin[i];
	 thmax[i] = Kthetamax[i];
      }
   }
   else if ( secondary == "K0s" )
   {
      seclink = 310;
      NBins = NThetaBinsK0s;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = K0sthetamin[i];
	 thmax[i] = K0sthetamax[i];
      }
   }
   else if ( secondary == "Lambda" )
   {
      seclink = 3122;
      NBins = NThetaBinsLambda;
      if ( thmin ) delete [] thmin;
      thmin = new std::string[NBins];
      if ( thmax ) delete [] thmax;
      thmax = new std::string[NBins];
      for ( int i=0; i<NBins; ++i )
      {
         thmin[i] = Lthetamin[i];
	 thmax[i] = Lthetamax[i];
      }
   }

   for ( int ith=0; ith<NBins; ++ith )
   {

      readMomSpectrum( "proton", "C", secondary, thmin[ith], thmax[ith], "../na61-exp-data/Pub-2015/" );
      if ( NPointsNA61 <= 0 )
      {
         std::cout << " Invalid record for " << secondary << 
	              " at " << thmin[ith] << "<theta<" << thmax[ith] << std::endl;
	 continue;
      }

      jsonfile << "{" << std::endl;
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":" << RefPub2015 <<",\"mcdetaillnk\":1," << std::endl;
   
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
               << thmin[ith] << "<theta<" << thmax[ith] << " [mrad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsNA61 << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/rad/(GeV/c)]\"]," << std::endl;

      DumpNA61( jsonfile );
      
      jsonfile << "}," << std::endl; // end of datatable 

      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle [mrad]\"], \"parvalues\":[\"" << thmin[ith] 
               << "-" << thmax[ith] << "\"]" << std::endl;

      if ( ith == NBins-1 )
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

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << SecSigma[0];
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << SecSigma[i];
   }

   jsonfile << "]," << std::endl;

   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << (SecESigma[0])/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESigma[i])/2.;
   }
   jsonfile << "]," << std::endl;
   
   // stat err "minus" ("down")
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << (SecESigma[0])/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESigma[i])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err up
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << (SecESys[0])/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESys[0])/2.;
   }
   jsonfile << "]," << std::endl;

   // sys err down
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << (SecESys[0])/2.;
   for ( int i=1; i<NPointsNA61; ++i )
   {
      jsonfile << ", " << (SecESys[0])/2.;
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
