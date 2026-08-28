#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include "NA61.h"

/* this part has moved to NA61.h

int    NPointsNA61    = 0;

float* SecMom     = 0;

float* SecSigma   = 0; 
float* SecESigma  = 0;
float* SecESys    = 0;

float* K2PiRatio  = 0;
float* K2PiERatio = 0;
*/

#ifndef G4VAL_READNA61DATA_C
#define G4VAL_READNA61DATA_C


void readMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max,
		      std::string dirname="../test23/na61-exp-data/",
                      bool useStatEr=true, bool useSysEr=true )
{

//   std::string dirname = "../test23/na61-exp-data/";
   std::string sec = secondary;
   if ( secondary == "lambda" ) sec = "Lambda";
   if ( secondary == "k0s" ) sec = "K0s";
   
   std::string filename = beam + "_" + target + "_" + sec;
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPointsNA61;
   
   // std::cout << "NPoints = " << NPointsNA61 << std::endl;
   
   if ( PMin ) delete [] PMin;
   PMin = new float[NPointsNA61];
   if ( PMax ) delete [] PMax;
   PMax = new float[NPointsNA61];   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPointsNA61];
   
   if ( SecSigma ) delete [] SecSigma;
   SecSigma = new float[NPointsNA61];
   
   if ( SecESigma ) delete [] SecESigma;
   SecESigma = new float[NPointsNA61];
   
   if ( SecESys ) delete [] SecESys;
   SecESys = new float[NPointsNA61];
   
   // float pmin, pmax;
   for ( int i=0; i<NPointsNA61; i++ )
   {  
      PMin[i] = 0.;
      PMax[i] = 0.;
      SecMom[i]    = 0.;
      SecSigma[i]  = 0.;
      SecESigma[i] = 0.;
      SecESys[i]   = 0.;
      
      if ( secondary == "proton" && dirname.find("Pub-2015") == std::string::npos ) // FIXME !!! 
                                                                                    // TMP STUFF untill I frigure out 
										    // about publication of the NA61's 
										    // proton data from the 2007 run
      {
         infile >> PMin[i] >> PMax[i] >> SecSigma[i] >> SecESigma[i] ;
      }
      else
      {
         infile >> PMin[i] >> PMax[i] >> SecSigma[i] >> SecESigma[i] >> SecESys[i];
      }

      SecMom[i] = 0.5 * ( PMin[i] + PMax[i] );
      // std::cout << SecMom[i] << " " << SecSigma[i] << " " << SecESigma[i] << std::endl;
      if ( secondary != "proton" && dirname.find("Pub-2015") == std::string::npos )
      {
         // for the 2007 run data, scale to the total xsec since proton data come normalized
	 SecSigma[i]  /= 229.3 ; // 229.3 is production xsec, while inelastic xsec = 257.2 
	                         // which one should be used for normalization ???        
         SecESigma[i] /= 229.3;
         SecESys[i]   /= 229.3 ;
      }
   }
      
   infile.close();
   
   return;

}

void readKPlus2PiPlusRatio( std::string beam, std::string target, 
		            std::string sec_angle_min, std::string sec_angle_max )
{

   std::string dirname = std::string( gSystem->Getenv("G4INSTALL") );
   dirname += "/tests/ctests_integration/test23/na61-exp-data/";
   std::string filename = beam + "_" + target + "_kplus2piplus";
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPointsNA61;
   
   if ( PMin ) delete [] PMin;
   PMin = new float[NPointsNA61];
   if ( PMax ) delete [] PMax;
   PMax = new float[NPointsNA61];   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPointsNA61];

   if ( K2PiRatio ) delete [] K2PiRatio;
   K2PiRatio = new float[NPointsNA61];
   
   if ( K2PiERatio ) delete [] K2PiERatio;
   K2PiERatio = new float[NPointsNA61];
   
   // float pmin, pmax;
   for ( int i=0; i<NPointsNA61; i++ )
   {  
      PMin[i] = 0.;
      PMax[i] = 0.;
      SecMom[i] = 0.;
      K2PiRatio[i] = 0.;
      K2PiERatio[i] = 0.;
      infile >> PMin[i] >> PMax[i] >> K2PiRatio[i] >> K2PiERatio[i] ;
      SecMom[i] = 0.5 * ( PMin[i] + PMax[i] );
      // std::cout << SecMom[i] << " " << K2PiRatio[i] << " " << K2PiERatio[i] << std::endl;
   }
   
   infile.close();

   return;

}

TGraphErrors* getNA61AsGraph()
{

   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     /* if ( useStatEr ) */ tmp += SecESigma[ii]*SecESigma[ii];
     /* if ( useSysEr ) */ tmp += SecESys[ii]*SecESys[ii];
      SumErr[ii] = std::sqrt( tmp );
   }

   double ymin = 10000.;
   double ymax = -1.;
   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (SecSigma[ip]+SumErr[ip]) > ymax ) ymax = SecSigma[ip]+SumErr[ip];
      if ( (SecSigma[ip]-SumErr[ip]) < ymin ) ymin = SecSigma[ip]-SumErr[ip];
      if ( ymin < 0. ) ymin = 0.;
   }

   TGraphErrors* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );
   gr->SetMinimum(ymin);
   gr->SetMaximum(ymax);

   return gr;   

}

#endif
