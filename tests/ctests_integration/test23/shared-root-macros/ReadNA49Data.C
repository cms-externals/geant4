#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"

// #include <boost/algorithm/string/split.hpp>
// #include <boost/algorithm/string/classification.hpp>

#include "NA49.h"

/* this part has moved to NA49.h

int    NPointsNA49         = 0;

float* xF              = 0;

float* dXSecdxF     = 0;
float* err_dXSecdxF = 0;
float* dNdxF           = 0;
float* err_dNdxF       = 0;
float* pT              = 0;
float* err_pT          = 0;
float* pT2             = 0;
float* err_pT2         = 0;

float* XSec            = 0;
float* err_XSec        = 0;

// for pions (pi+ & pi-)
//
float* Y =0;
float* dNdY[2]     = { 0, 0 };
float* err_dNdY[2] = { 0, 0 };

int NSubSetsNA49 = 0;
double*     DDiff_xF = 0;
TObjArray*  DDiffDataHolder = 0;
*/

#include "SplitString.C"

#ifndef G4VAL_READNA49DATA_C
#define G4VAL_READNA49DATA_C

void readDDiffSpectra( std::string beam, std::string target, std::string secondary,
                       bool useStatEr=true, bool useSysEr=true  )
{

// first of all, do the cleanups if necessary
   if ( DDiff_xF != 0 ) 
   {
      delete [] DDiff_xF;
      DDiff_xF = 0;
   }
   if ( DDiffDataHolder != 0 )
   {
      for ( int i=0; i<NSubSetsNA49; ++i )
      {
         DDiffDataHolder[i].Clear();
      }
      delete [] DDiffDataHolder;
      DDiffDataHolder = 0;
   }
   NSubSetsNA49 = 0;
   
   std::string dirname = std::string( gSystem->Getenv("G4ValHAD") );
   dirname += "/test23/na49-exp-data/";
      
   std::string filename = beam + "_" + target + "_" + secondary;
   
   std::string filename1;
   
   filename1 = filename + "_ddiffxsec.dat";
   
   std::string file1 = dirname + filename1;

//   std::cout << " opening datafile " << file1 << std::endl;
         
   ifstream infile1;
   infile1.open( file1.c_str() );
   
   // std::string line = "";
   char line[256];
   for ( int i=0; i<256; ++i ) line[i] = ' ';
   std::vector<std::string> tokens;
   std::string del = " ";
   int counter = -1;
   
   while ( !infile1.eof()  )
   {

      for ( int i=0; i<256; ++i ) line[i] = ' '; // cleanup the read-in buffer before next iteration

      infile1.getline( line, 256 );

      std::string line1 = std::string( line ); // , 256 );
      
      if ( line1.find("#") == 0 ) // comment line
      {	 
	 continue;
      }
      if ( line1.find_first_not_of(del) == std::string::npos ) 
      {
	 continue;
      }
      if ( line1.find( "<NSubSets>" ) != std::string::npos )
      {
         infile1 >> NSubSetsNA49;
	 // Here I have to do some memory management !!!... or maybe not ???
	 DDiff_xF = new double[NSubSetsNA49];
	 DDiffDataHolder = new TObjArray[NSubSetsNA49];
	 continue; 
      }
      if ( line1.find("<nextXF>") != std::string::npos )
      {
	 for ( int i=0; i<line1.size(); ++i ) line[i] = ' ';
	 infile1.getline( line, 256  );
	 line1 = std::string( line );
	 // first line containing data in a subset
	 // now need to split into tokens
	 tokens.clear();
	 SplitString( line1, del, tokens );
	 if ( tokens.size() != 4 )
	 {
	    std::cout << " EMERGENCY !!!" << std::endl;
	 }
	 int sz = tokens.size();
	 counter++;
	 DDiff_xF[counter] = atof(tokens[0].c_str());
	 double x  = atof(tokens[1].c_str());
	 double y  = atof(tokens[2].c_str());
	 double ey =  atof(tokens[3].c_str()); 
	 DDiffDataHolder[counter].Add( new TVector3( x, y, ey ) );
	 	 
	 continue;
      }
      // we get here for "normal" data line
      tokens.clear();
      SplitString( line1, del, tokens );
      DDiff_xF[counter] = atof(tokens[0].c_str());
      double x  = atof(tokens[1].c_str());
      double y  = atof(tokens[2].c_str());
      double ey =  atof(tokens[3].c_str()); 
      // now add 3.8% systematic error !!! (see Eur. Phys. J. C49 (2007) 897-917)
      double syst = y * 0.038;
      double err2 = 0.;
      if ( useStatEr ) err2 += ey*ey;
      if ( useSysEr )  err2 += syst*syst;
      ey = sqrt(err2);
      DDiffDataHolder[counter].Add( new TVector3( x, y, ey ) );
   }

//   for ( int i=0; i<NSubSetsNA49; ++i )
//   {
//      std::cout << " size of subset " << i << " is " << DDiffDataHolder[i].GetEntries() << std::endl;
//   }
   
   return;

}

TGraphErrors* get1DDiffXSecAsGraph( const int icount )
{

  if ( !DDiff_xF || !DDiffDataHolder ) return 0;

  // in principle, it's not safe as it explicitly assumes 
  //  that the data are in for a particular bea.target/secondary...

   int NPt = DDiffDataHolder[icount].GetEntries();
   double*    tmpPT    = new double[NPt];
   double*    tmpXSec  = new double[NPt];
   double*    tmpEXSec = new double[NPt];
   
   // std::cout << " icount = " << icount << " xF = " << DDiff_xF[icount] << std::endl;
   
   for ( int i=0; i<NPt; ++i )
   {
      
      TVector3* tmpVec = (TVector3*)(DDiffDataHolder[icount].At(i));
      tmpPT[i]    = tmpVec->x();
      tmpXSec[i]  = tmpVec->y();
      tmpEXSec[i] = tmpVec->z();
      // std::cout << " pt, xsec, exsec = " << tmpPT[i] << " " << tmpXSec[i] << " " << tmpEXSec[i] << std::endl;
   }

   //   delete [] tmpPt;
   // delete [] tmpXSec;
   // delete [] tmpEXSec;
   
   return new TGraphErrors( NPt, tmpPT, tmpXSec, 0, tmpEXSec );


}

/*
void SplitString( const std::string& str, const std::string del, std::vector<std::string>& tokens )
{
         
   string::size_type last = str.find_first_not_of(del, 0);
   string::size_type pos  = str.find_first_of(del, last);
   while (std::string::npos != pos || std::string::npos != last)
   {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(last, pos-last));
      // Skip delimiters.  Note the "not_of"
      last = str.find_first_not_of(del, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(del, last);
   }

   return;

}
*/

void readIntegratedSpectra( std::string beam, std::string target, std::string secondary, 
                            bool useStatEr=true, bool useSysEr=true )
{

   std::string dirname = std::string( gSystem->Getenv("G4ValHAD") );
   dirname += "/test23/na49-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
      
   std::string filename1;
   
   if ( secondary == "pion" )
   { 
      filename1 = filename + "_dndy.dat";
   }
   else
   {
      filename1 = filename + "_integr.dat";
   }
   
   std::string file1 = dirname + filename1;
   
//   std::cout << " input file = " << file1 << std::endl;
   
   ifstream infile1;
   infile1.open( file1.c_str() );

   infile1 >> NPointsNA49;
   
   if ( secondary == "pion" )
   {
      if (Y) delete Y;
      Y = new float[NPointsNA49];
      if ( dNdY[0] ) delete dNdY[0];
      if ( dNdY[1] ) delete dNdY[1];
      dNdY[0] = new float[NPointsNA49]; // pi+
      dNdY[1] = new float[NPointsNA49]; // pi-
      if ( err_dNdY[0] ) delete err_dNdY[0];
      if ( err_dNdY[1] ) delete err_dNdY[1];      
      err_dNdY[0] = new float[NPointsNA49];
      err_dNdY[1] = new float[NPointsNA49];
      for ( int i=0; i<NPointsNA49; i++ )
      {
         infile1 >> Y[i] >> dNdY[0][i] >> dNdY[1][i];
	 err_dNdY[0][i] = dNdY[0][i]*0.02; // common systematic errors are said to be 2%
	 err_dNdY[1][i] = dNdY[1][i]*0.02; // but that's for the dN/dy spectra only !!!
      } 
      return;     
   }
   

   if (xF) delete [] xF;
   xF = new float[NPointsNA49];
   
  
   if (dNdxF) delete dNdxF;
   dNdxF = new float[NPointsNA49];
   
   if (err_dNdxF) delete err_dNdxF;
   err_dNdxF = new float[NPointsNA49];
   
   if ( secondary != "neutron" )
   {
     if  (dXSecdxF) delete [] dXSecdxF;
     dXSecdxF = new float[NPointsNA49];
     if  (err_dXSecdxF) delete [] err_dXSecdxF;
     err_dXSecdxF = new float[NPointsNA49];
     if (pT) delete [] pT;
     pT = new float[NPointsNA49];
     if ( err_pT ) delete [] err_pT;
     err_pT = new float[NPointsNA49];
     if (pT2) delete [] pT2;
     pT2 = new float[NPointsNA49];
     if ( err_pT2 ) delete [] err_pT2;
     err_pT2 = new float[NPointsNA49];
   }

   float syserr = 0.;
   if ( beam == "H" )
   {
      syserr = 2; // in %
   }
   else
   {
      syserr = 2.5;
   }
   float err2 = 0.;
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( secondary == "neutron" )
      {
         infile1 >> xF[i] >> dNdxF[i] >> err_dNdxF[i];
	 err2 = 0.; 
         if ( useStatEr ) err2 += err_dNdxF[i]*err_dNdxF[i];
         if ( useSysEr )  err2 += syserr*syserr;
	 err_dNdxF[i]  = std::sqrt(err2)*0.01;
	 err_dNdxF[i] *= dNdxF[i];
      }
      else
      {
         infile1 >> xF[i] >> dXSecdxF[i] >> err_dXSecdxF[i]  >> dNdxF[i] >> err_dNdxF[i]  
	         >> pT[i] >> err_pT[i] >> pT2[i] >> err_pT2[i];
	 //
	 // here I also should account for the 2% (p+p) or 2.5% (p+C) systematics
	 //
	 err2 = err_dXSecdxF[i]*err_dXSecdxF[i] + syserr*syserr;
	 err_dXSecdxF[i]     = std::sqrt(err2)*0.01; 
	 err_dXSecdxF[i]    *= dXSecdxF[i];
         err2 =0.;
         if ( useStatEr ) err2 += err_dNdxF[i]*err_dNdxF[i];
         if ( useSysEr )  err2 += syserr*syserr;	 
	 err_dNdxF[i]        = std::sqrt(err2) * 0.01;
	 err_dNdxF[i]       *= dNdxF[i];
	 err2 = err_pT[i]*err_pT[i] + syserr*syserr;
	 err_pT[i]           = std::sqrt(err2)*0.01;
	 err_pT[i]          *= pT[i];
	 err2 = 0.;
         if ( useStatEr ) err2 += err_pT2[i]*err_pT2[i];
         if ( useSysEr )  err2 += syserr*syserr;
	 err_pT2[i]          = std::sqrt(err2)*0.01;
	 err_pT2[i]         *= pT2[i];
      }
   }
   
//   std::cout << " dNdxF" << std::endl;
//   for ( int np=0; np<NPointsNA49; ++np )
//   {
//      std::cout << dNdxF[np] << ", ";
//   }
//   std::cout << std::endl;
//   std::cout << " err_dNdxF" << std::endl;
//   
//   for ( int np=0; np<NPointsNA49; ++np )
//   {
//      std::cout << err_dNdxF[np] << ", ";
//   }
//   std::cout << std::endl;

   return;
   

}

TGraphErrors* getIntegratedSpectrumAsGraph( std::string histo )
{

   double ymin = 1000000.;
   double ymax = -1.;

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }   

   //   delete [] Value;
   //   delete [] Error;

   TGraphErrors* gr = new TGraphErrors( NPointsNA49, xF, Value, 0, Error );
   gr->SetMinimum(ymin);
   gr->SetMaximum(ymax);
   return gr;

}

#endif
