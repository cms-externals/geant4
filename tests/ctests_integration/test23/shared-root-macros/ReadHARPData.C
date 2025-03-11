/* this part has moved into HARP.h
// read exp.data
//
// FW data (forward = small angle with respect to projectile)
//
const int NSetsFW = 4;
//
// LA data (lateral = larger angle with respect to projectile)
//
const int NSetsLA = 9;

// general purpose counter
//
static int NSetsHARP = 0;

float** AngleBinHARP = 0;
int*    NPointsHARP = 0; 
float** XMinHARP = 0;
float** XMaxHARP = 0;
float** YHARP = 0;
float** EYHARP = 0;
float   NORM_UNCERTANTY = 0.02;

// NOTE: Based on the info from the ref. below, errors on the HARP data are published
//       as a total of stat&sys (except 2% overal normalization uncertanty which is not included)
//       See here:
//       http://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2008-26.pdf
//       Journal publications are not very explicit on this matter
*/

#include "HARP.h"

#include "TGraphErrors.h"

#ifndef G4VAL_READHARPDATA_C
#define G4VAL_READHARPDATA_C

void ReadHARPData( std::string beam, std::string target, std::string energy, 
                   std::string secondary, 
		   std::string region ) // should be either "FW" or "LA"
{

/* ---> old Wilson cluster
   std::string dirname = std::string( gSystem->Getenv("G4INSTALL") );
   if ( dirname.empty() )
   {
      std::cout << " Please set G4INSTALL env.var. to point to the area where your Geant4 source code resides" << std::endl;
      std::cout << " Example: export G4INSTALL=/path/to/geant4-XX-YY-ref-ZZ" << std::endl;
      exit(1);
   }
   dirname += "/tests/test23/harp-exp-data/";
*/
   
   std::string dirname = std::string( gSystem->Getenv("G4ValHAD") );
   dirname += "/test23/harp-exp-data/";
   
   dirname += ( target + "/" );
   
   std::string filename = beam + "_" + target + "_" + energy + "GeV_" + secondary + "_" + region + ".dat";
   
   std::string file = dirname + filename;
   
//   std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   if ( !infile.is_open() )
   {
   
      std::cout << " can NOT open data file " << file << std::endl;
   
   }
      
   if (AngleBinHARP) 
   {
      delete [] AngleBinHARP[0];
      delete [] AngleBinHARP[1];
      delete [] AngleBinHARP;
   }    
   
   if ( NPointsHARP ) delete [] NPointsHARP;   
   
   int i = 0;
   
   if ( NSetsHARP > 0 )
   {
      if ( XMinHARP )
      {
         for ( i=0; i<NSetsHARP; i++ )
         {
            delete [] XMinHARP[i];
         }
         delete [] XMinHARP;
      }   
      if ( XMaxHARP )
      {
         for ( i=0; i<NSetsHARP; i++ )
         {
            delete [] XMaxHARP[i];
         }
         delete [] XMaxHARP;
      }   
      if ( YHARP )
      {
         for ( i=0; i<NSetsHARP; i++ )
         {
            delete [] YHARP[i];
         }
         delete [] YHARP;
      }      
      if ( EYHARP )
      {
         for ( i=0; i<NSetsHARP; i++ )
         {
            delete [] EYHARP[i];
         }
         delete [] EYHARP;
      }
      NSetsHARP = 0;
   }   

   if ( region == "FW" )
   {
      NSetsHARP = NSetsFW;
   }
   else if ( region == "LA" )
   {
      NSetsHARP = NSetsLA;
   }
   else
   {
      return;
   }
   
   AngleBinHARP = new float*[2];
   AngleBinHARP[0] = new float[NSetsHARP];
   AngleBinHARP[1] = new float[NSetsHARP]; 

   NPointsHARP = new int[NSetsHARP];

   XMinHARP = new float*[NSetsHARP];
   XMaxHARP = new float*[NSetsHARP];
   YHARP    = new float*[NSetsHARP];
   EYHARP   = new float*[NSetsHARP]; 

   for ( i=0; i<NSetsHARP; i++ )
   {
      infile >> AngleBinHARP[0][i] >> AngleBinHARP[1][i];
      infile >> NPointsHARP[i];
      
      // std::cout << "Angle Bin: " << AngleBin[0][i] << " " << AngleBin[1][i] << std::endl;
      // std::cout << "NPoints: " << NPoints[i] << std::endl;
      
      XMinHARP[i] = new float[NPointsHARP[i]];
      XMaxHARP[i] = new float[NPointsHARP[i]];
      YHARP[i]    = new float[NPointsHARP[i]];
      EYHARP[i]   = new float[NPointsHARP[i]];
      for ( int j=0; j<NPointsHARP[i]; j++ )
      {
         infile >> XMinHARP[i][j] >> XMaxHARP[i][j] >> YHARP[i][j] >> EYHARP[i][j];
	 YHARP[i][j] *= 1000.; // why do I miltiply by 1000. ? I think it has something with the unites (mb,...)
	 EYHARP[i][j] *= 1000.; // I think there's something in Valdimir's scripts about that...
      }
   }
   
   infile.close();

   return;

}

TGraphErrors* getHARPAsGraph( int ibin )
{

   float* X = new float[NPointsHARP[ibin]];
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      X[i] = 0.5 * (XMinHARP[ibin][i]+XMaxHARP[ibin][i]);
   }
   
   TGraphErrors* gr = new TGraphErrors( NPointsHARP[ibin], X, YHARP[ibin], 0, EYHARP[ibin] );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5); 

   delete [] X;

   return gr;
}

TGraphErrors* getHARPAsThetaGraph( float pmin, float pmax, std::string region="LA" )
{

   int NN = NSetsLA ;
   if ( region == "FW" ) NN = NSetsFW;
   
   float* X  = new float[NN];
   float* Y  = new float[NN];
   float* EY = new float[NN];
   
   for ( int i=0; i<NN; ++i )
   {
      X[i] = 0.;
      Y[i] = 0.;
      EY[i] = 0.;
   }
   
   float pmean = 0.5 * ( pmin + pmax );
   
   float ymax = -1.;
   
   int ifound=0;
  
   for ( int i=0; i<NN; ++i)
   {
//      X[i] = 0.5 * ( AngleBinHARP[0][i] + AngleBinHARP[1][i] );
//      Y[i] = 0.;
//      EY[i] = 0.;
      for ( int j=0; j<NPointsHARP[i]; ++j )
      {
         if ( pmean > XMinHARP[i][j] && pmean <= XMaxHARP[i][j] )
	 {
	    ifound++;
	    X[ifound] = 0.5 * ( AngleBinHARP[0][i] + AngleBinHARP[1][i] );
	    Y[ifound] = YHARP[i][j];
	    EY[ifound] = EYHARP[i][j];
	    if ( ymax < Y[ifound]+EY[ifound] ) ymax = Y[ifound]+EY[ifound];
	    break;
	 }
      }
   } 
   
   NN=ifound;

   TGraphErrors* gr = new TGraphErrors( NN, X, Y, 0, EY );
   gr->GetYaxis()->SetRangeUser( 0., ymax );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5); 

   delete [] X;
   delete [] Y;
   delete [] EY;

   return gr;

}

int findHARPMomBin( float pmin, float pmax, std::string region="LA" )
{

   int bin = -1;

   int NN = NSetsLA ;
   if ( region == "FW" ) NN = NSetsFW;
   
   // find largets subset
   //
   int NPtMax = 0;
   int is = -1;
   for ( int i=0; i<NN; ++i )
   {
      if ( NPointsHARP[i] > NPtMax )
      {
         NPtMax = NPointsHARP[i];
	 is = i;
      }
   }
            
   float pmean = 0.5 * ( pmin + pmax );

//   for ( int i=0; i<NN; ++ i)
//   {
//      for ( int j=0; j<NPointsHARP[i]; ++j )
      for ( int j=0; j<NPointsHARP[is]; ++j )
      {
//         if ( pmean > XMinHARP[i][j] && pmean <= XMaxHARP[i][j] )
         if ( pmean > XMinHARP[is][j] && pmean <= XMaxHARP[is][j] )
	 {
	    bin = j;
	    return bin;
	 }
      }
//   } 
   

   return bin;

}

#endif
