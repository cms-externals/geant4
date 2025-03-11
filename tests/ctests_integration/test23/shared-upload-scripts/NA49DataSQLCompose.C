#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

// #include <list>

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"
// #include "TROOT.h"
// #include "TRint.h"
// #include "TObject.h"
// #include "TFile.h"
// #include "TH1F.h"
// #include "TCanvas.h"
// #include "TStyle.h"
// #include "TGraph.h"

#include "../shared-root-macros/NA49.h"
#include "../shared-root-macros/ReadNA49Data.C"

double xFmin_pi[] = { -0.055, -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035,
                       0.045,  0.0675, 0.0875, 0.1125, 0.1375, 0.175, 0.225, 0.275, 0.359, 0.450 } ;  
double xFmax_pi[] = { -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045,
                       0.0675, 0.0875, 0.1125, 0.1375, 0.175, 0.225, 0.275, 0.359, 0.450, 0.550 } ;  

void NA49DataSQLCompose()
{
     
   std::string fsqlname = "NA49-integrated-spectra.sql";
   std::ofstream sqlfile;
   sqlfile.open( fsqlname.c_str(), ios::ate|ios::app ); 
   

   // 158 GeV/c p+C -> pi+

   readIntegratedSpectra( "proton", "C", "piplus", true, false ); // last arg=false means "no sys errors added !"
                                                                  // will do it explicitly here: sys.err=2%
  
//   std::cout << "NPointsNA49=" << NPointsNA49 << std::endl;
   
   std::string title = "'Production of pi+ in proton-Carbon interactions at 158GeV/c (NA49)'";
  
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;
   
   DumpSQL_xF( sqlfile );
      
   // now the METADATA !!!

   int beam = 7; // this is an internal ID which actually means beamtype+energy/momentum entry in the table BEAMTABLE
   int target = 6 ; // Carbon... for H it'll be 1
   int observable = 1; // diff xsec
   int secondary = 211;
//   int reaction = 1; // partile production
   int dtcount = 181;
//   int score = 2; // N/A - because we only score MC results
//   int access = 1; // 1=public
   sqlfile << "INSERT INTO PUBLIC.METADATA(TID,MCDTID,BEAM,TARGET,OBSERVABLE,SECONDARY,REACTION,DATATABLEID,SCORE,ACCESS,PARNAMES,PARVALUES) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "10001,1," << std::endl; // the 2nd parameters means exp.data (=1) or MC detailed ID (Geant4, GENIE,...)
   sqlfile << beam << "," << target << "," << observable << "," << secondary << ",1," << dtcount << ",2,1,null,null);" << std::endl;


   // <pT> vs xF  spectrum for p+C -> pi+
   //
   // data table
   //
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;
   DumpSQL_pT( sqlfile );
   //
   // METADATA
   //
   dtcount++; // 182 in this case
   sqlfile << "INSERT INTO PUBLIC.METADATA(TID,MCDTID,BEAM,TARGET,OBSERVABLE,SECONDARY,REACTION,DATATABLEID,SCORE,ACCESS,PARNAMES,PARVALUES) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "10001,1," << std::endl; // the 2nd parameters means exp.data (=1) or MC detailed ID (Geant4, GENIE,...)
   sqlfile << beam << "," << target << "," << observable << "," << secondary << ",1," << dtcount << ",2,1,null,null);" << std::endl;
   
   
   // 158 GeV/c p+C -> pi-
   
   readIntegratedSpectra( "proton", "C", "piminus", true, false ); // last arg=false means "no sys errors added !"
                                                                   // will do it explicitly here: sys.err=2%
   std::string title = "'Production of pi- in proton-Carbon interactions at 158GeV/c (NA49)'";
  
   // dN/dxF spectrum
   //
   
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;   
   DumpSQL_xF( sqlfile );

   // now the METADATA again !!!
   // METADATA
   //
   secondary = -211;
   dtcount++; // 183 in this case
   sqlfile << "INSERT INTO PUBLIC.METADATA(TID,MCDTID,BEAM,TARGET,OBSERVABLE,SECONDARY,REACTION,DATATABLEID,SCORE,ACCESS,PARNAMES,PARVALUES) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "10001,1," << std::endl; // the 2nd parameters means exp.data (=1) or MC detailed ID (Geant4, GENIE,...)
   sqlfile << beam << "," << target << "," << observable << "," << secondary << ",1," << dtcount << ",2,1,null,null);" << std::endl;


   // <pT> vs xF spectrum
   //
   
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;   
   DumpSQL_pT( sqlfile );

   // now the METADATA again !!!
   // METADATA
   //
   dtcount++; // 184 in this case
   sqlfile << "INSERT INTO PUBLIC.METADATA(TID,MCDTID,BEAM,TARGET,OBSERVABLE,SECONDARY,REACTION,DATATABLEID,SCORE,ACCESS,PARNAMES,PARVALUES) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "10001,1," << std::endl; // the 2nd parameters means exp.data (=1) or MC detailed ID (Geant4, GENIE,...)
   sqlfile << beam << "," << target << "," << observable << "," << secondary << ",1," << dtcount << ",2,1,null,null);" << std::endl;


   return; 

}  
  
void DumpSQL_xF(   std::ofstream& sqlfile )
{

/*
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;
*/
   // dN/dxF value
   //
   sqlfile << "'{" << dNdxF[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << dNdxF[i];
   }
   sqlfile << "}'," << std::endl;
   
   // stat error "up" - half of the total value
   //
   sqlfile << "'{" << err_dNdxF[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << err_dNdxF[i]/2.;
   }
   sqlfile << "}'," << std::endl;

   // stat error "down" - half of the total value
   //
   sqlfile << "'{" << err_dNdxF[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << err_dNdxF[i]/2.;
   }
   sqlfile << "}'," << std::endl;

   // sys error "up" - half of the 2.5% of the total xsec value
   //
   sqlfile << "'{" << (0.025*dNdxF[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << (0.025*dNdxF[i])/2.;
   }
   sqlfile << "}'," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   sqlfile << "'{" << (0.025*dNdxF[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << (0.025*dNdxF[i])/2.;
   }
   sqlfile << "}'," << std::endl;

   // xf bin(s) - min edges
   //
   sqlfile << "'{" << xFmin_pi[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << xFmin_pi[i];
   }
   sqlfile << "}'," << std::endl;

   // xf bin(s) - max edges
   //
   sqlfile << "'{" << xFmax_pi[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << xFmax_pi[i];
   }
   sqlfile << "}');" << endl;

   return;

}

void DumpSQL_pT(   std::ofstream& sqlfile )
{

/*
   sqlfile << "INSERT INTO PUBLIC.DATATABLE(DTYPE,TITLE,NPOINTS,NBINS,AXIS_TITLE,VAL,ERR_STAT_PLUS,ERR_STAT_MINUS,ERR_SYS_PLUS,ERR_SYS_MINUS,BIN_MIN,BIN_MAX) " << std::endl;
   sqlfile << "VALUES" << std::endl;
   sqlfile << "(" << std::endl;
   sqlfile << "1," << title << ",null,'{" << NPointsNA49 << "}'," << std::endl;
   sqlfile << "'{\" xF \"}'," << std::endl;
*/
   // dN/dxF value
   //
   sqlfile << "'{" << pT[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << pT[i];
   }
   sqlfile << "}'," << std::endl;
   
   // stat error "up" - half of the total value
   //
   sqlfile << "'{" << err_pT[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << err_pT[i]/2.;
   }
   sqlfile << "}'," << std::endl;

   // stat error "down" - half of the total value
   //
   sqlfile << "'{" << err_pT[0]/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << err_pT[i]/2.;
   }
   sqlfile << "}'," << std::endl;

   // sys error "up" - half of the 2.5% of the total xsec value
   //
   sqlfile << "'{" << (0.025*pT[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << (0.025*pT[i])/2.;
   }
   sqlfile << "}'," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   sqlfile << "'{" << (0.025*pT[0])/2.;
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << (0.025*pT[i])/2.;
   }
   sqlfile << "}'," << std::endl;

   // xf bin(s) - min edges
   //
   sqlfile << "'{" << xFmin_pi[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << xFmin_pi[i];
   }
   sqlfile << "}'," << std::endl;

   // xf bin(s) - max edges
   //
   sqlfile << "'{" << xFmax_pi[0];
   for ( int i=1; i<NPointsNA49; ++i )
   {
      sqlfile << ", " << xFmax_pi[i];
   }
   sqlfile << "}');" << endl;

   return;

}
