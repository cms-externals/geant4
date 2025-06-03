#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>
#include <algorithm>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

std::string TEST_NAME="test19";
#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/shared-root-macros/Chi2Calc.C"

const int NV = 3;
std::string version[3] = { "geant4-10-04-patch-02", "geant4-10-04-patch-01", "geant4-10-03-patch-03" };
// ---> int colorV[5] = { kRed, kGreen, kBlack, kMagenta, 7 };
int colorV[5] = { kBlack, kGreen, kRed, kMagenta, 7 };
const int NM=2;
std::string model[2] = { "qgsp", "ftfp" };
int colorM[2] = { kBlack, 7 };
int symbolM[4]     = { 25, 8, 21, 23 };

double chi2V[NV];
double chi2M = 0.;
int NDFV[NV];
int NDFM = 0;

const int NPt = 3;
std::string PT[3] = { "0.18", "0.3", "0.5" };
std::string hlabel[3] = { "pt180", "pt300", "pt500" };
std::string htitle[3] = { "pt_{sec} = 0.18 GeV/c", "pt_{sec} = 0.3 GeV/c", "pt_{sec} = 0.5 GeV/c"};

double* Data[3] = { NULL, NULL, NULL };
double* StatErr[3] = { NULL, NULL, NULL };
double* SysErr[3] = { NULL, NULL, NULL };
double* TotalErr[3] = { NULL, NULL, NULL };
int     NData[3] = { 0, 7, 3 };
double* XX[3] = { NULL, NULL, NULL };

TGraphErrors* gr[3] = { NULL, NULL, NULL };

// fwd declarations
//
void DrawSASM6ERegre( std::string, std::string, std::string, std::string, std::string, int, bool, double*, int* );
//
void DumpJSON( std::ofstream&, std::string, std::string, int, std::string );
void DumpSpectrum( std::ofstream&, int );

void GetExData( std::string beam, std::string target, std::string secondary )
{

   for ( int ipt=0; ipt<NPt; ++ipt )
   {
   
      if ( Data[ipt] )     delete [] Data[ipt];     Data[ipt] = NULL;
      if ( StatErr[ipt] )  delete [] StatErr[ipt];  StatErr[ipt] = NULL;
      if ( SysErr[ipt] )   delete [] SysErr[ipt];   SysErr[ipt] = NULL;
      if ( TotalErr[ipt] ) delete [] TotalErr[ipt]; TotalErr[ipt] = NULL;
      if ( XX[ipt] )       delete [] XX[ipt];       XX[ipt] = NULL;
   } 
      
   if ( secondary == "piplus" || secondary == "kplus" || secondary == "proton" )
   {
      NData[0] = 0;
   }
   else
   {
      NData[0] = 2;
   }
   
   for ( int ipt=0; ipt<NPt; ++ipt )
   {  
      if ( NData[ipt] > 0 ) 
      {
         Data[ipt] = new double[NData[ipt]];
         StatErr[ipt] = new double[NData[ipt]];
         SysErr[ipt] = new double[NData[ipt]];
	 TotalErr[ipt] = new double[NData[ipt]]; 
	 XX[ipt] = new double[NData[ipt]];
      }
   }
   
   if ( XX[0] )
   {
      XX[0][0] = 50.; XX[0][1] = 80.;
   }
   XX[1][0] = 30.; XX[1][1] = 40.; XX[1][2] = 50., XX[1][3] = 60.; XX[1][4] = 70.; XX[1][5] = 80.; XX[1][6] = 88.;
   XX[2][0] = 30.; XX[2][1] = 50.; XX[2][2] = 80.;
      
   if ( beam == "piplus" )
   {
      if ( secondary == "piplus" )
      {
         if ( target == "C" )
	 {
	    Data[1][0]=54.15;    Data[1][1]=41.08;    Data[1][2]=37.03;    Data[1][3]=32.15;    Data[1][4]=30.01;    Data[1][5]=28.17;    Data[1][6]=29.76; 
	    StatErr[1][0]=2.191; StatErr[1][1]=0.845; StatErr[1][2]=1.011; StatErr[1][3]=1.243; StatErr[1][4]=0.708; StatErr[1][5]=0.606; StatErr[1][6]=0.456;
	    Data[2][0]=30.10;    Data[2][1]=15.72;    Data[2][2]=11.11;
	    StatErr[2][0]=1.386; StatErr[2][1]=0.827; StatErr[2][2]=0.425;
	 }
	 else if ( target == "Cu" )
	 {
	    Data[1][0]=171.4;    Data[1][1]=127.2;    Data[1][2]=101.8;    Data[1][3]=91.68;    Data[1][4]=86.67;    Data[1][5]=75.84;    Data[1][6]=75.73;
	    StatErr[1][0]=4.746; StatErr[1][1]=3.453; StatErr[1][2]=2.990; StatErr[1][3]=3.265; StatErr[1][4]=2.347; StatErr[1][5]=1.798; StatErr[1][6]=1.480;
	    Data[2][0]=99.50;    Data[2][1]=46.44;    Data[2][2]=25.86;
	    StatErr[2][0]=4.416; StatErr[2][1]=2.919; StatErr[2][2]=1.081;
	 }
	 else if ( target == "Pb" )
	 {
	    Data[1][0]=350.3;    Data[1][1]=265.9;    Data[1][2]=184.9;    Data[1][3]=172.3;    Data[1][4]=172.5;    Data[1][5]=137.9;    Data[1][6]=119.0;
	    StatErr[1][0]=11.50; StatErr[1][1]=8.093; StatErr[1][2]=6.682; StatErr[1][3]=8.836; StatErr[1][4]=10.25; StatErr[1][4]=3.989; StatErr[1][5]=3.463;
	    Data[2][0]=191.9;    Data[2][1]=91.60;    Data[2][2]=47.24;
	    StatErr[2][0]=9.118; StatErr[2][1]=6.510; StatErr[2][2]=2.616;
	 }
      }
      else if ( secondary == "piminus" )
      {
         if ( target == "C" )
	 {
	    Data[0][0]=16.04;    Data[0][1]=10.93;
	    StatErr[0][0]=0.919; StatErr[0][1]=0.641;
	    Data[1][0]=28.12;    Data[1][1]=17.70;    Data[1][2]=8.071;    Data[1][3]=6.939;    Data[1][4]=6.282;    Data[1][5]=5.242;    Data[1][6]=4.035;
	    StatErr[1][0]=1.414; StatErr[1][1]=0.708; StatErr[1][2]=0.518; StatErr[1][3]=0.355; StatErr[1][4]=0.325; StatErr[1][5]=0.327; StatErr[1][6]=0.240;
	    Data[2][0]=16.41;    Data[2][1]=4.765;    Data[2][2]=1.085;
	    StatErr[2][0]=0.863; StatErr[2][1]=0.369; StatErr[2][2]=0.127;
	 }
	 else if ( target == "Cu" )
	 {
	    Data[0][0]=46.15;    Data[0][1]=32.91;
	    StatErr[0][0]=2.570; StatErr[0][1]=1.683;
	    Data[1][0]=99.06;    Data[1][1]=52.93;    Data[1][2]=21.14;    Data[1][3]=23.03;    Data[1][4]=21.35;    Data[1][5]=15.32;    Data[1][6]=11.18;
	    StatErr[1][0]=4.875; StatErr[1][1]=1.848; StatErr[1][2]=1.392; StatErr[1][3]=1.177; StatErr[1][4]=1.066; StatErr[1][5]=0.894; StatErr[1][6]=0.692;
	    Data[2][0]=48.47;    Data[2][1]=15.05;    Data[2][2]=3.088;
	    StatErr[2][0]=2.697; StatErr[2][1]=1.057; StatErr[2][2]=0.350;
	 }
	 else if ( target == "Pb" )
	 {
	    Data[0][0]=100.1;    Data[0][1]=62.87;
	    StatErr[0][0]=6.964; StatErr[0][1]=4.161;
	    Data[1][0]=208.9;    Data[1][1]=106.4;    Data[1][2]=44.18;    Data[1][3]=48.83;    Data[1][4]=40.42;    Data[1][5]=30.90;    Data[1][6]=16.45;
	    StatErr[1][0]=10.30; StatErr[1][1]=6.643; StatErr[1][2]=3.708; StatErr[1][3]=2.750; StatErr[1][4]=2.541; StatErr[1][5]=2.594; StatErr[1][6]=1.556;
	    Data[2][0]=98.59;    Data[2][1]=29.02;    Data[2][2]=5.038;
	    StatErr[2][0]=5.644; StatErr[2][1]=2.551; StatErr[2][2]=1.061;
	 }
      }
      else if ( secondary == "kplus" )
      {
         if ( target == "C" )
	 {
	    Data[1][0]=5.200;    Data[1][1]=3.624;    Data[1][2]=2.792;    Data[1][3]=2.052;    Data[1][4]=1.087;    Data[1][5]=1.453;    Data[1][6]=0.987;
	    StatErr[1][0]=0.868; StatErr[1][1]=0.304; StatErr[1][2]=0.309; StatErr[1][3]=0.303; StatErr[1][4]=0.185; StatErr[1][5]=0.144; StatErr[1][6]=0.087;
	    Data[2][0]=3.222;    Data[2][1]=1.615;    Data[2][2]=0.854;
	    StatErr[2][0]=0.563; StatErr[2][1]=0.303; StatErr[2][2]=0.115;
	 }
	 else if ( target == "Cu" )
	 {
	    Data[1][0]=18.23;    Data[1][1]=12.31;    Data[1][2]=7.853;    Data[1][3]=6.635;    Data[1][4]=5.049;    Data[1][5]=5.962;    Data[1][6]=3.244;
	    StatErr[1][0]=2.019; StatErr[1][1]=1.246; StatErr[1][2]=0.945; StatErr[1][3]=0.902; StatErr[1][4]=0.609; StatErr[1][5]=0.490; StatErr[1][6]=0.317;
	    Data[2][0]=14.66;    Data[2][1]=5.676;    Data[2][2]=2.979;
	    StatErr[2][0]=2.100; StatErr[2][1]=1.183; StatErr[2][2]=0.361;
	 }
	 else if ( target == "Pb" )
	 {
	    Data[1][0]=27.98;    Data[1][1]=23.08;    Data[1][2]=18.50;    Data[1][3]=11.82;    Data[1][4]=7.969;    Data[1][5]=8.074;    Data[1][6]=7.676;
	    StatErr[1][0]=4.870; StatErr[1][1]=2.999; StatErr[1][2]=2.355; StatErr[1][3]=2.007; StatErr[1][4]=2.612; StatErr[1][5]=1.045; StatErr[1][6]=0.884;
	    Data[2][0]=22.72;    Data[2][1]=13.09;    Data[2][2]=3.920;
	    StatErr[2][0]=4.232; StatErr[2][1]=2.925; StatErr[2][2]=0.852;
	 }
      }
      else if ( secondary == "kminus" )
      {
          if ( target == "C" )
	 {
	    Data[0][0]=1.391;    Data[0][1]=0.417;
	    StatErr[0][0]=0.369; StatErr[0][1]=0.109;
	    Data[1][0]=2.8858;   Data[1][1]=2.127;    Data[1][2]=1.139;    Data[1][3]=0.633;    Data[1][4]=0.674;    Data[1][5]=0.337;    Data[1][6]=0.073;
	    StatErr[1][0]=0.570; StatErr[1][1]=0.314; StatErr[1][2]=0.231; StatErr[1][3]=0.112; StatErr[1][4]=0.119; StatErr[1][5]=0.075; StatErr[1][6]=0.032;
	    Data[2][0]=2.505;    Data[2][1]=0.842;    Data[2][2]=0.1167;
	    StatErr[2][0]=0.429; StatErr[2][1]=0.178; StatErr[2][2]=0.058;
	 }
	 else if ( target == "Cu" )
	 {
	    Data[0][0]=3.714;    Data[0][1]=1.057;
	    StatErr[0][0]=1.024; StatErr[0][1]=0.292;
	    Data[1][0]=8.791;    Data[1][1]=6.773;    Data[1][2]=2.727;    Data[1][3]=2.732;    Data[1][4]=2.410;    Data[1][5]=0.761;    Data[1][6]=0.319;
	    StatErr[1][0]=1.840; StatErr[1][1]=0.820; StatErr[1][2]=0.604; StatErr[1][3]=0.434; StatErr[1][4]=0.355; StatErr[1][5]=0.181; Data[1][6]=0.113;
	    Data[2][0]=7.349;    Data[2][1]=2.524;    Data[2][2]=0.403;
	    StatErr[2][0]=1.327; StatErr[2][1]=0.503; StatErr[2][2]=0.188;
	 }
	 else if ( target == "Pb" )
	 {
	    Data[0][0]=11.02;    Data[0][1]=1.665;
	    StatErr[0][0]=3.457; StatErr[0][1]=0.812;
	    Data[1][0]=21.22;    Data[1][1]=11.40;    Data[1][2]=7.353;    Data[1][3]=5.521;    Data[1][4]=3.988;    Data[1][5]=2.886;    Data[1][6]=-1.;
	    StatErr[1][0]=4.196; StatErr[1][1]=2.967; StatErr[1][2]=1.863; StatErr[1][3]=1.004; StatErr[1][4]=0.793; StatErr[1][5]=0.573; StatErr[1][6]=-1.;
	    Data[2][0]=14.62;    Data[2][1]=4.872;    Data[2][2]=1.761;
	    StatErr[2][0]=3.028; StatErr[2][1]=1.249; StatErr[2][2]=0.572;
	 }
     }
     else if ( secondary == "proton" )
     {
        if ( target == "C" )
	{
	   Data[1][0]=6.540;    Data[1][1]=4.340;    Data[1][2]=3.005;    Data[1][3]=1.540;    Data[1][4]=0.829;    Data[1][5]=0.386;    Data[1][6]=0.827;
	   StatErr[1][0]=0.757; StatErr[1][1]=0.269; StatErr[1][2]=0.269; StatErr[1][3]=0.228; StatErr[1][4]=0.109; StatErr[1][5]=0.078; StatErr[1][6]=0.074;
	   Data[2][0]=3.708;    Data[2][1]=2.074;    Data[2][2]=0.234;
	   StatErr[2][0]=0.507; StatErr[2][1]=0.280; StatErr[2][2]=0.060;
	}
	else if ( target == "Cu" )
	{
	   Data[1][0]=21.85;    Data[1][1]=14.40;    Data[1][2]=7.397;    Data[1][3]=5.097;    Data[1][4]=2.380;    Data[1][5]=1.095;    Data[1][6]=1.241;
	   StatErr[1][0]=1.707; StatErr[1][1]=1.133; StatErr[1][2]=0.784; StatErr[1][3]=0.700; StatErr[1][4]=0.350; StatErr[1][5]=0.243; StatErr[1][6]=0.210;
	   Data[2][0]=13.40;    Data[2][1]=5.476;    Data[2][2]=0.589;
	   StatErr[2][0]=1.673; StatErr[2][1]=0.970; StatErr[2][2]=0.166;
	}
	else if ( target == "Pb" )
	{
	   Data[1][0]=44.45;    Data[1][1]=27.26;    Data[1][2]=16.08;    Data[1][3]=7.424;    Data[1][4]=4.627;    Data[1][5]=2.758;    Data[1][6]=0.559;
	   StatErr[1][0]=4.355; StatErr[1][1]=2.693; StatErr[1][2]=1.929; StatErr[1][3]=1.448; StatErr[1][4]=1.534; StatErr[1][5]=0.638; StatErr[1][6]=0.643;
	   Data[2][0]=31.27;    Data[2][1]=13.94;    Data[2][2]=0.881;
	   StatErr[2][0]=4.203; StatErr[2][1]=2.374; StatErr[2][2]=0.437;
	}
     }
     else if ( secondary == "antiproton" )
     {
        if ( target == "C" )
	{
	   Data[0][0]=1.281;    Data[0][1]=0.128;
	   StatErr[0][0]=0.263; StatErr[0][1]=0.079;
	   Data[1][0]=2.812;    Data[1][1]=0.868;    Data[1][2]=0.474;    Data[1][3]=0.171;    Data[1][4]=0.658;    Data[1][5]=0.062;    Data[1][6]=0.019;
	   StatErr[1][0]=0.433; StatErr[1][1]=0.191; StatErr[1][2]=0.166; StatErr[1][3]=0.114; StatErr[1][4]=0.098; StatErr[1][5]=0.034; StatErr[1][6]=0.024;
	   Data[2][0]=1.719;    Data[2][1]=0.272;    Data[2][2]=0.227;
	   StatErr[2][0]=0.281; StatErr[2][1]=0.083; StatErr[2][2]=0.115;
	}
	else if ( target == "Cu" )
	{
	   Data[0][1]=2.774;    Data[0][1]=0.138;
	   StatErr[0][0]=0.696; StatErr[0][1]=0.202;
	   Data[1][0]=7.883;    Data[1][1]=3.239;    Data[1][2]=1.240;    Data[1][3]=0.648;    Data[1][4]=0.873;    Data[1][5]=0.557;    Data[1][6]=0.026;
	   StatErr[1][0]=1.353; StatErr[1][1]=0.527; StatErr[1][2]=0.481; StatErr[1][3]=0.378; StatErr[1][4]=0.361; StatErr[1][5]=0.207; StatErr[1][6]=0.048;
	   Data[2][0]=5.289;    Data[2][1]=1.267;    Data[2][2]=0.288;
	   StatErr[2][0]=0.896; StatErr[2][1]=0.290; StatErr[2][2]=0.294;
	}
	else if ( target == "Pb" )
	{
	   Data[0][0]=5.401;    Data[0][1]=1.016;
	   StatErr[0][0]=2.169; StatErr[0][1]=0.773;
	   Data[1][0]=23.00;    Data[1][1]=9.992;    Data[1][2]=4.469;    Data[1][3]=1.384;    Data[1][4]=3.048;    Data[1][5]=1.004;    Data[1][6]=1.224;
	   StatErr[1][0]=3.351; StatErr[1][1]=2.304; StatErr[1][2]=2.183; StatErr[1][3]=1.375; StatErr[1][4]=1.034; StatErr[1][5]=0.538; StatErr[1][6]=0.261;
	   Data[2][0]=14.52;    Data[2][1]=1.738;    Data[2][2]=0.708;
	   StatErr[2][0]=2.240; StatErr[2][1]=0.616; StatErr[2][2]=1.956;
	}
     }
   }
   else if ( beam == "kplus" )
   {
     if ( secondary == "kplus" )
     {
        if ( target == "C" )
	{
	}
	else if ( target == "Cu" )
	{
	}
	else if ( target == "Pb" )
	{
	}
     }
     else if ( secondary == "kminus" )
     {
        if ( target == "C" )
	{
	}
	else if ( target == "Cu" )
	{
	}
	else if ( target == "Pb" )
	{
	}
     }
     else if ( secondary == "proton" )
     {
        if ( target == "C" )
	{
	}
	else if ( target == "Cu" )
	{
	}
	else if ( target == "Pb" )
	{
	}
     }
     else if ( secondary == "antiproton" )
     {
        if ( target == "C" )
	{
	}
	else if ( target == "Cu" )
	{
	}
	else if ( target == "Pb" )
	{
	}
     }      
   } 
   else if ( beam == "proton" )
   {
     if ( secondary == "piplus" )
     {
        if ( target == "C" )
	{
	   Data[1][0]=36.70;    Data[1][1]=23.21;    Data[1][2]=11.91;    Data[1][3]=5.506;    Data[1][4]=2.082;    Data[1][5]=0.740;    Data[1][6]=0.332;
	   StatErr[1][0]=2.188; StatErr[1][1]=0.781; StatErr[1][2]=0.691; StatErr[1][3]=0.576; StatErr[1][4]=0.219; StatErr[1][5]=0.131; StatErr[1][6]=0.068;
	   Data[2][0]=18.45;    Data[2][1]=5.564;    Data[2][2]=0.175;
	   StatErr[2][0]=1.380; StatErr[2][1]=0.634; StatErr[2][2]=0.075;
	}
	else if ( target == "Cu" )
	{
	   Data[1][0]=120.7;    Data[1][1]=58.08;    Data[1][2]=30.97;    Data[1][3]=14.55;    Data[1][4]=5.227;    Data[1][5]=1.166;    Data[1][6]=0.719;
	   StatErr[1][0]=5.003; StatErr[1][1]=2.982; StatErr[1][2]=2.043; StatErr[1][3]=1.621; StatErr[1][4]=0.648; StatErr[1][5]=0.340; StatErr[1][6]=0.224;
	   Data[2][0]=51.59;    Data[2][1]=16.50;    Data[2][2]=0.846;
	   StatErr[2][0]=4.138; StatErr[2][1]=2.250; StatErr[2][2]=0.276; 
	}
	else if ( target == "Pb" )
	{
	   Data[1][0]=195.9;    Data[1][1]=113.3;    Data[1][2]=54.70;    Data[1][3]=28.31;    Data[1][4]=11.81;    Data[1][5]=3.689;    Data[1][6]=2.104;
	   StatErr[1][0]=11.53; StatErr[1][1]=6.957; StatErr[1][2]=4.948; StatErr[1][3]=3.761; StatErr[1][4]=2.861; StatErr[1][5]=0.825; StatErr[1][6]=0.610;
	   Data[2][0]=108.3;    Data[2][1]=20.81;    Data[2][2]=1.530;
	   StatErr[2][0]=9.629; StatErr[2][1]=4.655; StatErr[2][2]=0.679;
	}
     }
     else if ( secondary == "piminus" )
     {
        if ( target == "C" )
	{
	   Data[0][0]=7.374;    Data[0][1]=0.181;
	   StatErr[0][0]=0.726; StatErr[0][1]=0.080;
	   Data[1][0]=20.13;    Data[1][1]=9.396;    Data[1][2]=2.602;    Data[1][3]=1.506;    Data[1][4]=0.333;    Data[1][5]=0.094;    Data[1][6]=0.031;
	   StatErr[1][0]=1.502; StatErr[1][1]=0.659; StatErr[1][2]=0.361; StatErr[1][3]=0.198; StatErr[1][4]=0.080; StatErr[1][5]=0.067; StatErr[1][6]=0.021;
	   Data[2][0]=11.45;    Data[2][1]=2.447;    Data[2][2]=0.092;
	   StatErr[2][0]=0.920; StatErr[2][1]=0.316; StatErr[2][2]=0.048;
	}
	else if ( target == "Cu" )
	{
	   Data[0][0]=18.98;    Data[0][1]=0.236;
	   StatErr[0][0]=2.036; StatErr[0][1]=0.216;
	   Data[1][0]=50.99;    Data[1][1]=24.85;    Data[1][2]=9.145;    Data[1][3]=4.070;    Data[1][4]=0.936;    Data[1][5]=0.289;    Data[1][6]=0.114;
	   StatErr[1][0]=4.265; StatErr[1][1]=1.666; StatErr[1][2]=1.100; StatErr[1][3]=0.630; StatErr[1][4]=0.253; StatErr[1][5]=0.138; StatErr[1][6]=0.085;
	   Data[2][0]=30.38;    Data[2][1]=5.299;    Data[2][2]=0.154;
	   StatErr[2][0]=2.649; StatErr[2][1]=0.762; StatErr[2][2]=0.080;
	}
	else if ( target == "Pb" )
	{
	   Data[0][0]=28.48;    Data[0][1]=0.554;
	   StatErr[0][0]=5.053; StatErr[0][1]=0.654;
	   Data[1][0]=124.7;    Data[1][1]=37.66;    Data[1][2]=10.55;    Data[1][3]=5.755;    Data[1][4]=1.273;    Data[1][5]=0.377;    Data[1][6]=0.304;
	   StatErr[1][0]=9.885; StatErr[1][1]=5.761; StatErr[1][2]=2.743; StatErr[1][3]=1.435; StatErr[1][4]=0.602; StatErr[1][5]=0.396; StatErr[1][6]=0.209;
	   Data[2][0]=53.21;    Data[2][1]=8.025;    Data[2][2]=0.363;
	   StatErr[2][0]=5.132; StatErr[2][1]=1.613; StatErr[2][2]=0.279;
	}
     }
     else if ( secondary == "kplus" )
     {
        if ( target == "C" )
	{
	   Data[1][0]=5.279;    Data[1][1]=3.386;    Data[1][2]=1.650;    Data[1][3]=1.360;    Data[1][4]=0.592;    Data[1][5]=0.160;    Data[1][6]=0.085;
	   StatErr[1][0]=0.978; StatErr[1][1]=0.346; StatErr[1][2]=0.289; StatErr[1][3]=0.292; StatErr[1][4]=0.122; StatErr[1][5]=0.063; StatErr[1][6]=0.028;
	   Data[2][0]=3.230;    Data[2][1]=1.206;    Data[2][2]=0.158;
	   StatErr[2][0]=0.729; StatErr[2][1]=0.354; StatErr[2][2]=0.095;
	}
	else if ( target == "Cu" )
	{
	   Data[1][0]=14.79;    Data[1][1]=9.613;    Data[1][2]=4.575;    Data[1][3]=1.905;    Data[1][4]=1.008;    Data[1][5]=0.216;    Data[1][6]=0.135;
	   StatErr[1][0]=2.097; StatErr[1][1]=1.394; StatErr[1][2]=0.858; StatErr[1][3]=0.682; StatErr[1][4]=0.311; StatErr[1][5]=0.162; StatErr[1][6]=0.050;
	   Data[2][0]=7.315;    Data[2][1]=1.925;    Data[2][2]=0.360;
	   StatErr[2][0]=1.991; StatErr[2][1]=1.286; StatErr[2][2]=0.168;
	}
	else if ( target == "Pb" )
	{
	   Data[1][0]=37.32;    Data[1][1]=21.47;    Data[1][2]=8.101;    Data[1][3]=4.300;    Data[1][4]=1.397;    Data[1][5]=1.006;    Data[1][6]=0.329;
	   StatErr[1][0]=5.564; StatErr[1][1]=3.267; StatErr[1][2]=2.040; StatErr[1][3]=1.747; StatErr[1][4]=1.199; StatErr[1][5]=0.487; StatErr[1][6]=0.394;
	   Data[2][0]=20.67;    Data[2][1]=7.836;    Data[2][2]=0.432;
	   StatErr[2][0]=5.743; StatErr[2][1]=3.165; StatErr[2][2]=0.426; 
	}
     }
     else if ( secondary == "kminus" )
     {
        if ( target == "C" )
	{
	   Data[0][0]=0.485;    Data[0][1]=-1.;
	   StatErr[0][0]=0.212; StatErr[0][1]=-1.;
	   Data[1][0]=1.144;    Data[1][1]=0.715;    Data[1][2]=0.241;    Data[1][3]=0.157;    Data[1][4]=-1.,    Data[1][5]=-1.;    Data[1][6]=-1.;
	   StatErr[1][0]=0.563; StatErr[1][1]=0.212; StatErr[1][2]=0.124; StatErr[1][3]=0.072; StatErr[1][4]=-1.; StatErr[1][5]=-1.; StatErr[1][6]=-1.;
	   Data[2][0]=1.033;    Data[2][1]=0.073;    Data[2][2]=-1.;
	   StatErr[2][0]=1.105; StatErr[2][1]=0.149; StatErr[2][2]=-1.; 
	}
	else if ( target == "Cu" )
	{
	   Data[0][0]=1.236;    Data[0][1]=-1.;
	   StatErr[0][0]=0.597; StatErr[0][1]=-1.;
	   Data[1][0]=2.356;    Data[1][1]=2.353;    Data[1][2]=0.639;    Data[1][3]=0.142;    Data[1][4]=-1.;    Data[1][5]=-1.;    Data[1][6]=-1.;
	   StatErr[1][0]=1.599; StatErr[1][1]=0.572; StatErr[1][2]=0.341; StatErr[1][3]=0.254; StatErr[1][4]=-1.; StatErr[1][5]=-1.; StatErr[1][6]=-1.;
	   Data[2][0]=2.254;    Data[2][1]=0.439;    Data[2][2]=-1.;
	   StatErr[2][0]=1.105; StatErr[2][1]=0.502; StatErr[2][2]=-1.;
	}
	else if ( target == "Pb" )
	{
	   Data[0][0]=4.208;    Data[0][1]=-1.;
	   StatErr[0][0]=2.468; StatErr[0][1]=-1.;
	   Data[1][0]=5.944;    Data[1][1]=2.850;    Data[1][2]=1.084;    Data[1][3]=0.432;    Data[1][4]=0.204;    Data[1][5]=-1.;    Data[1][6]=-1.;
	   StatErr[1][0]=5.669; StatErr[1][1]=1.954; StatErr[1][2]=1.015; StatErr[1][3]=0.426; StatErr[1][4]=0.353; StatErr[1][5]=-1.; StatErr[1][6]=-1.;
	   Data[2][0]=5.215;    Data[2][1]=0.988;    Data[2][2]=-1.;
	   StatErr[2][0]=3.324; StatErr[2][1]=0.974; StatErr[2][2]=-1.;
	}
     }
     else if ( secondary == "proton" )
     {
        if ( target == "C" )
	{
	   Data[1][0]=42.25;    Data[1][1]=49.24;    Data[1][2]=49.87;    Data[1][3]=54.44;    Data[1][4]=60.95;    Data[1][5]=58.83;    Data[1][6]=61.66;
	   StatErr[1][0]=2.391; StatErr[1][1]=1.184; StatErr[1][2]=1.463; StatErr[1][3]=2.305; StatErr[1][4]=1.412; StatErr[1][5]=1.246; StatErr[1][6]=0.952;
	   Data[2][0]=25.74;    Data[2][1]=29.59;    Data[2][2]=24.14;
	   StatErr[2][0]=1.646; StatErr[2][1]=1.591; StatErr[2][2]=0.997; 
	}
	else if ( target == "Cu" )
	{
	   Data[1][0]=127.4;    Data[1][1]=133.2;    Data[1][2]=131.8;    Data[1][3]=129.0;    Data[1][4]=132.3;    Data[1][5]=137.2;    Data[1][6]=137.6;
	   StatErr[1][0]=5.190; StatErr[1][1]=4.615; StatErr[1][2]=4.267; StatErr[1][3]=5.352; StatErr[1][4]=3.927; StatErr[1][5]=3.588; StatErr[1][6]=2.879;
	   Data[2][0]=74.63;    Data[2][1]=81.22;    Data[2][2]=53.94;
	   StatErr[2][0]=4.977; StatErr[2][1]=5.056; StatErr[2][2]=2.541;
	}
	else if ( target == "Pb" )
	{
	   Data[1][0]=242.1;    Data[1][1]=247.5;    Data[1][2]=228.2;    Data[1][3]=224.7;    Data[1][4]=249.2;    Data[1][5]=247.7;    Data[1][6]=217.7;
	   StatErr[1][0]=12.56; StatErr[1][1]=10.53; StatErr[1][2]=9.564; StatErr[1][3]=10.70; StatErr[1][4]=16.83; StatErr[1][5]=7.393; StatErr[1][6]=6.767;
	   Data[2][0]=138.4;    Data[2][1]=115.1;    Data[2][2]=89.44;
	   StatErr[2][0]=10.90; StatErr[2][1]=9.570; StatErr[2][2]=5.390;
	}
     }
     else if ( secondary == "antiproton" )
     {
        if ( target == "C" )
	{
	   Data[0][0]=0.079;    Data[0][1]=0.0;
	   StatErr[0][0]=0.098; StatErr[0][1]=0.043;
	   Data[1][0]=1.022;    Data[1][1]=0.122;    Data[1][2]=0.113;    Data[1][3]=0.018;    Data[1][4]=0.;       Data[1][5]=0.;        Data[1][6]=0.;
	   StatErr[1][0]=0.352; StatErr[1][1]=0.100; StatErr[1][2]=0.086; StatErr[1][3]=0.034; StatErr[1][4]=0.016; StatErr[1][5]=0.027; StatErr[1][6]=0.012;
	   Data[2][0]=0.098;    Data[2][1]=-1.;    Data[2][2]=0.;
	   StatErr[2][0]=0.161; StatErr[2][1]=-1.; StatErr[2][2]=0.025;
 	}
	else if ( target == "Cu" )
	{
	   Data[0][0]=0.352;    Data[0][1]=0.;
	   StatErr[0][0]=0.305; StatErr[0][1]=0.133;
	   Data[1][0]=1.566;    Data[1][1]=0.482;    Data[1][2]=0.301;    Data[1][3]=0.036;    Data[1][4]=0.;       Data[1][5]=0.;       Data[1][6]=0.;
	   StatErr[1][0]=0.856; StatErr[1][1]=0.296; StatErr[1][2]=0.342; StatErr[1][3]=0.154; StatErr[1][4]=0.048; StatErr[1][5]=0.083; StatErr[1][6]=0.036;
	   Data[2][0]=1.870;    Data[2][1]=-1.;    Data[2][2]=0.;
	   StatErr[2][0]=0.775; StatErr[2][1]=-1.; StatErr[2][2]=0.077; 
	}
	else if ( target == "Pb" )
	{
	   Data[0][0]=0.534;    Data[0][1]=0.;
	   StatErr[0][0]=1.759; StatErr[0][1]=0.587;
	   Data[1][0]=8.948;    Data[1][1]=0.342;    Data[1][2]=0.241;    Data[1][3]=0.492;    Data[1][4]=-1.;    Data[1][5]=0.;       Data[1][6]=0.;
	   StatErr[1][0]=3.055; StatErr[1][1]=1.177; StatErr[1][2]=1.155; StatErr[1][3]=0.721; StatErr[1][4]=-1.; StatErr[1][5]=0.372; StatErr[1][6]=0.161;
	   Data[2][0]=0.547;    Data[2][1]=-1.;    Data[2][2]=0.;
	   StatErr[2][0]=3.584; StatErr[2][1]=-1.; StatErr[2][2]=0.343;
	}
     }      
   }
   
   for ( int ipt=0; ipt<NPt; ++ipt )
   {
      for ( int id=0; id<NData[ipt]; ++id )
      {
         // systematic uncertainty is 10% (due to normalization)
	 //
	 SysErr[ipt][id] = Data[ipt][id] * 0.1;
	 TotalErr[ipt][id] = std::sqrt( StatErr[ipt][id]*StatErr[ipt][id] + SysErr[ipt][id]*SysErr[ipt][id] );
      }
   }
         
   return;

}

void SASM6E( std::string beam, std::string target, std::string energy,
             std::string secondary )
{

   
   GetExData( beam, target, secondary );
   
   for ( int ip=0; ip<NPt; ++ip )
   {
      if ( gr[ip] ) delete gr[ip];
      gr[ip] = 0;
      if ( NData[ip] > 0 )
      {      
         gr[ip] = new TGraphErrors( NData[ip], XX[ip], Data[ip], 0, TotalErr[ip] );
         gr[ip]->SetMarkerColor(kBlue);
         gr[ip]->SetMarkerStyle(22);
         gr[ip]->SetMarkerSize(1.5);   
      }      
   }
   
   chi2M = 0.;
   NDFM = 0;
   for ( int i=0; i<NV; ++i )
   {
      chi2V[i] = 0.;
      NDFV[i] = 0;
   }
   
   TCanvas* cnv = new TCanvas( "cnv", "", 1100, 600 );
   TPad*    pdm = new TPad( "pdm", "", 0.01, 0.20, 0.99, 0.95 ); 

   std::string txt1 = energy + " GeV/c " + beam + " + " + target + " #rightarrow " + secondary;
   txt1 += " + X";
   TLatex* ltxt1 = new TLatex( 0.35, 0.96, txt1.c_str() );
   ltxt1->SetTextSize(0.035);
   cnv->cd();
   ltxt1->Draw();
   cnv->cd();
   pdm->Draw();
   pdm->Divide(3.,1.,0.,0.);
   
   TH1F* hm[NM];
   
   double ymin = 10000.;
   double ymax = -1.;
      
   for ( int im=0; im<NM; ++im )
   {
   
      std::string hfile = "sasm6e-histo/" + beam + target + energy + "GeV" + model[im] + ".root";
      TFile* f = new TFile( hfile.c_str() );

      TCanvas* myc = new TCanvas( "myc", "", 1100, 600 ); 
      TPad*    pad = new TPad("pad", "", 0.01, 0.20, 0.99, 0.95 );
      myc->cd();
      ltxt1->Draw();
      std::string mtxt = ", model: " + model[im];
      TText* txt2 = new TText( 0.655, 0.96, mtxt.c_str() );
      txt2->SetTextSize(0.035);
      txt2->Draw();
      myc->cd();
      pad->Draw();
      pad->Divide(3.,1.,0.,0.);
   
      chi2M = 0.;
      NDFM = 0;
      for ( int iv=0; iv<NV; ++iv )
      {
         chi2V[iv] = 0.;
         NDFV[iv] = 0;
      }
      
      for ( int ipt=0; ipt<NPt; ++ipt )
      {  
         std::string hname = secondary + "_" + hlabel[ipt]; 
         hm[im] = (TH1F*)f->Get( hname.c_str() );
         if ( hm[im] == NULL )
         {
            std::cout << " NULL histo " << hname << " from file " << hfile << std::endl;
	    continue;
         }      
         hm[im]->SetStats(0);
         hm[im]->SetLineColor(colorM[im]);
         hm[im]->SetLineWidth(2);
         int nx = hm[im]->GetNbinsX();
         for (int k=1; k <= nx; k++) 
         {
	    double yy  = hm[im]->GetBinContent(k);
	    double eyy = hm[im]->GetBinError(k);
	    if ( (yy+eyy) > ymax ) ymax = yy+eyy;
	    if ( (yy-eyy) < ymin && (yy-eyy) > 0. ) ymin = yy-eyy;
         }
	 cnv->cd();
	 pdm->cd(ipt+1);
         gPad->SetRightMargin(0.1);
         gPad->SetLeftMargin(0.15);
         ymin = std::min( 0., ymin );
	 for ( int j=0; j<=im; ++j )
	 {
	    hm[j]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); 
	 }
	 gPad->Update();
         if ( im == 0 )
         {
	    hm[im]->SetTitle( htitle[ipt].c_str() );
	    hm[im]->GetXaxis()->SetTitle("momentum_{sec} (GeV/c)");
	    hm[im]->GetXaxis()->SetTitleSize(0.05);
	    hm[im]->GetXaxis()->SetTitleFont(62);
	    hm[im]->GetXaxis()->SetLabelFont(62);
	    hm[im]->GetXaxis()->CenterTitle();
	    hm[im]->GetYaxis()->SetTitle("E d^{3}#sigma / dp^{3} [mb/(GeV/c)^{2}]");
	    hm[im]->GetYaxis()->SetTitleOffset(1.0);
	    hm[im]->GetYaxis()->SetTitleSize(0.05);
	    hm[im]->GetYaxis()->SetTitleFont(62);
	    hm[im]->GetYaxis()->SetLabelFont(62);
	    hm[im]->GetYaxis()->CenterTitle();
            hm[im]->Draw("histE1");
         }
         else hm[im]->Draw("histE1same");
//	 gPad->Update();
//	 cnv->Update();
	 if ( gr[ipt] ) chi2M += Chi2( gr[ipt], hm[im], NDFM ); 	 
	 //
	 //
         myc->cd();
         pad->cd(ipt+1);
         gPad->SetRightMargin(0.1);
         gPad->SetLeftMargin(0.15);
         DrawSASM6ERegre( beam, target, energy, secondary, model[im], ipt, true, chi2V, NDFV );
      }
      
      // std::cout << " for " << model[im] << " chi2/NDF = " << chi2M << "/" << NDFM << std::endl;
      std::ostringstream os;
      os << chi2M/(double)NDFM;
      std::string chi2txt = model[im] + ": chi2/NDF = " + os.str();
      TText* ttxt = new TText( 0.5, 0.15-0.035*im, chi2txt.c_str() );
      ttxt->SetTextFont(62);
      ttxt->SetTextSize(0.035);
      cnv->cd();
      ttxt->Draw();

      TLegend* leg = new TLegend(0.05,0.01,0.34,0.19);
      TLegendEntry* entry = 0;
      for ( int m=0; m<NV; ++m )
      {
         std::string info = model[im] +", " + version[m];
         entry = leg->AddEntry( "", info.c_str(),"p" );
	 entry->SetMarkerStyle(kFullCross); // font=34
	 entry->SetMarkerSize(1.5);
	 entry->SetMarkerColor( colorV[m] );
         entry->SetTextFont(62);
      }
      TLegendEntry* dentry = leg->AddEntry("","exp.data","p");
      dentry->SetMarkerStyle(22);
      dentry->SetMarkerColor(kBlue);
      dentry->SetMarkerSize(1.5);
      dentry->SetTextFont(62);
      leg->SetFillColor(kWhite);
      myc->cd();
      leg->Draw();
      myc->Update();
   
      for ( int iv=0; iv<NV; ++iv )
      {
         //std::cout << " chi2=" << chi2V[iv] << " ndf=" << NDFV[iv] << " for version " << version[iv] << std::endl; 
         std::ostringstream os1;
         os1 << chi2V[iv]/(double)NDFV[iv];
         std::string chi2txt1 = version[iv] + ": chi2/NDF = " + os1.str();
	 //std::cout << chi2txt1 << std::endl;
         TText* ttxt1 = new TText( 0.5, 0.15-0.035*iv, chi2txt1.c_str() );
         ttxt1->SetTextFont(62);
         ttxt1->SetTextSize(0.035);
         myc->cd();
         ttxt1->Draw();
	 myc->Update();
      }

      std::string output = beam + "-" + target + "-" + energy + "GeV-" + secondary + "-" ;
      output += model[im];
      output += "-regre.gif";
      myc->Print( output.c_str() );
      
      
      delete leg;
      delete pad;
      delete myc;
   }

   
   if ( beam == "kplus" && ( secondary == "piplus" || secondary == "piminus " ) ) return; // no data for kaon + A --> pion


   cnv->cd();
   for ( int ip=0; ip<NPt; ++ip )
   {
      pdm->cd(ip+1);
      if ( NData[ip] > 0 )
      {      
         gr[ip]->Draw("p");     
      }
      else
      {
         TText* txt = new TText( 60., ymax*0.8, "No data" );
	 txt->SetTextFont(62);
	 txt->SetTextSize(0.05);
	 txt->Draw();
      }
   }
 
   TLegend* leg1 = new TLegend(0.05,0.01,0.34,0.19);
   TLegendEntry* entry1 = 0;
   for ( int m=0; m<NM; ++m )
   {
         // std::string info = model[m] + ", geant4-10-02-ref-10";
         std::string info = model[m] + ", " + version[0];
	 entry1 = leg1->AddEntry( hm[m], info.c_str(),"L" );
         entry1->SetLineColor( colorM[m] );
         entry1->SetLineWidth(3);
         entry1->SetTextFont(62);
   }
   entry1 = leg1->AddEntry("","exp.data","p");
   entry1->SetMarkerStyle(22);
   entry1->SetMarkerColor(kBlue);
   entry1->SetMarkerSize(1.5);
   entry1->SetTextFont(62);
   leg1->SetFillColor(kWhite);
   cnv->cd();
   leg1->Draw();
       
   std::string output1 = beam + "-" + target + "-" + energy + "GeV-" + secondary + "-models.gif";
   cnv->Print( output1.c_str() );
   
   return;

}

void DrawSASM6ERegre( std::string beam, std::string target, std::string energy, std::string secondary, 
                      std::string model, 
		      int ipt, 
		      bool incdata,
		      double* chi2, int* ndf )
{

   if ( NV <= 0 ) return;
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   // this assumes that GetExData has been called already...
   //
   for ( int i=0; i<NData[ipt]; ++i )
   {
      if ( (Data[ipt][i]+TotalErr[ipt][i]) > ymax ) ymax = Data[ipt][i]+TotalErr[ipt][i];
      if ( (Data[ipt][i]-TotalErr[ipt][i]) < ymin && (Data[ipt][i]-TotalErr[ipt][i]) > 0. ) ymin = (Data[ipt][i]-TotalErr[ipt][i]);
   }
   
   TH1F* hi[NV]; // ftf for 10.2.ref07, etc.
   
   std::string location = ".";
   for ( int iv=0; iv<NV; ++iv )
   {

      if ( iv > 0 )
      {
         location = regre_test_dir + "/" + TEST_NAME + "/" + version[iv];
      }

      std::string histofile = location + "/sasm6e-histo/" + beam + target + energy + "GeV" + model + ".root";
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + hlabel[ipt]; 
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );

      if ( hi[iv] == NULL )
      {
         std::cout << " NULL histo " << histoname << " from file " << histofile << std::endl;
	 continue;
      }
      
      hi[iv]->SetStats(0);
      hi[iv]->SetLineColor(colorV[iv]);
      hi[iv]->SetLineWidth(2);

      int nx = hi[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy  = hi[iv]->GetBinContent(k);
	double eyy = hi[iv]->GetBinError(k);
	if ( (yy+eyy) > ymax ) ymax = yy+eyy;
	if ( (yy-eyy) < ymin && (yy-eyy) > 0. ) ymin = yy-eyy;
      }
      
      if ( iv == 0 )
      {
	 hi[iv]->SetTitle( htitle[ipt].c_str() );
	 hi[iv]->GetXaxis()->SetTitle("momentum_{sec} (GeV/c)");
	 hi[iv]->GetXaxis()->SetTitleSize(0.05);
	 hi[iv]->GetXaxis()->SetTitleFont(62);
	 hi[iv]->GetXaxis()->SetLabelFont(62);
	 hi[iv]->GetXaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetTitle("E d^{3}#sigma / dp^{3} [mb/(GeV/c)^{2}]");
	 hi[iv]->GetYaxis()->SetTitleOffset(1.0);
	 hi[iv]->GetYaxis()->SetTitleSize(0.05);
	 hi[iv]->GetYaxis()->SetTitleFont(62);
	 hi[iv]->GetYaxis()->SetLabelFont(62);
	 hi[iv]->GetYaxis()->CenterTitle();
         hi[iv]->Draw("histE1");
      }
      else hi[iv]->Draw("histE1same"); 
      
   }
   
   ymin = std::min( 0., ymin );

   for ( int iv=0; iv<NV; ++iv )
   {
      hi[iv]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
   }
   
   if ( beam == "kplus" && ( secondary == "piplus" || secondary == "piminus " ) ) return; // no data for kaon + A --> pion

   if ( incdata )
   {
      if ( NData[ipt] > 0 )
      {      
         gr[ipt]->Draw("p");  
	 for ( int iv=0; iv<NV; ++iv )
	 {   
	    chi2[iv] += Chi2( gr[ipt], hi[iv], ndf[iv] ); 
	 }	 
      }
      else
      {
         TText* txt = new TText( 60., ymax*0.8, "No data" );
	 txt->SetTextFont(62);
	 txt->SetTextSize(0.05);
	 txt->Draw();
      }
   }
     
   return;

}

void DumpJSONMaster( std::string beam, std::string target, int tglnk )
{

   std::string fjsonname = "SASM6E-spectra-" + beam + "-" + target + ".json";
   
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;

   DumpJSON( jsonfile, beam, target, tglnk, "piplus" );
   DumpJSON( jsonfile, beam, target, tglnk, "piminus" );
   DumpJSON( jsonfile, beam, target, tglnk, "kplus" );
   DumpJSON( jsonfile, beam, target, tglnk, "kminus" );
   DumpJSON( jsonfile, beam, target, tglnk, "proton" );
   DumpJSON( jsonfile, beam, target, tglnk, "antiproton" );

   jsonfile << "]}" << std::endl;
   
   return;
  
}

void DumpJSON( std::ofstream& jsonfile, std::string beam, std::string target, int tglnk, std::string secondary )
{

   GetExData( beam, target, secondary );
   
   int reflnk = 56; // Barton paper

   // FIXME !!!
   //
   int bmlnk = -1;
   if ( beam == "piplus" )
   {
      bmlnk = 10;
   }
   else if ( beam == "kplus" )
   {
      bmlnk = 11;
   }
   else if ( beam == "proton" )
   {
      bmlnk = 12;
   }
   
   std::string bm = beam;
   if ( beam == "piplus" )
   {
      bm = "pi+";
   }
   else if ( beam == "kplus" )
   {
      bm = "K+";
   }
   
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
   }
   else if ( secondary == "antiproton" )
   {
      seclink = -2212;
   }
   else if ( secondary == "neutron" )
   {
      seclink = 2112;
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

   for ( int ipt=0; ipt<NPt; ++ipt )
   {
      
      if ( NData[ipt] <= 0 ) continue;
      
      jsonfile << "{" << std::endl;
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":" << reflnk <<",\"mcdetaillnk\":1," << std::endl;

      jsonfile << "\"beamlnk\":" << bmlnk << ",\"targetlnk\":" << tglnk 
               << ",\"observablelnk\":9,\"secondarylnk\":"  << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":1000," << std::endl;

      jsonfile << "\"title\":\"Production of " << sec << " in " << bm << " on " << target << " interactions at 100GeV/c, " 
               << "pt=" << PT[ipt] << " [GeV/c]\"," << std::endl;

      jsonfile << "\"npoints\":" << NData[ipt] << "," << std::endl;
      jsonfile << "\"nbins\":[]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"E(d3sigma/dp3) [mb/(GeV/c)**2]\"]," << std::endl;
      
      DumpSpectrum( jsonfile, ipt );
      
      jsonfile << "}," << std::endl; // end of datatable 

      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"transerse mometum of secondary particle [GeV/c]\"], \"parvalues\":[\"" << PT[ipt] << "\"]" << std::endl;

      jsonfile << "}," << std::endl; 
      
   }
  
   return;

}

void DumpSpectrum( std::ofstream& jsonfile, int ipt )
{

   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << XX[ipt][0];
   for ( int i=1; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << XX[ipt][i];
   }
   for ( int i=0; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << Data[ipt][i];
   }
   jsonfile << "]," << std::endl;

   // stat errors "up/down"
   //
   
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NData[ipt]; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << StatErr[ipt][i] * 0.5;
   }
   jsonfile << "]," << std::endl;

   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NData[ipt]; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << StatErr[ipt][i] * 0.5;
   }
   jsonfile << "]," << std::endl;

   // sys errors "up/down"
   //
   
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NData[ipt]; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << SysErr[ipt][i] * 0.5;
   }
   jsonfile << "]," << std::endl;
   
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NData[ipt]; ++i )
   {
      jsonfile << ", 0.";
   }
   for ( int i=0; i<NData[ipt]; ++i )
   {
      jsonfile << ", " << SysErr[ipt][i] * 0.5;
   }
   jsonfile << "]," << std::endl;
   
   jsonfile << "\"binMin\":[], \"binmax\":[]" << std::endl;
   
   return;

}
