#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

const int NSetsFW = 4;
const int NDataPointsFW[4] = { 10, 10, 10, 10 };
const int NSetsLA = 9;
const int NDataPointsLA[9] = { 10, 11, 10, 9, 8, 8, 8, 8, 8 };

void read_write( std::string input, std::string output, std::string region )
{

   std::string location = "/home/yarba_j/work/g4tests/verification/hadronic/test35/data/";
   location += input;
   
   ifstream infile;
   infile.open( location.c_str() );
   
   ofstream outfile;
   outfile.open ( output.c_str(),ios::ate|ios::app);
   
   float the1, the2, x1, x2, val, err;
   int NN1, NN2;
   
   if ( region == "FW" )
   {
      for ( int i=0; i<NSetsFW; ++i )
      {
         for ( int j=0; j<NDataPointsFW[i]; ++j )
	 {
	    infile >> the1 >> the2 >> x1 >> x2 >> val >> err ;
//	    std::cout << the1 << " " << the1 << " " << x1 << " " << x2 << " " << val << " " << err << std::endl;
	    if ( j == 0 )
	    {
	       outfile << the1 << " " << the2 << std::endl;
	       outfile << NDataPointsFW[i] << std::endl;
	    }
	    outfile << x1 << " " << x2 << " " << val << " " << err << std::endl;
	 }
      }  
   }
   else if ( region == "LA" )
   {
      for ( int i=0; i<NSetsLA; ++i )
      {
         for ( int j=0; j<NDataPointsLA[i]; ++j )
	 {
	    infile >> the1 >> the2 >> x1 >> x2 >> val >> err ;
	    if ( j == 0 )
	    {
	       outfile << the1 << " " << the2 << std::endl;
	       outfile << NDataPointsLA[i] << std::endl;
	    }
	    outfile << x1 << " " << x2 << " " << val << " " << err << std::endl;
	 }
      }  
   }
   
   return;

}
