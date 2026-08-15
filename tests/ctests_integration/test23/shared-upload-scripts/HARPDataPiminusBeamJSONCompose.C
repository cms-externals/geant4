#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>
#include <vector>

#include "TSystem.h"
#include "Rtypes.h"

#include "../shared-root-macros/HARP.h"
#include "../shared-root-macros/ReadHARPData.C"

#include "DumpHARP.C"

void HARPDataJSONCompose( std::ofstream&, std::string, int, std::string );

void HARPDataJSONMaster( std::string tgt, int tgtlnk )
{

   std::string fjsonname = "HARP-spectra-piminus-";
   fjsonname += tgt;
   fjsonname += ".json";
   
   std::ofstream jsonfile;
   jsonfile.open( fjsonname.c_str(), ios::ate|ios::app ); 
   
   jsonfile << "{\"ResultList\":[" << std::endl;
   
   HARPDataJSONCompose( jsonfile, tgt, tgtlnk, "piplus" );
      
   HARPDataJSONCompose( jsonfile, tgt, tgtlnk, "piminus" );
   
// save proton production business for later...
// the data aren't fantastic anyway - the errors are huge !
//
//   HARPDataJSONCompose( jsonfile, "proton" );   
   
   jsonfile << "]}" << std::endl;
    
   return;

}

void HARPDataJSONCompose( std::ofstream& jsonfile, std::string tgt, int tgtlnk, std::string secondary )
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

   //
   // 3GeV/c pi- on A
   //
   // FW production

   ReadHARPData( "piminus", tgt, "3.0", secondary, "FW" ); 
   
       
   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
      
      // NOTE: referencelnk=44 is FWD pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":44,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE: beamlnk=46 is 3 GeV/c pi- 
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":46,\"targetlnk\":"<< tgtlnk <<",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on "<< tgt << " interactions at 3GeV/c, " 
               << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }
   
   //
   // LA production
   //
   
   ReadHARPData( "piminus",tgt, "3.0", secondary, "LA" );

   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
   
      // NOTE: referencelnk=45 is LA pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":45,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE-1: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE-2: beamlnk=46 is 3 GeV/c pi- 
      //
      // NOTE-3: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":46,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 3GeV/c, " 
               << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }

   //
   // 5GeV/c pi- on A
   //
   // FW production

   ReadHARPData( "piminus", tgt, "5.0", secondary, "FW" ); 
   
       
   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
      
      // NOTE: referencelnk=44 is FWD pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":44,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE: beamlnk=47 is 5 GeV/c pi- 
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":47,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 5GeV/c, " 
               << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }
   
   //
   // LA production
   //
   
   ReadHARPData( "piminus",tgt, "5.0", secondary, "LA" );

   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
   
      // NOTE: referencelnk=45 is LA pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":45,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE-1: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE-2: beamlnk=47 is 5 GeV/c pi-
      //
      // NOTE-3: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":47,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 5GeV/c, " 
      << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }

   //
   // 8GeV/c pi- on A
   //
   // FW production

   ReadHARPData( "piminus", tgt, "8.0", secondary, "FW" ); 
   
       
   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
      
      // NOTE: referencelnk=44 is FWD pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":44,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE: beamlnk=48 is 8 GeV/c pi- 
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":48,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 8GeV/c, " 
      << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }
   
   //
   // LA production
   //
   
   ReadHARPData( "piminus",tgt, "8.0", secondary, "LA" );

   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
   
      // NOTE: referencelnk=45 is LA pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":45,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE-1: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE-2: beamlnk=48 is 8 GeV/c pi- 
      //
      // NOTE-3: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":48,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 8GeV/c, " 
      << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }

   //
   // 12GeV/c pi- on A
   //
   // FW production

   ReadHARPData( "piminus", tgt, "12.0", secondary, "FW" ); 
   
       
   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
      
      // NOTE: referencelnk=44 is FWD pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":44,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE: beamlnk=49 is 12 GeV/c pi- 
      //
      // NOTE: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":49,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 12GeV/c, " 
      << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }
   
   //
   // LA production
   //
   
   ReadHARPData( "piminus",tgt, "12.0", secondary, "LA" );

   if ( NSetsHARP <= 0 )
   {
      std::cout << " Invalid record, NSetsHARP = " << NSetsHARP << std::endl;
      return;
   }
   
   for ( int i=0; i<NSetsHARP; ++i )
   {
      if ( NPointsHARP[i] <= 0 )
      {
         std::cout << " Invalid record, NPointsHARP[" << i << "] = " << NPointsHARP[i] << std::endl;
	 return;
      }
   
      jsonfile << "{" << std::endl;
   
      // NOTE: referencelnk=45 is LA pion production in pi+/- on A from HARP
      //
      jsonfile << "\"trid\":1,\"testlnk\":0,\"referencelnk\":45,\"mcdetaillnk\":1," << std::endl;
   
      // NOTE-1: observablelnk=8 means diff. xsec dsigma / dtheta dp
      //
      // NOTE-2: beamlnk=49 is 12 GeV/c pi- 
      //
      // NOTE-3: reactionlnk=1 means "particle production"
      //
      jsonfile << "\"beamlnk\":49,\"targetlnk\":" << tgtlnk << ",\"observablelnk\":8,\"secondarylnk\":" 
               << seclink << ",\"reactionlnk\":1," << std::endl;

      jsonfile << "\"datatable\":{" << std::endl; // start datatable
      
      jsonfile << "\"dtid\":1,\"datatypeslnk\":2," << std::endl;
      jsonfile << "\"title\":\"Production of " << sec << " in pi- on " << tgt << " interactions at 12GeV/c, " 
      << AngleBinHARP[0][i] << "<theta<" << AngleBinHARP[1][i] << " [rad]\"," << std::endl;
      jsonfile << "\"npoints\":0," << std::endl;
      jsonfile << "\"nbins\":[" << NPointsHARP[i] << "]," << std::endl;
      jsonfile << "\"axisTitle\":[\"momentum of " << sec << " [GeV/c]\",\"Cross section [mb/(GeV/c/rad)]\"]," << std::endl;
      
      DumpHARP( jsonfile, i );

      jsonfile << "}," << std::endl; // end of datatable 

// -->      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,\"parvalue\":{}" << std::endl;
      jsonfile << "\"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1," << std::endl;
      jsonfile << "\"parnames\":[\"theta interval for secondary particle\"], \"parvalues\":[\"" << AngleBinHARP[0][i] 
               << "-" << AngleBinHARP[1][i] << "\"]" << std::endl;

      jsonfile << "}," << std::endl;  

   }
         
   return;

}
