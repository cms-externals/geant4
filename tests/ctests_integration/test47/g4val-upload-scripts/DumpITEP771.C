#ifndef DumpITEP771_C
#define DumpITEP771_C

void FillNeuKEBins();

int NumNeuKEBins=16;
double NeuKEBins[17];

void DumpITEP771(   std::ofstream& jsonfile, std::string& secondary )
{

   // dN/dxF
   //
   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << ExpValue[0];
   for ( int i=1; i<NPtKE; ++i )
   {
      jsonfile << ", " << ExpValue[i];
   }
   jsonfile << "]," << std::endl;
   
   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << ErrStat[0]/2.;
   for ( int i=1; i<NPtKE; ++i )
   {
      jsonfile << ", " << ErrStat[i]/2.;
   }
   jsonfile << "]," << std::endl;
   
   // stat err "minus" ("down")
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << ErrStat[0]/2.;
   for ( int i=1; i<NPtKE; ++i )
   {
      jsonfile << ", " << ErrStat[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "up" 
   //
   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[" << ErrSys[0]/2.;
   for ( int i=1; i<NPtKE; ++i )
   {
      jsonfile << ", " << ErrSys[i]/2.;
   }
   jsonfile << "]," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[" << ErrSys[0]/2.;
   for ( int i=1; i<NPtKE; ++i )
   {
      jsonfile << ", " << ErrSys[i]/2.;
   }
   jsonfile << "]," << std::endl;
 
   if ( secondary == "proton" )
   {
      // KE bin(s) - min edge
      //
      jsonfile << "\"binMin\":" << std::endl;
      jsonfile << "[" << KE[0]-0.01;
      for ( int i=1; i<NPtKE; ++i )
      {
         jsonfile << ", " << KE[i]-0.01;
      }
      jsonfile << "]," << std::endl;

      // KE bin(s) - max edges
      //
      jsonfile << "\"binMax\":" << std::endl;
      jsonfile << "[" << KE[0]+0.01;
      for ( int i=1; i<NPtKE; ++i )
      {
         jsonfile << ", " << KE[i]+0.01;
      }
      jsonfile << "]" << std::endl; // NOTE: NO comma after this last square bracket 
                                    // since it's the last one in the list of data arrays
   }
   else if ( secondary == "neutron" )
   {
      FillNeuKEBins();

      // KE bin(s) - min edge
      //
      jsonfile << "\"binMin\":" << std::endl;
      jsonfile << "[" << NeuKEBins[0];
      for ( int i=1; i<NumNeuKEBins; ++i )
      {
         jsonfile << ", " << NeuKEBins[i];
      }
      jsonfile << "]," << std::endl;

      // KE bin(s) - max edges
      //
      jsonfile << "\"binMax\":" << std::endl;
      jsonfile << "[" << NeuKEBins[1];
      for ( int i=1; i<NumNeuKEBins; ++i )
      {
         jsonfile << ", " << NeuKEBins[i+1];
      }
      jsonfile << "]" << std::endl; // NOTE: NO comma after this last square bracket 
                                    // since it's the last one in the list of data arrays
   }
   
   return;

}

void FillNeuKEBins()
{

   for ( int ike=0; ike<=NumNeuKEBins; ++ike )
   {
      NeuKEBins[ike] = 0.;
   }
   
   NeuKEBins[0] = 0.012;
   for ( int ike=1; ike<5; ++ike )
   {
      NeuKEBins[ike] = NeuKEBins[ike-1] + 0.002;
   }
   NeuKEBins[5] = 0.024;
   NeuKEBins[6] = 0.030;
   for ( int ike=7; ike<=9; ++ike )
   {
      NeuKEBins[ike] = NeuKEBins[ike-1] + 0.010;
   }
   for ( int ike=10; ike<=NumNeuKEBins; ++ike )
   {
      NeuKEBins[ike] = NeuKEBins[ike-1] + 0.020;
   }

   return;   

}


#endif
