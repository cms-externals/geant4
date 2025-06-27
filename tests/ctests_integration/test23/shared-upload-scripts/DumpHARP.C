#ifndef DumpHARP_C
#define DumpHARP_C

void DumpHARP(   std::ofstream& jsonfile, int cnt )
{

   // dN/dxF
   //
   jsonfile << "\"val\":" << std::endl;
   jsonfile << "[" << YHARP[cnt][0];
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", " << YHARP[cnt][i];
   }
   jsonfile << "]," << std::endl;
   
   // stat err "plus" ("up") - half of the total value
   //
   jsonfile << "\"errStatPlus\":" << std::endl;
   jsonfile << "[" << EYHARP[cnt][0]/2.;
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", " << EYHARP[cnt][i]/2.;
   }
   jsonfile << "]," << std::endl;
   
   // stat err "minus" ("down")
   //
   jsonfile << "\"errStatMinus\":" << std::endl;
   jsonfile << "[" << EYHARP[cnt][0]/2.;
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", " << EYHARP[cnt][i]/2.;
   }
   jsonfile << "]," << std::endl;

   // NOTE: NO explicit SYS errors 
   //       HARP data are believe to be published with TOTAL errors,
   //       i.e. the errors incl. both stat and sys. errors:
   //       http://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2008-26.pdf
   // sys error "up" - half of the 2.5% of the total xsec value
   //

   jsonfile << "\"errSysPlus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", 0.";
   }
   jsonfile << "]," << std::endl;

   // sys error "down" - half of the 2.5% of the total xsec value
   //
   jsonfile << "\"errSysMinus\":" << std::endl;
   jsonfile << "[0.";
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", 0.";
   }
   jsonfile << "]," << std::endl;
   
   // xf bin(s) - min edges
   //
   jsonfile << "\"binMin\":" << std::endl;
   jsonfile << "[" << XMinHARP[cnt][0];
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", " << XMinHARP[cnt][i];
   }
   jsonfile << "]," << std::endl;

   // xf bin(s) - max edges
   //
   jsonfile << "\"binMax\":" << std::endl;
   jsonfile << "[" << XMaxHARP[cnt][0];
   for ( int i=1; i<NPointsHARP[cnt]; ++i )
   {
      jsonfile << ", " << XMaxHARP[cnt][i];
   }
   jsonfile << "]" << std::endl; // NOTE: NO comma after this last square bracket 
                                 // since it's the last one in the list of data arrays

   return;

}

#endif
