#include "../test48/PlotPiMinusStopping.C"

void PiMinus_chi2_table_regre( /* std::string model */ )
{

   std::vector< std::pair<double,int> > C;
   std::vector< std::pair<double,int> > N;
   std::vector< std::pair<double,int> > O;
   std::vector< std::pair<double,int> > Al;
   std::vector< std::pair<double,int> > Cu;
   std::vector< std::pair<double,int> > Pb;
   std::vector< std::pair<double,int> > Ta;
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 = calcChi2Madey( "C", "BertiniPreCo", NDF, Versions[iv] );
      C.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "N", "BertiniPreCo", NDF, Versions[iv] );
      N.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "O", "BertiniPreCo", NDF, Versions[iv] );
      O.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "Al", "BertiniPreCo", NDF, Versions[iv] );
      Al.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "Cu", "BertiniPreCo", NDF, Versions[iv] );
      Cu.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "Pb", "BertiniPreCo", NDF, Versions[iv] );
      Pb.push_back( std::pair<double,int>(chi2,NDF) );
      chi2 = 0.;
      NDF = 0;
      chi2 = calcChi2Madey( "Ta", "BertiniPreCo", NDF, Versions[iv] );
      Ta.push_back( std::pair<double,int>(chi2,NDF) );
   }

   std::ofstream out_table;
   std::string table_name = "BertiniPreCo_piminus_capture.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );

   out_table << std::right << std::setw(75) << "Bertini(PreCo): caputure of pi- on nucleus --> neutron + X" << std::endl;
   out_table << std::right << std::setw(85) << " chi2/NDF calculated vs data from R.Madey et al., Phys.Rev.C25, 3050 (1982)" << std::endl;
   out_table << " " << std::endl;
   
   out_table << std::right << std::setw(30) << "Carbon (C)" << std::right << std::setw(40) << "Nitrogen (N)" << std::endl;
   
   for ( int ic=0; ic<2; ++ic )
   {
      out_table << "        ";
      for ( size_t iv=0; iv<NVersions; ++iv )
      {
         std::string ver = Versions[iv];
         if ( ver.find("geant4-") != std::string::npos )
         {
            ver.erase( 0, std::string("geant4-").length() );
         }
         if ( ver.find( "patch-") != std::string::npos )
         {
            ver.replace( ver.find("patch-"), std::string( "patch-").length(), "p" ); 
         }
         if ( ver.find( "beta-") != std::string::npos )
         {
            ver.replace( ver.find("beta-"), std::string("beta-").length(), "b" );
         }
         while ( ver.find("-") != std::string::npos )
         {
            ver.replace( ver.find("-"), std::string("-").length(), "." );
         }
         out_table << std::setw(12) << ver;
      }
   }

   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << C[iv].first/C[iv].second ;
   }   
   out_table << "        ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << N[iv].first/N[iv].second ;
   }
   
   out_table << " " << std::endl;

   out_table << std::right << std::setw(30) << "Oxygen (O)" << std::right << std::setw(40) << "Aluminum (Al)" << std::endl;
   
   for ( int ic=0; ic<2; ++ic )
   {
      out_table << "        ";
      for ( size_t iv=0; iv<NVersions; ++iv )
      {
         std::string ver = Versions[iv];
         if ( ver.find("geant4-") != std::string::npos )
         {
            ver.erase( 0, std::string("geant4-").length() );
         }
         if ( ver.find( "patch-") != std::string::npos )
         {
            ver.replace( ver.find("patch-"), std::string( "patch-").length(), "p" ); 
         }
         if ( ver.find( "beta-") != std::string::npos )
         {
            ver.replace( ver.find("beta-"), std::string("beta-").length(), "b" );
         }
         while ( ver.find("-") != std::string::npos )
         {
            ver.replace( ver.find("-"), std::string("-").length(), "." );
         }
         out_table << std::setw(12) << ver;
      }
   }

   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << O[iv].first/O[iv].second ;
   }   
   out_table << "        ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << Al[iv].first/Al[iv].second ;
   }
   
   out_table << " " << std::endl;
   
   out_table << std::right << std::setw(30) << "Copper (Cu)" << std::right << std::setw(40) << "Lead (Pb)" << std::endl;
   
   for ( int ic=0; ic<2; ++ic )
   {
      out_table << "        ";
      for ( size_t iv=0; iv<NVersions; ++iv )
      {
         std::string ver = Versions[iv];
         if ( ver.find("geant4-") != std::string::npos )
         {
            ver.erase( 0, std::string("geant4-").length() );
         }
         if ( ver.find( "patch-") != std::string::npos )
         {
            ver.replace( ver.find("patch-"), std::string( "patch-").length(), "p" ); 
         }
         if ( ver.find( "beta-") != std::string::npos )
         {
            ver.replace( ver.find("beta-"), std::string("beta-").length(), "b" );
         }
         while ( ver.find("-") != std::string::npos )
         {
            ver.replace( ver.find("-"), std::string("-").length(), "." );
         }
         out_table << std::setw(12) << ver;
      }
   }

   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << Cu[iv].first/Cu[iv].second ;
   }   
   out_table << "        ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << Pb[iv].first/Pb[iv].second ;
   }
   
   out_table << " " << std::endl;
   
   out_table << std::right << std::setw(30) << "Tantalum (Ta)" << std::endl;
   
   for ( int ic=0; ic<1; ++ic )
   {
      out_table << "        ";
      for ( size_t iv=0; iv<NVersions; ++iv )
      {
         std::string ver = Versions[iv];
         if ( ver.find("geant4-") != std::string::npos )
         {
            ver.erase( 0, std::string("geant4-").length() );
         }
         if ( ver.find( "patch-") != std::string::npos )
         {
            ver.replace( ver.find("patch-"), std::string( "patch-").length(), "p" ); 
         }
         if ( ver.find( "beta-") != std::string::npos )
         {
            ver.replace( ver.find("beta-"), std::string("beta-").length(), "b" );
         }
         while ( ver.find("-") != std::string::npos )
         {
            ver.replace( ver.find("-"), std::string("-").length(), "." );
         }
         out_table << std::setw(12) << ver;
      }
   }

   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << Ta[iv].first/Ta[iv].second ;
   }   

   out_table.close();

   return;

}
