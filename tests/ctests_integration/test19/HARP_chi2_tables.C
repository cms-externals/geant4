std::string TEST_NAME="test19";
#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_IntermediateEnergy.h"

// #include "../test23/shared-root-macros/HARP.h"
#include "../test23/shared-root-macros/ReadHARPData.C"

#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"

// #include <vector>
// #include <map>

#include <iostream>
#include <fstream>

void HARP_chi2_table_regre( std::string beam, std::string target, /* std::string energy, */
                            /* std::string secondary, */ 
		            std::string model )
{

   std::vector< std::pair<double,int> > FW_3GeV_piplus;
   std::vector< std::pair<double,int> > FW_5GeV_piplus;
   std::vector< std::pair<double,int> > FW_8GeV_piplus;
   std::vector< std::pair<double,int> > FW_12GeV_piplus;

   std::vector< std::pair<double,int> > LA_3GeV_piplus;
   std::vector< std::pair<double,int> > LA_5GeV_piplus;
   std::vector< std::pair<double,int> > LA_8GeV_piplus;
   std::vector< std::pair<double,int> > LA_12GeV_piplus;
   
   std::vector< std::pair<double,int> > FW_3GeV_piminus;
   std::vector< std::pair<double,int> > FW_5GeV_piminus;
   std::vector< std::pair<double,int> > FW_8GeV_piminus;
   std::vector< std::pair<double,int> > FW_12GeV_piminus;

   std::vector< std::pair<double,int> > LA_3GeV_piminus;
   std::vector< std::pair<double,int> > LA_5GeV_piminus;
   std::vector< std::pair<double,int> > LA_8GeV_piminus;
   std::vector< std::pair<double,int> > LA_12GeV_piminus;

   // FW production
   //
   
   // piplus
   //
   
   ReadHARPData( beam, target, "3.0", "piplus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "3.0", "piplus", "FW", model, NDF_FW, Versions[iv] );
      FW_3GeV_piplus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }
   
   ReadHARPData( beam, target, "5.0", "piplus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "5.0", "piplus", "FW", model, NDF_FW, Versions[iv] );
      FW_5GeV_piplus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }

   ReadHARPData( beam, target, "8.0", "piplus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "8.0", "piplus", "FW", model, NDF_FW, Versions[iv] );
      FW_8GeV_piplus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }
   
   ReadHARPData( beam, target, "12.0", "piplus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "12.0", "piplus", "FW", model, NDF_FW, Versions[iv] );
      FW_12GeV_piplus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }
   
   // piminus
   //

   ReadHARPData( beam, target, "3.0", "piminus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "3.0", "piminus", "FW", model, NDF_FW, Versions[iv] );
      FW_3GeV_piminus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }
   
   ReadHARPData( beam, target, "5.0", "piminus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "5.0", "piminus", "FW", model, NDF_FW, Versions[iv] );
      FW_5GeV_piminus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }

   ReadHARPData( beam, target, "8.0", "piminus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "8.0", "piminus", "FW", model, NDF_FW, Versions[iv] );
      FW_8GeV_piminus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }

   ReadHARPData( beam, target, "12.0", "piminus", "FW" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_FW = 0;
      double chi2_FW = Chi2MomSpectrumHARP( beam, target, "12.0", "piminus", "FW", model, NDF_FW, Versions[iv] );
      FW_12GeV_piminus.push_back( std::pair<double,int>(chi2_FW,NDF_FW) );
   }

   // LA production
   //
   
   // piplus
   //
   
   ReadHARPData( beam, target, "3.0", "piplus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "3.0", "piplus", "LA", model, NDF_LA, Versions[iv] );
      LA_3GeV_piplus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "5.0", "piplus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "5.0", "piplus", "LA", model, NDF_LA, Versions[iv] );
      LA_5GeV_piplus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "8.0", "piplus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "8.0", "piplus", "LA", model, NDF_LA, Versions[iv] );
      LA_8GeV_piplus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "12.0", "piplus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "12.0", "piplus", "LA", model, NDF_LA, Versions[iv] );
      LA_12GeV_piplus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }
   
   // piminus
   //
   ReadHARPData( beam, target, "3.0", "piminus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "3.0", "piminus", "LA", model, NDF_LA, Versions[iv] );
      LA_3GeV_piminus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "5.0", "piminus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "5.0", "piminus", "LA", model, NDF_LA, Versions[iv] );
      LA_5GeV_piminus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "8.0", "piminus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "8.0", "piminus", "LA", model, NDF_LA, Versions[iv] );
      LA_8GeV_piminus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

   ReadHARPData( beam, target, "12.0", "piminus", "LA" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF_LA = 0;
      double chi2_LA = Chi2MomSpectrumHARP( beam, target, "12.0", "piminus", "LA", model, NDF_LA, Versions[iv] );
      LA_12GeV_piminus.push_back( std::pair<double,int>(chi2_LA,NDF_LA) );
   }

    // PRINT THE FORMATTED MAP INTO 4x4 TABLE !!!
   
   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_HARP.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;
   
   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(40) << model << ":  " << beam << " on " << target << " --> pions" << std::endl;
   out_table << std::right << std::setw(66) << " chi2/NDF calculated vs HARP data " << std::endl;
   out_table << " " << std::endl;

   out_table << std::right << std::setw(30) << "3 GeV/c" << std::right << std::setw(40) << "5 GeV/c" << std::endl;
   out_table << "     ";
   for ( int ic=0; ic<2; ++ic )
   {
   out_table << "     ";
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
   
   out_table << "\n pi+" << std::endl;
   out_table << "      FW  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_3GeV_piplus[iv].first/FW_3GeV_piplus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_5GeV_piplus[iv].first/FW_5GeV_piplus[iv].second;  
   }
   out_table << "\n      LA  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_3GeV_piplus[iv].first/LA_3GeV_piplus[iv].second; 
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_5GeV_piplus[iv].first/LA_5GeV_piplus[iv].second;  
   }
   out_table << "\n     Total";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_3GeV_piplus[iv].first+LA_3GeV_piplus[iv].first)/(FW_3GeV_piplus[iv].second+LA_3GeV_piplus[iv].second);  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_5GeV_piplus[iv].first+LA_5GeV_piplus[iv].first)/(FW_5GeV_piplus[iv].second+LA_5GeV_piplus[iv].second);
   }
   
   out_table << "\n pi-" << std::endl;
   out_table << "      FW  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_3GeV_piminus[iv].first/FW_3GeV_piminus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_5GeV_piminus[iv].first/FW_5GeV_piminus[iv].second;  
   }
   out_table << "\n      LA  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_3GeV_piminus[iv].first/LA_3GeV_piminus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_5GeV_piminus[iv].first/LA_5GeV_piminus[iv].second;  
   }
   out_table << "\n     Total";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_3GeV_piminus[iv].first+LA_3GeV_piminus[iv].first)/(FW_3GeV_piminus[iv].second+LA_3GeV_piminus[iv].second);  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_5GeV_piminus[iv].first+LA_5GeV_piminus[iv].first)/(FW_5GeV_piminus[iv].second+LA_5GeV_piminus[iv].second);
   }
   
   out_table << " " << std::endl;
   out_table << "\n" << std::right << std::setw(30) << "8 GeV/c" << std::right << std::setw(40) << "12 GeV/c" << std::endl;
   out_table << "     ";
   for ( int ic=0; ic<2; ++ic )
   {
   out_table << "     ";
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
   
   out_table << "\n pi+" << std::endl;
   out_table << "      FW  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_8GeV_piplus[iv].first/FW_8GeV_piplus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_12GeV_piplus[iv].first/FW_12GeV_piplus[iv].second;  
   }
   out_table << "\n      LA  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_8GeV_piplus[iv].first/LA_8GeV_piplus[iv].second; 
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_12GeV_piplus[iv].first/LA_12GeV_piplus[iv].second;  
   }
   out_table << "\n     Total";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_8GeV_piplus[iv].first+LA_8GeV_piplus[iv].first)/(FW_8GeV_piplus[iv].second+LA_8GeV_piplus[iv].second);  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_12GeV_piplus[iv].first+LA_12GeV_piplus[iv].first)/(FW_12GeV_piplus[iv].second+LA_12GeV_piplus[iv].second);
   }
   
   out_table << "\n pi-" << std::endl;
   out_table << "      FW  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_8GeV_piminus[iv].first/FW_8GeV_piminus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << FW_12GeV_piminus[iv].first/FW_12GeV_piminus[iv].second; 
   }
   out_table << "\n      LA  ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_8GeV_piminus[iv].first/LA_8GeV_piminus[iv].second;  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << LA_12GeV_piminus[iv].first/LA_12GeV_piminus[iv].second;  
   }
   out_table << "\n     Total";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_8GeV_piminus[iv].first+LA_8GeV_piminus[iv].first)/(FW_8GeV_piminus[iv].second+LA_8GeV_piminus[iv].second);  
   }
   out_table << "     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::setw(12) << (FW_12GeV_piminus[iv].first+LA_12GeV_piminus[iv].first)/(FW_12GeV_piminus[iv].second+LA_12GeV_piminus[iv].second);
   }

   out_table.close();
      
   return;
   
}
