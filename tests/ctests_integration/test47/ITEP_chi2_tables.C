#include "../test47/PlotITEPAnalysis.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include <iostream>
#include <fstream>

void ITEP_chi2_table_regre( std::string beam,  std::string en, std::string model )
{

   std::vector< std::pair<double,int> > C_proton;
   std::vector< std::pair<double,int> > C_neutron;
   
   std::vector< std::pair<double,int> > Cu_proton;
   std::vector< std::pair<double,int> > Cu_neutron;

   std::vector< std::pair<double,int> > Pb_proton;
   std::vector< std::pair<double,int> > Pb_neutron;
      
   std::vector< std::pair<double,int> > U_proton;
   std::vector< std::pair<double,int> > U_neutron;

   // production of proton/neutron on C
   //
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "proton", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "proton", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "proton", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "proton", "159.6",  model, NDF, Versions[iv] );
      C_proton.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "neutron", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "neutron", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "neutron", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "C", en, "neutron", "159.6",  model, NDF, Versions[iv] );
      C_neutron.push_back( std::pair<double,int>(chi2,NDF) );
   }   

   // production of proton/neutron on Cu
   //
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "proton", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "proton", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "proton", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "proton", "159.6",  model, NDF, Versions[iv] );
      Cu_proton.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "neutron", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "neutron", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "neutron", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Cu", en, "neutron", "159.6",  model, NDF, Versions[iv] );
      Cu_neutron.push_back( std::pair<double,int>(chi2,NDF) );
   }   

   // production of proton/neutron on Pb
   //
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "proton", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "proton", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "proton", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "proton", "159.6",  model, NDF, Versions[iv] );
      Pb_proton.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "neutron", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "neutron", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "neutron", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "Pb", en, "neutron", "159.6",  model, NDF, Versions[iv] );
      Pb_neutron.push_back( std::pair<double,int>(chi2,NDF) );
   }   

   // production of proton/neutron on U
   //
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "proton", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "proton", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "proton", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "proton", "159.6",  model, NDF, Versions[iv] );
      U_proton.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "neutron", "59.1",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "neutron", "89.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "neutron", "119.0",  model, NDF, Versions[iv] );
      chi2 += Chi2KESpectrumITEP( beam, "U", en, "neutron", "159.6",  model, NDF, Versions[iv] );
      U_neutron.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   
    // PRINT THE FORMATTED MAP INTO 4x4 TABLE !!!
   
   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + en + "GeV_ITEP771.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;
   
   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(20) << model << ":  " << en << " GeV/c " << beam << " on nucleus --> proton/neutron" << std::endl;
   out_table << std::right << std::setw(60) << " chi2/NDF calculated vs ITEP771 data " << std::endl;
   out_table << " " << std::endl;

   out_table << std::right << std::setw(30) << "Carbon (C)" << std::right << std::setw(40) << "Copper (Cu)" << std::endl;
   
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

   out_table << "\n proton "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << C_proton[iv].first/C_proton[iv].second ;
   }   
   out_table << "       ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << Cu_proton[iv].first/Cu_proton[iv].second ;
   }
   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << C_neutron[iv].first/C_neutron[iv].second ;
   }   
   out_table << "       ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << Cu_neutron[iv].first/Cu_neutron[iv].second ;
   }
   
   out_table << "\n " << std::endl;

   out_table << std::right << std::setw(30) << "Lead (Pb) " << std::right << std::setw(40) << "Uranium (U)" << std::endl;
   
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

   out_table << "\n proton "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << Pb_proton[iv].first/Pb_proton[iv].second ;
   }   
   out_table << "       ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << U_proton[iv].first/U_proton[iv].second ;
   }
   out_table << "\n neutron"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << Pb_neutron[iv].first/Pb_neutron[iv].second ;
   }   
   out_table << "       ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << std::right << setw(12) << U_neutron[iv].first/U_neutron[iv].second ;
   }

   out_table.close();

   return;

}
