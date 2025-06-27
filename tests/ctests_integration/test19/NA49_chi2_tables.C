std::string TEST_NAME="test19";
#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test19/G4MODELS_HighEnergy.h"

#include "../test23/shared-root-macros/ReadNA49Data.C"

#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA49.C"

#include <iostream>
#include <fstream>

void NA49_chi2_table_regre( std::string beam, std::string target, std::string model )
{

   // first of, check if G4INSTALL env.var is setup;
   // if not, bail out (because the Read<...>Data.C will bomb anyway)
   //
   std::string g4install = std::string( gSystem->Getenv("G4INSTALL") );
   if ( g4install.empty() ) 
   {
      std::cout << " G4INSTALL env.var. is NOT setup; bail out" << std::endl;
      return;
   } 

   std::vector< std::pair<double,int> > dN_dxF_piplus;
   std::vector< std::pair<double,int> > pT_dxF_piplus;
   std::vector< std::pair<double,int> > d2sigma_dxFdpT_piplus;
   
   std::vector< std::pair<double,int> > dN_dxF_piminus;
   std::vector< std::pair<double,int> > pT_dxF_piminus;
   std::vector< std::pair<double,int> > d2sigma_dxFdpT_piminus;

   std::vector< std::pair<double,int> > dN_dxF_proton;
   std::vector< std::pair<double,int> > pT_dxF_proton;
   std::vector< std::pair<double,int> > d2sigma_dxFdpT_proton;
   
   std::vector< std::pair<double,int> > dN_dxF_antiproton;
   std::vector< std::pair<double,int> > pT_dxF_antiproton;
   std::vector< std::pair<double,int> > d2sigma_dxFdpT_antiproton;

   std::vector< std::pair<double,int> > dN_dxF_neutron;
   
   // piplus
   //
   
   readIntegratedSpectra( beam, target, "piplus" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "dNdxF", model, NDF, Versions[iv] );
      dN_dxF_piplus.push_back( std::pair<double,int>(chi2,NDF) );
      NDF = 0;
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "pT", model, NDF, Versions[iv] );
      pT_dxF_piplus.push_back( std::pair<double,int>(chi2,NDF) );
   }
   
   readDDiffSpectra( beam, target, "piplus" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, "piplus", model, NDF, Versions[iv] );
      d2sigma_dxFdpT_piplus.push_back( std::pair<double,int>(chi2,NDF) );
   } 
   
   // piminus
   //

   readIntegratedSpectra( beam, target, "piminus" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "dNdxF", model, NDF, Versions[iv] );
      dN_dxF_piminus.push_back( std::pair<double,int>(chi2,NDF) );
      NDF = 0;
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "pT", model, NDF, Versions[iv] );
      pT_dxF_piminus.push_back( std::pair<double,int>(chi2,NDF) );
   }
     
   readDDiffSpectra( beam, target, "piminus" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, "piminus", model, NDF, Versions[iv] );
      d2sigma_dxFdpT_piminus.push_back( std::pair<double,int>(chi2,NDF) );
   } 
   
   // proton
   //

   readIntegratedSpectra( beam, target, "proton" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "proton", "dNdxF", model, NDF, Versions[iv] );
      dN_dxF_proton.push_back( std::pair<double,int>(chi2,NDF) );
      NDF = 0;
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, "proton", "pT", model, NDF, Versions[iv] );
      pT_dxF_proton.push_back( std::pair<double,int>(chi2,NDF) );
   }

   readDDiffSpectra( beam, target, "proton" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, "proton", model, NDF, Versions[iv] );
      d2sigma_dxFdpT_proton.push_back( std::pair<double,int>(chi2,NDF) );
   } 
   
   // antiproton
   //

   readIntegratedSpectra( beam, target, "antiproton" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "dNdxF", model, NDF, Versions[iv] );
      dN_dxF_antiproton.push_back( std::pair<double,int>(chi2,NDF) );
      NDF = 0;
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "pT", model, NDF, Versions[iv] );
      pT_dxF_antiproton.push_back( std::pair<double,int>(chi2,NDF) );
   }

   std::cout << "Now processing DDiff for Pbar" << std::endl;

/* FIXME !!! NEEDS TO BE ADDED LATER !!! */
   readDDiffSpectra( beam, target, "antiproton" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, "antiproton", model, NDF, Versions[iv] );
      d2sigma_dxFdpT_antiproton.push_back( std::pair<double,int>(chi2,NDF) );
   } 
   std::cout << "Read-in DDiff for Pbar, calculated chi2" << std::endl;
/* */   
   // neutron
   //
   
   readIntegratedSpectra( beam, target, "neutron" );
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "neutron", "dNdxF", model, NDF, Versions[iv] );
      dN_dxF_neutron.push_back( std::pair<double,int>(chi2,NDF) );
   }

   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_NA49.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(30) << model << ":  158GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(62) << " chi2/NDF calculated vs NA49 data " << std::endl;
   out_table << " " << std::endl;

   out_table << std::right << std::setw(30) << "dN / dxF" << std::right << std::setw(45) << "<pt> vs xF" << std::endl;
   for ( int ic=0; ic<2; ++ic )
   {
   out_table << "           ";
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

   out_table << "\n pi+       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << dN_dxF_piplus[iv].first/dN_dxF_piplus[iv].second ;
   }   
   out_table << "           ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << setw(12) << pT_dxF_piplus[iv].first/pT_dxF_piplus[iv].second ;
   }
   out_table << "\n pi-       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << dN_dxF_piminus[iv].first/dN_dxF_piminus[iv].second ;
   }   
   out_table << "           ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << setw(12) << pT_dxF_piminus[iv].first/pT_dxF_piminus[iv].second ;
   }
   out_table << "\n proton    "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << dN_dxF_proton[iv].first/dN_dxF_proton[iv].second ;
   }   
   out_table << "           ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << setw(12) << pT_dxF_proton[iv].first/pT_dxF_proton[iv].second ;
   }
   out_table << "\n antiproton"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << dN_dxF_antiproton[iv].first/dN_dxF_antiproton[iv].second ;
   }   
   out_table << "           ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      out_table << setw(12) << pT_dxF_antiproton[iv].first/pT_dxF_antiproton[iv].second ;
   }
   out_table << "\n neutron   "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << dN_dxF_neutron[iv].first/dN_dxF_neutron[iv].second ;
   }   

   out_table << "\n " << std::endl;

   out_table << std::right << std::setw(40) << "d2sigma / dxF dpT" << std::endl;
   out_table << "           ";
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

   out_table << "\n pi+       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << d2sigma_dxFdpT_piplus[iv].first/d2sigma_dxFdpT_piplus[iv].second ;
   }   
   out_table << "\n pi-       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << d2sigma_dxFdpT_piminus[iv].first/d2sigma_dxFdpT_piminus[iv].second ;
   }   
   out_table << "\n proton    ";
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << d2sigma_dxFdpT_proton[iv].first/d2sigma_dxFdpT_proton[iv].second ;
   }   
   out_table << "\n antiproton"; //          will be added shortly !!!";
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << d2sigma_dxFdpT_antiproton[iv].first/d2sigma_dxFdpT_antiproton[iv].second ;
   }   

   out_table.close();

   return;

}

void NA49_chi2_table_models( std::string beam, std::string target, std::string cur_ver="." )
{

   // first of, check if G4INSTALL env.var is setup;
   // if not, bail out (because the Read<...>Data.C will bomb anyway)
   //
   std::string g4install = std::string( gSystem->Getenv("G4INSTALL") );
   if ( g4install.empty() ) 
   {
      std::cout << " G4INSTALL env.var. is NOT setup; bail out" << std::endl;
      return;
   } 
   
   std::ofstream out_table;
   std::string table_name = "models_" + beam + "_" + target + "_NA49.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   std::string ver = cur_ver;
   if ( cur_ver == "." ) ver = std::string( gSystem->Getenv("G4RELEASE") );
   out_table << std::right << std::setw(45) << "  FTFP vs QGSP as of " << ver << std::endl;
   out_table << std::right << std::setw(35) << "158 GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(55) << " chi2/NDF calculated vs NA49 data " << std::endl;
   out_table << " " << std::endl;

   out_table << " " << std::endl;

   out_table << std::right << std::setw(25) << "dN / dxF" 
             << std::right << std::setw(25) << "<pt> vs xF"
	     << std::right << std::setw(25) << "d2sigma / dxF dpt" << std::endl;

   out_table << "           " << std::setw(10) << "QGSP" << std::setw(10) << "FTFP" 
                              << std::setw(10) << "QGSP" << std::setw(10) << "FTFP"
			      << std::setw(10) << "QGSP" << std::setw(10) << "FTFP"; 

   double chi2_qgsp_dNdxF = 0.;
   int NDF_qgsp_dNdxF = 0;
   double chi2_ftfp_dNdxF = 0.;
   int NDF_ftfp_dNdxF = 0;
   
   double chi2_qgsp_pt = 0.;
   int NDF_qgsp_pt = 0;
   double chi2_ftfp_pt = 0.;
   int NDF_ftfp_pt = 0;

   double chi2_qgsp_ddiff = 0.;
   int NDF_qgsp_ddiff = 0;
   double chi2_ftfp_ddiff = 0.;
   int NDF_ftfp_ddiff = 0;


   readIntegratedSpectra( beam, target, "piplus" );
   chi2_qgsp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "dNdxF", "qgsp", NDF_qgsp_dNdxF, cur_ver );
   chi2_qgsp_pt = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "pT", "qgsp", NDF_qgsp_pt, cur_ver );
   chi2_ftfp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "dNdxF", "ftfp", NDF_ftfp_dNdxF, cur_ver );
   chi2_ftfp_pt = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "pT", "ftfp", NDF_ftfp_pt, cur_ver );

   readDDiffSpectra( beam, target, "piplus" );
   chi2_qgsp_ddiff = Chi2DDiffXSecNA49( beam, target, "piplus", "qgsp", NDF_qgsp_ddiff, cur_ver );
   chi2_ftfp_ddiff = Chi2DDiffXSecNA49( beam, target, "piplus", "ftfp", NDF_ftfp_ddiff, cur_ver );
   
   out_table << "\n pi+       "; 
   out_table << std::right << std::setw(10) << chi2_qgsp_dNdxF/NDF_qgsp_dNdxF << std::right << std::setw(10) << chi2_ftfp_dNdxF/NDF_ftfp_dNdxF
             << std::right << std::setw(10) << chi2_qgsp_pt/NDF_qgsp_pt << std::setw(10) << chi2_ftfp_pt/NDF_ftfp_pt
             << std::right << std::setw(10) << chi2_qgsp_ddiff/NDF_qgsp_ddiff << std::setw(10) << chi2_ftfp_ddiff/NDF_ftfp_ddiff; // << std::endl;
   
   chi2_qgsp_dNdxF = 0.;
   NDF_qgsp_dNdxF = 0;
   chi2_ftfp_dNdxF = 0.;
   NDF_ftfp_dNdxF = 0;
   
   chi2_qgsp_pt = 0.;
   NDF_qgsp_pt = 0;
   chi2_ftfp_pt = 0.;
   NDF_ftfp_pt = 0;

   chi2_qgsp_ddiff = 0.;
   NDF_qgsp_ddiff = 0;
   chi2_ftfp_ddiff = 0.;
   NDF_ftfp_ddiff = 0;

   readIntegratedSpectra( beam, target, "piminus" );
   chi2_qgsp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "dNdxF", "qgsp", NDF_qgsp_dNdxF, cur_ver );
   chi2_qgsp_pt = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "pT", "qgsp", NDF_qgsp_pt, cur_ver );
   chi2_ftfp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "dNdxF", "ftfp", NDF_ftfp_dNdxF, cur_ver );
   chi2_ftfp_pt = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "pT", "ftfp", NDF_ftfp_pt, cur_ver );

   readDDiffSpectra( beam, target, "piminus" );
   chi2_qgsp_ddiff = Chi2DDiffXSecNA49( beam, target, "piminus", "qgsp", NDF_qgsp_ddiff, cur_ver );
   chi2_ftfp_ddiff = Chi2DDiffXSecNA49( beam, target, "piminus", "ftfp", NDF_ftfp_ddiff, cur_ver );

   out_table << "\n pi-       "; 
   out_table << std::right << std::setw(10) << chi2_qgsp_dNdxF/NDF_qgsp_dNdxF << std::right << std::setw(10) << chi2_ftfp_dNdxF/NDF_ftfp_dNdxF
             << std::right << std::setw(10) << chi2_qgsp_pt/NDF_qgsp_pt << std::setw(10) << chi2_ftfp_pt/NDF_ftfp_pt
             << std::right << std::setw(10) << chi2_qgsp_ddiff/NDF_qgsp_ddiff << std::setw(10) << chi2_ftfp_ddiff/NDF_ftfp_ddiff; // << std::endl;

   chi2_qgsp_dNdxF = 0.;
   NDF_qgsp_dNdxF = 0;
   chi2_ftfp_dNdxF = 0.;
   NDF_ftfp_dNdxF = 0;
   
   chi2_qgsp_pt = 0.;
   NDF_qgsp_pt = 0;
   chi2_ftfp_pt = 0.;
   NDF_ftfp_pt = 0;

   chi2_qgsp_ddiff = 0.;
   NDF_qgsp_ddiff = 0;
   chi2_ftfp_ddiff = 0.;
   NDF_ftfp_ddiff = 0;

   readIntegratedSpectra( beam, target, "proton" );
   chi2_qgsp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "proton", "dNdxF", "qgsp", NDF_qgsp_dNdxF, cur_ver );
   chi2_qgsp_pt = Chi2IntegratedSpectrumNA49( beam, target, "proton", "pT", "qgsp", NDF_qgsp_pt, cur_ver );
   chi2_ftfp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "proton", "dNdxF", "ftfp", NDF_ftfp_dNdxF, cur_ver );
   chi2_ftfp_pt = Chi2IntegratedSpectrumNA49( beam, target, "proton", "pT", "ftfp", NDF_ftfp_pt, cur_ver );

   readDDiffSpectra( beam, target, "proton" );
   chi2_qgsp_ddiff = Chi2DDiffXSecNA49( beam, target, "proton", "qgsp", NDF_qgsp_ddiff, cur_ver );
   chi2_ftfp_ddiff = Chi2DDiffXSecNA49( beam, target, "proton", "ftfp", NDF_ftfp_ddiff, cur_ver );

   out_table << "\n proton    "; 
   out_table << std::right << std::setw(10) << chi2_qgsp_dNdxF/NDF_qgsp_dNdxF << std::right << std::setw(10) << chi2_ftfp_dNdxF/NDF_ftfp_dNdxF
             << std::right << std::setw(10) << chi2_qgsp_pt/NDF_qgsp_pt << std::setw(10) << chi2_ftfp_pt/NDF_ftfp_pt
             << std::right << std::setw(10) << chi2_qgsp_ddiff/NDF_qgsp_ddiff << std::setw(10) << chi2_ftfp_ddiff/NDF_ftfp_ddiff; // << std::endl;

   chi2_qgsp_dNdxF = 0.;
   NDF_qgsp_dNdxF = 0;
   chi2_ftfp_dNdxF = 0.;
   NDF_ftfp_dNdxF = 0;
   
   chi2_qgsp_pt = 0.;
   NDF_qgsp_pt = 0;
   chi2_ftfp_pt = 0.;
   NDF_ftfp_pt = 0;

   chi2_qgsp_ddiff = 0.;
   NDF_qgsp_ddiff = 0;
   chi2_ftfp_ddiff = 0.;
   NDF_ftfp_ddiff = 0;

   // FIXME !!!
   // DDiff. spectra for antiproton need to be added !!!
   //
   
   readIntegratedSpectra( beam, target, "antiproton" );
   chi2_qgsp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "dNdxF", "qgsp", NDF_qgsp_dNdxF, cur_ver );
   chi2_qgsp_pt = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "pT", "qgsp", NDF_qgsp_pt, cur_ver );
   chi2_ftfp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "dNdxF", "ftfp", NDF_ftfp_dNdxF, cur_ver );
   chi2_ftfp_pt = Chi2IntegratedSpectrumNA49( beam, target, "antiproton", "pT", "ftfp", NDF_ftfp_pt, cur_ver );

   readDDiffSpectra( beam, target, "antiproton" );
   chi2_qgsp_ddiff = Chi2DDiffXSecNA49( beam, target, "antiproton", "qgsp", NDF_qgsp_ddiff, cur_ver );
   chi2_ftfp_ddiff = Chi2DDiffXSecNA49( beam, target, "antiproton", "ftfp", NDF_ftfp_ddiff, cur_ver );

   out_table << "\n antiproton"; 
   out_table << std::right << std::setw(10) << chi2_qgsp_dNdxF/NDF_qgsp_dNdxF << std::right << std::setw(10) << chi2_ftfp_dNdxF/NDF_ftfp_dNdxF
             << std::right << std::setw(10) << chi2_qgsp_pt/NDF_qgsp_pt << std::setw(10) << chi2_ftfp_pt/NDF_ftfp_pt
             << std::right << std::setw(10) << chi2_qgsp_ddiff/NDF_qgsp_ddiff << std::setw(10) << chi2_ftfp_ddiff/NDF_ftfp_ddiff; 

   chi2_qgsp_dNdxF = 0.;
   NDF_qgsp_dNdxF = 0;
   chi2_ftfp_dNdxF = 0.;
   NDF_ftfp_dNdxF = 0;
   
   chi2_qgsp_pt = 0.;
   NDF_qgsp_pt = 0;
   chi2_ftfp_pt = 0.;
   NDF_ftfp_pt = 0;

   chi2_qgsp_ddiff = 0.;
   NDF_qgsp_ddiff = 0;
   chi2_ftfp_ddiff = 0.;
   NDF_ftfp_ddiff = 0;

   readIntegratedSpectra( beam, target, "neutron" );
   chi2_qgsp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "neutron", "dNdxF", "qgsp", NDF_qgsp_dNdxF, cur_ver );
   chi2_ftfp_dNdxF = Chi2IntegratedSpectrumNA49( beam, target, "neutron", "dNdxF", "ftfp", NDF_ftfp_dNdxF, cur_ver );

   out_table << "\n neutron   "; 
   out_table << std::right << std::setw(10) << chi2_qgsp_dNdxF/NDF_qgsp_dNdxF << std::right << std::setw(10) << chi2_ftfp_dNdxF/NDF_ftfp_dNdxF;

   out_table.close();

   return;

}

void NA49_chi2_table_per_spectrum( std::string beam, std::string target, std::string model )
{
   
   std::vector< std::vector<double> > piplus;
   std::vector< std::vector<double> > piminus;
   std::vector< std::vector<double> > proton;
   std::vector< std::vector<double> > antiproton;
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      
      std::string location = "";
      if ( Versions[iv] == CurrentVersion  )
      {
         location = ".";
      }
      else
      {
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[iv];
      }
      std::string histofile = location + "/na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 
      
      TFile* f = TFile::Open( histofile.c_str() );
      
      readDDiffSpectra( beam, target, "piplus" );
      piplus.push_back( std::vector<double>() );
      for ( int ic=0; ic<NSubSetsNA49; ++ic )
      {
         std::string hname = "pTpip" + std::to_string(ic+7);
         TH1F* h = (TH1F*)f->Get( hname.c_str() );
         h->Scale(226.); // or should it be xsec=251. as in Geant4 ???   
         TGraphErrors* gdata = get1DDiffXSecAsGraph( ic );
         int NDF = 0;
	 double chi2 = Chi2( gdata, h, NDF ); 
	 // std::cout << " htitle = " << h->GetTitle() << std::endl;
	 // std::cout << " xF=" << DDiff_xF[ic] << "  chi2=" << chi2 << "  NDF=" << NDF << std::endl; 
	 piplus.back().push_back( (chi2/((double)NDF)) ); 
      }
      
      readDDiffSpectra( beam, target, "piminus" );
      piminus.push_back( std::vector<double>() );
      for ( int ic=0; ic<NSubSetsNA49; ++ic )
      {
         std::string hname = "pTpim" + std::to_string(ic+7);
         TH1F* h = (TH1F*)f->Get( hname.c_str() );
         h->Scale(226.); // or should it be xsec=251. as in Geant4 ???   
         TGraphErrors* gdata = get1DDiffXSecAsGraph( ic );
         int NDF = 0;
	 double chi2 = Chi2( gdata, h, NDF ); 
	 piminus.back().push_back( (chi2/((double)NDF)) ); 
      }
      
      readDDiffSpectra( beam, target, "proton" );
      proton.push_back( std::vector<double>() );
      for ( int ic=0; ic<NSubSetsNA49; ++ic )
      {
         std::string hname = "pTpro" + std::to_string(ic);
         TH1F* h = (TH1F*)f->Get( hname.c_str() );
         h->Scale(226.); // or should it be xsec=251. as in Geant4 ???   
         TGraphErrors* gdata = get1DDiffXSecAsGraph( ic );
         int NDF = 0;
	 double chi2 = Chi2( gdata, h, NDF ); 
	 proton.back().push_back( (chi2/((double)NDF)) ); 
      }
      
      // FIXME !!! Add later !!!
      //
      // histoname will be pTpbar<...>
      //
      // readDDiffSpectra( beam, target, "antiproton" );
      // antiproton.push_back( std::vector<double>() );
            
      f->Close();
   }

   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_per_spectrum_NA49.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(45) << model << ":  158GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(55) << " chi2/NDF calculated vs NA49 data " << std::endl;
   out_table << " " << std::endl;

   out_table << "              ";
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
         out_table << std::right << std::setw(14) << ver;
   }

   out_table << "\n pi+       "; 
   
   readIntegratedSpectra( beam, target, "piplus" );
   
   out_table << "\n dN/dxF     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "dNdxF", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }

   out_table << "\n <pT> vs xF ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piplus", "pT", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }
   
   readDDiffSpectra( beam, target, "piplus" );   
   for ( int i=0; i<NSubSetsNA49; ++i )
   {
      out_table << "\n xF=" << std::right << std::setw(6) << DDiff_xF[i] << " ";      
      for ( int iv=0; iv<NVersions; ++iv )
      {
	 out_table << std::right << std::setw(13) << (piplus[iv])[i];
      }
   } 

   out_table << "\n pi-       "; 
   
   readIntegratedSpectra( beam, target, "piminus" );
   
   out_table << "\n dN/dxF     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "dNdxF", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }

   out_table << "\n <pT> vs xF ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "piminus", "pT", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }

   readDDiffSpectra( beam, target, "piminus" );   
   for ( int i=0; i<NSubSetsNA49; ++i )
   {
      out_table << "\n xF=" << std::right << std::setw(6) << DDiff_xF[i] << " ";      
      for ( int iv=0; iv<NVersions; ++iv )
      {
	 out_table << std::right << std::setw(13) << (piminus[iv])[i];
      }
   } 

   out_table << "\n proton    "; 
   
   readIntegratedSpectra( beam, target, "proton" );
   
   out_table << "\n dN/dxF     ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "proton", "dNdxF", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }

   out_table << "\n <pT> vs xF ";
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, "proton", "pT", model, NDF, Versions[iv] );
      out_table << std::right << std::setw(13) << chi2/NDF ;      
   }

   readDDiffSpectra( beam, target, "proton" );   
   for ( int i=0; i<NSubSetsNA49; ++i )
   {
      out_table << "\n xF=" << std::right << std::setw(6) << DDiff_xF[i] << " ";      
      for ( int iv=0; iv<NVersions; ++iv )
      {
	 out_table << std::right << std::setw(13) << (proton[iv])[i];
      }
   } 

   return;

}
