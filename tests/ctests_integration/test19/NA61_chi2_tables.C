std::string TEST_NAME="test19";
#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test19/G4MODELS_HighEnergy.h"

#include "../test23/shared-root-macros/ReadNA61Data.C"

#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"

#include <iostream>
#include <fstream>

void NA61_chi2_table_regre( std::string beam, std::string target, std::string model )
{

   std::vector< std::pair<double,int> > piplus;
   std::vector< std::pair<double,int> > piminus;
   std::vector< std::pair<double,int> > kplus;
   std::vector< std::pair<double,int> > kminus;
   std::vector< std::pair<double,int> > k0s;
   std::vector< std::pair<double,int> > lambda;
   std::vector< std::pair<double,int> > proton;
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "piplus", model, NDF, Versions[iv] );
      piplus.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "piminus", model, NDF, Versions[iv] );
      piminus.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "kplus", model, NDF, Versions[iv] );
      kplus.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "kminus", model, NDF, Versions[iv] );
      kminus.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "k0s", model, NDF, Versions[iv] );
      k0s.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "lambda", model, NDF, Versions[iv] );
      lambda.push_back( std::pair<double,int>(chi2,NDF) );
   }   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2NA61IntegralData2015( beam, target, "proton", model, NDF, Versions[iv] );
      proton.push_back( std::pair<double,int>(chi2,NDF) );
   }   

   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_NA61.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(30) << model << ":  31GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(62) << " chi2/NDF calculated vs NA61 data " << std::endl;
   out_table << " " << std::endl;

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
      out_table << std::right << setw(12) << piplus[iv].first/piplus[iv].second ;
   }   
   out_table << "\n pi-       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << piminus[iv].first/piminus[iv].second ;
   }   
   out_table << "\n K+        "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << kplus[iv].first/kplus[iv].second ;
   }   
   out_table << "\n K-        "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << kminus[iv].first/kminus[iv].second ;
   }   
   out_table << "\n K0s       "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << k0s[iv].first/k0s[iv].second ;
   }   
   out_table << "\n Lambda    "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << lambda[iv].first/lambda[iv].second ;
   }   
   out_table << "\n proton    "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << proton[iv].first/proton[iv].second ;
   }   

   out_table.close();
   
   return;

}

void NA61_chi2_table_models( std::string beam, std::string target, std::string cur_ver="." )
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
   std::string table_name = "models_" + beam + "_" + target + "_NA61.txt"; 

   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );

   std::string ver = cur_ver;
   if ( cur_ver == "." ) ver = std::string( gSystem->Getenv("G4RELEASE") );
   out_table << std::right << std::setw(45) << "  QGSP vs FTFP as of " << ver << std::endl;
   out_table << std::right << std::setw(35) << "31 GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(55) << " chi2/NDF calculated vs NA61 data " << std::endl;
   out_table << " " << std::endl;

   out_table << " " << std::endl;

   out_table << "           " << std::right << std::setw(12) << "QGSP" << std::right << std::setw(12) << "FTFP";
   
   int NDF = 0;
   double chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "piplus", "qgsp", NDF, cur_ver );
   NDF = 0;
   double chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "piplus", "ftfp", NDF, cur_ver );
   
   out_table << "\n pi+       "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 
   
   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "piminus", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "piminus", "ftfp", NDF, cur_ver );
   
   out_table << "\n pi-       "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "kplus", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "kplus", "ftfp", NDF, cur_ver );
   
   out_table << "\n K+        "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "kminus", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "kminus", "ftfp", NDF, cur_ver );
   
   out_table << "\n K-        "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "k0s", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "k0s", "ftfp", NDF, cur_ver );
   
   out_table << "\n K0s       "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "lambda", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "lambda", "ftfp", NDF, cur_ver );
   
   out_table << "\n Lambda    "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   NDF = 0;
   chi2_qgsp = Chi2NA61IntegralData2015( beam, target, "proton", "qgsp", NDF, cur_ver );
   NDF = 0;
   chi2_ftfp = Chi2NA61IntegralData2015( beam, target, "proton", "ftfp", NDF, cur_ver );
   
   out_table << "\n proton    "; 
   out_table << std::right << std::setw(12) << chi2_qgsp/NDF << std::right << std::setw(12) << chi2_ftfp/NDF; 

   return;

}

void NA61_chi2_table_per_bin_regre( std::string beam, std::string target, std::string model )
{

   // the idea:
   // 
   //               10.x1.y1    10.x2.y2 ...........
   // pi+
   //      0.-10.   chi2/ndf    chi2/ndf ...........
   //     10.-20.   chi2/ndf .......................
   
   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_per_bin_NA61.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(30) << model << ":  31GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(62) << " chi2/NDF calculated vs NA61 data " << std::endl;
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
   
   std::vector< std::pair<std::string,std::string> > theta;
   
   // pi+/pi-
   //
   theta.push_back( std::pair<std::string,std::string>("0","10")  );
   theta.push_back( std::pair<std::string,std::string>("10","20")  );
   theta.push_back( std::pair<std::string,std::string>("20","40")  );
   theta.push_back( std::pair<std::string,std::string>("40","60")  );
   theta.push_back( std::pair<std::string,std::string>("60","100")  );
   theta.push_back( std::pair<std::string,std::string>("100","140")  );
   theta.push_back( std::pair<std::string,std::string>("140","180")  );
   theta.push_back( std::pair<std::string,std::string>("180","240")  );
   theta.push_back( std::pair<std::string,std::string>("240","300")  );
   theta.push_back( std::pair<std::string,std::string>("300","360")  );
   theta.push_back( std::pair<std::string,std::string>("360","420")  );
   
   out_table << "\n pi+       "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "piplus", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "piplus", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }
   
   out_table << "\n pi-       "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "piminus", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "piminus", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }

   // proton 
   //
   // remove just one bin at the end (360-420)
   //   
   theta.pop_back();  

   out_table << "\n proton    "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "proton", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "proton", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }
   
   // K+/K-
   //
   
   theta.clear();
   theta.push_back( std::pair<std::string,std::string>("0","20") );
   theta.push_back( std::pair<std::string,std::string>("20","40") );
   theta.push_back( std::pair<std::string,std::string>("40","60") );
   theta.push_back( std::pair<std::string,std::string>("60","100") );
   theta.push_back( std::pair<std::string,std::string>("100","140") );
   theta.push_back( std::pair<std::string,std::string>("140","180") );
   theta.push_back( std::pair<std::string,std::string>("180","240") );
   theta.push_back( std::pair<std::string,std::string>("240","300") );

   out_table << "\n K+        "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "kplus", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "kplus", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }

   out_table << "\n K-        "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "kminus", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "kminus", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }

   // K0s/Lambda
   //
   
   theta.clear();
   theta.push_back( std::pair<std::string,std::string>("0","40") );
   theta.push_back( std::pair<std::string,std::string>("40","60") );
   theta.push_back( std::pair<std::string,std::string>("60","100") );
   theta.push_back( std::pair<std::string,std::string>("100","140") );
   theta.push_back( std::pair<std::string,std::string>("140","180") );
   theta.push_back( std::pair<std::string,std::string>("180","240") );
   theta.push_back( std::pair<std::string,std::string>("240","300") );

   out_table << "\n K0s       "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "k0s", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "k0s", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }

   out_table << "\n Lambda    "; 
   for ( size_t i=0; i<theta.size(); ++i )
   {
      out_table << "\n      " << std::right << std::setw(3) << theta[i].first << " " << std::right << std::setw(3) << theta[i].second << " ";
      for ( int iv=0; iv<NVersions; iv++ )
      {
         readMomSpectrum( "proton", "C", "lambda", theta[i].first, theta[i].second, "../test23/na61-exp-data/Pub-2015/" );
         int NDF = 0;
         double chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "lambda", 
                                                    theta[i].first, theta[i].second,
				                    model,
					            NDF, Versions[iv], true );
         out_table << std::right << std::setw(13) << chi2/NDF ;
      } 
      out_table << " " << std::endl;
   }

   return;

}
