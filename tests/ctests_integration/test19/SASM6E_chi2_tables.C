
#include "../test19/SASM6E.C"

double Chi2SASM6E( std::string beam, std::string target, std::string secondary,
                   std::string model,
		   int& NDF,
		   std::string version = "." )
{

   double chi2 = 0.;

   std::string location = "";
   if ( version == CurrentVersion || version == "." )
   {
         location = ".";
   }
   else
   {
//         location = regre_test_dir + "/test19/" + version;
         location = regre_test_dir + "/" + TEST_NAME + "/" + version;
   }
   std::string histofile = location + "/sasm6e-histo/" + beam + target + "100.0GeV" + model + ".root"; 

   TFile* f = TFile::Open( histofile.c_str() );
   
   GetExData( beam, target, secondary );
   
   for ( int ipt=0; ipt<NPt; ++ipt )
   {
      std::string hname = secondary + "_" + hlabel[ipt];
      TH1F* h = (TH1F*)f->Get( hname.c_str() );
      TGraphErrors* gr = new TGraphErrors( NData[ipt], XX[ipt], Data[ipt], 0, TotalErr[ipt] ); 
      chi2 += Chi2( gr, h, NDF );
      delete gr;
   }
   
   return chi2;

}

void SASM6E_chi2_table_regre( std::string beam, std::string target, std::string model )
{

   std::vector< std::pair<double,int> > piplus;
   std::vector< std::pair<double,int> > piminus;
   std::vector< std::pair<double,int> > kplus;
   std::vector< std::pair<double,int> > kminus;
   std::vector< std::pair<double,int> > proton;
   std::vector< std::pair<double,int> > antiproton;
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "piplus", model, NDF, Versions[iv] ); 
      piplus.push_back( std::pair<double,int>(chi2,NDF) );
   }
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "piminus", model, NDF, Versions[iv] ); 
      piminus.push_back( std::pair<double,int>(chi2,NDF) );
   }
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "kplus", model, NDF, Versions[iv] ); 
      kplus.push_back( std::pair<double,int>(chi2,NDF) );
   }
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "kminus", model, NDF, Versions[iv] ); 
      kminus.push_back( std::pair<double,int>(chi2,NDF) );
   }
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "proton", model, NDF, Versions[iv] ); 
      proton.push_back( std::pair<double,int>(chi2,NDF) );
   }
   for ( int iv=0; iv<NVersions; ++iv )
   {
      int NDF = 0;
      double chi2 = Chi2SASM6E( beam, target, "antiproton", model, NDF, Versions[iv] ); 
      antiproton.push_back( std::pair<double,int>(chi2,NDF) );
   }

   std::ofstream out_table;
   std::string table_name = model +"_" + beam + "_" + target + "_SASM6E.txt"; 
   
   std::cout << " Opening ascii file: " << table_name << std::endl;

   out_table.open( table_name.c_str() );
     
   out_table << std::right << std::setw(30) << model << ":  100GeV/c " << beam << " on " << target << " --> hadrons" << std::endl;
   out_table << std::right << std::setw(62) << " chi2/NDF calculated vs SASM6E data " << std::endl;
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
   out_table << "\n proton    "; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << proton[iv].first/proton[iv].second ;
   }   
   out_table << "\n antiproton"; 
   for ( int iv=0; iv<NVersions; iv++ )
   {
      out_table << std::right << setw(12) << antiproton[iv].first/antiproton[iv].second ;
   }   

   out_table.close();   
   
   return;

}
