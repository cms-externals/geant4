#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
// #include <list>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
// #include "TStyle.h"

void ReScaleNA49DDiff( std::string file_in )
{

   TFile* finput = TFile::Open( file_in.c_str(), "READ" );
   
   // std::string file_out = file_in + "-ReScaleNA49DDiff";
   size_t pos = file_in.find(".root");
   std::string file_out = file_in.substr( 0, pos );
   file_out += "-ReScaleNA49DDiff.root" ;
   
   TFile* foutput = TFile::Open( file_out.c_str(), "RECREATE" );
   
   TIter     next( finput->GetListOfKeys() );
   TKey*     key = (TKey*)next();
   TH1*      h     = 0;

   while ( key )
   {   
         if ( !(TClass::GetClass(key->GetClassName())->InheritsFrom(TH1::Class())) ) continue;
         const char* kName = key->GetName();
	 h = (TH1*)key->ReadObj();
         // const char* hName = h->GetName();
	 std::string hName( h->GetName() );
	 TH1D* h1 = (TH1D*)h->Clone();
	 if ( h1->GetSumw2() == 0 ) h1->Sumw2();
	 if ( hName.find("pTp") != std::string::npos )
	 {
	    h1->Scale(226.);
	 }
	 foutput->cd();
	 h1->Write();
	 key = (TKey*)next();
   }
   
   foutput->Close();
   finput->Close();
   
   return;

}
