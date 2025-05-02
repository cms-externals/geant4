#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
// #include "TList.h"

void crudeMerge( std::string beam, std::string target, std::string energy, std::string model )
{

// Note: if one merges weighted/normilized histograms, in this case 
//       the resulting histogram(s) will be 32 times the statistics;
//       at the analysis stage they'd need to be scaled by 1/32 

   TFileMerger* fm = new TFileMerger();
   
   int nfiles = 32;
   
   std::string output = beam + target + model + energy + "GeV.root";
   
   fm->OutputFile( output.c_str() );
   
   for ( int id=1; id<=nfiles; id++ )
   {
      std::string filename = beam + target + model + energy + "GeV-";
      char buf[5];
      sprintf( buf, "%i", id );
      filename.append( buf ); 
      filename += ".root";    
      // std::cout << " file name = " << file << std::endl;           
      fm->AddFile( filename.c_str() );
   }
   
   fm->Merge();
   
   return;

}

void fancyMerge( std::string beam, std::string target, std::string energy, std::string model, std::string dir, bool doScale=true )
{
      
   TH1::SetDefaultSumw2();
   
   std::string output = beam + target + energy + "GeV" + model + ".root" ;
   
   TFile* targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   int nfiles = 32;
//   int nfiles = 2;
   double scale = 1./((double)nfiles);
      
   //
   // NOTE (JVY): PBS job/task/array numbering starts at 1, while SLURM starts counting at 0
   //
   // std::string input = dir + "/" + beam + target + model + energy + "GeV-1.root";
   std::string input = dir + "/" + beam + target + model + energy + "GeV-0.root";
   
   TFile*    iFile1 = TFile::Open( input.c_str() );
   TIter     next( iFile1->GetListOfKeys() );
   TKey*     key = (TKey*)next();
   TH1*      h     = 0;
   TProfile* hprof = 0;
//    TList* list = new TList;
//   list->Clear();
   while ( key )
   {   
         if ( !(TClass::GetClass(key->GetClassName())->InheritsFrom(TH1::Class())) ) continue;
         const char* kName = key->GetName();
	 h = (TH1*)key->ReadObj();
         const char* hName = h->GetName();
         std::cout << " histoname = " << hName << std::endl;
	 TH1F* h1 = (TH1F*)h->Clone();
	 if ( h1->GetSumw2() == 0 ) h1->Sumw2();
// NOTE (JVY): PBS job/task/array count from 1 to NN, while SLURM counts from 0 to NN-1
//	 for ( int id=2; id<=nfiles; id++ )
	 for ( int id=2; id<nfiles; id++ )
	 {
	    std::string input_t = dir + "/" + beam + target + model + energy + "GeV-" ;
            char buf[5];
            sprintf( buf, "%i", id );
            input_t.append( buf ); 
            input_t += ".root"; 
	    TFile* iFile_t = TFile::Open( input_t.c_str() );
	    TH1F* h_t = (TH1F*)iFile_t->Get( h->GetName() );
	    if ( h_t->GetEntries() <= 0 ) continue;
	    h1->Add( h_t );  
	    iFile_t->Close("R");
	 }
	 if ( doScale )
	 {
	    if ( (strcmp(key->GetClassName(),"TProfile")) != 0 ) // 0 means that the 2 strings are equal (no diff)
	    {
	       h1->Scale( scale ); // do I need to Scale(scale,"width") ???
	    }
	 }
	 targetFile->cd();
	 h1->Write();
         key = (TKey*)next();
   }
   
   targetFile->Close();
     
   return;

}
