//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TRandom3.h"

#include <iostream>
   
void test15_plots_root6()
{

// qbbc_nokiller_10k_ref07.root
// qbbc_10k_ref07.root
// bic_10k_ref07.root
// bert_10k_ref07.root

// qbbc_nokiller100_2.5.root
// qbbc_nokiller100_3.5.root
// qbbc100_2.5.root
// qbbc100_3.5.root
// bic100_2.5.root
// bic100_3.5.root
// bert100_2.5.root
// bert100_3.5.root

  TFile * histo_file[1];

  TString names[1] = {"Test15_output"};

  for(Long_t j=0; j<1; ++j) {
    TString file_name = names[j] + ".root";
    std::cout << " Rootfile is : " << file_name  << std::endl;
    histo_file[j] = TFile::Open(file_name,"r");

  // // histo_file[0] = TFile::Open("Test15_2.5_100.root","r");
  // // histo_file[0] = TFile::Open("Test15_new.root","r");
  // // histo_file[0] = TFile::Open("tarc_working.root","r");
  // // histo_file[1] = TFile::Open("Test15_3.5_100.root","r");
  // histo_file[0] = TFile::Open("Test15_ouput.root","r");
  // // histo_file[0] = TFile::Open("bic_10k_ref07.root","r");
  // histo_file[1] = TFile::Open("bert_10k_ref07.root","r");
  // histo_file[2] = TFile::Open("qbbc_10k_ref07.root","r");
  // histo_file[3] = TFile::Open("qbbc_nokiller_10k_ref07.root","r");

    TH1F * exiting_spectrum;
    TH2F *radial_h1;
    TH1F *tarc_data_fluence_high;
    TH1F *tarc_data_fluence_high_error;
    TH1F *tarc_g4_fluence_high;
    TH1F *tarc_g4_ratio_high;
    TH1F *tarc_data_fluence_he3;
    TH1F *tarc_data_fluence_he3_error;
    TH1F *tarc_g4_fluence_he3;
    TH1F *tarc_g4_ratio_he3;
    TH1F *tarc_data_fluence_li;
    TH1F *tarc_data_fluence_li_error;
    TH1F *tarc_g4_fluence_li;
    TH1F *tarc_g4_ratio_li;
    
  // for(int j=0; j<2; ++j) {

    TNtuple* exiting_tuple = (TNtuple*) histo_file[j]->Get("Test15 Exiting");
    
    TH2 * h1 = (TH2*)histo_file[j]->Get("1");
    TH2 * h2 = (TH2*)histo_file[j]->Get("2");
    TNtuple* my_tupleD = (TNtuple*) histo_file[j]->Get("Test15 Flux 4002");
    Int_t count3 = my_tupleD->GetEntries();
    std::cout << " Count: " << count3 << std::endl;

    TNtuple* my_tupleD2 = (TNtuple*) histo_file[j]->Get("Test15 Flux 4004");
    Int_t count4 = my_tupleD2->GetEntries();
    std::cout << " Count4: " << count4 << std::endl;

    TNtuple* my_tupleD3 = (TNtuple*) histo_file[j]->Get("Test15 Flux 4005");
    Int_t count5 = my_tupleD3->GetEntries();
    std::cout << " Count5: " << count5 << std::endl;

    TNtuple* my_tupleD_rad = (TNtuple*) histo_file[j]->Get("Test15 Radial Shell Fluence");
    TNtuple* my_tupleD_rad_data = (TNtuple*) histo_file[j]->Get("Test15 Radial Fluence Data");
    TNtuple* my_tupleD_rad_he3 = (TNtuple*) histo_file[j]->Get("Test15 Radial Fluence He3");
    
    // std::cout << " double ntuple: " << my_tupleD->GetEntries() << std::endl;

     //  data xbins /59500.,109000.,158500.,208000.,257500.,307000.,356500.
     // +     ,406000.,455500.,505000.,554500.,604000.,653500.,703000.
     // +     ,752500.,802000.,901000.,1000000.,1162308.,1350960.,1570232.
     // +     ,1825092./

    TTree* my_tuple9 = (TTree*) histo_file[j]->Get("Test15 Radial Shell Fluence");

    radial_h1 = new TH2F(TString("radial_h1_")+j,"Radial in Test15",1000,-200.,200.,1000,1.,25000000);

    Float_t xbins_high[22] = {59500, 109000, 158500, 208000, 257500, 307000, 356500, 406000, 455500, 505000, 554500, 604000, 653500, 703000, 752500, 802000, 901000, 1000000, 1162308, 1350960, 1570232, 1825092};

    Float_t xbins_low[102];

    xbins_low[0] = 0.01;

    Float_t bin_width = (std::log(pow(10.0,5.0))-std::log(0.01))/100.0;
    // bin_width = 0.1611809
    std::cout << " bin width is: " << bin_width << " should be 0.1611809 " << std::endl;
    std::cout << " now press enter " << std::endl;
    // getchar();

    int index = 0;
    for(int i=1; i<102; ++i) {
      xbins_low[i] = std::exp(bin_width+std::log(xbins_low[i-1]));
    }

    std::cout << " xbins_low max: " << xbins_low[101] << " xbins_high min: " << xbins_high[0] << std::endl;
     // 	call hbookb(20000,' Test15 Data Fluence ',21,xbins,0.)
     // 	call hbookb(20001,' 4pi Geant4 Shell Fluence ',21,xbins,0.)

    // tarc_data_fluence[j] = new TH1F(TString("fluence_data_")+j,"Fluence Data in Test15",123,xbins);
    tarc_data_fluence_high = new TH1F(TString("fluence_data_")+j,"Fluence Data in Test15 High",21,xbins_high);
    tarc_data_fluence_he3 = new TH1F(TString("fluence_data_")+j,"Fluence Data in Test15 He3",101,xbins_low);
    tarc_data_fluence_li = new TH1F(TString("fluence_data_li_")+j,"Fluence Data in Test15 Li",101,xbins_low);
    tarc_data_fluence_high_error = new TH1F(TString("fluence_data_error_")+j,"Fluence Data in Test15 High Combined Error",21,xbins_high);
    tarc_data_fluence_he3_error = new TH1F(TString("fluence_data_error_")+j,"Fluence Data in Test15 He3 Combined Error",101,xbins_low);
    tarc_data_fluence_li_error = new TH1F(TString("fluence_data_li_error_")+j,"Fluence Data in Test15 Li Combined Error",101,xbins_low);

    tarc_g4_fluence_high = new TH1F(TString("fluence_g4_high_")+j,"Fluence G4 in Test15",21,xbins_high);
    tarc_g4_fluence_he3 = new TH1F(TString("fluence_g4_he3_")+j,"Fluence G4 in Test15",101,xbins_low);
    tarc_g4_fluence_li = new TH1F(TString("fluence_g4_li_")+j,"Fluence G4 in Test15",101,xbins_low);

    tarc_g4_ratio_high = new TH1F(TString("ratio_g4_high_")+j,"Fluence Ratio G4 in Test15",21,xbins_high);
    tarc_g4_ratio_he3 = new TH1F(TString("ratio_g4_he3_")+j,"Fluence Ratio G4 in Test15",101,xbins_low);
    tarc_g4_ratio_li = new TH1F(TString("ratio_g4_li_")+j,"Fluence Ratio G4 in Test15",101,xbins_low);

  // analysisManager->CreateNtupleDColumn("radius");
  // analysisManager->CreateNtupleDColumn("energy");
  // analysisManager->CreateNtupleDColumn("fluence");
  // analysisManager->CreateNtupleDColumn("true_e");
  // analysisManager->CreateNtupleDColumn("true_f");
    // double radius, energy, fluence, true_e, true_f;

  // analysisManager->CreateNtuple("Test15 Flux 4002", "Neutrons Test15 flux"); 
  // analysisManager->CreateNtupleDColumn("energy");
  // analysisManager->CreateNtupleDColumn("tarcflux");
  // analysisManager->CreateNtupleDColumn("errstat");
  // analysisManager->CreateNtupleDColumn("errsyst");
  // analysisManager->CreateNtupleDColumn("g4flux");
  // analysisManager->CreateNtupleDColumn("g4perp");
  // analysisManager->CreateNtupleDColumn("gfluence");
  // analysisManager->CreateNtupleDColumn("g4err");
  // analysisManager->CreateNtupleDColumn("rawflux");
  // analysisManager->CreateNtupleDColumn("trceflux");
  // analysisManager->CreateNtupleDColumn("g4eflux");
  // analysisManager->CreateNtupleDColumn("gstep");
  // analysisManager->CreateNtupleDColumn("gfl_cyl");
  // analysisManager->CreateNtupleDColumn("g4front");
  // analysisManager->CreateNtupleDColumn("g4_shell");

    double* row_content;
    std::cout << " got here - before my_tuple access " << std::endl;
    // getchar();
    std::cout << " Entries: " << my_tuple9->GetEntries() << std::endl;
    std::cout << " Entries: " << my_tupleD->GetEntries() << std::endl;

    double energy, tarcflux, errstat, errsyst, g4flux, g4perp, fluence, g4err, rawflux, trceflux, g4eflux, gstep, gfl_cyl, g4front, g4_shell;

    exiting_spectrum = new TH1F(TString("exiting_spectrum_")+j,"Exiting Neutron Spectrum Test15",100,0.01,2.);
    for (int irow=0;irow<exiting_tuple->GetEntries();++irow){
      exiting_tuple->SetBranchAddress("energy",&energy);
      exiting_tuple->GetEntry(irow);
      exiting_spectrum->Fill(energy);
    }

    for (int irow=0;irow<my_tupleD->GetEntries();++irow){
      // std::cout << " got here 999999 " << std::endl;
      // getchar();
      // my_tupleD->SetBranchAddress("energy",&row_content); 
      double * output;
      my_tupleD->SetBranchAddress("energy",&energy);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD->GetEntry(irow);
      my_tupleD->SetBranchAddress("tarcflux",&tarcflux);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      std::cout << " tarcflux: " << tarcflux << " for row: " << irow << std::endl;
      my_tupleD->GetEntry(irow);
      my_tupleD->SetBranchAddress("g4_shell",&g4_shell);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD->SetBranchAddress("g4err",&g4err);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
	// my_tupleD->SetBranchAddress("g4perp",&g4perp);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD->GetEntry(irow);
      my_tupleD->SetBranchAddress("errstat",&errstat);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD->GetEntry(irow);
      my_tupleD->SetBranchAddress("errsyst",&errsyst);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD->GetEntry(irow);
	// radial_h1[j]->Fill(radius/10,fluence/energy,1.);
	// tarc_data_fluence[j]->Fill(energy,tarcflux);
	// // radial_h1[j]->Fill(radius/10,fluence/energy,1.);

      tarc_data_fluence_high->Fill(energy,tarcflux);
      TAxis *xaxis = tarc_data_fluence_high->GetXaxis();
      Int_t binx = xaxis->FindBin(energy);
      tarc_data_fluence_high->SetBinError(binx,errstat);
      std::cout << " energy: " << energy << " tarcflux: " << tarcflux << " error: " << errstat << " errsyst: " << errsyst << std::endl; 
      
      tarc_data_fluence_high_error->Fill(energy,tarcflux);
      TAxis *xaxis2 = tarc_data_fluence_high_error->GetXaxis();
      Int_t binx2 = xaxis2->FindBin(energy);
      tarc_data_fluence_high_error->SetBinError(binx2,errstat+errsyst);
      
      
      double corr_g4perp = 1.0*g4_shell;
      tarc_g4_fluence_high->Fill(energy,corr_g4perp);
      xaxis = tarc_g4_fluence_high->GetXaxis();
      binx = xaxis->FindBin(energy);
      tarc_g4_fluence_high->SetBinError(binx,g4err);
      if(tarcflux != 0.) {
	double ratio = corr_g4perp/tarcflux;
	tarc_g4_ratio_high->Fill(energy,ratio);
      }
    }

    for (int irow=0;irow<my_tupleD2->GetEntries();++irow){
      // std::cout << " got here 999999 " << std::endl;
      // getchar();
      // my_tupleD->SetBranchAddress("energy",&row_content); 
      double * output;
      my_tupleD2->SetBranchAddress("energy",&energy);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD2->GetEntry(irow);
      my_tupleD2->SetBranchAddress("tarcflux",&tarcflux);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      std::cout << " tarcflux: " << tarcflux << " for row: " << irow << std::endl;
      my_tupleD2->GetEntry(irow);
      my_tupleD2->SetBranchAddress("g4_shell",&g4_shell);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD2->GetEntry(irow);
      my_tupleD2->SetBranchAddress("errstat",&errstat);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD2->GetEntry(irow);
      my_tupleD2->SetBranchAddress("errsyst",&errsyst);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD2->GetEntry(irow);
      // radial_h1->Fill(radius/10,fluence/energy,1.);
      // tarc_data_fluence->Fill(energy,tarcflux);
      // // radial_h1->Fill(radius/10,fluence/energy,1.);
      
      tarc_data_fluence_he3->Fill(energy,tarcflux);
      TAxis *xaxis = tarc_data_fluence_he3->GetXaxis();
      Int_t binx = xaxis->FindBin(energy);
      tarc_data_fluence_he3->SetBinError(binx,errstat);
      
      tarc_data_fluence_he3_error->Fill(energy,tarcflux);
      TAxis *xaxis2 = tarc_data_fluence_he3_error->GetXaxis();
      Int_t binx2 = xaxis2->FindBin(energy);
      tarc_data_fluence_he3_error->SetBinError(binx2,errsyst+errstat);
      
      double corr_g4perp = 1.0*g4_shell;
      tarc_g4_fluence_he3->Fill(energy,corr_g4perp);
      xaxis = tarc_g4_fluence_he3->GetXaxis();
      binx = xaxis->FindBin(energy);
      tarc_g4_fluence_he3->SetBinError(binx,g4err);
      if(tarcflux != 0.) {
	double ratio = corr_g4perp/tarcflux;
	tarc_g4_ratio_he3->Fill(energy,ratio);
      }
    }
    
    for (int irow=0;irow<my_tupleD3->GetEntries();++irow){
      // std::cout << " got here 999999 " << std::endl;
      // getchar();
      // my_tupleD->SetBranchAddress("energy",&row_content); 
      double * output;
      my_tupleD3->SetBranchAddress("energy",&energy);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD3->GetEntry(irow);
      my_tupleD3->SetBranchAddress("tarcflux",&tarcflux);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      std::cout << " tarcflux: " << tarcflux << " for row: " << irow << std::endl;
      my_tupleD3->GetEntry(irow);
      my_tupleD3->SetBranchAddress("g4_shell",&g4_shell);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD3->GetEntry(irow);
      my_tupleD3->SetBranchAddress("errstat",&errstat);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD3->GetEntry(irow);
      my_tupleD3->SetBranchAddress("errsyst",&errsyst);//,"tarcflux",&tarcflux);//,"g4_shell",&g4_shell); 
      my_tupleD3->GetEntry(irow);
      // radial_h1->Fill(radius/10,fluence/energy,1.);
      // tarc_data_fluence->Fill(energy,tarcflux);
      // // radial_h1->Fill(radius/10,fluence/energy,1.);
      
      tarc_data_fluence_li->Fill(energy,tarcflux);
      TAxis *xaxis = tarc_data_fluence_li->GetXaxis();
      Int_t binx = xaxis->FindBin(energy);
      tarc_data_fluence_li->SetBinError(binx,errstat);
      
      tarc_data_fluence_li_error->Fill(energy,tarcflux);
      TAxis *xaxis2 = tarc_data_fluence_li_error->GetXaxis();
      Int_t binx2 = xaxis2->FindBin(energy);
      tarc_data_fluence_li_error->SetBinError(binx2,errsyst+errstat);
      std::cout << " lithium fluence: " << tarcflux << std::endl;
      
      double corr_g4perp = 1.0*g4_shell;
      tarc_g4_fluence_li->Fill(energy,corr_g4perp);
      xaxis = tarc_g4_fluence_li->GetXaxis();
      binx = xaxis->FindBin(energy);
      tarc_g4_fluence_li->SetBinError(binx,g4err);
      if(tarcflux != 0.) {
	double ratio = corr_g4perp/tarcflux;
	tarc_g4_ratio_li->Fill(energy,ratio);
      }
    }
    
    TCanvas* c2 = new TCanvas("c2","Test15 Fluence summary",700,500);
    gStyle->SetHistLineWidth(3);  // Default:1(GetHistLineWidth) --- Line width of histogram
    gStyle->SetLineWidth(0.3);        // Default:1 --- Line width of axis.
    gStyle->SetTitleX(0.2);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->DrawFrame(0.001,1000,10000000,25000000,"; Energy/eV; EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
    TString my_title2 = names[j] + " Gev/c protons 100 events ";
    // TLatex* tlx=new TLatex(0.23, 0.93, "100 events binary 2.5 GeV/c proton");
    TLatex* tlx=new TLatex(0.23, 0.93, my_title2);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.05);
    tlx->Draw();
    
    tarc_data_fluence_high->SetLineColor(kRed);
    tarc_data_fluence_high->SetLineWidth(0.85);
    tarc_data_fluence_he3->SetLineColor(kBlue);
    tarc_data_fluence_he3->SetLineWidth(0.85);
    tarc_data_fluence_li->SetLineColor(kCyan);
    tarc_data_fluence_li->SetLineWidth(0.85);
    tarc_data_fluence_high->Draw("SAME E1");
    tarc_data_fluence_he3->Draw("SAME E1");
    tarc_data_fluence_li->Draw("SAME E1");
    
    tarc_data_fluence_high_error->SetLineColor(kGreen);
    tarc_data_fluence_high_error->SetLineWidth(0.7);
    tarc_data_fluence_he3_error->SetLineColor(kYellow);
    tarc_data_fluence_he3_error->SetLineWidth(0.7);
    tarc_data_fluence_li_error->SetLineColor(kMagenta);
    tarc_data_fluence_li_error->SetLineWidth(0.7);
    tarc_data_fluence_high_error->Draw("SAME E");
    tarc_data_fluence_he3_error->Draw("SAME E");
    tarc_data_fluence_li_error->Draw("SAME E");
    
    // gStyle->OptStat(0);
    tarc_g4_fluence_high->SetLineColor(kBlack);
    // tarc_g4_fluence_high->SetLineWidth(2);
    tarc_g4_fluence_high->SetStats(kFALSE);
    tarc_g4_fluence_high->Draw("SAME E1");
    tarc_g4_fluence_he3->SetLineColor(kBlue);
    // tarc_g4_fluence_he3->SetLineWidth(2);
    tarc_g4_fluence_he3->SetStats(kFALSE);
    tarc_g4_fluence_he3->Draw("SAME E1");
    tarc_g4_fluence_li->SetLineColor(kCyan);
    // tarc_g4_fluence_li->SetLineWidth(2);
    tarc_g4_fluence_li->SetStats(kFALSE);
    tarc_g4_fluence_li->Draw("SAME E1");
    
    TString file2 = names[j] + "fluence1.pdf";
    c2->Print(file2);
    TString file2g = names[j] + "fluence1.png";
    c2->Print(file2g);
    c2->Close();
    
    TCanvas* c3 = new TCanvas("c3","Test15 Fluence Ratio G4/Data summary",700,500);
    gStyle->SetHistLineWidth(3);  // Default:1(GetHistLineWidth) --- Line width of histogram
    gStyle->SetLineWidth(0.3);        // Default:1 --- Line width of axis.
    gStyle->SetTitleX(0.2);
    gPad->SetLogx();
    gPad->DrawFrame(0.001,0,10000000,2,"; Energy/eV; G4/Data")->GetXaxis()->SetTitleOffset(1.2);
    TString my_title3 = names[j] + " Gev/c protons 100 events ";
    // TLatex* tlx=new TLatex(0.23, 0.93, "100 events binary 2.5 GeV/c proton");
    tlx=new TLatex(0.23, 0.93, my_title3);
    // TLatex* tlx=new TLatex(0.23, 0.93, "10000 events binary 2.5 GeV/c proton");
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.05);
    tlx->Draw();
    
    tarc_g4_ratio_high->SetMarkerColor(kRed);
    tarc_g4_ratio_high->SetMarkerStyle(21);
    tarc_g4_ratio_high->SetLineColor(kRed);
    tarc_g4_ratio_high->SetLineWidth(0.85);
    tarc_g4_ratio_he3->SetLineColor(kBlue);
    tarc_g4_ratio_he3->SetMarkerStyle(22);
    tarc_g4_ratio_he3->SetMarkerColor(kBlue);
    tarc_g4_ratio_he3->SetLineWidth(0.85);
    tarc_g4_ratio_li->SetLineColor(kCyan);
    tarc_g4_ratio_li->SetMarkerStyle(23);
    tarc_g4_ratio_li->SetMarkerColor(kCyan);
    tarc_g4_ratio_li->SetLineWidth(0.85);
    tarc_g4_ratio_high->Draw("SAME P");
    tarc_g4_ratio_he3->Draw("SAME P");
    tarc_g4_ratio_li->Draw("SAME P");
    
    TString file3 = names[j] + "_ratio3.pdf";
    c3->Print(file3);
    TString file3g = names[j] + "_ratio3.png";
    c3->Print(file3g);
    c3->Close();
    
    // radial:
    TCanvas* c5 = new TCanvas("c2","Test15 Radial Fluence",700,500);
    gStyle->SetHistLineWidth(3);  // Default:1(GetHistLineWidth) --- Line width of histogram
    gStyle->SetLineWidth(0.3);        // Default:1 --- Line width of axis.
    gStyle->SetTitleX(0.2);
    gPad->SetLogy();
    gPad->DrawFrame(-200,1,200,25000000,"; Radial Distance/cm; dF/dE n/cm^{2}/eV/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
    TString my_title4 = names[j] + " Gev/c protons 100 events ";
    // TLatex* tlx=new TLatex(0.23, 0.93, "100 events binary 2.5 GeV/c proton");
    tlx=new TLatex(0.23, 0.93, my_title4);
    // TLatex* tlx=new TLatex(0.23, 0.93, "10000 events binary 2.5 GeV/c proton");
    // TLatex* tlx=new TLatex(0.23, 0.93, "10000 events binary 3.5 GeV/c proton");
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.05);
    tlx->Draw();
    my_tupleD_rad->SetMarkerStyle(21);
    my_tupleD_rad->SetMarkerColor(kRed);
    my_tupleD_rad->SetMarkerSize(0.8);
    my_tupleD_rad->Draw("fluence/energy : radius/10","","SAME");
    my_tupleD_rad_data->SetMarkerStyle(28);
    my_tupleD_rad_data->SetMarkerSize(0.8);
    my_tupleD_rad_data->SetMarkerColor(kBlue);
    my_tupleD_rad_data->Draw("data : radius","energy<1001","SAME");
    
    my_tupleD_rad_he3->SetMarkerStyle(29);
    my_tupleD_rad_he3->SetMarkerSize(0.8);
    my_tupleD_rad_he3->SetMarkerColor(kBlue);
    my_tupleD_rad_he3->Draw("data : radius","energy>1000","SAME");
    
      /*
set mtyp 21
set pmci 2
set mscf 0.7
nt/pl 4999.(fluence/energy)%radius/10 ! ! ! ! s

set mtyp 28
set pmci 4
set mscf 0.9
nt/pl 4998.data%radius energy.lt.1001 ! ! ! s

set mtyp 29
set pmci 4
set mscf 0.9
nt/pl 4997.data%radius energy.gt.1000 ! ! ! s

mytext 0. 500 ' 50 keV ' 0.3 ! c
mytext 0. 1000 ' 10 keV ' 0.3 ! c
mytext 0. 2000 ' 1 keV ' 0.3 ! c
mytext 0. 5000 ' 480 eV ' 0.3 ! c
mytext 0. 11000 ' 100 eV ' 0.3 ! c
mytext 0. 80000 ' 18 eV ' 0.3 ! c
mytext 0. 200000 ' 5 eV ' 0.3 ! c
mytext 0. 600000 ' 1.5 eV ' 0.3 ! c
mytext 0. 1400000 ' 0.1 eV ' 0.3 ! c
      */
      // tarc_data_fluence_high->SetLineColor(kRed);
      // tarc_data_fluence_high->SetLineWidth(0.85);
      // tarc_data_fluence_he3->SetLineColor(kBlue);
      // tarc_data_fluence_he3->SetLineWidth(0.85);
      // tarc_data_fluence_li->SetLineColor(kCyan);
      // tarc_data_fluence_li->SetLineWidth(0.85);
      // tarc_data_fluence_high->Draw("SAME E1");
      // tarc_data_fluence_he3->Draw("SAME E1");
      // tarc_data_fluence_li->Draw("SAME E1");

    TString file5 = names[j] + "_radial1.pdf";
    c5->Print(file5);
    TString file5g = names[j] + "_radial1.png";
    c5->Print(file5g);
    
    c5->Close();
    
  }    

//  exit();

}

