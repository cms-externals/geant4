void plotG()
{
  gROOT->Reset();
  TFile f = TFile("jeff2.root");    
  TCanvas* c1 = new TCanvas("c1", "  ");
  c1->SetLogy(1);
  c1->Update();   
  
  TH1D* hist2 = (TH1D*)f.Get("2");
  hist2->Draw("HIST");

  c1->Print("jeff2.png");
  
}  
