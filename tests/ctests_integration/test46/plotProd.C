
void plotProd()
{ 
  gROOT->Reset();
  gROOT->SetStyle("Plain"); 
  //gStyle->SetOptStat(0);
  gStyle->SetNdivisions(210, "x");
  gStyle->SetNdivisions(10, "y");
  gStyle->SetTextFont(2);
  gStyle->SetLabelOffset(0.005, "x");
  gStyle->SetLabelOffset(0.005, "y");
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetLabelSize(0.03, "y");
  gStyle->SetTickLength(0.05, "x");
  gStyle->SetTickLength(0.03, "y");
  gStyle->SetPadBorderMode(0);

  gStyle->SetMarkerStyle(21); 
  gStyle->SetMarkerColor(1); 
  gStyle->SetMarkerSize(.6);
  gStyle->SetTitleOffset(.95, "x");
  gStyle->SetTitleOffset(.95, "y");
  gStyle->SetTitleSize(.06, "x");
  gStyle->SetTitleSize(.06, "y");
  gStyle->SetPadBottomMargin(.15);          
  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.03);
  gStyle->SetMarkerSize(1.5);

  double x1[4] = {-5.,-5.,-5.,-5.}; 
  double x2[4] = { 5., 5., 5., 5.}; 
  double y1[4] = {0.001,0.001,0.00001,0.00001}; 
  double y2[4] = {10000.,1000.,100.,100.}; 

  TH1F*  hh[4] = {0, 0, 0, 0};

  TCanvas c1("c1","",0, 5, 800, 600);
  c1.Divide(2,2);

  TString hist[4] = {"histo/h16","histo/h17","histo/h18","histo/h19"};

  Int_t col[6] = {4, 3, 2, 2, 13, 13};
  Int_t mar[6] = {20, 21, 22, 25, 26, 26};
  Int_t lin[6] = {2, 1, 2, 1, 1, 1};

  TString tit[3] = {"log_{10} (E/MeV)","","#pi- 50 GeV"};
  TString tit1[4] = {"Gamma","Electron","Proton","Neutron"};

  TLegend* leg[5];

  TString fileR = "test46_50gev.root";
  TFile* ff = new TFile(fileR);

  if(!ff) { 
    cout << "### File <" << fileR << "> is not opened!" << endl;
    exit(1);
  }

  for(int i=0; i<4; ++i) {
    c1.cd(i+1);
    gPad->SetGrid();
    gPad->SetLogy(); 
    hh[i] = gPad->DrawFrame(x1[i],y1[i],x2[i],y2[i]);
    if(!hh[i]) {
      cout << "Fail draw frame " << i << endl;
      exit(1);
    }
    hh[i]->GetXaxis()->SetTitle(tit[0]);
    hh[i]->GetYaxis()->SetTitle(tit[1]);

    hh[i]->Draw("AXIS SAME 9");
  
    TH1F* hhi = (TH1F*)ff->Get(hist[i]);
    if(!hhi) {
      cout << "### Fail access " << hist[i] << " for pad " << i << endl;
      exit(1);
    }
    hhi->SetLineColor(col[0]);
    hhi->SetLineStyle(lin[1]);
    //hhi->SetLineWidth(2);
    hhi->Draw("HISTO SAME 9");

    leg[i] = new TLegend(0.22, 0.70, 0.42, 0.86);
    leg[i]->SetHeader(tit1[i]); 
    leg[i]->Draw("SAME 9");
    if(i == 0) {
      leg[4] = new TLegend(0.68, 0.70, 0.92, 0.86);
      leg[4]->SetHeader(tit[2]); 
      leg[4]->Draw("SAME 9");
    }
  }

  c1.Print("prod_50gev.png");

}
