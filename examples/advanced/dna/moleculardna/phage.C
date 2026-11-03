//-------------------------------------------------------------------------------//
// This macrofile was developed by Konstantinos Chatzipapas at LP2iB (ex. CENBG) //
// in collaboration with the whole team of molecularDNA Geant4-DNA example       //
// For any question please contact through:                                      //
// k.chatzipapas@yahoo.com                                                       //
// dousatsu.sakata@aist.go.jp                                                    //
//-------------------------------------------------------------------------------//
//
// This macro requires the molecular-dna.root file generated from molecularDNA example
// To run this file just insert this command to the terminal:
// root .X phage.C
// ROOT6.x should be installed
//
//************************************//
// Please define the parameters below //
//************************************//

{
//*******************************************************************************//
// Define these parameters of the simulation
char ifile[256] = "molecular-dna.root";  // input filepath
Double_t r3 = 7000e-9 * 3500e-9 * 3500e-9; // a * b * c   // Chromosome size, as defined in the mac file, but in meters. If sphere, a=b=c
Double_t Nbp = 0.141158; // Mbp // Length of the DNA chain in Mbp
Double_t mass = 997 * 3.141592 * r3;  // waterDensity * pi * r3 in kg // cylinder
//*******************************************************************************//

typedef std::pair <int64_t, int64_t> ipair;
bool greaterPair(const ipair &l, const ipair &r);
bool smallerPair(const ipair &l, const ipair &r);

void BinLogX(TH1 *h);

// Standard error of the mean of a sample of size n:
//   s2 = sum of squared values, s = sum of values.
// Two overloads preserve the exact integer vs. floating division of the caller.
Double_t SEM(Float_t s2, Int_t s,   Long64_t n);
Double_t SEM(Double_t s2, Double_t s, Long64_t n);

// Relative error contribution of a yield: sqrt( (n*SD/total)^2 + (SD_dose/dose)^2 ).
Double_t RelErr(Long64_t n, Double_t SD, Double_t total, Double_t SD_dose, Double_t dose);


gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(00000);

// Initialize output histograms
TCanvas *cfragment = new TCanvas("cfragment","DNA Fragments Distribution", 900, 120, 600,400);
cfragment->SetLogx();
cfragment->SetLogy();
TH1F *h1fragments = new TH1F("h1fragments","h1fragments",40,0,5);
BinLogX(h1fragments);

TCanvas *c1 = new TCanvas("c1", "Molecular DNA - Damage Quantification", 60, 120, 800, 800);
c1->SetBorderSize(0);
c1->SetFillColor(0);
c1->SetFillStyle(4000);
gPad->SetLeftMargin(0.13);

TPad* pad1 = new TPad("pad1","Species", 0, 0.51, 0.49, 1);
pad1->SetBorderSize(0);
pad1->SetFillColor(0);
pad1->SetFillStyle(4000);
pad1->SetLeftMargin(0.15);
pad1->SetRightMargin(0.01);
pad1->SetBottomMargin(0.2);

TPad* pad2 = new TPad("pad2","Damage Yield", 0.51, 0.5, 1, 1);
pad2->SetBorderSize(0);
pad2->SetFillColor(0);
pad2->SetFillStyle(4000);
pad2->SetLeftMargin(0.15);
pad2->SetRightMargin(0.05);
pad2->SetBottomMargin(0.2);

TPad* pad3 = new TPad("pad3","Breaks Yield SSB", 0, 0, 0.49, 0.49);
pad3->SetBorderSize(0);
pad3->SetFillColor(0);
pad3->SetFillStyle(4000);
pad3->SetLeftMargin(0.15);
pad3->SetRightMargin(0.01);
//pad3->SetTopMargin(0.2);
pad3->SetBottomMargin(0.2);

TPad* pad4 = new TPad("pad4","Breaks Yield DSB", 0.51, 0, 1, 0.49);
pad4->SetBorderSize(0);
pad4->SetFillColor(0);
pad4->SetFillStyle(4000);
pad4->SetLeftMargin(0.15);
pad4->SetRightMargin(0.05);
//pad3->SetTopMargin(0.2);
pad4->SetBottomMargin(0.2);

pad1->Draw();
pad2->Draw();
pad3->Draw();
pad4->Draw();

// Open root file
TFile *f = TFile::Open(ifile);

// Initialize Variables
Int_t EB, ES, OHB, OHS, HB, HS;
Int_t total_EB, total_ES, total_OHB, total_OHS, total_HB, total_HS;
Float_t total_EB2, total_ES2, total_OHB2, total_OHS2, total_HB2, total_HS2;
Float_t SD_EB, SD_ES, SD_OHB, SD_OHS, SD_HB, SD_HS;
Float_t SD_SSB, SD_SSBp, SD_SSB2p, SD_sSSB, SD_SSBd, SD_SSBi, SD_SSBm;
Float_t SD_DSB, SD_DSBp, SD_DSBpp, SD_sDSB, SD_DSBd, SD_DSBi, SD_DSBm, SD_DSBh;

Int_t SSB, SSBp, SSB2p;
Int_t total_SSB, total_SSBp, total_SSB2p;
Float_t total_SSB2, total_SSBp2, total_SSB2p2;
Int_t DSB, DSBp, DSBpp;
Int_t total_DSB, total_DSBp, total_DSBpp;
Float_t total_DSB2, total_DSBp2, total_DSBpp2;

Int_t SSBd, SSBi, SSBm;
Int_t total_sSSB, total_SSBd, total_SSBi, total_SSBm;
Float_t total_sSSB2, total_SSBd2, total_SSBi2, total_SSBm2;
Int_t DSBd, DSBi, DSBm, DSBh;
Int_t total_sDSB, total_DSBd, total_DSBi, total_DSBm, total_DSBh;
Float_t total_sDSB2, total_DSBd2, total_DSBi2, total_DSBm2, total_DSBh2;

Double_t dose = 0;
Double_t SD_dose = 0;

Double_t EB_yield = 0; Double_t ES_yield = 0; Double_t OHB_yield = 0; Double_t OHS_yield = 0; Double_t HB_yield = 0; Double_t HS_yield = 0;
Double_t SD_EB_yield = 0; Double_t SD_ES_yield = 0; Double_t SD_OHB_yield = 0; Double_t SD_OHS_yield = 0; Double_t SD_HB_yield = 0; Double_t SD_HS_yield = 0;

Double_t SSB_yield = 0; Double_t SSBp_yield = 0; Double_t SSB2p_yield = 0;
Double_t SD_SSB_yield = 0; Double_t SD_SSBp_yield = 0; Double_t SD_SSB2p_yield = 0;
Double_t DSB_yield = 0; Double_t DSBp_yield = 0; Double_t DSBpp_yield = 0;
Double_t SD_DSB_yield = 0; Double_t SD_DSBp_yield = 0; Double_t SD_DSBpp_yield = 0;

Double_t sSSB_yield = 0; Double_t SSBi_yield = 0; Double_t SSBd_yield = 0; Double_t SSBm_yield = 0;
Double_t SD_sSSB_yield = 0; Double_t SD_SSBi_yield = 0; Double_t SD_SSBd_yield = 0; Double_t SD_SSBm_yield = 0;
Double_t sDSB_yield = 0; Double_t DSBi_yield = 0; Double_t DSBd_yield = 0; Double_t DSBm_yield = 0; Double_t DSBh_yield = 0;
Double_t SD_sDSB_yield = 0; Double_t SD_DSBi_yield = 0; Double_t SD_DSBd_yield = 0; Double_t SD_DSBm_yield = 0; Double_t SD_DSBh_yield = 0;

total_EB  = 0; total_ES = 0; total_OHB = 0; total_OHS = 0; total_HB = 0; total_HS = 0;
total_EB2 = 0; total_ES2 = 0; total_OHB2 = 0; total_OHS2 = 0; total_HB2 = 0; total_HS2 = 0;

total_SSB  = 0; total_SSBp = 0; total_SSB2p = 0;
total_SSB2  = 0; total_SSBp2 = 0; total_SSB2p2 = 0;
total_DSB  = 0; total_DSBp = 0; total_DSBpp = 0;
total_DSB2  = 0; total_DSBp2 = 0; total_DSBpp2 = 0;

total_sSSB = 0; total_SSBd = 0; total_SSBi = 0; total_SSBm = 0;
total_sSSB2 = 0; total_SSBd2 = 0; total_SSBi2 = 0; total_SSBm2 = 0;
total_sDSB = 0; total_DSBd = 0; total_DSBi = 0; total_DSBm = 0; total_DSBh = 0;
total_sDSB2 = 0; total_DSBd2 = 0; total_DSBi2 = 0; total_DSBm2 = 0;  total_DSBh2 = 0;

Double_t eVtoJ = 1.60218e-19;
Double_t EnergyDeposited_eV = 0;
Double_t acc_edep = 0;
Double_t acc_edep2 = 0;

Double_t Energy;
Double_t BPID;
Char_t Primary;
char *primaryName = new char[32];
char *type= new char[256];

// Read trees and leaves from root file, and give values to variables
TTree* tree = (TTree*) f->Get("tuples/primary_source");
Float_t number = (Float_t) tree->GetEntries();

if (number<2) {
    std::cout << "Not enough entries in the \"primary_source\" TTree (" << (long)number << " entries)\n";
    gApplication->Terminate(0);
  }

vector<pair<int,int64_t>> DSBBPID;

// For reading species production
tree = (TTree*) f->Get("tuples/damage");
tree->SetBranchAddress("Primary",           &Primary);
tree->SetBranchAddress("Energy",            &Energy);
tree->SetBranchAddress("EaqBaseHits",       &EB);
tree->SetBranchAddress("EaqStrandHits",     &ES);
tree->SetBranchAddress("OHBaseHits",        &OHB);
tree->SetBranchAddress("OHStrandHits",      &OHS);
tree->SetBranchAddress("HBaseHits",         &HB);
tree->SetBranchAddress("HStrandHits",       &HS);
tree->SetBranchAddress("TypeClassification", type);
tree->SetBranchAddress("BasePair",          &BPID);


Long64_t nentries = tree->GetEntries();
for(int i = 0;i<nentries;i++){
  tree->GetEntry(i);

  total_EB   += EB;
  total_EB2  += pow(EB,2);
  total_ES   += ES;
  total_ES2  += pow(ES,2);
  total_OHB  += OHB;
  total_OHB2 += pow(OHB,2);
  total_OHS  += OHS;
  total_OHS2 += pow(OHS,2);
  total_HB   += HB;
  total_HB2  += pow(HB,2);
  total_HS   += HS;
  total_HS2  += pow(HS,2);

  if((string)type=="DSB"||(string)type=="DSB+"||(string)type=="DSB++"){
    //cout << "DSB:"<<type<<endl;
    DSBBPID.push_back(make_pair(i,(int64_t)BPID));
    }

  }

// Sort DSBs from the one with lower ID value to the one with higher ID value
// Then find the number of fragments that have been produced
sort(DSBBPID.begin(),  DSBBPID.end(), smallerPair);
for(int ie = 0;ie<DSBBPID.size()-1;ie++){
  int64_t dsbfragment = DSBBPID[ie+1].second-DSBBPID[ie].second;

  double val       = (double)dsbfragment/1000.;
  double meanw     = h1fragments->GetBinCenter(h1fragments->FindBin(val));
  double binw      = h1fragments->GetBinWidth (h1fragments->FindBin(val));
  h1fragments->Fill(val,1./binw/1000);//bp-1
  //cout <<"val:"<<val<<endl;
  }

// Calculate the SEM of species
SD_EB  = SEM(total_EB2, total_EB, nentries);
SD_ES  = SEM(total_ES2, total_ES, nentries);
SD_OHB = SEM(total_OHB2, total_OHB, nentries);
SD_OHS = SEM(total_OHS2, total_OHS, nentries);
SD_HB  = SEM(total_HB2, total_HB, nentries);
SD_HS  = SEM(total_HS2, total_HS, nentries);

// Read damage classification SSB, SSB+, 2SSB, DSB, DSB+, DSB++
// As they have been defined in: Nikjoo, H., O’Neill, O., Goodhead, T., & Terrissol, M. 1997,
// Computational modelling of low-energy electron-induced DNA damage by early physical
// and chemical events, International Journal of Radiation Biology, 71, 467.
tree = (TTree *) f->Get("tuples/classification");
tree->SetBranchAddress("Primary",&Primary);
tree->SetBranchAddress("Energy", &Energy);
tree->SetBranchAddress("SSB",    &SSB);
tree->SetBranchAddress("SSBp",   &SSBp);
tree->SetBranchAddress("2SSB",   &SSB2p);
tree->SetBranchAddress("DSB",    &DSB);
tree->SetBranchAddress("DSBp",   &DSBp);
tree->SetBranchAddress("DSBpp",  &DSBpp);


Long64_t nentriesC = tree->GetEntries();
for(int i = 0;i<nentriesC;i++){
  tree->GetEntry(i);

  total_SSBp   += SSBp;
  total_SSBp2  += pow(SSBp,2);
  total_SSB2p  += SSB2p;
  total_SSB2p2 += pow(SSB2p,2);
  total_SSB    += SSB;
  total_SSB2   += pow(SSB,2);

  total_DSBp   += DSBp;
  total_DSBp2  += pow(DSBp,2);
  total_DSBpp  += DSBpp;
  total_DSBpp2 += pow(DSBpp,2);
  total_DSB    += DSB;
  total_DSB2   += pow(DSB,2);

  }

// Calculate the SEM
SD_SSB   = SEM(total_SSB2, total_SSB, nentriesC);
SD_SSBp  = SEM(total_SSBp2, total_SSBp, nentriesC);
SD_SSB2p = SEM(total_SSB2p2, total_SSB2p, nentriesC);
SD_DSB   = SEM(total_DSB2, total_DSB, nentriesC);
SD_DSBp  = SEM(total_DSBp2, total_DSBp, nentriesC);
SD_DSBpp = SEM(total_DSBpp2, total_DSBpp, nentriesC);

// Read damage classification SSBd, SSBi, SSBm, DSBd, DSBi, DSBm, DSBh
// As they have been defined in: Nikjoo, H., O’Neill, O., Goodhead, T., & Terrissol, M. 1997,
// Computational modelling of low-energy electron-induced DNA damage by early physical
// and chemical events, International Journal of Radiation Biology, 71, 467.
tree = (TTree *) f->Get("tuples/source");
tree->SetBranchAddress("Primary",primaryName);
tree->SetBranchAddress("Energy", &Energy);
tree->SetBranchAddress("SSBd",   &SSBd);
tree->SetBranchAddress("SSBi",   &SSBi);
tree->SetBranchAddress("SSBm",   &SSBm);
tree->SetBranchAddress("DSBd",   &DSBd);
tree->SetBranchAddress("DSBi",   &DSBi);
tree->SetBranchAddress("DSBm",   &DSBm);
tree->SetBranchAddress("DSBh",   &DSBh);

Long64_t nentriesS = tree->GetEntries();
for(int i = 0;i<nentriesS;i++){
  tree->GetEntry(i);

  total_SSBd += SSBd;
  total_SSBd2 += pow((SSBd),2);
  total_SSBi += SSBi;
  total_SSBi2 += pow((SSBi),2);
  total_SSBm += SSBm;
  total_SSBm2 += pow((SSBm),2);
  total_sSSB += SSBd + SSBi + SSBm;
  total_sSSB2 += pow((SSBd+SSBi+SSBm),2);

  total_DSBd  += DSBd;
  total_DSBd2 += pow(DSBd,2);
  total_DSBi  += DSBi;
  total_DSBi2 += pow(DSBi,2);
  total_DSBm  += DSBm;
  total_DSBm2 += pow(DSBm,2);
  total_DSBh  += DSBh;
  total_DSBh2 += pow(DSBh,2);
  total_sDSB  += DSBd + DSBi + DSBm + DSBh;
  total_sDSB2 += pow((DSBd+DSBi+DSBm+DSBh),2);

  }

// Calculate the SEM
SD_sSSB = SEM(total_sSSB2, total_sSSB, nentriesS);
SD_SSBd = SEM(total_SSBd2, total_SSBd, nentriesS);
SD_SSBi = SEM(total_SSBi2, total_SSBi, nentriesS);
SD_SSBm = SEM(total_SSBm2, total_SSBm, nentriesS);
SD_sDSB = SEM(total_sDSB2, total_sDSB, nentriesS);
SD_DSBd = SEM(total_DSBd2, total_DSBd, nentriesS);
SD_DSBi = SEM(total_DSBi2, total_DSBi, nentriesS);
SD_DSBm = SEM(total_DSBm2, total_DSBm, nentriesS);
SD_DSBh = SEM(total_DSBh2, total_DSBh, nentriesS);

// Measure the Deposited Energy in the whole volume that includes DNA chain
Long64_t totnentries =0;
tree = (TTree *) f->Get("tuples/chromosome_hits");
tree->SetBranchAddress("e_chromosome_kev",&EnergyDeposited_eV);
Long64_t nentriesD = tree->GetEntries();
totnentries+=nentriesD;
for(int i = 0;i<nentriesD;i++){
  tree->GetEntry(i);
  acc_edep += EnergyDeposited_eV *1e3;
  acc_edep2 += EnergyDeposited_eV *EnergyDeposited_eV *1e6;
}

// Close the root file to free space
f->Close();

// Calculate the absorbed dose
dose = acc_edep * eVtoJ / mass;
SD_dose = totnentries*SEM(acc_edep2, acc_edep, totnentries)* eVtoJ / mass;

// This is a normalization factor to produce the output in Gy-1 Gbp-1, or else.
// Default value is 1 to produce the result in Gy-1 Mbp-1
// It changes Mbp to Gbp. Some other changes may be needed in graphs section (name of axes)
double norm = 1;

// Calculate the yields, together with their error
EB_yield  = (Double_t) total_EB  / dose / Nbp;
ES_yield  = (Double_t) total_ES  / dose / Nbp;
OHB_yield = (Double_t) total_OHB / dose / Nbp;
OHS_yield = (Double_t) total_OHS / dose / Nbp;
HB_yield  = (Double_t) total_HB  / dose / Nbp;
HS_yield  = (Double_t) total_HS  / dose / Nbp;

if(0<total_EB )SD_EB_yield  = RelErr(nentries, SD_EB, total_EB, SD_dose, dose) ;
if(0<total_ES )SD_ES_yield  = RelErr(nentries, SD_ES, total_ES, SD_dose, dose) ;
if(0<total_OHB)SD_OHB_yield = RelErr(nentries, SD_OHB, total_OHB, SD_dose, dose);
if(0<total_OHS)SD_OHS_yield = RelErr(nentries, SD_OHS, total_OHS, SD_dose, dose);
if(0<total_HB )SD_HB_yield  = RelErr(nentries, SD_HB, total_HB, SD_dose, dose) ;
if(0<total_HS )SD_HS_yield  = RelErr(nentries, SD_HS, total_HS, SD_dose, dose) ;
//Previous method
//SD_EB_yield  = SD_EB  / dose / Nbp;
//SD_ES_yield  = SD_ES  / dose / Nbp;
//SD_OHB_yield = SD_OHB / dose / Nbp;
//SD_OHS_yield = SD_OHS / dose / Nbp;
//SD_HB_yield  = SD_HB  / dose / Nbp;
//SD_HS_yield  = SD_HS  / dose / Nbp;

SSB_yield   = (Double_t) norm * total_SSB   / dose / Nbp;
SSBp_yield  = (Double_t) norm * total_SSBp  / dose / Nbp;
SSB2p_yield = (Double_t) norm * total_SSB2p / dose / Nbp;

DSB_yield   = (Double_t) norm * total_DSB   / dose / Nbp;
DSBp_yield  = (Double_t) norm * total_DSBp  / dose / Nbp;
DSBpp_yield = (Double_t) norm * total_DSBpp / dose / Nbp;

if(0<total_SSB  )SD_SSB_yield   = RelErr(nentriesC, SD_SSB, total_SSB, SD_dose, dose)  ;
if(0<total_SSBp )SD_SSBp_yield  = RelErr(nentriesC, SD_SSBp, total_SSBp, SD_dose, dose) ;
if(0<total_SSB2p)SD_SSB2p_yield = RelErr(nentriesC, SD_SSB2p, total_SSB2p, SD_dose, dose);
if(0<total_DSB  )SD_DSB_yield   = RelErr(nentriesC, SD_DSB, total_DSB, SD_dose, dose)  ;
if(0<total_DSBp )SD_DSBp_yield  = RelErr(nentriesC, SD_DSBp, total_DSBp, SD_dose, dose) ;
if(0<total_DSBpp)SD_DSBpp_yield = RelErr(nentriesC, SD_DSBpp, total_DSBpp, SD_dose, dose);
//Previous method
//SD_SSB_yield   = norm * SD_SSB   / dose / Nbp;
//SD_SSBp_yield  = norm * SD_SSBp  / dose / Nbp;
//SD_SSB2p_yield = norm * SD_SSB2p / dose / Nbp;
//SD_DSB_yield   = norm * SD_DSB   / dose / Nbp;
//SD_DSBp_yield  = norm * SD_DSBp  / dose / Nbp;
//SD_DSBpp_yield = norm * SD_DSBpp / dose / Nbp;

sSSB_yield = (Double_t) norm * total_sSSB / dose / Nbp;
SSBi_yield = (Double_t) norm * total_SSBi / dose / Nbp;
SSBd_yield = (Double_t) norm * total_SSBd / dose / Nbp;
SSBm_yield = (Double_t) norm * total_SSBm / dose / Nbp;

sDSB_yield = (Double_t) norm * total_sDSB / dose / Nbp;
DSBi_yield = (Double_t) norm * total_DSBi / dose / Nbp;
DSBd_yield = (Double_t) norm * total_DSBd / dose / Nbp;
DSBm_yield = (Double_t) norm * total_DSBm / dose / Nbp;
DSBh_yield = (Double_t) norm * total_DSBh / dose / Nbp;

if(0<total_sSSB)SD_sSSB_yield = RelErr(nentriesS, SD_sSSB, total_sSSB, SD_dose, dose);
if(0<total_SSBi)SD_SSBi_yield = RelErr(nentriesS, SD_SSBi, total_SSBi, SD_dose, dose);
if(0<total_SSBd)SD_SSBd_yield = RelErr(nentriesS, SD_SSBd, total_SSBd, SD_dose, dose);
if(0<total_SSBm)SD_SSBm_yield = RelErr(nentriesS, SD_SSBm, total_SSBm, SD_dose, dose);
if(0<total_sDSB)SD_sDSB_yield = RelErr(nentriesS, SD_sDSB, total_sDSB, SD_dose, dose);
if(0<total_DSBi)SD_DSBi_yield = RelErr(nentriesS, SD_DSBi, total_DSBi, SD_dose, dose);
if(0<total_DSBd)SD_DSBd_yield = RelErr(nentriesS, SD_DSBd, total_DSBd, SD_dose, dose);
if(0<total_DSBm)SD_DSBm_yield = RelErr(nentriesS, SD_DSBm, total_DSBm, SD_dose, dose);
if(0<total_DSBh)SD_DSBh_yield = RelErr(nentriesS, SD_DSBh, total_DSBh, SD_dose, dose);
//Previous method
//SD_sSSB_yield = norm * SD_sSSB / dose / Nbp;
//SD_SSBi_yield = norm * SD_SSBi / dose / Nbp;
//SD_SSBd_yield = norm * SD_SSBd / dose / Nbp;
//SD_SSBm_yield = norm * SD_SSBm / dose / Nbp;
//SD_sDSB_yield = norm * SD_sDSB / dose / Nbp;
//SD_DSBi_yield = norm * SD_DSBi / dose / Nbp;
//SD_DSBd_yield = norm * SD_DSBd / dose / Nbp;
//SD_DSBm_yield = norm * SD_DSBm / dose / Nbp;
//SD_DSBh_yield = norm * SD_DSBh / dose / Nbp;


// Print output in terminal

float total_SSB_totalYield = SSB_yield + SSBp_yield + SSB2p_yield;
float total_DSB_totalYield = DSB_yield + DSBp_yield + DSBpp_yield;


cout<<"\n"                      <<ifile         <<'\n'
    <<"\nDose Absorbed (Gy): "  <<dose          <<" relative error %: "  <<100*SD_dose/dose <<'\n'
    <<"Particle : "             <<primaryName   <<'\t'
    <<"Energy (MeV) : "         <<Energy        <<'\t'
    <<"Number of Primaries : "  <<number        <<'\n'
    <<"  Output Damage : "                      <<'\n'<<'\t'
    <<"     Species Hits (Gy-1 Mbp-1) "      <<'\n'<<'\t'
    <<"EaqBaseHits   : "     <<EB_yield     <<"   \t"      <<" relative error %: "   <<100*SD_EB_yield       <<'\n'<<'\t'
    <<"EaqStrandHits : "     <<ES_yield     <<"   \t"      <<" relative error %: "   <<100*SD_ES_yield       <<'\n'<<'\t'
    <<"OHBaseHits    : "     <<OHB_yield    <<"   \t"      <<" relative error %: "   <<100*SD_OHB_yield      <<'\n'<<'\t'
    <<"OHStrandHits  : "     <<OHS_yield    <<"   \t"      <<" relative error %: "   <<100*SD_OHS_yield      <<'\n'<<'\t'
    <<"HBaseHits     : "     <<HB_yield     <<"   \t"      <<" relative error %: "   <<100*SD_HB_yield       <<'\n'<<'\t'
    <<"HStrandHits   : "     <<HS_yield     <<"   \t"      <<" relative error %: "   <<100*SD_HS_yield       <<'\n'<<'\n'<<'\t'
    <<"     Damage yield (Gy-1 Mbp-1) "      <<'\n'<<'\t'
    <<"SSB           : "     <<SSB_yield    <<"   \t"      <<" relative error %: "   <<100*SD_SSB_yield    <<'\n'<<'\t'
    <<"SSB+          : "     <<SSBp_yield   <<"   \t"      <<" relative error %: "   <<100*SD_SSBp_yield   <<'\n'<<'\t'
    <<"2SSB          : "     <<SSB2p_yield  <<"   \t"      <<" relative error %: "   <<100*SD_SSB2p_yield  <<'\n'<<'\t'
    <<"SSB total     : "     <<total_SSB_totalYield        <<'\n'<<'\t'
    <<"DSB           : "     <<DSB_yield    <<"   \t"      <<" relative error %: "   <<100*SD_DSB_yield    <<'\n'<<'\t'
    <<"DSB+          : "     <<DSBp_yield   <<"   \t"      <<" relative error %: "   <<100*SD_DSBp_yield   <<'\n'<<'\t'
    <<"DSB++         : "     <<DSBpp_yield  <<"   \t"      <<" relative error %: "   <<100*SD_DSBpp_yield  <<'\n'<<'\t'
    <<"DSB total     : "     <<total_DSB_totalYield        <<'\n'<<'\n'<<'\t'
    <<"     Breaks yield (Gy-1 Mbp-1) "      <<'\n'<<'\t'
    <<"SSB direct    : "     <<SSBd_yield   <<"   \t"      <<" relative error %: "   <<100*SD_SSBd_yield   <<'\n'<<'\t'
    <<"SSB indirect  : "     <<SSBi_yield   <<"   \t"      <<" relative error %: "   <<100*SD_SSBi_yield   <<'\n'<<'\t'
    <<"SSB mixed     : "     <<SSBm_yield   <<"   \t"      <<" relative error %: "   <<100*SD_SSBm_yield   <<'\n'<<'\t'
    <<"SSB total     : "     <<sSSB_yield   <<"   \t"      <<" relative error %: "   <<100*SD_sSSB_yield   <<'\n'<<'\t'
    <<"DSB direct    : "     <<DSBd_yield   <<"   \t"      <<" relative error %: "   <<100*SD_DSBd_yield   <<'\n'<<'\t'
    <<"DSB indirect  : "     <<DSBi_yield   <<"   \t"      <<" relative error %: "   <<100*SD_DSBi_yield   <<'\n'<<'\t'
    <<"DSB mixed     : "     <<DSBm_yield   <<"   \t"      <<" relative error %: "   <<100*SD_DSBm_yield   <<'\n'<<'\t'
    <<"DSB hybrid    : "     <<DSBh_yield   <<"   \t"      <<" relative error %: "   <<100*SD_DSBh_yield   <<'\n'<<'\t'
    <<"DSB total     : "     <<sDSB_yield   <<"   \t"      <<" relative error %: "   <<100*SD_sDSB_yield   <<'\n'<<'\n'<<'\t'
    <<"SSB total /DSB total : "     <<sSSB_yield/sDSB_yield       <<'\n'<<'\n';


// Plot Histograms

cfragment->GetCanvas()->cd();
h1fragments->SetStats(false);
h1fragments->SetMarkerSize(0.1);
h1fragments->SetMarkerColor(kRed);
h1fragments->SetLineColor  (kRed);
h1fragments->Scale(1./(Nbp*1e6)); //bp^-1
h1fragments->SetTitle("");
h1fragments->SetYTitle("Number of Fragments (bp^{-2})");
h1fragments->SetXTitle("Fragment Length (kbp)");
//h1fragments->SetAxisRange(10,1e4);
h1fragments->SetAxisRange(1,10);
//h1fragments->SetMaximum(3e-11);
//h1fragments->SetMinimum(1e-15);
h1fragments->Draw();


c1->GetCanvas()->cd();
pad1->cd();
const Int_t n = 6;
Double_t x[n] = {1,2,3,4,5,6};
Double_t y[n] = {EB_yield,ES_yield,OHB_yield,OHS_yield,HB_yield,HS_yield};
Double_t err_y[n] = {EB_yield*SD_EB_yield,ES_yield*SD_ES_yield,OHB_yield*SD_OHB_yield,OHS_yield*SD_OHS_yield,HB_yield*SD_HB_yield,HS_yield*SD_HS_yield};
TGraph* gr = new TGraphErrors(n,x,y,0,err_y);
gr->SetTitle("Species");
gr->GetXaxis()->SetBinLabel(9, "EaqBaseHits");
gr->GetXaxis()->SetBinLabel(25,"EaqStrandHits");
gr->GetXaxis()->SetBinLabel(42,"OHBaseHits");
gr->GetXaxis()->SetBinLabel(58,"OHStrandHits");
gr->GetXaxis()->SetBinLabel(75,"HBaseHits");
gr->GetXaxis()->SetBinLabel(92,"HStrandHits");
gr->GetYaxis()->SetTitle("Species Hits (Gy^{-1} Mbp^{-1})");
gr->GetYaxis()->SetTitleOffset(2);

gr->SetFillColor(49);
gr->Draw("ba");


pad2->cd();
Double_t x2[n] = {1,2,3,4,5,6};
Double_t y2[n] = {SSBp_yield,SSB2p_yield,SSB_yield,DSBp_yield,DSBpp_yield,DSB_yield};
Double_t err_y2[n] = {SSBp_yield*SD_SSBp_yield,SSB2p_yield*SD_SSB2p_yield,SSB_yield*SD_SSB_yield,DSBp_yield*SD_DSBp_yield,DSBpp_yield*SD_DSBpp_yield,DSB_yield*SD_DSB_yield};
TGraph* gr2 = new TGraphErrors(n,x2,y2,0,err_y2);
gr2->SetTitle("Damage Yield");
gr2->GetXaxis()->SetBinLabel(9, "SSB+");
gr2->GetXaxis()->SetBinLabel(25,"2SSB");
gr2->GetXaxis()->SetBinLabel(42,"SSB");
gr2->GetXaxis()->SetBinLabel(58,"DSB+");
gr2->GetXaxis()->SetBinLabel(75,"DSB++");
gr2->GetXaxis()->SetBinLabel(92,"DSB");
gr2->GetYaxis()->SetTitle("Damage yield (Gy^{-1} Mbp^{-1})");
//gr2->GetYaxis()->SetTitle("Damage yield (particle^{-1})");
gr2->GetYaxis()->SetTitleOffset(2);

gr2->SetFillColor(8);
gr2->Draw("ba");


pad3->cd();
const Int_t m = 4;
Double_t x3[m] = {1,2,3,4};
Double_t y3[m] = {SSBd_yield,SSBi_yield,SSBm_yield,sSSB_yield};
Double_t err_y3[m] = {SSBd_yield*SD_SSBd_yield,SSBi_yield*SD_SSBi_yield,SSBm_yield*SD_SSBm_yield,sSSB_yield*SD_sSSB_yield};
TGraph* gr3 = new TGraphErrors(m,x3,y3,0,err_y3);
gr3->SetTitle("Breaks Yield");
gr3->GetXaxis()->SetBinLabel(8, "SSB direct");
gr3->GetXaxis()->SetBinLabel(35,"SSB indirect");
gr3->GetXaxis()->SetBinLabel(64,"SSB mixed");
gr3->GetXaxis()->SetBinLabel(92,"SSB all");
gr3->GetYaxis()->SetTitle("Breaks yield (Gy^{-1} Mbp^{-1})");
//gr3->GetYaxis()->SetTitle("SSB yield (particle^{-1}");
gr3->GetYaxis()->SetTitleOffset(2);

gr3->SetFillColor(7);
gr3->Draw("ba");


pad4->cd();
const Int_t k = 5;
Double_t x4[k] = {1,2,3,4,5};
Double_t y4[k] = {DSBd_yield,DSBi_yield,DSBm_yield,DSBh_yield,sDSB_yield};
Double_t err_y4[k] = {DSBd_yield*SD_DSBd_yield,DSBi_yield*SD_DSBi_yield,DSBm_yield*SD_DSBm_yield,DSBh_yield*SD_DSBh_yield,sDSB_yield*SD_sDSB_yield};
TGraph* gr4 = new TGraphErrors(k,x4,y4,0,err_y4);
gr4->SetTitle("Breaks Yield");
gr4->GetXaxis()->SetBinLabel(8,"DSB direct");
gr4->GetXaxis()->SetBinLabel(29,"DSB indirect");
gr4->GetXaxis()->SetBinLabel(50,"DSB mixed");
gr4->GetXaxis()->SetBinLabel(71,"DSB hybrid");
gr4->GetXaxis()->SetBinLabel(92,"DSB all");
gr4->GetYaxis()->SetTitle("Breaks yield (Gy^{-1} Mbp^{-1})");
//gr4->GetYaxis()->SetTitle("DSB yield (particle^{-1})");
gr4->GetYaxis()->SetTitleOffset(2);

gr4->SetFillColor(4);
gr4->Draw("ba");

}

// Some important bools that are needed to run the root macro file
bool greaterPair(const ipair& l, const ipair& r){return l.second > r.second;}
bool smallerPair(const ipair& l, const ipair& r){return l.second < r.second;}

void BinLogX(TH1 *h) {
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;

}

//==============================================================================//
//                           HELPER FUNCTION DEFINITIONS                         //
//==============================================================================//

// Standard error of the mean: sqrt( abs( (s2/n - ((double)s/n)^2) ) / (n-1) )
// The sum s is cast to Double_t before dividing (as in the original macro),
// while s2 retains its original type so the s2/n division keeps the same
// float vs. double behaviour of the caller.
Double_t SEM(Float_t s2, Int_t s, Long64_t n) {
  return sqrt(abs(((s2 / n) - pow((Double_t) s / n, 2))) / (n - 1));
}
Double_t SEM(Double_t s2, Double_t s, Long64_t n) {
  return sqrt(abs(((s2 / n) - pow(s / n, 2))) / (n - 1));
}

// Relative error of a yield: sqrt( (n*SD/total)^2 + (SD_dose/dose)^2 )
// The factor n reproduces the original formula (SDs are SEMs already scaled by 1/sqrt(n),
// so multiplying by n reconstructs the standard deviation).
// The n*SD/total term is computed in Float_t to reproduce the original macro, which
// multiplies/divides the Float_t SD and Int_t total variables in single precision.
Double_t RelErr(Long64_t n, Double_t SD, Double_t total, Double_t SD_dose, Double_t dose) {
  Float_t term = n * (Float_t) SD / (Int_t) total;
  return sqrt(pow(term, 2) + pow(SD_dose / dose, 2));
}
