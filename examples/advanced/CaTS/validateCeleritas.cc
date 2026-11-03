//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// ********************************************************************
//
/// \file validateCeleritas.cc
/// \brief Compare CaTS photon hits from a Geant4-only run against a
///        Celeritas distribution-offload run, producing overlay and
///        ratio plots for physics validation.
///
/// Usage:
///   validateCeleritas g4file.root celefile.root outbase [collection]
///
/// The collection defaults to "Det_PhotonDetector_HC".
/// Output: <outbase>.root and <outbase>.pdf
//
// Root headers
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
// Project headers
#include "Event.hh"
#include "PhotonHit.hh"

namespace
{
constexpr unsigned kBinsWl = 300;
constexpr double kWlMax = 300.0;  // nm
constexpr unsigned kBinsTime = 500;
constexpr double kTimeMax = 1000.0;  // ns
constexpr unsigned kBinsNph = 300;
constexpr double kNphMax = 600.0;  // photons / event
constexpr unsigned kBinsPos = 200;
constexpr double kYRange = 600.0;  // mm
constexpr double kZRange = 1100.0;  // mm

struct PhotonStats
{
    Long64_t events = 0;
    Long64_t nphotons = 0;
    double mean_nph = 0;
    double rms_nph = 0;
    double mean_wl = 0;
    double rms_wl = 0;
    double mean_time = 0;
    double rms_time = 0;
    double mean_z = 0;
    double mean_y = 0;
    Long64_t pid0 = 0;
    Long64_t pid1 = 0;
};

struct Plots
{
    TH1F* nph = nullptr;
    TH1F* wl = nullptr;
    TH1F* time = nullptr;
    TH1F* z = nullptr;
    TH1F* y = nullptr;
    TH2F* yz = nullptr;
    TH1F* pid = nullptr;
};

std::string g_collection = "Det_PhotonDetector_HC";

Plots MakePlots(const char* nameprefix)
{
  Plots p;
  p.nph = new TH1F(Form("%s_nph", nameprefix), "Detected photons per event;N_{photons};Events",
                   kBinsNph, 0, kNphMax);
  p.wl = new TH1F(Form("%s_wl", nameprefix), "Photon wavelength;#lambda [nm];Photons", kBinsWl, 0,
                  kWlMax);
  p.time = new TH1F(Form("%s_time", nameprefix), "Photon detection time;t [ns];Photons", kBinsTime,
                    0, kTimeMax);
  p.z = new TH1F(Form("%s_z", nameprefix), "Hit position Z;Z [mm];Photons", kBinsPos, -kZRange,
                 kZRange);
  p.y = new TH1F(Form("%s_y", nameprefix), "Hit position Y;Y [mm];Photons", kBinsPos, -kYRange,
                 kYRange);
  p.yz = new TH2F(Form("%s_yz", nameprefix), "Hit position;Z [mm];Y [mm]", kBinsPos, -kZRange,
                  kZRange, kBinsPos, -kYRange, kYRange);
  p.pid = new TH1F(Form("%s_pid", nameprefix), "Photon creation process;Process;Photons", 2, 0, 2);
  p.pid->GetXaxis()->SetBinLabel(1, "Cerenkov");
  p.pid->GetXaxis()->SetBinLabel(2, "Scintillation");
  p.nph->Sumw2();
  p.wl->Sumw2();
  p.time->Sumw2();
  p.z->Sumw2();
  p.y->Sumw2();
  p.yz->Sumw2();
  p.pid->Sumw2();
  return p;
}

PhotonStats FillPlots(const char* filename, Plots* p)
{
  PhotonStats s;
  TFile fo(filename);
  if (fo.IsZombie())
  {
    G4cerr << "ERROR: cannot open " << filename << G4endl;
    exit(1);
  }
  Event* event = new Event();
  TTree* Tevt = (TTree*)fo.Get("Events");
  if (!Tevt)
  {
    G4cerr << "ERROR: no Events tree in " << filename << G4endl;
    exit(1);
  }
  Tevt->SetBranchAddress("event.", &event);
  TBranch* fevtbranch = Tevt->GetBranch("event.");
  Long64_t nevent = fevtbranch->GetEntries();
  s.events = nevent;
  for (Long64_t i = 0; i < nevent; ++i)
  {
    fevtbranch->GetEntry(i);
    auto* hcmap = event->GetHCMap();
    for (const auto& ele : *hcmap)
    {
      if (ele.first.compare(g_collection) != 0)
      {
        continue;
      }
      const auto& hits = ele.second;
      Long64_t n = hits.size();
      s.nphotons += n;
      p->nph->Fill(n);
      for (auto hit : hits)
      {
        PhotonHit* ph = dynamic_cast<PhotonHit*>(hit);
        if (!ph)
        {
          continue;
        }
        p->wl->Fill(ph->GetWavelength());
        p->time->Fill(ph->GetTime());
        p->z->Fill(ph->GetPosition().getZ());
        p->y->Fill(ph->GetPosition().getY());
        p->yz->Fill(ph->GetPosition().getZ(), ph->GetPosition().getY());
        p->pid->Fill(ph->GetPid());
        if (ph->GetPid() == 0)
        {
          ++s.pid0;
        }
        else if (ph->GetPid() == 1)
        {
          ++s.pid1;
        }
      }
      break;
    }
  }
  s.mean_nph = p->nph->GetMean();
  s.rms_nph = p->nph->GetRMS();
  s.mean_wl = p->wl->GetMean();
  s.rms_wl = p->wl->GetRMS();
  s.mean_time = p->time->GetMean();
  s.rms_time = p->time->GetRMS();
  s.mean_z = p->z->GetMean();
  s.mean_y = p->y->GetMean();
  return s;
}

double Chi2PerNdf(TH1F* a, TH1F* b)
{
  // Both normalized to unit integral.  Only bins where either has content
  // contribute.
  double chi2 = 0;
  int ndf = 0;
  for (int i = 1; i <= a->GetNbinsX(); ++i)
  {
    double ca = a->GetBinContent(i);
    double cb = b->GetBinContent(i);
    if (ca <= 0 && cb <= 0)
    {
      continue;
    }
    double ea2 = a->GetBinError(i) * a->GetBinError(i);
    double eb2 = b->GetBinError(i) * b->GetBinError(i);
    if (ea2 + eb2 <= 0)
    {
      continue;
    }
    chi2 += (ca - cb) * (ca - cb) / (ea2 + eb2);
    ++ndf;
  }
  return ndf > 0 ? chi2 / ndf : -1;
}

TH1F* MakeRatio(TH1F* ce, TH1F* g4, const char* name)
{
  TH1F* r = (TH1F*)ce->Clone(name);
  r->Divide(g4);
  r->SetTitle("Ratio Celeritas/Geant4");
  return r;
}

void StylePlot(TH1F* h, Color_t color, int style)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(0.8);
  h->SetLineWidth(2);
  h->SetStats(false);
}

void Overlay(TCanvas* c, int pad, TH1F* g4, TH1F* ce, const char* xtitle, const char* ytitle)
{
  c->cd(pad);
  g4->SetTitle(Form("%s;%s;%s", g4->GetTitle(), xtitle, ytitle));
  if (pad == 3) gPad->SetLogy();
  else g4->SetMinimum(0);
  g4->Draw("HIST E");
  ce->Draw("HIST E SAME");
  TLegend* leg = new TLegend(0.62, 0.75, 0.89, 0.88);
  leg->AddEntry(g4, "Geant4", "L");
  leg->AddEntry(ce, "Celeritas", "L");
  leg->SetBorderSize(0);
  leg->Draw();
}
}  // namespace

int main(int argc, char** argv)
{
  TSystem ts;
  gSystem->Load("libCaTSClassesDict");
  if (argc < 4)
  {
    G4cerr << "Usage: validateCeleritas <g4file.root> <celefile.root> "
              "<outbase> [collection]"
           << G4endl;
    return 1;
  }
  const char* g4file = argv[1];
  const char* ceFile = argv[2];
  G4String outbase = argv[3];
  if (argc >= 5)
  {
    g_collection = argv[4];
  }

  Plots pg4 = MakePlots("g4");
  Plots pce = MakePlots("ce");
  PhotonStats sg4 = FillPlots(g4file, &pg4);
  PhotonStats sce = FillPlots(ceFile, &pce);

  // ---- Summary ----
  G4cout << "\n=== Validation summary ===" << G4endl;
  G4cout << "Collection: " << g_collection << G4endl;
  G4cout << "  " << g4file << ": " << sg4.events << " events, " << sg4.nphotons << " photon hits"
         << G4endl;
  G4cout << "  " << ceFile << ": " << sce.events << " events, " << sce.nphotons << " photon hits"
         << G4endl;
  G4cout << "  hits/event G4=" << sg4.mean_nph << " +/- " << sg4.rms_nph
         << "  Cele=" << sce.mean_nph << " +/- " << sce.rms_nph
         << "  ratio=" << sce.nphotons / double(sg4.nphotons) << G4endl;
  G4cout << "  mean wavelength G4=" << sg4.mean_wl << " nm  Cele=" << sce.mean_wl
         << " nm  diff=" << (sce.mean_wl - sg4.mean_wl) / sg4.mean_wl * 100 << "%" << G4endl;
  G4cout << "  mean time G4=" << sg4.mean_time << " ns  Cele=" << sce.mean_time
         << " ns  diff=" << (sce.mean_time - sg4.mean_time) / sg4.mean_time * 100 << "%" << G4endl;
  G4cout << "  mean Z G4=" << sg4.mean_z << " mm (RMS " << pg4.z->GetRMS()
         << ")  Cele=" << sce.mean_z << " mm (RMS " << pce.z->GetRMS()
         << ")  |diff|=" << TMath::Abs(sce.mean_z - sg4.mean_z) << " mm" << G4endl;
  G4cout << "  mean Y G4=" << sg4.mean_y << " mm (RMS " << pg4.y->GetRMS()
         << ")  Cele=" << sce.mean_y << " mm (RMS " << pce.y->GetRMS()
         << ")  |diff|=" << TMath::Abs(sce.mean_y - sg4.mean_y) << " mm" << G4endl;
  G4cout << "  PID G4: Cerenkov=" << sg4.pid0 << " Scint=" << sg4.pid1
         << "  Cele: Cerenkov=" << sce.pid0 << " Scint=" << sce.pid1 << G4endl;

  // Style
  StylePlot(pg4.nph, kBlack, 20);
  StylePlot(pg4.wl, kBlack, 20);
  StylePlot(pg4.time, kBlack, 20);
  StylePlot(pg4.z, kBlack, 20);
  StylePlot(pg4.y, kBlack, 20);
  StylePlot(pg4.pid, kBlack, 20);
  StylePlot(pce.nph, kRed, 24);
  StylePlot(pce.wl, kRed, 24);
  StylePlot(pce.time, kRed, 24);
  StylePlot(pce.z, kRed, 24);
  StylePlot(pce.y, kRed, 24);
  StylePlot(pce.pid, kRed, 24);

  // ---- Ratios ----
  TH1F* r_wl = MakeRatio(pce.wl, pg4.wl, "ratio_wl");
  TH1F* r_time = MakeRatio(pce.time, pg4.time, "ratio_time");
  TH1F* r_z = MakeRatio(pce.z, pg4.z, "ratio_z");
  TH1F* r_y = MakeRatio(pce.y, pg4.y, "ratio_y");
  for (TH1F* r : {r_wl, r_time, r_z, r_y})
  {
    r->SetMarkerStyle(20);
    r->SetMarkerSize(0.7);
    r->SetLineColor(kBlue + 2);
    r->SetStats(false);
    r->GetYaxis()->SetRangeUser(0.0, 2.0);
  }

  // ---- PDF ----
  TCanvas* c1 = new TCanvas("c1", "Validation 1D", 1600, 1200);
  c1->Divide(2, 2);
  Overlay(c1, 1, pg4.nph, pce.nph, "N_{photons}/event", "Events (norm)");
  Overlay(c1, 2, pg4.wl, pce.wl, "#lambda [nm]", "Photons (norm)");
  Overlay(c1, 3, pg4.time, pce.time, "t [ns]", "Photons (norm)");
  c1->cd(4)->SetLogy();
  pg4.pid->GetXaxis()->SetRangeUser(0.1, 2);
  pg4.pid->Draw("HIST E");
  pce.pid->Draw("HIST E SAME");
  TLegend* leg4 = new TLegend(0.62, 0.75, 0.89, 0.88);
  leg4->AddEntry(pg4.pid, "Geant4", "L");
  leg4->AddEntry(pce.pid, "Celeritas", "L");
  leg4->SetBorderSize(0);
  leg4->Draw();

  TCanvas* c2 = new TCanvas("c2", "Validation position", 1600, 1200);
  c2->Divide(2, 2);
  Overlay(c2, 1, pg4.z, pce.z, "Z [mm]", "Photons (norm)");
  Overlay(c2, 2, pg4.y, pce.y, "Y [mm]", "Photons (norm)");
  c2->cd(3);
  pg4.yz->Draw("COLZ");
  c2->cd(4);
  pce.yz->Draw("COLZ");

  TCanvas* c3 = new TCanvas("c3", "Ratios", 1600, 1200);
  c3->Divide(2, 2);
  c3->cd(1);
  r_wl->Draw("HIST E");
  TLine* l1 = new TLine(0, 1, kWlMax, 1);
  l1->SetLineColor(kBlack);
  l1->SetLineStyle(2);
  l1->Draw();
  c3->cd(2);
  r_time->Draw("HIST E");
  TLine* l2 = new TLine(0, 1, kTimeMax, 1);
  l2->SetLineColor(kBlack);
  l2->SetLineStyle(2);
  l2->Draw();
  c3->cd(3);
  r_z->Draw("HIST E");
  TLine* l3 = new TLine(-kZRange, 1, kZRange, 1);
  l3->SetLineColor(kBlack);
  l3->SetLineStyle(2);
  l3->Draw();
  c3->cd(4);
  r_y->Draw("HIST E");
  TLine* l4 = new TLine(-kYRange, 1, kYRange, 1);
  l4->SetLineColor(kBlack);
  l4->SetLineStyle(2);
  l4->Draw();

  // ---- Save ----
  G4String rootname = outbase + ".root";
  TFile out(rootname, "RECREATE");
  pg4.nph->Write();
  pg4.wl->Write();
  pg4.time->Write();
  pg4.z->Write();
  pg4.y->Write();
  pg4.yz->Write();
  pg4.pid->Write();
  pce.nph->Write();
  pce.wl->Write();
  pce.time->Write();
  pce.z->Write();
  pce.y->Write();
  pce.yz->Write();
  pce.pid->Write();
  r_wl->Write();
  r_time->Write();
  r_z->Write();
  r_y->Write();
  out.Close();

  G4String pdfname = outbase + ".pdf";
  c1->Print(Form("%s(", pdfname.data()));
  c2->Print(Form("%s", pdfname.data()));
  c3->Print(Form("%s)", pdfname.data()));

  G4cout << "\nWrote " << pdfname << " and " << rootname << G4endl;
  return 0;
}
