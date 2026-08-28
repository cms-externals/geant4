//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//---------------------------------------------------------------------------//
// Validation utility for G4MuonDecayChannel::DecayIt.
//---------------------------------------------------------------------------//

#include "G4MuonDecayChannel.hh"

#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4Positron.hh"

#include "CLHEP/Random/Random.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

//---------------------------------------------------------------------------//
// VALIDATION FIXTURES
//---------------------------------------------------------------------------//
/**
 * \brief Container for sampled and expected 1D distributions and their comparison metrics.
 *
 * This type owns four histograms sharing the same binning:
 * - sampled counts,
 * - expected counts,
 * - residuals (%),
 * - pulls.
 *
 * Expectation values are supplied per-bin through a callback, allowing either
 * analytic-model or data-driven expected distributions.
 */
struct ValidationHistogram1D
{
    /**
     * \brief Per-bin value returned by an expectation callback.
     */
    struct ExpectedBinValue
    {
        G4double count = 0.0;
        G4double sigma = 0.0;
        G4bool ok = true;
    };

    using ExpectedCallback =
      std::function<ExpectedBinValue(G4int bin_id, G4double x_low, G4double x_high)>;

    ValidationHistogram1D() = delete;
    /**
     * \brief Construct histograms with shared binning and axis titles.
     * \param prefix Name prefix used to label internal ROOT histograms.
     * \param x_axis_title X-axis title for all managed histograms.
     * \param y_axis_title Y-axis title used for sampled/expected histograms.
     * \param bins Number of bins.
     * \param low_edge Lower histogram edge.
     * \param high_edge Upper histogram edge.
     */
    ValidationHistogram1D(const G4String& prefix, const G4String& x_axis_title,
                          const G4String& y_axis_title, G4int bins, G4double low_edge,
                          G4double high_edge)
    {
      sample.SetNameTitle((prefix + "_sample").c_str(), "Sampled");
      sample.SetBins(bins, low_edge, high_edge);
      sample.GetXaxis()->SetTitle(x_axis_title);
      sample.GetYaxis()->SetTitle(y_axis_title);
      sample.SetStats(0);

      expect.SetNameTitle((prefix + "_expect").c_str(), "Expected");
      expect.SetBins(bins, low_edge, high_edge);
      expect.GetXaxis()->SetTitle(x_axis_title);
      expect.GetYaxis()->SetTitle(y_axis_title);
      expect.SetStats(0);

      residuals.SetNameTitle((prefix + "_residuals").c_str(), "Residuals");
      residuals.SetBins(bins, low_edge, high_edge);
      residuals.GetXaxis()->SetTitle(x_axis_title);
      residuals.GetYaxis()->SetTitle("Residual (%)");
      residuals.SetStats(0);

      pulls.SetNameTitle((prefix + "_pulls").c_str(), "Pulls");
      pulls.SetBins(bins, low_edge, high_edge);
      pulls.GetXaxis()->SetTitle(x_axis_title);
      pulls.GetYaxis()->SetTitle("Pull");
      pulls.SetStats(0);
    }

    /**
     * \brief Fill one sampled entry.
     * \param val Sample value.
     * \param weight Statistical weight applied to this entry.
     */
    void FillSampleHistogram(const G4double val, const G4double weight = 1.0)
    {
      sample.Fill(val, weight);
    }

    /**
     * \brief Fill the expected histogram bin-by-bin via callback.
     * \param cb Callback receiving bin id and bin boundaries.
     * \param scale_factor Global scale applied to expected counts (for example, N events).
     *
     * The callback contract is:
     * - input: bin id, lower edge, upper edge,
     * - output: expected count, expected uncertainty, success flag.
     */
    void FillExpectedHistogram(ExpectedCallback cb, G4double scale_factor = 1.0)
    {
      for (G4int i = 1; i <= expect.GetNbinsX(); ++i)
      {
        const G4double x_low = expect.GetBinLowEdge(i);
        const G4double x_high = x_low + expect.GetBinWidth(i);

        const ExpectedBinValue v = cb(i, x_low, x_high);
        if (!v.ok)
        {
          std::cout << expect.GetName() << ": expectation callback failed at bin " << i << "\n";
        }

        expect.SetBinContent(i, scale_factor * v.count);
        expect.SetBinError(i, v.sigma);
      }
    }

    /**
     * \brief Build residual and pull histograms from sampled and expected data.
     */
    void FillComparisonHistograms()
    {
      auto calculate_pull = [](const G4double observed, const G4double expected) {
        if (observed == 0.0) return std::pair{0.0, 0.0};
        return std::pair{(observed - expected) / std::sqrt(observed),
                         0.5 + 0.5 * expected / observed};
      };

      auto calculate_residual = [](const G4double observed, const G4double expected) {
        return std::pair{100.0 * (observed - expected) / expected,
                         100.0 * std::sqrt(observed) / expected};
      };

      for (G4int i = 1; i <= sample.GetNbinsX(); ++i)
      {
        const G4double data = sample.GetBinContent(i);
        const G4double expected = expect.GetBinContent(i);

        const auto [p, p_error] = calculate_pull(data, expected);
        pulls.SetBinContent(i, p);
        pulls.SetBinError(i, p_error);

        const auto [r, r_error] = calculate_residual(data, expected);
        residuals.SetBinContent(i, r);
        residuals.SetBinError(i, r_error);
      }
    }

    /**
     * \brief Summarize pull distribution with mean and standard deviation.
     * \return Pair of (mean, standard deviation) computed from pull bin contents.
     */
    std::pair<G4double, G4double> SummarizePullData()
    {
      G4double pull_mean = pulls.Integral() / pulls.GetNbinsX();
      G4double pull_stddev = 0.0;
      for (G4int i = 1; i <= pulls.GetNbinsX(); ++i)
      {
        pull_stddev += std::pow(pulls.GetBinContent(i) - pull_mean, 2.0);
      }
      pull_stddev /= pulls.GetNbinsX();
      return std::pair(pull_mean, std::sqrt(pull_stddev));
    }

    /**
     * \brief Persist all owned histograms into the current ROOT file directory.
     */
    void Write()
    {
      sample.Write();
      expect.Write();
      residuals.Write();
      pulls.Write();
    }

    TH1D sample;
    TH1D expect;
    TH1D residuals;
    TH1D pulls;
};

/**
 * \brief Composite ROOT canvas with pads for overlay, residuals, and pulls.
 *
 * The top pad shows expected and sampled distributions. The middle and
 * bottom pads show residuals and pulls with axis scaling adapted to pad height.
 */
struct ValidationHistogramPlotter
{
    TCanvas canvas;
    TPad overlay;
    TPad residual;
    TPad pull;

    ValidationHistogramPlotter() = delete;
    /**
     * \brief Construct a three-pad validation canvas.
     * \param id ROOT object identifier.
     * \param title Canvas title.
     * \param x Canvas width in pixels.
     * \param y Canvas height in pixels.
     */
    ValidationHistogramPlotter(const G4String& id, const G4String& title, G4int x, G4int y)
      : canvas(id, title, x, y),
        overlay((id + "_overlay").c_str(), "overlay", 0.0, 0.5, 1.0, 1.0),
        residual((id + "_residuals").c_str(), "residuals", 0.0, 0.25, 1.0, 0.5),
        pull((id + "_pulls").c_str(), "pulls", 0.0, 0.0, 1.0, 0.25)
    {
      canvas.SetMargin(0, 0, 0, 0);

      overlay.SetTopMargin(0.1);
      overlay.SetBottomMargin(0.1);
      overlay.SetLeftMargin(0.1);
      overlay.SetRightMargin(0.05);

      residual.SetTopMargin(0.15);
      residual.SetBottomMargin(0.0);
      residual.SetLeftMargin(0.1);
      residual.SetRightMargin(0.05);
      residual.SetGridy();

      pull.SetTopMargin(0.0);
      pull.SetBottomMargin(0.15);
      pull.SetLeftMargin(0.1);
      pull.SetRightMargin(0.05);
      pull.SetGridy();
    }

    /**
     * \brief Draw sampled, expected, residual, and pull histograms.
     * \param hist Histogram bundle to visualize.
     */
    void Draw(ValidationHistogram1D& hist)
    {
      canvas.cd();
      overlay.Draw();
      residual.Draw();
      pull.Draw();

      // Overlay - plot expectation then data so points of latter are on top.
      overlay.cd();
      hist.expect.Draw("HIST L");
      hist.sample.Draw("E1 SAME");

      // Because residual/pull pads are half-height, scale label/title/offset
      const G4double text_size = hist.sample.GetXaxis()->GetLabelSize();
      G4double tmp = hist.sample.GetYaxis()->GetTitleOffset();
      const G4double y_title_offset = (tmp == 0.0) ? 1.0 : tmp;

      // Residuals
      residual.cd();
      CenterYaxis(hist.residuals);
      hist.residuals.SetTitle(0);
      hist.residuals.GetYaxis()->SetTitleSize(2.0 * text_size);
      hist.residuals.GetYaxis()->SetTitleOffset(0.5 * y_title_offset);
      hist.residuals.GetYaxis()->SetLabelSize(2.0 * text_size);
      hist.residuals.GetYaxis()->CenterTitle(true);
      hist.residuals.Draw("E1");

      pull.cd();
      CenterYaxis(hist.pulls, 0.25);
      hist.pulls.SetTitle(0);
      hist.pulls.GetXaxis()->SetTitleSize(2.0 * text_size);
      hist.pulls.GetXaxis()->SetLabelSize(2.0 * text_size);
      hist.pulls.GetYaxis()->SetTitleSize(2.0 * text_size);
      hist.pulls.GetYaxis()->SetTitleOffset(0.5 * y_title_offset);
      hist.pulls.GetYaxis()->SetLabelSize(2.0 * text_size);
      hist.pulls.GetYaxis()->CenterTitle(true);
      hist.pulls.Draw("E1");
    }

    /**
     * \brief Save the canvas to a file format inferred from \p filename.
     * \param filename Output path, for example .svg or .pdf.
     */
    void SaveAs(const G4String& filename) { canvas.SaveAs(filename); }

    /**
     * \brief Persist the canvas into the current ROOT file directory.
     */
    void Write() { canvas.Write(); }

    /**
     * \brief Center Y-axis of a histogram on zero with symmetric limits
     * \param h Histogram to center.
     * \param padding_factor Percentage to add to Y-axis limit.
     */
    void CenterYaxis(TH1& h, double padding_factor = 0.1)
    {
      // 1. Find the maximum and minimum values currently in the histogram
      double y_max = h.GetMaximum();
      double y_min = h.GetMinimum();

      // 2. Determine which value is furthest from zero
      double max_deviation = std::max(std::abs(y_max), std::abs(y_min));

      // 3. Add a small padding so points don't clip the top/bottom frame lines
      double symmetric_limit = max_deviation * (1.0 + padding_factor);

      // 4. Force the symmetric range
      h.GetYaxis()->SetRangeUser(-symmetric_limit, symmetric_limit);
    }
};

//---------------------------------------------------------------------------//
// VALIDATION PROBLEM DATA
//---------------------------------------------------------------------------//
/**
 * \brief Kinematic constants and allowed energy ranges for muon decay at rest.
 *
 * The constructor ensures required particle definitions are available, then
 * computes masses, mass ratios, and energy boundaries used by
 * validation and theory helpers.
 *
 * \sa G4MuonDecayChannel
 */
struct MuonDecayParams
{
    /**
     * \brief Closed energy interval helper.
     */
    struct EnergyRange
    {
        G4double minimum_energy = 0.0;
        G4double maximum_energy = 0.0;

        /**
         * \brief Test whether a value lies outside this interval.
         * \param val Value to test.
         * \return true when val < minimum_energy or val > maximum_energy.
         */
        bool IsOutside(G4double val) const
        {
          return (val > maximum_energy || val < minimum_energy);
        }
    };

    /**
     * \brief Build and cache muon/electron masses and derived kinematic bounds.
     */
    MuonDecayParams()
    {
      // Make definitions available
      G4MuonMinus::Definition();
      G4MuonPlus::Definition();
      G4Electron::Definition();
      G4Positron::Definition();
      G4NeutrinoE::Definition();
      G4AntiNeutrinoE::Definition();
      G4NeutrinoMu::Definition();
      G4AntiNeutrinoMu::Definition();

      // Parameters
      muon_mass = G4MuonMinus::Definition()->GetPDGMass();
      electron_mass = G4Electron::Definition()->GetPDGMass();
      mass_ratio = electron_mass / muon_mass;
      mass_ratio_sq = mass_ratio * mass_ratio;
      electron = {electron_mass, 0.5 * muon_mass * (1.0 + mass_ratio_sq)};
      neutrino = {0.0, 0.5 * muon_mass * (1.0 - mass_ratio_sq)};
    }

    G4double muon_mass;
    G4double electron_mass;
    G4double mass_ratio;
    G4double mass_ratio_sq;
    EnergyRange electron;
    EnergyRange neutrino;
};

/**
 * \brief Analytic muon-decay energy probability density functions and callbacks.
 *
 * This helper builds normalized PDFs for electron and electron-antineutrino
 * energies in the muon-rest-frame parameterization and exposes callback
 * factories compatible with ValidationHistogram1D.
 *
 * \sa G4MuonDecayChannel
 */
class MuonDecayTheory
{
  public:

    MuonDecayTheory() = delete;
    /**
     * \brief Construct and normalize analytic PDFs from kinematic parameters.
     * \param p Precomputed muon-decay constants and energy ranges.
     */
    explicit MuonDecayTheory(const MuonDecayParams& p) : params(p)
    {
      // Change of variable/jacobian for x = (2/muon_mass)E and p(E) = |dx/dE|p(x)
      reduced_energy_jacobian = 2.0 / p.muon_mass;

      // Electron Energy PDF
      electron_reduced = {reduced_energy_jacobian * params.electron.minimum_energy,
                          reduced_energy_jacobian * params.electron.maximum_energy};

      electron_pdf_normalization = [this]() {
        const G4double r = this->params.mass_ratio;
        const G4double r2 = this->params.mass_ratio_sq;
        const G4double r4 = r2 * r2;
        const G4double r6 = r4 * r2;
        const G4double r8 = r6 * r2;

        return (1. / 12.) - (2. / 3.) * r2 + (2. / 3.) * r6 - (1. / 12.) * r8
               - 2. * r4 * std::log(r);
      }();

      electron_pdf = TF1(
        "neutrino_energy_pdf",
        [this](double* x, double*) { return this->electron_prob_density(x[0]); },
        params.electron.minimum_energy, params.electron.maximum_energy, 0);

      // Electron Antineutrino Energy PDF
      neutrino_reduced = {reduced_energy_jacobian * params.neutrino.minimum_energy,
                          reduced_energy_jacobian * params.neutrino.maximum_energy};

      neutrino_pdf_normalization = [this]() {
        const G4double a = this->neutrino_reduced.maximum_energy;
        const G4double a2 = a * a;
        const G4double a3 = a2 * a;
        const G4double a4 = a3 * a;
        return -(1. / 12.) * a4 - (1. / 3.) * a3 + (3. / 2.) * a2 - a
               - (1 - a) * (1 - a) * std::log(1 - a);
      }();

      neutrino_pdf = TF1(
        "neutrino_energy_pdf",
        [this](double* x, double*) { return this->neutrino_prob_density(x[0]); },
        params.neutrino.minimum_energy, params.neutrino.maximum_energy, 0);
    }

    /**
     * \brief Create expectation callback for electron-antineutrino energy bins.
     * \return Callback yielding per-bin expected count from the neutrino PDF integral.
     */
    ValidationHistogram1D::ExpectedCallback get_neutrino_pdf_callback()
    {
      return [this](G4int, G4double x_low, G4double x_high) {
        const double expected = this->neutrino_pdf.Integral(x_low, x_high);
        return ValidationHistogram1D::ExpectedBinValue{expected, 0.0, true};
      };
    }

    /**
     * \brief Create expectation callback for electron energy bins.
     * \return Callback yielding per-bin expected count from the electron PDF integral.
     */
    ValidationHistogram1D::ExpectedCallback get_electron_pdf_callback()
    {
      return [this](G4int, G4double x_low, G4double x_high) {
        const double expected = this->electron_pdf.Integral(x_low, x_high);
        return ValidationHistogram1D::ExpectedBinValue{expected, 0.0, true};
      };
    }

  private:

    /**
     * \brief Electron-antineutrino probability density in physical energy.
     * \param energy Physical energy (MeV).
     * \return Normalized \f$p(E_{\bar\nu_e})\f$.
     */
    G4double neutrino_prob_density(const double energy) const
    {
      const G4double y = reduced_energy_jacobian * energy;
      if (neutrino_reduced.IsOutside(y)) return 0.0;
      const G4double a = neutrino_reduced.maximum_energy;
      const G4double shape = (y * y * (a - y) * (a - y)) / (1.0 - y);
      return reduced_energy_jacobian * shape / neutrino_pdf_normalization;
    }

    /**
     * \brief Electron probability density in physical energy.
     * \param energy Physical energy (MeV).
     * \return Normalized \f$p(E_e)\f$.
     */
    G4double electron_prob_density(const double energy) const
    {
      const G4double x = reduced_energy_jacobian * energy;
      if (electron_reduced.IsOutside(x)) return 0.0;
      const G4double r2 = params.mass_ratio_sq;

      const G4double shape =
        std::sqrt(x * x - 4.0 * r2) * (0.5 * (1.0 + r2) * x - (1. / 3.) * x * x - (2. / 3.) * r2);

      return reduced_energy_jacobian * shape / electron_pdf_normalization;
    }

  private:

    const MuonDecayParams params;
    G4double reduced_energy_jacobian;

    MuonDecayParams::EnergyRange electron_reduced;
    G4double electron_pdf_normalization;
    TF1 electron_pdf;

    MuonDecayParams::EnergyRange neutrino_reduced;
    G4double neutrino_pdf_normalization;
    TF1 neutrino_pdf;
};

const G4DynamicParticle* FindDaughter(const G4DecayProducts* products,
                                      const G4ParticleDefinition* p)
{
  if (products == nullptr) return nullptr;
  for (G4int i = 0; i < products->entries(); ++i)
  {
    const G4DynamicParticle* daughter = (*products)[i];
    if (daughter == nullptr) continue;
    if (daughter->GetParticleDefinition() == p) return daughter;
  }
  return nullptr;
}

// Check under/overflows
bool HasSpillover(const TH1D& h, std::ostream& logger)
{
  const G4double underflow = h.GetBinContent(0);
  const G4double overflow = h.GetBinContent(h.GetNbinsX() + 1);

  bool kinematic_fail = false;
  if (underflow > 0.0)
  {
    kinematic_fail = true;
    logger << h.GetName() << ": underflow by " << underflow << "\n";
  }
  if (overflow > 0.0)
  {
    kinematic_fail = true;
    logger << h.GetName() << ": overflow by " << overflow << "\n";
  }
  return kinematic_fail;
};

//---------------------------------------------------------------------------//
// VALIDATION
//---------------------------------------------------------------------------//
int main(int argc, char** argv)
{
  long nEvents = 200000;
  if (argc > 1)
  {
    nEvents = std::strtol(argv[1], nullptr, 10);
  }

  if (nEvents <= 0)
  {
    std::cerr << "nEvents must be > 0" << std::endl;
    return 1;
  }

  //---------------------------------------------------------------------------//
  // Basic parameters
  CLHEP::HepRandom::setTheSeed(20659894);
  G4cout << "Input seed = " << CLHEP::HepRandom::getTheSeed() << G4endl;
  const std::size_t nBins = 120;

  const MuonDecayParams params;
  // Cannot be const because of callback chain and use of TF1::Integral.
  MuonDecayTheory theory{params};

  //---------------------------------------------------------------------------//
  // Create datasets
  ValidationHistogram1D electron_data("electron", "Energy (MeV)", "Counts", nBins,
                                      params.electron.minimum_energy,
                                      params.electron.maximum_energy);

  ValidationHistogram1D neutrino_data("neutrino", "Energy (MeV)", "Counts", nBins,
                                      params.neutrino.minimum_energy,
                                      params.neutrino.maximum_energy);

  // Sample the data
  G4MuonDecayChannel channel("mu-", 1.0);
  const auto electron_def = G4Electron::Definition();
  const auto neutrino_def = G4AntiNeutrinoE::Definition();

  for (long i = 0; i < nEvents; ++i)
  {
    G4DecayProducts* products = channel.DecayIt(0.0);
    const G4DynamicParticle* electron = FindDaughter(products, electron_def);
    const G4DynamicParticle* neutrino = FindDaughter(products, neutrino_def);
    if (electron != nullptr)
    {
      electron_data.FillSampleHistogram(electron->GetTotalEnergy());
    }
    if (neutrino != nullptr)
    {
      neutrino_data.FillSampleHistogram(neutrino->GetTotalEnergy());
    }
    delete products;
  }

  // Basic check that we don't have any samples outside kinematic boundaries
  if (HasSpillover(electron_data.sample, G4cerr)) return 1;
  if (HasSpillover(neutrino_data.sample, G4cerr)) return 1;

  // Create expected distributions
  neutrino_data.FillExpectedHistogram(theory.get_neutrino_pdf_callback(), nEvents);
  electron_data.FillExpectedHistogram(theory.get_electron_pdf_callback(), nEvents);

  // Create comparison plots
  neutrino_data.FillComparisonHistograms();
  electron_data.FillComparisonHistograms();

  // Summary data purely visual for now...
  // For now to confirm the Pull mean, sigma are ~0 and ~1 respectively.
  {
    auto const [mean, sigma] = neutrino_data.SummarizePullData();
    G4cout << "Neutrino pulls (mean, sigma) = (" << mean << ", " << sigma << ")" << G4endl;
  }
  {
    auto const [mean, sigma] = electron_data.SummarizePullData();
    G4cout << "Electron pulls (mean, sigma) = (" << mean << ", " << sigma << ")" << G4endl;
  }

  ValidationHistogramPlotter neutrino_plots("cMuonDecayValidationNuE",
                                            "G4MuonDecayChannel NuE Validation", 800, 900);
  neutrino_data.expect.SetLineColor(kGreen);
  neutrino_plots.Draw(neutrino_data);
  neutrino_plots.SaveAs("G4MuonDecayChannel_validation_neutrino.svg");

  ValidationHistogramPlotter electron_plots("cMuonDecayValidationElectron",
                                            "G4MuonDecayChannel E Validation", 800, 900);
  electron_data.expect.SetLineColor(kGreen);
  electron_plots.Draw(electron_data);
  electron_plots.SaveAs("G4MuonDecayChannel_validation_electron.svg");

  //---------------------------------------------------------------------------//
  // Persist data to ROOT file
  TFile outFile("G4MuonDecayChannel_validation.root", "RECREATE");
  neutrino_data.Write();
  neutrino_plots.Write();
  electron_data.Write();
  electron_plots.Write();
  outFile.Close();

  return 0;
}
