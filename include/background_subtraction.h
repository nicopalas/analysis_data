#pragma once
#include "types.h"
#include "config.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include "Math/PdfFuncMathCore.h"

// ---------------------------------------------------------------------------
//  Model
//
//    f(x) = flat  +  N_sig * w * Gauss(x; mu_s, sigma_s)
//                 +  N_U   * w * CrystalBall(x; alpha, n, sigma_U, mu_U)   [gold only]
//
//  The pdfs are unit-normalised in x, so multiplying by the bin width w makes
//  f(x) a counts-per-bin prediction and the yields N_sig, N_U literal counts.
//  (The old code did not do this: f->Integral() was only a count if w == 1.)
//
//  Fixed parameter layout for every sample, so index bugs cannot creep in.
// ---------------------------------------------------------------------------
enum ParIdx {
    P_FLAT   = 0,   // continuum level, counts / bin
    P_NSIG   = 1,   // signal yield, counts (integrated over -inf..inf)
    P_MUSIG  = 2,
    P_SIGSIG = 3,
    P_NU     = 4,   // uranium peak yield, counts
    P_MUU    = 5,
    P_SIGU   = 6,
    P_ALPHA  = 7,   // CB tail onset, in units of sigma_U; sign fixes the tail side
    P_NCB    = 8    // CB power-law index, must be > 1
};

// --- geometry of the uranium peak, used only for seeding --------------------
static constexpr double kMuU_seed    = -7.0;
static constexpr double kSigU_seed   =  1.5;
static constexpr double kAlphaU_fix  =  1.0;   // > 0  ->  tail toward negative x
static constexpr double kNCB_fix     =  3.0;

// ---------------------------------------------------------------------------
static std::string getModelFormula(Sample sample, double w)
{
    char buf[256];

    std::snprintf(buf, sizeof(buf),
        "[1]*%.10g*ROOT::Math::gaussian_pdf(x,[3],[2])", w);
    const std::string sig = buf;

    std::snprintf(buf, sizeof(buf),
        "[4]*%.10g*ROOT::Math::crystalball_pdf(x,[7],[8],[6],[5])", w);
    const std::string cb = buf;

    switch (sample) {
        case Sample::uranium:
            // no uranium peak in the spectrum: keep the parameters present but
            // frozen at zero yield so the layout stays identical
            return "[0] + " + sig + " + 0.0*([4]+[5]+[6]+[7]+[8])";
        case Sample::gold:
            return "[0] + " + sig + " + " + cb;
    }
    return "[0] + " + sig;
}

// robust seed for the flat level: median of the bins that belong to neither peak
static double seedFlatLevel(TH1D* h, double roi_min, double roi_max)
{
    std::vector<double> side;
    const double u_lo = kMuU_seed - 4.0 * kSigU_seed;
    const double u_hi = kMuU_seed + 4.0 * kSigU_seed;

    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double x = h->GetBinCenter(b);
        if (x > roi_min && x < roi_max) continue;
        if (x > u_lo    && x < u_hi   ) continue;
        side.push_back(h->GetBinContent(b));
    }
    if (side.empty()) return std::max(1e-3, h->GetBinContent(1));

    std::nth_element(side.begin(), side.begin() + side.size() / 2, side.end());
    return std::max(1e-3, side[side.size() / 2]);
}

// excess counts above a flat level in [lo, hi]
static double seedYield(TH1D* h, double lo, double hi, double flat)
{
    double s = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double x = h->GetBinCenter(b);
        if (x < lo || x > hi) continue;
        s += h->GetBinContent(b) - flat;
    }
    return std::max(1.0, s);
}

// ---------------------------------------------------------------------------
static void setModelParameters(TF1* f, Sample sample, TH1D* h,
                               double roi_min, double roi_max,
                               bool fix_cb_tails = true)
{
    const double w    = h->GetBinWidth(1);
    const double flat = seedFlatLevel(h, roi_min, roi_max);

    const double mu_s_seed  = 0.5 * (roi_min + roi_max);
    const double sig_s_seed = std::max(w, (roi_max - roi_min) / 6.0);
    const double n_s_seed   = seedYield(h, roi_min, roi_max, flat);

    f->SetParName(P_FLAT,   "flat");
    f->SetParName(P_NSIG,   "N_sig");
    f->SetParName(P_MUSIG,  "mu_sig");
    f->SetParName(P_SIGSIG, "sigma_sig");
    f->SetParName(P_NU,     "N_U");
    f->SetParName(P_MUU,    "mu_U");
    f->SetParName(P_SIGU,   "sigma_U");
    f->SetParName(P_ALPHA,  "alpha");
    f->SetParName(P_NCB,    "n");

    // --- continuum ---------------------------------------------------------
    f->SetParameter(P_FLAT, flat);
    f->SetParLimits(P_FLAT, 0.0, 1e6);

    // --- signal ------------------------------------------------------------
    // N_sig is allowed to go negative: a hard lower bound at zero biases the
    // estimator and makes the uncertainty meaningless when the signal is small.
    f->SetParameter(P_NSIG, n_s_seed);
    f->SetParLimits(P_NSIG, -10.0 * std::fabs(n_s_seed) - 100.0, 1e7);

    f->SetParameter(P_MUSIG, mu_s_seed);
    f->SetParLimits(P_MUSIG, roi_min, roi_max);

    f->SetParameter(P_SIGSIG, sig_s_seed);
    f->SetParLimits(P_SIGSIG, 0.5 * w, 0.5 * (roi_max - roi_min));

    // --- uranium peak ------------------------------------------------------
    if (sample == Sample::gold) {
        const double n_u_seed = seedYield(h,
                                          kMuU_seed - 3.0 * kSigU_seed,
                                          kMuU_seed + 3.0 * kSigU_seed,
                                          flat);
        f->SetParameter(P_NU, n_u_seed);
        f->SetParLimits(P_NU, 0.0, 1e7);

        f->SetParameter(P_MUU, kMuU_seed);
        f->SetParLimits(P_MUU, -9.0, -5.0);

        f->SetParameter(P_SIGU, kSigU_seed);
        f->SetParLimits(P_SIGU, 0.5, 2.5);

        // alpha and n are almost unconstrained by the data and are strongly
        // anticorrelated with each other and with the flat level. Fix them to
        // simulation values and vary them afterwards as a systematic.
        f->SetParameter(P_ALPHA, kAlphaU_fix);
        f->SetParameter(P_NCB,   kNCB_fix);
        if (fix_cb_tails) {
            f->FixParameter(P_ALPHA, kAlphaU_fix);
            f->FixParameter(P_NCB,   kNCB_fix);
        } else {
            // never let alpha cross zero: the sign is a discrete flip of the
            // tail side, not a direction MINUIT can walk along
            f->SetParLimits(P_ALPHA, 0.2, 5.0);
            f->SetParLimits(P_NCB,   1.01, 20.0);
        }
    } else {
        f->FixParameter(P_NU,    0.0);
        f->FixParameter(P_MUU,   kMuU_seed);
        f->FixParameter(P_SIGU,  kSigU_seed);
        f->FixParameter(P_ALPHA, kAlphaU_fix);
        f->FixParameter(P_NCB,   kNCB_fix);
    }
}

// ---------------------------------------------------------------------------
//  Component integrals inside the ROI, in counts.
//  Clone the full model and switch the other components off, so the parameter
//  indices are written down in exactly one place.
// ---------------------------------------------------------------------------
struct RoiCounts {
    double flat  = 0.0;
    double upeak = 0.0;
    double sig   = 0.0;
    double total = 0.0;
};

static RoiCounts decomposeRoi(TF1* f, double w, double roi_min, double roi_max)
{
    RoiCounts c;

    c.flat = f->GetParameter(P_FLAT) * (roi_max - roi_min) / w;

    std::unique_ptr<TF1> tmp(static_cast<TF1*>(f->Clone("f_component_tmp")));

    // signal only
    tmp->SetParameter(P_FLAT, 0.0);
    tmp->SetParameter(P_NU,   0.0);
    c.sig = tmp->Integral(roi_min, roi_max, 1e-6) / w;

    // uranium peak only
    tmp->SetParameter(P_NSIG, 0.0);
    tmp->SetParameter(P_NU,   f->GetParameter(P_NU));
    c.upeak = tmp->Integral(roi_min, roi_max, 1e-6) / w;

    c.total = c.flat + c.sig + c.upeak;
    return c;
}

// ---------------------------------------------------------------------------
struct BackgroundFit {
    // background under the ROI (flat + uranium peak), in counts
    double counts_subtract         = 0.0;
    double u_counts_subtract       = 0.0;
    double counts_subtract_bkg     = 0.0;
    double u_counts_subtract_bkg   = 0.0;
    double counts_subtract_upeak   = 0.0;
    double u_counts_subtract_upeak = 0.0;

    // signal
    double n_sig_val      = 0.0;   // total yield from the fit
    double n_sig_err      = 0.0;   // parabolic error
    double n_sig_err_lo   = 0.0;   // MINOS
    double n_sig_err_hi   = 0.0;
    double n_sig_roi      = 0.0;   // fraction of the yield inside the ROI
    double n_sig_boot     = 0.0;   // bootstrap cross-check of n_sig_err

    double chi2ndf        = 0.0;
    int    n_toys_ok      = 0;
    bool   converged      = false;

    TF1*  func            = nullptr;
    TH1D* hist_subtracted = nullptr;   // data - (flat + uranium peak)
    TFitResultPtr fit_result;
};

// ---------------------------------------------------------------------------
static BackgroundFit fitBackground(
    const AnalysisConfig& cfg,
    TH1D* h,
    double roi_min,
    double roi_max,
    int ebin)
{
    BackgroundFit result;

    const int    ntot = h->GetNbinsX();
    const double xmin = h->GetXaxis()->GetXmin();
    const double xmax = h->GetXaxis()->GetXmax();
    const double w    = h->GetBinWidth(1);

    const std::string formula = getModelFormula(cfg.sample, w);

    TF1* f = new TF1(Form("model_%s_%d", cfg.output_tag.c_str(), ebin),
                     formula.c_str(), xmin, xmax);
    setModelParameters(f, cfg.sample, h, roi_min, roi_max);

    // Binned extended maximum likelihood over the FULL range.
    //   L  Poisson likelihood, correct for low-occupancy bins
    //   S  return the fit result
    //   E  MINOS errors on the free parameters
    //   Q  quiet
    result.fit_result = h->Fit(f, "L S E R Q", "", xmin, xmax);
    result.converged  = (result.fit_result.Get() &&
                         result.fit_result->IsValid());

    // --- goodness of fit: Pearson chi2 over the full range ------------------
    double chi2 = 0.0;
    int    ndf  = 0;
    for (int b = 1; b <= ntot; ++b) {
        const double x   = h->GetBinCenter(b);
        const double obs = h->GetBinContent(b);
        const double exp = f->Eval(x);
        if (exp <= 0.0) continue;
        chi2 += (obs - exp) * (obs - exp) / exp;
        ndf++;
    }
    ndf -= f->GetNumberFreeParameters();
    result.chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    // --- signal -------------------------------------------------------------
    result.n_sig_val = f->GetParameter(P_NSIG);
    result.n_sig_err = f->GetParError(P_NSIG);
    if (result.fit_result.Get()) {
        result.n_sig_err_lo = result.fit_result->LowerError(P_NSIG);
        result.n_sig_err_hi = result.fit_result->UpperError(P_NSIG);
    }

    // --- background under the ROI ------------------------------------------
    const RoiCounts c = decomposeRoi(f, w, roi_min, roi_max);
    result.counts_subtract_bkg   = c.flat;
    result.counts_subtract_upeak = c.upeak;
    result.counts_subtract       = c.flat + c.upeak;
    result.n_sig_roi             = (std::fabs(result.n_sig_val) > 0.0)
                                 ? c.sig / result.n_sig_val : 0.0;

    // -----------------------------------------------------------------------
    //  Bootstrap. With a proper likelihood fit the errors already come from
    //  MINOS; the toys are a cross-check and the source of the bin-by-bin
    //  errors on the subtracted histogram. A large disagreement between
    //  n_sig_err and n_sig_boot means the likelihood is not parabolic.
    // -----------------------------------------------------------------------
    const int ntoys       = cfg.n_toys;
    const int max_attempts = 4 * ntoys;   // hard cap: never spin forever

    std::vector<double> toy_sub (ntot, 0.0);
    std::vector<double> toy_sub2(ntot, 0.0);
    std::vector<double> toy_bkg, toy_upeak, toy_total, toy_nsig;
    toy_bkg.reserve(ntoys); toy_upeak.reserve(ntoys);
    toy_total.reserve(ntoys); toy_nsig.reserve(ntoys);

    TRandom3 rng(42 + ebin);

    for (int attempt = 0; attempt < max_attempts &&
                          (int)toy_nsig.size() < ntoys; ++attempt) {

        std::unique_ptr<TH1D> h_toy(
            static_cast<TH1D*>(h->Clone(Form("h_toy_tmp_%d", attempt))));
        h_toy->SetDirectory(nullptr);

        for (int b = 1; b <= ntot; ++b) {
            const double fluct = rng.Poisson(h->GetBinContent(b));
            h_toy->SetBinContent(b, fluct);
            h_toy->SetBinError(b, fluct > 0 ? std::sqrt(fluct) : 1.0);
        }

        std::unique_ptr<TF1> f_toy(new TF1(Form("f_toy_tmp_%d", attempt),
                                           formula.c_str(), xmin, xmax));
        setModelParameters(f_toy.get(), cfg.sample, h_toy.get(),
                           roi_min, roi_max);
        for (int p = 0; p < f->GetNpar(); ++p)
            f_toy->SetParameter(p, f->GetParameter(p));

        // identical estimator to the nominal fit, minus MINOS
        TFitResultPtr r = h_toy->Fit(f_toy.get(), "L S R Q N 0", "", xmin, xmax);
        if (!r.Get() || !r->IsValid()) continue;

        const RoiCounts ct = decomposeRoi(f_toy.get(), w, roi_min, roi_max);
        toy_bkg  .push_back(ct.flat);
        toy_upeak.push_back(ct.upeak);
        toy_total.push_back(ct.flat + ct.upeak);
        toy_nsig .push_back(f_toy->GetParameter(P_NSIG));

        // subtract only the background components, keep the signal visible
        std::unique_ptr<TF1> f_bkg(
            static_cast<TF1*>(f_toy->Clone(Form("f_bkg_tmp_%d", attempt))));
        f_bkg->SetParameter(P_NSIG, 0.0);

        for (int b = 1; b <= ntot; ++b) {
            const double x   = h_toy->GetBinCenter(b);
            const double sub = h_toy->GetBinContent(b) - f_bkg->Eval(x);
            toy_sub [b-1] += sub;
            toy_sub2[b-1] += sub * sub;
        }
    }

    result.n_toys_ok = (int)toy_nsig.size();

    auto boot_std = [](const std::vector<double>& v) -> double {
        if (v.size() < 2) return 0.0;
        double mean = 0.0;
        for (double x : v) mean += x;
        mean /= v.size();
        double var = 0.0;
        for (double x : v) var += (x - mean) * (x - mean);
        return std::sqrt(var / (v.size() - 1));
    };

    result.u_counts_subtract       = boot_std(toy_total);
    result.u_counts_subtract_bkg   = boot_std(toy_bkg);
    result.u_counts_subtract_upeak = boot_std(toy_upeak);
    result.n_sig_boot              = boot_std(toy_nsig);

    // --- background-subtracted histogram ------------------------------------
    result.hist_subtracted = static_cast<TH1D*>(
        h->Clone(Form("hsub_%s_%d", cfg.output_tag.c_str(), ebin)));

    const int n_ok = std::max(1, result.n_toys_ok);
    for (int b = 0; b < ntot; ++b) {
        const double mean_s = toy_sub[b] / n_ok;
        const double var    = toy_sub2[b] / n_ok - mean_s * mean_s;
        result.hist_subtracted->SetBinContent(b + 1, mean_s);
        result.hist_subtracted->SetBinError(b + 1, var > 0 ? std::sqrt(var) : 0.0);
    }

    result.func = f;
    return result;
}