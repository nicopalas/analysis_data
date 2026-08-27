#pragma once
#include "types.h" 
#include "config.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <vector>
#include <string>
#include <cmath>
#include "Math/PdfFuncMathCore.h"  



struct BackgroundFit {
    double counts_subtract        = 0.0;
    double u_counts_subtract      = 0.0;
    double counts_subtract_bkg    = 0.0;
    double u_counts_subtract_bkg  = 0.0;
    double counts_subtract_upeak  = 0.0;
    double u_counts_subtract_upeak= 0.0;
    double chi2ndf                = 0.0;
    TF1*   func                   = nullptr;  // conservado para compatibilidad
    TH1D*  hist_subtracted        = nullptr;
    // nuevo: resultado RooFit
    RooFitResult* roo_result      = nullptr;
    double n_sig_val              = 0.0;
    double n_sig_err              = 0.0;
};

static std::string getBackgroundFormula(Sample sample){
    switch(sample){
        case Sample::uranium: return "[0]";
        case Sample::gold:    return "[0] + [1]*ROOT::Math::crystalball_pdf(x,[4],[5],[3],[2])";
    }
    return "[0]+[1]*x";
}

static void setBackgroundParameters(TF1* f, Sample sample, TH1D* h){
    switch(sample){
        case Sample::uranium:
            f->SetParLimits(0, 0.0, 1e6);
            f->SetParameter(0, std::max(1e-3, h->GetBinContent(h->FindBin(0.0))));
            break;

        case Sample::gold: {
            // rough seed for the flat level: average of two sidebands far from the peak
            const double lvl = 0.5*( h->GetBinContent(h->FindBin(-15.0))
                                   + h->GetBinContent(h->FindBin( 12.0)) );

            // rough seed for the yield: peak height x sigma x sqrt(2pi)
            const double pk   = std::max(1.0, h->GetBinContent(h->FindBin(-6.0)) - lvl);
            const double sig0 = 1.5;
            const double yld  = pk * sig0 * TMath::Sqrt(TMath::TwoPi()) * h->GetBinWidth(1);

            f->SetParLimits(0,  0.0,   1e6);   // flat level
            f->SetParLimits(1,  0.0,   1e7);   // yield >= 0
            f->SetParLimits(2, -7.0,  -5.0);   // mean
            f->SetParLimits(3,  0.5,   2.0);   // sigma
            f->SetParLimits(4,  0.2,   5.0);   // alpha != 0
            f->SetParLimits(5,  1.01, 20.0);   // n > 1 mandatory

            f->SetParameters(std::max(1e-3, lvl), yld, -6.0, sig0, 1.0, 3.0);

            f->SetParName(0, "flat");
            f->SetParName(1, "N_U");
            f->SetParName(2, "mu_U");
            f->SetParName(3, "sigma_U");
            f->SetParName(4, "alpha");
            f->SetParName(5, "n");
            break;
        }
    }
}

// decompose gold fit into continuum and uranium peak integrals
static void decomposeGoldIntegrals(
    TF1* f,
    double roi_min, double roi_max,
    double& counts_bkg,
    double& counts_upeak)
{
    // flat continuum: [0]
    counts_bkg = f->GetParameter(0) * (roi_max - roi_min);

    // uranium peak: clone the full model and switch off the continuum,
    // so the parameter indices never have to be repeated here
    std::unique_ptr<TF1> f_upeak(static_cast<TF1*>(f->Clone("f_upeak_decomp")));
    f_upeak->SetParameter(0, 0.0);
    counts_upeak = f_upeak->Integral(roi_min, roi_max);
}


static BackgroundFit fitBackground(
    const AnalysisConfig& cfg,
    TH1D* h,
    double roi_min,
    double roi_max,
    int ebin)
{
    BackgroundFit result;
    int    ntot = h->GetNbinsX();
    double xmin = h->GetBinCenter(1);
    double xmax = h->GetBinCenter(ntot);

    std::string formula = getBackgroundFormula(cfg.sample);
    TF1* f = new TF1(Form("bkg_%s_%d", cfg.output_tag.c_str(), ebin),
                     formula.c_str(), xmin, xmax);
    setBackgroundParameters(f, cfg.sample, h);

    TH1D* h_tails = (TH1D*)h->Clone(Form("h_tails_%s_%d", cfg.output_tag.c_str(), ebin));
    for(int b = 1; b <= ntot; ++b){
        double x = h_tails->GetBinCenter(b);
        if(x > roi_min && x < roi_max){
            h_tails->SetBinContent(b, 0);
            h_tails->SetBinError(b, 1e10);
        }
    }
    h_tails->Fit(f, "IRS", "", xmin, xmax);
    delete h_tails;

    // chi2
    double chi2 = 0.0;
    int    ndf  = 0;
    for(int b = 1; b <= ntot; ++b){
        double x   = h->GetBinCenter(b);
        double obs = h->GetBinContent(b);
        double err = h->GetBinError(b);
        if(err <= 0) continue;
        if(x > roi_min && x < roi_max) continue;
        chi2 += std::pow((obs - f->Eval(x)) / err, 2);
        ndf++;
    }
    ndf -= f->GetNumberFreeParameters();
    result.chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    // total nominal integral
    result.counts_subtract = f->Integral(roi_min, roi_max);

    // decompose for gold
    if(cfg.sample == Sample::gold){
        decomposeGoldIntegrals(f, roi_min, roi_max,
            result.counts_subtract_bkg,
            result.counts_subtract_upeak);
    } else {
        result.counts_subtract_bkg   = result.counts_subtract;
        result.u_counts_subtract_bkg = 0.0;  // filled by bootstrap below
    }

    // bootstrap
    int ntoys = cfg.n_toys;
    std::vector<double> toy_sub(ntot, 0.0);
    std::vector<double> toy_sub2(ntot, 0.0);
    std::vector<double> toy_integrals(ntoys, 0.0);
    std::vector<double> toy_integrals_bkg(ntoys, 0.0);
    std::vector<double> toy_integrals_upeak(ntoys, 0.0);

    TRandom3 rng(42 + ebin);

    for(int itoy = 0; itoy < ntoys; ++itoy){
        TH1D* h_toy = (TH1D*)h->Clone(Form("h_toy_tmp_%d", itoy));
        for(int b = 1; b <= ntot; ++b){
            double fluct = rng.Poisson(h->GetBinContent(b));
            h_toy->SetBinContent(b, fluct);
            h_toy->SetBinError(b, fluct > 0 ? std::sqrt(fluct) : 1.0);
        }

        TH1D* h_toy_tails = (TH1D*)h_toy->Clone(Form("h_toy_tails_tmp_%d", itoy));
        for(int b = 1; b <= ntot; ++b){
            double x = h_toy_tails->GetBinCenter(b);
            if(x > roi_min && x < roi_max){
                h_toy_tails->SetBinContent(b, 0.0);
                h_toy_tails->SetBinError(b, 1e10);
            }
        }

        TF1* f_toy = new TF1(Form("f_toy_tmp_%d", itoy),
                             formula.c_str(), xmin, xmax);
        setBackgroundParameters(f_toy, cfg.sample, h_toy);
        for(int p = 0; p < f->GetNpar(); ++p)
            f_toy->SetParameter(p, f->GetParameter(p));

        auto fit_status = h_toy_tails->Fit(f_toy, "RS Q 0 N", "", xmin, xmax);
        delete h_toy_tails;

        // Saltar toy si el fit no convergió
        if (fit_status != 0 && fit_status != 4000) {
            delete h_toy;
            delete f_toy;
            --itoy;  // reintentar, o simplemente continuar
            continue;
        }

        toy_integrals[itoy] = f_toy->Integral(roi_min, roi_max);

        if(cfg.sample == Sample::gold){
            double cb, ucb, cup, ucup;
            decomposeGoldIntegrals(f_toy, roi_min, roi_max, cb, cup);
            toy_integrals_bkg[itoy]   = cb;
            toy_integrals_upeak[itoy] = cup;
        } else {
            toy_integrals_bkg[itoy] = toy_integrals[itoy];
        }

        for(int b = 1; b <= ntot; ++b){
            double x   = h_toy->GetBinCenter(b);
            double sub = h_toy->GetBinContent(b) - f_toy->Eval(x);
            toy_sub[b-1]  += sub;
            toy_sub2[b-1] += sub * sub;
        }

        delete h_toy;
        delete f_toy;
    }

    // uncertainty from bootstrap std
    auto bootstrap_std = [&](const std::vector<double>& vals) -> double {
        double mean = 0.0;
        for(double v : vals) mean += v;
        mean /= ntoys;
        double var = 0.0;
        for(double v : vals) var += (v - mean) * (v - mean);
        return std::sqrt(var / (ntoys - 1));
    };

    result.u_counts_subtract     = bootstrap_std(toy_integrals);
    result.u_counts_subtract_bkg = bootstrap_std(toy_integrals_bkg);
    if(cfg.sample == Sample::gold)
        result.u_counts_subtract_upeak = bootstrap_std(toy_integrals_upeak);

    // subtracted histogram
    result.hist_subtracted = (TH1D*)h->Clone(
        Form("hclone_%s_%d", cfg.output_tag.c_str(), ebin));
    for(int b = 0; b < ntot; ++b){
        double mean_s = toy_sub[b]  / ntoys;
        double var    = toy_sub2[b] / ntoys - mean_s * mean_s;
        result.hist_subtracted->SetBinContent(b+1, mean_s);
        result.hist_subtracted->SetBinError(b+1, var > 0 ? std::sqrt(var) : 0.0);
    }

    result.func = f;
    return result;
}

