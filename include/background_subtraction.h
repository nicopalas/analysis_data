#pragma once
#include "config.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <vector>
#include <string>
#include <cmath>

struct BackgroundFit {
    double counts_subtract  = 0.0;
    double u_counts_subtract = 0.0;
    TF1*   func             = nullptr;
    TH1D*  hist_subtracted  = nullptr;
}; 

// -------------------------------------------------------
// background model depends on sample and energy
// uranium: simple linear (one peak, flat background)
// gold:    two gaussians + polynomial (U peak + Au peak + continuum)
// -------------------------------------------------------
static std::string getBackgroundFormula(Sample sample){
    switch(sample){
        case Sample::Uranium: return "[0]+[1]*x";
        case Sample::Gold:    return "[0]+[1]*x+[2]*x*x"
                                     " + [3]*TMath::Gaus(x,[4],[5],1)"
                                     " + [6]*TMath::Gaus(x,[7],[8],1)";
    }
    return "[0]+[1]*x";
}

static void setBackgroundParameters(TF1* f, Sample sample, TH1D* h){
    switch(sample){
        case Sample::Uranium:
            f->SetParLimits(0, 0, 1e5);
            break;
        case Sample::Gold:
            f->SetParLimits(0,  0.0,  1e6);
            f->SetParLimits(1, -5.0,  5.0);
            f->SetParLimits(2, -1.0,  1.0);
            f->SetParLimits(3,  0.0,  1e6);  // A_U
            f->SetParLimits(4, -15.0, -5.0); // mu_U
            f->SetParLimits(5,  0.5,  4.0);  // sigma_U
            f->SetParLimits(6,  0.0,  1e6);  // A_Au
            f->SetParLimits(7, -2.0,  2.0);  // mu_Au
            f->SetParLimits(8,  0.5,  4.0);  // sigma_Au
            f->SetParameters(
                1.0, 0.0, 1e-3,
                std::max(1.0, h->GetBinContent(h->FindBin(-10.0))), -10.0, 1.5,
                std::max(1.0, h->GetBinContent(h->FindBin(  0.0))),   0.0, 1.5
            );
            break;
    }
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

    // --- build and fit background model ---
    std::string formula = getBackgroundFormula(cfg.sample);
    TF1* f = new TF1(Form("bkg_%s_%d", cfg.output_tag.c_str(), ebin),
                     formula.c_str(), xmin, xmax);
    setBackgroundParameters(f, cfg.sample, h);

    TH1D* h_tails = (TH1D*)h->Clone("h_tails_tmp");  
    for(int b = 1; b <= ntot; ++b){                
        double x = h_tails->GetBinCenter(b);    
        if(x > roi_min && x < roi_max){
            h_tails->SetBinContent(b, 0);
            h_tails->SetBinError(b, 1e10);
        }
    }
    h_tails->Fit(f, "IRS Q", "", xmin, xmax);
    delete h_tails;

    // --- nominal background integral ---
    result.counts_subtract = f->Integral(roi_min, roi_max);

    // --- bootstrap ---
    int ntoys = cfg.n_toys;
    std::vector<double> toy_sub(ntot, 0.0);
    std::vector<double> toy_sub2(ntot, 0.0);
    std::vector<double> toy_integrals(ntoys, 0.0);

    TRandom3 rng(42 + ebin);

    for(int itoy = 0; itoy < ntoys; ++itoy){
        TH1D* h_toy = (TH1D*)h->Clone(Form("h_toy_tmp_%d", itoy));
        for(int b = 1; b <= ntot; ++b){
            double fluct = rng.Poisson(h->GetBinContent(b));
            h_toy->SetBinContent(b, fluct);            
            h_toy->SetBinError(b, std::sqrt(fluct));  
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
        // warm start from nominal fit
        for(int p = 0; p < f->GetNpar(); ++p)
            f_toy->SetParameter(p, f->GetParameter(p));

        h_toy_tails->Fit(f_toy, "IRS Q", "", xmin, xmax);
        delete h_toy_tails;

        // integral of background in ROI for this toy
        toy_integrals[itoy] = f_toy->Integral(roi_min, roi_max); 

        // bin-by-bin subtracted signal
        for(int b = 1; b <= ntot; ++b){
            double x       = h_toy->GetBinCenter(b);
            double sub     = h_toy->GetBinContent(b) - f_toy->Eval(x);  
            toy_sub[b-1]  += sub;
            toy_sub2[b-1] += sub * sub; 
        }

        delete h_toy;
        delete f_toy;
    }

    // --- uncertainty from std of toy integrals ---
    double mean_int = 0.0;
    for(double v : toy_integrals) mean_int += v;
    mean_int /= ntoys;

    double var_int = 0.0;
    for(double v : toy_integrals) var_int += (v - mean_int) * (v - mean_int);
    result.u_counts_subtract = std::sqrt(var_int / (ntoys - 1)); 

    // --- subtracted histogram ---
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