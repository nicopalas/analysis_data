#pragma once
#include "types.h"
#include "config.h"
#include "constants.h"
#include "TRandom3.h"
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TLine.h"

struct CrossSection {
    double sigma;
    double u_sigma;
};

static CrossSection cross_section_ratio(
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal_gold,
    const Vec3D& u_counts_signal_gold,
    const Vec3D& counts_signal_uranium,
    const Vec3D& u_counts_signal_uranium,
    const std::vector<double> eps_gold,
    const std::vector<double> u_eps_gold,
    const std::vector<double> eps_uranium,
    const std::vector<double> u_eps_uranium,
    const Vec2D& dOmega,
    int ebin,
    const AnalysisConfig& cfg)
{
    CrossSection result{0.0, 0.0};

    // --- valor nominal ---
    double sum_gold    = 0.0;
    double sum_uranium = 0.0;

    for (int j = 0; j < nbins_beam; j++) {
        for (int ii = 0; ii < nbins_det; ii++) {
            double dO = dOmega[j][ii];
            if (dO <= 0.0 || ii*dcos_det<cfg.cos_det_cut) continue;

            double eg = eps_gold[ii];
            double eu = eps_uranium[ii];
            if (eg <= 0.0 || eu <= 0.0 || counts_signal_gold[ebin][j][ii]<0 || counts_signal_uranium[ebin][j][ii]<0 ) continue;

            sum_gold    += counts_signal_gold[ebin][j][ii]    / (eg*dO);
            sum_uranium += counts_signal_uranium[ebin][j][ii] / (eu*dO);
        }
    }

    if (sum_uranium <= 0.0 || sum_gold <= 0.0) return result;

    result.sigma = sum_gold / sum_uranium;

    // --- bootstrapping ---
    TRandom3 rng(0);
    double n_bootstrap = cfg.n_toys;
    std::vector<double> bootstrap_ratios(n_bootstrap);

    for (int b = 0; b < n_bootstrap; b++) {
        double bs_gold    = 0.0;
        double bs_uranium = 0.0;

        for (int j = 0; j < nbins_beam; j++) {
            for (int ii = 0; ii < nbins_det; ii++) {
                double dO = dOmega[j][ii];
                if (dO <= 0.0) continue;

                double eg = eps_gold[ii];
                double eu = eps_uranium[ii];
                if (eg <= 0.0 || eu <= 0.0) continue;

                double ng = rng.Gaus(
                    counts_signal_gold[ebin][j][ii],
                    u_counts_signal_gold[ebin][j][ii]);

                double nu = rng.Gaus(
                    counts_signal_uranium[ebin][j][ii],
                    u_counts_signal_uranium[ebin][j][ii]);

                if (nu <0 || ng<0) continue;

                bs_gold    += ng / (eg * dO);
                bs_uranium += nu / (eu * dO);
            }
        }

        if (bs_uranium <= 0.0 || bs_gold <= 0.0) {
            bootstrap_ratios[b] = result.sigma;
            continue;
        }

        bootstrap_ratios[b] = bs_gold / bs_uranium;
    }

    // bootstrapping-
    double mean = 0.0;
    for (double v : bootstrap_ratios) mean += v;
    mean /= n_bootstrap;

    double var = 0.0;
    for (double v : bootstrap_ratios) var += (v - mean) * (v - mean);
    result.u_sigma = std::sqrt(var / (n_bootstrap - 1));

    return result;
}

bool readFlux(const std::string& fname,
              const std::string& hname,
              int nbins,
              const std::vector<double>& E_low,
              const std::vector<double>& E_high,
              std::vector<double>& flux,
              std::vector<double>& u_flux)
{
    flux.assign(nbins, 0.0);
    u_flux.assign(nbins, 0.0);

    TFile *f = TFile::Open(fname.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Cannot open flux file " << fname << "\n";
        return false;
    }

    TH1D *h = (TH1D*) f->Get(hname.c_str());
    if (!h) {
        std::cerr << "[ERROR] Cannot find hist '" << hname
                  << "' in " << fname << "\n";
        f->Close();
        return false;
    }

    h = (TH1D*) h->Clone("h_flux_clone");
    h->SetDirectory(nullptr);
    f->Close();

    for (int e = 0; e < nbins; e++) {

        // MeV->eV
        Int_t bin1 = h->FindBin(E_low[e]  * 1e6 + 1e-9);
        Int_t bin2 = h->FindBin(E_high[e] * 1e6 - 1e-9);

        bin1 = std::max(bin1, 1);
        bin2 = std::min(bin2, h->GetNbinsX());

        if (bin1 > bin2) {
            std::cerr << "[WARN] Ebin " << e
                      << " [" << E_low[e] << ", " << E_high[e]
                      << "] out of range.\n";
            continue;
        }

        double integral = 0.0;
        double err2_sum = 0.0;

        for (int b = bin1; b <= bin2; b++) {
            double E_center_b = h->GetBinCenter(b);   // eV
            double width      = h->GetBinWidth(b);    // eV
            double content    = h->GetBinContent(b);  // E·dΦ/dE
            double error      = h->GetBinError(b);

            if (E_center_b <= 0.0) continue;

            // from lethargic to normal
            integral  += (content / E_center_b) * width;
            err2_sum  += TMath::Power((error / E_center_b) * width, 2);
        }

        flux[e]   = integral;
        u_flux[e] = std::sqrt(err2_sum);
    }

    delete h;
    return true;
}

static CrossSection cross_section_absolute(
    const std::string&          flux_file,
    const std::string&          flux_hist,
    int                         nbins_beam,
    int                         nbins_det,
    const Vec3D&                counts_signal,
    const Vec3D&                u_counts_signal,
    const std::vector<double>&  eps,
    const std::vector<double>&  u_eps,
    const Vec2D&                dOmega,
    int                         ebin,
    const std::vector<double>&  E_low,
    const std::vector<double>&  E_high,
    double                      N_atoms,
    const AnalysisConfig&       cfg)
{
    CrossSection result{0.0, 0.0};

    std::vector<double> flux, u_flux;
    int nbins_e = (int)E_low.size();
    if (!readFlux(flux_file, flux_hist, nbins_e, E_low, E_high, flux, u_flux))
        return result;

    double Phi   = flux[ebin];
    double u_Phi = u_flux[ebin];
    if (Phi <= 0.0) return result;

    double sum_counts = 0.0;
    double u2_counts  = 0.0;

    for (int j = 0; j < nbins_beam; j++) {
        for (int ii = 0; ii < nbins_det; ii++) {
            double dO = dOmega[j][ii];
            if (dO <= 0.0 || ii * dcos_det < cfg.cos_det_cut) continue;

            double e  = eps[ii];
            double ue = u_eps[ii];
            if (e <= 0.0 || counts_signal[ebin][j][ii] < 0) continue;

            double c  = counts_signal[ebin][j][ii];
            double uc = u_counts_signal[ebin][j][ii];
            if (c <= 0.0) continue;

            double contrib = c / (e);
            sum_counts += contrib;

            // dc/c)^2 + (de/e)^2
            double rel2 = (uc * uc) / (c * c) + (ue * ue) / (e * e);
            u2_counts += contrib * contrib * rel2;

            printf("  ebin=%d j=%d ii=%d | c=%.4f+/-%.4f (rel=%.4f) | e=%.6f+/-%.6f (rel=%.4f) | contrib=%.4e +/- %.4e\n",
                ebin, j, ii,
                c, uc, uc/c,
                e, ue, ue/e,
                contrib, contrib * std::sqrt(rel2));
        }
    }

    if (sum_counts <= 0.0) return result;

    const double barn = 1.0e-24;
    result.sigma = sum_counts / (Phi * N_atoms * barn);

    double rel2_total = u2_counts / (sum_counts * sum_counts)
                      + (u_Phi * u_Phi) / (Phi * Phi);
    result.u_sigma = result.sigma * std::sqrt(rel2_total);

    printf("ebin=%d | sum_counts=%.6e +/- %.6e (rel=%.4f) | Phi=%.6e +/- %.6e (rel=%.4f) | sigma=%.6e +/- %.6e (rel=%.4f)\n",
        ebin,
        sum_counts, std::sqrt(u2_counts), std::sqrt(u2_counts) / sum_counts,
        Phi, u_Phi, u_Phi / Phi,
        result.sigma, result.u_sigma, result.u_sigma / result.sigma);

    return result;
}

void normalised_xs(TH1D* cross_section, double xmin, double xmax, double integral, TH1D* normalised_histo){
    int bin_min = cross_section->FindBin(xmin);
    int bin_max = cross_section->FindBin(xmax);
    double integral_exp = 0.;
    for (int i = bin_min; i<=bin_max; i++){
        integral_exp += cross_section->GetBinContent(i) * cross_section->GetBinWidth(i);
    }
    std::cout << "Exp integral = " << integral_exp << std::endl;
    std::cout << "Integral (barn*MeV) = " << integral << std::endl;
    std::cout << "Conversion factor = " << integral_exp/integral << std::endl;
    double conversion_factor = integral_exp / integral;

    for (int i = 1; i <= cross_section->GetNbinsX(); i++){
        double bin             = cross_section->GetBinContent(i);
        double new_bin_content = bin / conversion_factor;
        normalised_histo->SetBinContent(i, new_bin_content);
        normalised_histo->SetBinError  (i, cross_section->GetBinError(i) / conversion_factor);
    }

    // ── fix value at 9 MeV to 1.007 barn ────────────────────────────
    int bin_9MeV      = normalised_histo->FindBin(9.0);
    double val_9MeV   = normalised_histo->GetBinContent(bin_9MeV);
    double scale_9MeV = 1.007 / val_9MeV;

    std::cout << "Value at 9 MeV before = " << val_9MeV   << " barn\n";
    std::cout << "Scale factor at 9 MeV = " << scale_9MeV << "\n";

    for (int i = 1; i <= normalised_histo->GetNbinsX(); i++){
    normalised_histo->SetBinContent(i, normalised_histo->GetBinContent(i) * scale_9MeV);
    normalised_histo->SetBinError  (i, normalised_histo->GetBinError(i)   * scale_9MeV);
    std::cout << "E=" << normalised_histo->GetBinCenter(i) << " MeV"
              << "  sigma=" << normalised_histo->GetBinContent(i) << " barn"
              << " +/- "    << normalised_histo->GetBinError(i)   << " barn\n";
}

    std::cout << "Value at 9 MeV after  = "
              << normalised_histo->GetBinContent(bin_9MeV) << " barn\n";
}