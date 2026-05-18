#pragma once
#include "types.h"
#include "config.h"
#include "constants.h"
#include "TRandom3.h"
#include <cmath>
#include <vector>

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
            if (dO <= 0.0) continue;

            double eg = eps_gold[ii];
            double eu = eps_uranium[ii];
            if (eg <= 0.0 || eu <= 0.0 || counts_signal_gold[ebin][j][ii]<0 || counts_signal_uranium[ebin][j][ii]<0 ) continue;

            sum_gold    += counts_signal_gold[ebin][j][ii]    / (eg * dO);
            sum_uranium += counts_signal_uranium[ebin][j][ii] / (eu * dO);
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

    // --- std del bootstrap como incertidumbre ---
    double mean = 0.0;
    for (double v : bootstrap_ratios) mean += v;
    mean /= n_bootstrap;

    double var = 0.0;
    for (double v : bootstrap_ratios) var += (v - mean) * (v - mean);
    result.u_sigma = std::sqrt(var / (n_bootstrap - 1));

    return result;
}