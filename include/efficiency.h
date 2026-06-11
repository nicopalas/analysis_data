#pragma once 
#include "types.h" 
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "TRandom3.h"


struct EfficiencyResult{
    std::vector<double> eps;
    std::vector<double>u_eps;
};


static EfficiencyResult computeEfficiency(
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal,
    const Vec3D& u_counts_signal,
    int ebin,
    double plateau_tol = 0.05,
    int n_bootstrap = 1000)
{
    EfficiencyResult efficiency;
    efficiency.eps.resize(nbins_det, 0.0);
    efficiency.u_eps.resize(nbins_det, 0.0);

    TRandom3 rng(0);

    // ── central value: suma bruta sobre bin_beam ──────────────────────────────
    for (int i = 0; i < nbins_det; i++) {
        double raw = 0.0;
        for (int j = 0; j < nbins_beam; j++)
            raw += counts_signal[ebin][j][i];
        efficiency.eps[i] = raw;
    }

    // ── plateau detection ─────────────────────────────────────────────────────
    double plateau_mean = 1.0;
    {
        std::vector<double> smooth(nbins_det, 0.0);
        for (int i = 1; i < nbins_det - 1; i++)
            smooth[i] = (efficiency.eps[i-1] + efficiency.eps[i] + efficiency.eps[i+1]) / 3.0;
        smooth[0]           = efficiency.eps[0];
        smooth[nbins_det-1] = efficiency.eps[nbins_det-1];

        int i_max = std::max_element(smooth.begin(), smooth.end()) - smooth.begin();
        double val_max = smooth[i_max];

        int plateau_start = i_max;
        for (int i = i_max - 1; i >= 0; i--) {
            if (smooth[i] < val_max * (1.0 - plateau_tol)) break;
            plateau_start = i;
        }

        double sum = 0.0; int cnt = 0;
        for (int i = plateau_start; i <= i_max; i++) {
            if (efficiency.eps[i] > 0) { sum += efficiency.eps[i]; cnt++; }
        }
        plateau_mean = (cnt > 0) ? sum / cnt : 1.0;

        std::cout << "[INFO] ebin=" << ebin
                  << "  plateau bins=[" << plateau_start << "," << i_max << "]"
                  << "  plateau_mean=" << plateau_mean << "\n";
    }

    for (int i = 0; i < nbins_det; i++)
        efficiency.eps[i] = (plateau_mean > 0) ? efficiency.eps[i] / plateau_mean : 0.0;

    // ── bootstrap ─────────────────────────────────────────────────────────────
    std::vector<std::vector<double>> bootstrap_eps(n_bootstrap,
                                                    std::vector<double>(nbins_det, 0.0));
    for (int b = 0; b < n_bootstrap; b++) {
        for (int i = 0; i < nbins_det; i++) {
            double raw = 0.0;
            for (int j = 0; j < nbins_beam; j++)
                raw += rng.Gaus(counts_signal[ebin][j][i],
                                u_counts_signal[ebin][j][i]);
            bootstrap_eps[b][i] = (plateau_mean > 0) ? raw / plateau_mean : 0.0;
        }
    }

    for (int i = 0; i < nbins_det; i++) {
        double mean = 0.0;
        for (int b = 0; b < n_bootstrap; b++) mean += bootstrap_eps[b][i];
        mean /= n_bootstrap;
        double var = 0.0;
        for (int b = 0; b < n_bootstrap; b++)
            var += std::pow(bootstrap_eps[b][i] - mean, 2);
        efficiency.u_eps[i] = std::sqrt(var / n_bootstrap);
    }

    return efficiency;
}