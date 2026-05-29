#pragma once
#include "types.h"
#include "config.h"
#include <cmath>

// ------------------signal_selection.h---------------------
// subtract the background in ROI
// for uranium: single background component, angular weight from counts_bkg
// for gold:    two components — continuum (counts_bkg) and uranium peak
//              (counts_upeak) — each subtracted with its own angular weight

static void computeSignal(
    const AnalysisConfig& cfg,
    const Vec3D& counts_roi,
    const Vec3D& counts_bkg,
    const Vec3D& counts_upeak,                      // uranium peak window (gold only)
    const std::vector<double>& counts_subtract_bkg,
    const std::vector<double>& u_counts_subtract_bkg,
    const std::vector<double>& counts_subtract_upeak,
    const std::vector<double>& u_counts_subtract_upeak,
    Vec3D& counts_signal,
    Vec3D& u_counts_signal)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    for(int i = 0; i < nbins; i++){

        double total_bkg   = 0.0;
        double total_upeak = 0.0;
        for(int j = 0; j < nbins_beam; j++)
            for(int ii = 0; ii < nbins_det; ii++){
                total_bkg   += counts_bkg[i][j][ii];
                total_upeak += counts_upeak[i][j][ii];
            }

        for(int j = 0; j < nbins_beam; j++){
            for(int ii = 0; ii < nbins_det; ii++){

                double w_bkg = (total_bkg > 0) ?
                    counts_bkg[i][j][ii] / total_bkg : 0.0;

                double sub   = w_bkg * counts_subtract_bkg[i];
                double u_sub2 = std::pow(w_bkg * u_counts_subtract_bkg[i], 2);

                if(cfg.sample == Sample::gold && total_upeak > 0){
                    double w_upeak = counts_upeak[i][j][ii] / total_upeak;
                    sub    += w_upeak * counts_subtract_upeak[i];
                    u_sub2 += std::pow(w_upeak * u_counts_subtract_upeak[i], 2);
                }

                counts_signal[i][j][ii]   = counts_roi[i][j][ii] - sub;
                u_counts_signal[i][j][ii] = std::sqrt(
                    counts_roi[i][j][ii] + u_sub2);
            }
        }
    }
}