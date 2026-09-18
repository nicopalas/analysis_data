#pragma once
#include "types.h"
#include "config.h"
#include "accidental_coincidences.h"
#include <cmath>

// ------------------signal_selection.h---------------------
// subtract the background in ROI
// for uranium: single background component, angular weight from counts_bkg
// for gold:    two components — acc coincidences (counts_bkg) and uranium peak
//              (counts_upeak) — each subtracted with its own angular weight

static void computeSignal(
    const AnalysisConfig& cfg,
    const Vec3D& counts_roi,
    const Vec3D& counts_acc,      // accidentals file, ROI dt window, ROI cuts
    const Vec3D& counts_upeak,    // gold only: real correlated component
    const std::vector<double>& counts_subtract_bkg,
    const std::vector<double>& u_counts_subtract_bkg,
    const std::vector<double>& counts_subtract_upeak,
    const std::vector<double>& u_counts_subtract_upeak,
    Vec3D& counts_signal,
    Vec3D& u_counts_signal)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    for (int i = 0; i < nbins; ++i) {

        double total_acc = 0.0, total_upeak = 0.0;
        for (int j = 0; j < nbins_beam; ++j)
            for (int ii = 0; ii < nbins_det; ++ii) {
                total_acc   += counts_acc[i][j][ii];
                total_upeak += counts_upeak[i][j][ii];
            }

        const double B  = counts_subtract_bkg[i];
        const double uB = u_counts_subtract_bkg[i];

        for (int j = 0; j < nbins_beam; ++j) {
            for (int ii = 0; ii < nbins_det; ++ii) {

                double w = (total_acc > 0) ? counts_acc[i][j][ii] / total_acc : 0.0;

                double sub    = w * B;
                double u_sub2 = w * w * uB * uB;
                if (total_acc > 0)
                    u_sub2 += B * B * w * (1.0 - w) / total_acc;   // -> 0 as stats grow

                if (cfg.sample == Sample::gold && total_upeak > 0) {
                    double wu = counts_upeak[i][j][ii] / total_upeak;
                    double U  = counts_subtract_upeak[i];
                    sub    += wu * U;
                    u_sub2 += wu * wu * u_counts_subtract_upeak[i]
                                     * u_counts_subtract_upeak[i];
                    u_sub2 += U * U * wu * (1.0 - wu) / total_upeak;
                }

                counts_signal[i][j][ii]   = counts_roi[i][j][ii] - sub;
                u_counts_signal[i][j][ii] = std::sqrt(counts_roi[i][j][ii] + u_sub2);
            }
        }
    }
}