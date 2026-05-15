#pragma once 
#include <cmath>

// ------------------signal_selection.h---------------------
// subtract the background in ROI
// take into account the angular distribution


using Vec3D = std::vector<std::vector<std::vector<double>>>;

static void signalSelection(
    const AnalysisConfig& cfg,
    const Vec3D& counts_roi,
    const Vec3D& counts_bkg,
    const std::vector<double>& counts_subtract,
    const std::vector<double>& u_counts_subtract,
    Vec3D& counts_signal,
    Vec3D& u_counts_signal)
    {
    int nbins = (int)cfg.energy_bins.size()-1;
    for (int i = 0; i<nbins; i++){
        double totals = 0.0;
        for (int j = 0; j<nbins_beam; j++){
            for (int ii = 0; ii<nbins_det; ii++){
                totals+=counts_bkg[i][j][ii];
            }
        }
        for (int j = 0; j<nbins_beam; ++j){
            for (int ii = 0; ii<nbins_det; ii+){
                double w = counts_bkg[i][j][ii]/totals;
                counts_signal[i][j][ii] = counts_roi[i][j][ii]-w*counts_subtract[i][j][ii];
                u_counts_signal[i][j][ii] = std::sqrt(counts_roi[e][j][ii]+std::pow(w*u_counts_subtract,2));
            }
        }
    }
}