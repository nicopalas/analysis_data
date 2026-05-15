#pragma once 
#include "types.h" 
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "TRandom3.h"


struct AnisotropyResult{
    std::vector<double> w;
    std::vector<double> u_w;
};


static AnisotropyResult anisotropy (
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal,
    const Vec3D& u_counts_signal,
    const Vec2D& dOmega,
    int ebin,
    std::vector<double>& eps,
    std::vector<double>& u_eps,
    const AnalysisConfig& cfg)
{   
    AnisotropyResult result;
    int bin_90 = nbins_beam-1; // last bin
    std::vector<double> counts_theta;
    for (int i = 0; i<nbins_beam; i++){
        double counts = 0.;
        for (int j = 0; j<nbins_det; j++){
            if (j*dcos_det>= cfg.cos_det_cut && eps[j]>0 && dOmega[i][j]>0){
                counts+= counts_signal[ebin][i][j]/(eps[j]*dOmega[i][j]);
            }
        }
        counts_theta.push_back(counts);
    }
    for (auto& x:counts_theta) x/=counts_theta[bin_90];

    // bootstrapping
    TRandom3 rng(42 + ebin);
    Vec2D anisotropy_toys;
    for (int ntoy = 0; ntoy<cfg.n_toys; ntoy++){
        std::vector<double> eps_b(nbins_det,0.0);
        std::vector<double> counts_theta_toy;
        for (int j = 0; j<nbins_det;j++){
            eps_b[j] = rng.Gaus(eps[j],u_eps[j]);
            while (eps_b[j]<=0) eps_b[j]= rng.Gaus(eps[j],u_eps[j]);
        }
        for (int i = 0; i<nbins_beam; i++){
            double counts_b = 0.0;
                for (int j = 0 ; j < nbins_det; j++){
                    if (j*dcos_det>cfg.cos_det_cut && dOmega[i][j]>0){
                        double counts_signal_b = rng.Poisson (int(counts_signal[ebin][i][j]));
                        counts_b+=counts_signal_b/(eps_b[j]*dOmega[i][j]);
                    }
                }
                counts_theta_toy.push_back(counts_b);
            }
            for (auto& x: counts_theta_toy) x/=counts_theta_toy[bin_90];
        anisotropy_toys.push_back(counts_theta_toy); 
    }
    std::vector<double> mean_theta_toy(nbins_beam,0.0);
    std::vector<double> mean2_theta_toy(nbins_beam,0.0);
    for (int i = 0; i<cfg.n_toys; i++){
        for (int k = 0; k<nbins_beam; k++){
            mean_theta_toy[k]+=anisotropy_toys[i][k];
            mean2_theta_toy[k]+= anisotropy_toys[i][k]*anisotropy_toys[i][k];
            }
        }

    for (int k = 0; k<nbins_beam;k++){
        mean_theta_toy[k]/=cfg.n_toys;
        mean2_theta_toy[k]/=cfg.n_toys;
        }

    result.w = counts_theta;
    std::vector<double> u_anisotropy(nbins_beam,0.0);
    for (int i = 0; i<nbins_beam; i++){
        double var = mean2_theta_toy[i]-mean_theta_toy[i]*mean_theta_toy[i];
        u_anisotropy[i]= std::sqrt(var);
        }
    result.u_w = u_anisotropy;

    return result;
}