#pragma once 
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "TRandom3.h"


struct EfficiencyResult{
    std::vector<double> eps;
    std::vector<double>u_eps;
}

using Vec3D = std::vector<std::vector<std::vector<double>>>;

static EfficiencyResult computeEfficiency(
    int ref_bin,
    Vec3D& counts_signal,
    Vec3D& u_counts_signal,
    Vec3D& dOmega_eff,
    int ebin
){
    EfficiencyResult efficiency;
    efficiency.eps[ref_bin] = 1.0;
    efficiency.u_eps[ref_bin] = 0.0;
    for (int i = ref_bin-1; i>=0; i--){
        double num = 0.0, u_num2 = 0.0;
        double denom = 0.0, u_denom2 = 0.0;
        double eff = 0.0;
        for (int j = 0; j<nbins_beam; j++){
            
            double Ni =counts_signal[ebin][j][i];
            double Nii = counts_signal[ebin][j][i+1];
            double wii = dOmega_eff[ebin][j][i];
            double wi = dOmega_eff[ebin][j][i+1];

            if(Ni>0 && Nii>0 && wii>0 && wi>0){
                num+=Nii/wii;
                denom+=Ni/wi;
                u_num2+=std::pow(u_counts_signal[ebin][j][i]/wii,2);
                u_denom2+=std::pow(u_counts_signal[ebin][j][i+1]/wi,2);
            }
        }
        if (num>0 && denom>0){
            efficiency.eps[i] = efficiency.eps[i+1]*num/denom;
            double rel2 = 0.0;
            if (efficiency.eps[i]>0 && efficiency.eps[i+1]>0){
                rel2+=std::pow(efficiency.u_eps[i]/efficiency.eps[i+1],2);
            }
            rel2+=u_num2/(num*num);
            rel2+=u_denom2/(denom*denom);
            efficiency.u_eps[i] = efficiency.eps[i+1] * std::sqrt(2);
        }
    }
    return efficiency;
}

