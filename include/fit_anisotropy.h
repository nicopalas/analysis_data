#pragma once
#include "types.h"
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

struct LegendreResult {
    double a2;
    double u_a2;
    double anisotropy;    // W(0°)/W(90°) = (1 + a2) / (1 - a2/2)
    double u_anisotropy;
    double chi2ndf;
};

// W(cos_theta) = norm * (1 + a2*P2(cos_theta))
// P2(x) = (3x^2 - 1) / 2
static double legendre_W2(double* x, double* par) {
    double cos_theta = x[0];
    double P2 = (3.*cos_theta*cos_theta - 1.) / 2.;
    return par[0] * (1. + par[1] * P2);
}

static LegendreResult legendre_fit(
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal,
    const Vec2D& dOmega,
    int ebin,
    std::vector<double>& eps,
    std::vector<double>& u_eps,
    const AnalysisConfig& cfg)
{
    LegendreResult result{};

    // --- 1. Build points ---
    std::vector<double> cos_theta_vals, w_vals, u_w_vals;

    for (int i = 0; i < nbins_beam; i++) {
        double cos_theta = (i + 0.5) * dcos_beam; // adjust to your convention
        double counts = 0., denom = 0.;

        for (int j = 0; j < nbins_det; j++) {
            if (j * dcos_det > cfg.cos_det_cut
                && eps[j] > 0.01
                && dOmega[i][j] > 1e-6) {
                double sig = counts_signal[ebin][i][j];
                if (sig < 0) continue;
                counts += sig;
                denom  += eps[j] * dOmega[i][j];
            }
        }

        if (denom < 1e-10 || counts <= 0) continue;

        cos_theta_vals.push_back(cos_theta);
        w_vals.push_back(counts / denom);
        u_w_vals.push_back(std::sqrt(counts) / denom); // Poisson
    }

    int npts = cos_theta_vals.size();
    if (npts < 2) {
        std::cerr << "legendre_fit: too few points at ebin=" << ebin << "\n";
        return result;
    }

    // --- 2. Fit ---
    TGraphErrors gr(npts,
                    cos_theta_vals.data(),
                    w_vals.data(),
                    nullptr,
                    u_w_vals.data());

    TF1 f("f_leg2", legendre_W2, 0., 1., 2);
    f.SetParNames("norm", "a2");
    f.SetParameter(0, w_vals[npts / 2]);
    f.SetParameter(1, 0.1);

    TFitResultPtr fit_res = gr.Fit(&f, "S Q R");

    if (!fit_res.Get() || !fit_res->IsValid()) {
        std::cerr << "legendre_fit: did not converge at ebin=" << ebin << "\n";
        return result;
    }

    result.a2    = fit_res->Parameter(1);
    result.u_a2  = fit_res->Error(1);
    result.chi2ndf = fit_res->Chi2() / fit_res->Ndf();

    // --- 3. Anisotropy W(0°)/W(90°), norm cancels ---
    // W(0°)  = norm*(1 + a2*P2(1)) = norm*(1 + a2)
    // W(90°) = norm*(1 + a2*P2(0)) = norm*(1 - a2/2)
    double num = 1. + result.a2;
    double den = 1. - result.a2 / 2.;
    result.anisotropy = num / den;

    // dR/da2 = (3/2) / den^2
    double dR_da2 = 1.5 / (den * den);
    result.u_anisotropy = std::abs(dR_da2) * result.u_a2;

    std::cout << "ebin=" << ebin
              << "  a2=" << result.a2 << " +/- " << result.u_a2
              << "  W(0)/W(90)=" << result.anisotropy
              << " +/- " << result.u_anisotropy
              << "  chi2/ndf=" << result.chi2ndf << "\n";

    return result;
}