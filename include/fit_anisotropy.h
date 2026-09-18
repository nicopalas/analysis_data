#pragma once
#include "types.h"
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include <cmath>
#include <iostream>
#include <vector>

struct LegendreResult {
    double A0   = 0., u_A0 = 0.;
    double a2   = 0., u_a2 = 0.;
    double a4   = 0., u_a4 = 0.;
    double W0   = 0.;             // ajuste evaluado en theta = 0°  (cos = 1)
    double W90  = 0.;             // ajuste evaluado en theta = 90° (cos = 0)
    double anisotropy   = 0.;     // W(0°)/W(90°)
    double u_anisotropy = 0.;
    std::vector<double> cos_theta, w, u_w;  
    double chi2ndf = 0.;
    int    ndf     = 0;
    bool   valid   = false;
};

// Polinomios de Legendre pares
inline double legP2(double c) { return (3.*c*c - 1.) / 2.; }
inline double legP4(double c) { double c2 = c*c; return (35.*c2*c2 - 30.*c2 + 3.) / 8.; }

// W(cos_theta) = A0 * (1 + a2*P2(cos_theta) + a4*P4(cos_theta))
static double legendre_W(double* x, double* par) {
    double c = x[0];
    return par[0] * (1. + par[1]*legP2(c) + par[2]*legP4(c));
}

static LegendreResult legendre_fit(
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal,
    const Vec2D& dOmega,
    int ebin,
    std::vector<double>& eps,
    std::vector<double>& u_eps,   
    const AnalysisConfig& cfg,
    bool fit_a4 = false)
{
    LegendreResult result{};

    // --- 1. Construir los puntos W(cos_theta) ---
    std::vector<double> cos_theta_vals, w_vals, u_w_vals;

    for (int i = 0; i < nbins_beam; i++) {
        // Abscisa: sqrt(<x^2>) en el bin, exacta para una función lineal en x^2
        double lo = i * dcos_beam;
        double hi = (i + 1) * dcos_beam;
        double cos_theta = std::sqrt((hi*hi*hi - lo*lo*lo) / (3. * dcos_beam));

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
        u_w_vals.push_back(std::sqrt(counts) / denom); // Poisson (si hay fondo restado, usar su varianza)
    }

    const int npar_free = fit_a4 ? 3 : 2;
    const int npts = static_cast<int>(cos_theta_vals.size());
    result.cos_theta = cos_theta_vals;
    result.w         = w_vals;
    result.u_w       = u_w_vals;
    if (npts <= npar_free) {
        std::cerr << "legendre_fit: too few points (" << npts
                  << ") for " << npar_free << " parameters at ebin=" << ebin << "\n";
        return result;
    }

    // --- 2. Ajuste a la serie de Legendre ---
    TGraphErrors gr(npts,
                    cos_theta_vals.data(),
                    w_vals.data(),
                    nullptr,
                    u_w_vals.data());

    TF1 f("f_leg", legendre_W, 0., 1., 3);
    f.SetParNames("A0", "a2", "a4");
    f.SetParameters(w_vals.front(), 0.1, 0.);
    if (!fit_a4) f.FixParameter(2, 0.);

    TFitResultPtr fit_res = gr.Fit(&f, "S Q R N");
    if (!fit_res.Get() || !fit_res->IsValid()) {
        std::cerr << "legendre_fit: did not converge at ebin=" << ebin << "\n";
        return result;
    }

    result.A0 = fit_res->Parameter(0);  result.u_A0 = fit_res->Error(0);
    result.a2 = fit_res->Parameter(1);  result.u_a2 = fit_res->Error(1);
    result.a4 = fit_res->Parameter(2);  result.u_a4 = fit_res->Error(2);
    result.ndf     = fit_res->Ndf();
    result.chi2ndf = (result.ndf > 0) ? fit_res->Chi2() / result.ndf : -1.;

    // --- 3. Evaluar el ajuste en 0° y 90° ---
    result.W0  = f.Eval(1.);   // cos(0°)  = 1
    result.W90 = f.Eval(0.);   // cos(90°) = 0

    if (result.W90 <= 0.) {
        std::cerr << "legendre_fit: W(90) <= 0 at ebin=" << ebin << "\n";
        return result;
    }
    result.anisotropy = result.W0 / result.W90;

    // --- 4. Error de R = W(0°)/W(90°) con la covarianza completa ---
    // R = N/D,  N = 1 + a2*P2(1) + a4*P4(1),  D = 1 + a2*P2(0) + a4*P4(0)
    // A0 se cancela -> dR/dA0 = 0
    // dR/dak = (Pk(1)*D - N*Pk(0)) / D^2
    double N = 1. + result.a2*legP2(1.) + result.a4*legP4(1.);
    double D = 1. + result.a2*legP2(0.) + result.a4*legP4(0.);

    double grad[3] = {
        0.,
        (legP2(1.)*D - N*legP2(0.)) / (D*D),
        (legP4(1.)*D - N*legP4(0.)) / (D*D)
    };

    TMatrixDSym cov = fit_res->GetCovarianceMatrix();  // parámetros fijos tienen fila/col = 0
    double var = 0.;
    for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            var += grad[k] * cov(k, l) * grad[l];

    result.u_anisotropy = std::sqrt(std::max(var, 0.));
    result.valid = true;

    std::cout << "ebin=" << ebin
              << "  a2=" << result.a2 << " +/- " << result.u_a2;
    if (fit_a4)
        std::cout << "  a4=" << result.a4 << " +/- " << result.u_a4;
    std::cout << "  W(0)=" << result.W0
              << "  W(90)=" << result.W90
              << "  W(0)/W(90)=" << result.anisotropy
              << " +/- " << result.u_anisotropy
              << "  chi2/ndf=" << result.chi2ndf << "\n";

    return result;
}