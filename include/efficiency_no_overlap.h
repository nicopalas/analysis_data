#pragma once
#include "types.h"
#include "config.h"
#include "utils.h"
#include "constants.h"
#include "efficiency.h"   // reusa EfficiencyResult
#include <cmath>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────
// computeEfficiencyNoOverlap
//
// Misma lógica que computeEfficiency() (cadena recursiva bin i vs i+1,
// normalizada por dOmega_eff, anclada en ref_bin=1), PERO sin exigir que
// haya solape físico de aceptancia entre bins vecinos: se usa el bin i y el
// bin i+1 de forma independiente (detector homogéneo asumido), sin la
// condición de continuidad de solape que sí impone la versión "overlap".
//
// En la práctica esto significa: se elimina el requisito de que wii y wi
// (dOmega de bins vecinos) representen una región de aceptancia compartida
// -- se combinan directamente aunque los bins sean disjuntos y no compartan
// campo de visión de detector.
// ─────────────────────────────────────────────────────────────────────────
static EfficiencyResult computeEfficiencyNoOverlap(
    int ref_bin,
    int nbins_beam,
    int nbins_det,
    const Vec3D& counts_signal,
    const Vec3D& u_counts_signal,
    const Vec2D& dOmega_eff,
    int ebin)
{
    EfficiencyResult efficiency;
    efficiency.eps.resize(nbins_det, 0.0);
    efficiency.u_eps.resize(nbins_det, 0.0);
    efficiency.eps[ref_bin] = 1.0;
    efficiency.u_eps[ref_bin] = 0.0;

    for (int i = ref_bin - 1; i >= 0; i--) {
        double num = 0.0, u_num2 = 0.0;
        double denom = 0.0, u_denom2 = 0.0;

        for (int j = 0; j < nbins_beam; j++) {

            double Ni  = counts_signal[ebin][j][i + 1];
            double Nii = counts_signal[ebin][j][i];
            double wii = dOmega_eff[j][i];
            double wi  = dOmega_eff[j][i + 1];

            // Único cambio respecto a computeEfficiency(): no se exige que
            // wii y wi correspondan a una región de solape compartida entre
            // los bins i e i+1 -- se toma cada bin de forma independiente
            // (detector homogéneo), siempre que haya señal y aceptancia > 0.
                num    += Nii;
                denom  += Ni;
                u_num2   += std::pow(u_counts_signal[ebin][j][i], 2);
                u_denom2 += std::pow(u_counts_signal[ebin][j][i + 1] ,  2);
        }

        if (num > 0 && denom > 0) {
            efficiency.eps[i] = efficiency.eps[i + 1] * num / denom;

            double rel2 = 0.0;
            if (efficiency.eps[i + 1] > 0) {
                rel2 += std::pow(efficiency.u_eps[i + 1] / efficiency.eps[i + 1], 2);
            }
            rel2 += u_num2   / (num   * num);
            rel2 += u_denom2 / (denom * denom);

            efficiency.u_eps[i] = efficiency.eps[i] * std::sqrt(rel2);
        }
    }
    return efficiency;
}