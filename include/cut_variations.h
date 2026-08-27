#pragma once
#include <string>
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Knobs applied on top of the nominal getCuts(). One global instance is read
// by getCuts()/passAmplitudeCut(), so a variation automatically propagates to
// EVERY call site (fillHistograms, the efficiency pass, the anisotropy pass).
// That is the point: a selection change must move the efficiency and the
// signal together, otherwise the systematic is grossly overestimated.
// ---------------------------------------------------------------------------
struct CutVariation {
    std::string name = "nominal";

    // multiplicative, applied first (1.0 = nominal)
    double f_amp_min   = 1.0;
    double f_amp_max   = 1.0;
    double f_roi_width = 1.0;   // scales the ROI half-width about its centre

    // additive, applied second (same units as the cut)
    double d_roi_min   = 0.0;   // ns
    double d_roi_max   = 0.0;   // ns
    double d_amp_min   = 0.0;   // ADC
    double d_amp_max   = 0.0;   // ADC
    double d_ratio_max = 0.0;
    double d_upeak_min = 0.0;   // ns (gold only)
    double d_upeak_max = 0.0;   // ns (gold only)

    // structural switches
    bool   symmetric_ratio   = false; // use |ratio| < ratio_max instead of one-sided
    bool   smooth_energy     = false; // interpolate cuts in log(E) instead of steps
    double e_threshold_shift = 1.0;   // scales the energy thresholds inside getCuts

    // which sample the variation applies to (-1 = all). Lets you decouple the
    // U-238 threshold from the Au reference threshold, which do NOT cancel in
    // the ratio.
    int    only_sample = -1;          // cast of Sample, or -1
};

inline CutVariation& cutVar(){
    static CutVariation v;
    return v;
}

// Convenience builders -------------------------------------------------------
inline CutVariation varNominal(){ CutVariation v; v.name = "nominal"; return v; }

inline CutVariation varAmpMin(double frac, int sample = -1, const char* tag = ""){
    CutVariation v;
    v.name = std::string("amp_min_x") + std::to_string(frac) + tag;
    v.f_amp_min = frac;
    v.only_sample = sample;
    return v;
}

inline CutVariation varAmpMax(double frac, int sample = -1){
    CutVariation v;
    v.name = std::string("amp_max_x") + std::to_string(frac);
    v.f_amp_max = frac;
    v.only_sample = sample;
    return v;
}

inline CutVariation varRatio(double dr){
    CutVariation v;
    v.name = std::string("ratio_max_") + std::to_string(dr);
    v.d_ratio_max = dr;
    return v;
}



inline CutVariation varRoiWidth(double frac){
    CutVariation v;
    v.name = std::string("roi_width_x") + std::to_string(frac);
    v.f_roi_width = frac;
    return v;
}

inline CutVariation varRoiShift(double dt){
    CutVariation v;
    v.name = std::string("roi_shift_") + std::to_string(dt);
    v.d_roi_min = dt;
    v.d_roi_max = dt;
    return v;
}

inline CutVariation varEThreshold(double frac){
    CutVariation v;
    v.name = std::string("e_threshold_x") + std::to_string(frac);
    v.e_threshold_shift = frac;
    return v;
}

inline CutVariation varSmoothE(){
    CutVariation v;
    v.name = "smooth_cut_vs_E";
    v.smooth_energy = true;
    return v;
}

inline CutVariation varUPeak(double dlo, double dhi){
    CutVariation v;
    v.name = std::string("upeak_") + std::to_string(dlo) + "_" + std::to_string(dhi);
    v.d_upeak_min = dlo;
    v.d_upeak_max = dhi;
    return v;
}