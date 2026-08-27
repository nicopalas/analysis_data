#pragma once
#include <string>
#include <cmath>
#include <algorithm>
#include "config.h"
#include "cut_variations.h"

struct EventCuts {
    double roi_min, roi_max;
    double amp_min, amp_max;
    double ratio_max, ratio_min;
    double upeak_min = 0.0;   // uranium contamination window
    double upeak_max = 0.0;   // zero = not applicable
};

// ---------------------------------------------------------------------------
// NOMINAL cuts — byte-for-byte the original table. Do not edit this when
// studying systematics; edit the CutVariation instead, so the nominal result
// is always reproducible.
// ---------------------------------------------------------------------------
static EventCuts getCutsNominal(Sample sample, double neutron_energy){
    switch(sample){
        case Sample::uranium:
            if(neutron_energy >= 1000)              return {-5.,  5.,  10e3, 35e3, 0.1, -0.9, 0.0, 0.0};
            if(neutron_energy >= 500)               return {-4.,  4.,  13e3, 35e3, 0., -0.9, 0.0, 0.0};
            if(neutron_energy >= 100)               return {-4.,  5.,  13e3, 35e3, 0.1, -.9, 0.0, 0.0};
            if(neutron_energy >= 10)                return {-5.5, 6.,  10e3, 35e3, 0.1, -.9, 0.0, 0.0};
            return                                         {-5.,  6.,   9e3, 38e3, 0.1, -1, 0.0, 0.0};

        case Sample::gold:
            if(neutron_energy >= 1000)              return {-3.0, 3.2, 18e3, 37e3, 0.3, -0.7,  -14.5, -3.0};
            if(neutron_energy >= 600)               return {-2., 3., 18e3, 40e3, 0.4, -0.5,  -14.5, -3.0};
            if(neutron_energy >= 300)               return {-2., 2.5, 20e3, 40e3, 0.4, -0.5,  -14.5, -3.0};
            if(neutron_energy >= 150)               return {-2.,  2.3, 20e3, 40e3, 0.3, -0.5, -14., -3.0};
            return                                         {-2.,  2.2, 21e3, 40e3, 0.2, -0.4, -14., -2.5};
        case Sample::uranium_mc:
            if(neutron_energy >= 500)               return {-4.,  4.,  10e3, 35e3, 0.3, 0.8, 0.0, 0.0};
            if(neutron_energy >= 100)               return {-5.,  5.,  11e3, 35e3, 0.3, 0.9, 0.0, 0.0};
            if(neutron_energy >= 10)                return {-4.5, 9.,  10e3, 34e3, 0.3, 0.9, 0.0, 0.0};
            return                                         {-7.,  9.,   8e3, 38e3, 0.3, 0.9, 0.0, 0.0};
    }
    return {-8., 15., 6e3, 39e3, 0.3, 0.0, 0.0};
}

// ---------------------------------------------------------------------------
// Smooth (log-E interpolated) version of the same table. Used to quantify the
// artefact introduced by the step function at 10/100/500 MeV (U) and
// 150/300/600 MeV (Au), which your 50 log bins straddle.
// ---------------------------------------------------------------------------
static EventCuts getCutsSmooth(Sample sample, double E){
    // anchor energies = geometric centres of the nominal step regions
    static const double anchU[5] = {3.162, 31.62, 223.6, 707.1, 1414.};
    static const double anchA[5] = {12.25, 212.1, 424.3, 774.6, 1414.};
    const double* anchor = (sample == Sample::uranium) ? anchU : anchA;

    EventCuts c[5];
    for(int k = 0; k < 5; ++k) c[k] = getCutsNominal(sample, anchor[k]);

    double x = std::log(std::max(E, 1e-3));
    if(x <= std::log(anchor[0])) return c[0];
    if(x >= std::log(anchor[4])) return c[4];

    int k = 0;
    while(k < 3 && x > std::log(anchor[k+1])) ++k;
    double x0 = std::log(anchor[k]), x1 = std::log(anchor[k+1]);
    double t  = (x - x0) / (x1 - x0);

    EventCuts o;
    o.roi_min   = c[k].roi_min   + t * (c[k+1].roi_min   - c[k].roi_min);
    o.roi_max   = c[k].roi_max   + t * (c[k+1].roi_max   - c[k].roi_max);
    o.amp_min   = c[k].amp_min   + t * (c[k+1].amp_min   - c[k].amp_min);
    o.amp_max   = c[k].amp_max   + t * (c[k+1].amp_max   - c[k].amp_max);
    o.ratio_max = c[k].ratio_max + t * (c[k+1].ratio_max - c[k].ratio_max);
    o.upeak_min = c[k].upeak_min + t * (c[k+1].upeak_min - c[k].upeak_min);
    o.upeak_max = c[k].upeak_max + t * (c[k+1].upeak_max - c[k].upeak_max);
    return o;
}

// ---------------------------------------------------------------------------
// Public entry point: nominal cuts with the active variation applied.
// ---------------------------------------------------------------------------
static EventCuts getCuts(Sample sample, double neutron_energy){
    const CutVariation& v = cutVar();

    const bool apply = (v.only_sample < 0) ||
                       (v.only_sample == static_cast<int>(sample));

    // shifting the energy thresholds == evaluating the table at a rescaled E
    double Eeval = apply ? neutron_energy / v.e_threshold_shift : neutron_energy;

    EventCuts c = (apply && v.smooth_energy) ? getCutsSmooth(sample, Eeval)
                                             : getCutsNominal(sample, Eeval);
    if(!apply) return c;

    // ROI: scale the half-width about the centre, then shift the edges
    double roi_c = 0.5 * (c.roi_min + c.roi_max);
    double roi_h = 0.5 * (c.roi_max - c.roi_min) * v.f_roi_width;
    c.roi_min = roi_c - roi_h + v.d_roi_min;
    c.roi_max = roi_c + roi_h + v.d_roi_max;

    c.amp_min   = c.amp_min * v.f_amp_min + v.d_amp_min;
    c.amp_max   = c.amp_max * v.f_amp_max + v.d_amp_max;
    c.ratio_max = c.ratio_max + v.d_ratio_max;

    if(c.upeak_min != 0.0 || c.upeak_max != 0.0){
        c.upeak_min += v.d_upeak_min;
        c.upeak_max += v.d_upeak_max;
    }

    // guards
    if(c.roi_max <= c.roi_min) c.roi_max = c.roi_min + 1e-3;
    if(c.amp_max <= c.amp_min) c.amp_max = c.amp_min + 1.0;
    if(c.ratio_max < 0.0)      c.ratio_max = 0.0;
    return c;
}

static bool passAmplitudeCut(float a8, float a9, const EventCuts& c){
    double sum   = a8 + a9;
    double ratio = (a9 - a8) / (a9 + a8);
    if(sum < c.amp_min || sum > c.amp_max) return false;
    if(cutVar().symmetric_ratio){
        if(ratio < c.ratio_min || ratio > c.ratio_max) return false;
    } else {
        if(ratio < c.ratio_min || ratio > c.ratio_max) return false;
    }
    return true;
}

static bool inUraniumPeak(double dt, const EventCuts& c){
    if(c.upeak_min == 0.0 && c.upeak_max == 0.0) return false;
    return dt >= c.upeak_min && dt <= c.upeak_max;
}