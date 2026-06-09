#pragma once
#include <string>
#include "config.h"

struct EventCuts {
    double roi_min, roi_max;
    double amp_min, amp_max;
    double ratio_max;
    double upeak_min = 0.0;   // uranium contamination window
    double upeak_max = 0.0;   // zero = not applicable
};

static EventCuts getCuts(Sample sample, double neutron_energy){
    switch(sample){
        case Sample::uranium:
            if(neutron_energy > 1000 && neutron_energy < 1500) return {-6.,  6., 9e3, 50e3, 0.3, 0.0, 0.0};
            if(neutron_energy > 500 && neutron_energy < 1000) return {-10.,  8., 8e3, 50e3, 0.3, 0.0, 0.0};
            if(neutron_energy > 100)                          return {-8.,  10.,  8e3, 50e3, 0.3, 0.0, 0.0};
            if(neutron_energy > 10)                           return {-7.,  10.,  7e3, 50e3, 0.4, 0.0, 0.0};
            return                                                   {-8.,  10.,  6e3, 50e3, 0.3, 0.0, 0.0};

        case Sample::gold:
            if (neutron_energy>1000 && neutron_energy<1500) return {-2.0, 4.0, 17e3, 35e3, 0.1, -13., -4.0};
            if(neutron_energy > 600 && neutron_energy < 1000) return {-3.0, 4.0, 16e3, 38e3, 0.2, -15.0, -4.0};
            if(neutron_energy > 200)                          return {-2.5, 4.0, 17e3, 37e3, 0.2, -15.0, -4.0};
            return                                                   {-2.0, 3.0, 17e3, 36e3, 0.2, -15.0, -4.0};
    }
    return {-8., 15., 6e3, 39e3, 0.3, 0.0, 0.0};
}

static bool passAmplitudeCut(float a8, float a9, const EventCuts& c){
    double sum   = a8 + a9;
    double ratio = (a9 - a8) / (a9 + a8);
    if(sum < c.amp_min || sum > c.amp_max || ratio > c.ratio_max) return false;
    return true;
}

static bool inUraniumPeak(double dt, const EventCuts& c){
    if(c.upeak_min == 0.0 && c.upeak_max == 0.0) return false;
    return dt >= c.upeak_min && dt <= c.upeak_max;
}