#pragma once
#include <string>
#include "config.h"


// -----------------------cuts.h---------------------
// define the cuts to distinguish background and fission regions


struct EventCuts {
    double roi_min, roi_max;
    double amp_min, amp_max;
    double ratio_max;
};
static EventCuts getCuts(Sample sample, double neutron_energy){  // ← use energy not index
    switch(sample){
        case Sample::uranium:
            if(neutron_energy > 500 && neutron_energy<1000) return {-8.,  10.,  13e3, 33e3, 0.2};
            if(neutron_energy > 100) return {-8., 12.,  8e3, 35e3, 0.3};
            if(neutron_energy > 10)  return {-8.,  12.,  8e3, 36e3, 0.3};
            return                          {-8.,  16.,  6e3, 39e3, 0.3};

        case Sample::gold:
            if(neutron_energy > 600 && neutron_energy<1000) return {-3.0, 3.5, 18e3, 38e3, 0.3};
            if(neutron_energy > 300) return {-2.5, 3.0, 18e3, 34e3, 0.2};
            return                          {-2.0, 3.0, 20e3, 36e3, 0.2};
    }
    return {-8., 15., 6e3, 39e3, 0.3};  // fallback
}
static bool passAmplitudeCut(float a8, float a9, const EventCuts& c){
    double sum = a8+a9;
    double ratio = (a9-a8)/(a9+a8);
    if (sum<c.amp_min || sum>c.amp_max || ratio>c.ratio_max) return false;
    return true;
}