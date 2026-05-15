#pragma once
#include <string>


// -----------------------cuts.h---------------------
// define the cuts to distinguish background and fission regions


struct EventCuts {
    double roi_min, roi_max;
    double amp_min, amp_max;
    double ratio_max;
}
static EventCuts getCuts (Sample sample, int ecat){
    switch(sample){
        cas Sample::uranium:
    if (ecat == 0) return {-8.,15.6e3,39e3,0.3};
    if (ecat == 1) return { -7.0, 10., 8e3, 36e3, 0.3};
    if (ecat ==2) return {-10., 10., 8e3, 35e3, 0.3};
    return {-8.0,7.0,13e3, 33e3, 0.2};
    case Sample::gold:
        if (ecat==2) return {-3.0, 3.5, 18e3, 38e3, 0.3};
            if (ecat==1) return {-2.5, 3.0, 18e3, 34e3, 0.2};
            return {-2.0, 3.0, 20e3, 36e3, 0.2};
    }
}
static bool passAmplitudeCut(float a8, float a9, const EventCuts& c){
    double sum = a8+a9;
    double ratio = (a9-a8)/(a9+a8);
    if (sum<c.amp_min || sum>c.amp_max || ratio>ratio_max) return false;
    return true;
}