#pragma once
#include <vector>
#include <algorithm>
#include <cmath>


// ---------------------utils.h----------------------
// functions that are used frequently
static int findBin( const std::vector<double>& edges, double value){
    for (int i = 0; i<edges.size(); i++){
        if (value>edges[i] && value<edges[i+1]){
            return i;
        }
        if (edges.size()==0){
            std::cout << "WARNING: vector is empty. No index found." << std::endl;
            return -1;
        }
    }
    return -1;
}
static double mean(const std::vector<double>& v){
    if (v.size()==0) return 0;
    double mean = 0.0;
    for (int i = 0; i<v.size(); i++){
        mean+=v[i];
    }
    return mean/v.size();
}

static double computeStd(const std::vector<double>& v){
    if (v.size()) return 0;
    double mean = 0.;
    double mean2 = 0.;
    for (double x : v){
        mean+=x;
        mean2+=x*x;
    }
    double var = mean2-mean*mean;
    return std::sqrt(var);
}

static std::vector<double> buildLogBins(int nbins, double e_min, double e_max){
    std::vector<double> bins;
    double log_min = std::log10(e_min);
    double log_max = std::log10(e_max);
    for(int i = 0; i <= nbins; ++i){
        double log_e = log_min + (log_max - log_min) * std::pow(i / double(nbins), 0.7);
        bins.push_back(std::pow(10.0, log_e));
    }
    return bins;
}