#pragma once
#include <vector>
#include <string>

enum class Sample {uranium, gold};

struct AnalysisConfig {
    Sample sample;
    std::string tree_name;
    std::string input_file;
    std::string output_tag;

    std::vector<double> energy_bins;

    std::string efficiency_file;
    std::vector<double> energy_bins_eff;
    double cos_det_cut = 0.5;
    int    bins_beam  = 40;
    int    bins_det   = 40;
    int    n_toys      = 500;   
    double atoms = 0.0;
    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root";
    std::string flux_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/Energy_EAR1_FLUKA_500bdp.root";
    std::string flux_hist = "h_Flux";

};

static AnalysisConfig makeUraniumConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::uranium;
    c.tree_name       = "events_uranium";
    c.input_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "output_efficiency_uranium.root";
    c.energy_bins_eff = {1, 10, 100, 500, 1000};
    c.atoms = 5.04e17;
    return c;
}

static AnalysisConfig makeGoldConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::gold;
    c.tree_name       = "events_gold";
    c.input_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "output_efficiency_gold.root";
    c.energy_bins_eff = {100, 300, 600, 1000};
    c.atoms = 9.17e17;
    return c;
}