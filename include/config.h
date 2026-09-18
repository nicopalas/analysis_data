#pragma once
#include <vector>
#include <string>

enum class Sample {uranium, gold, uranium_mc};

struct AnalysisConfig {
    Sample sample;
    std::string tree_name;
    std::string input_file;
    std::string output_tag;
    std::string acc_file;        
    std::string acc_tree_name;

    std::vector<double> energy_bins;

    std::string efficiency_file;
    std::vector<double> energy_bins_eff;
    double cos_det_cut = 0.5;
    int    bins_beam  = 20;
    int    bins_det   = 20;
    int    n_toys      = 1000;   
    double atoms = 0.0;
    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root";
    std::string flux_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/evalFlux_prelim.root";
    std::string flux_hist = "hEval_Abs";

};

static AnalysisConfig makeUraniumConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::uranium;
    c.tree_name       = "events_uranium";
    c.input_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences_final.root";
    c.output_tag      = tag;
    c.acc_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acc_selection.root";
    c.acc_tree_name = "events_uranium";

    c.energy_bins     = energy_bins;
    c.efficiency_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/output_efficiency_uranium.root";
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
    c.input_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences_final.root";
    c.acc_file      = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acc_selection.root";
    c.acc_tree_name = "events_gold";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/Au-197/output_efficiency_gold.root";
    c.energy_bins_eff = {50, 500, 1000};
    c.atoms = 9.17e17;
    return c;
}

static AnalysisConfig makeUraniumMCConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::uranium_mc;
    c.tree_name       = "CoincTree";
    c.input_file      = "/Users/nico/Desktop/Tese/Analysis/montecarlo/data/output_mc_final.root";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "/Users/nico/Desktop/Tese/Analysis/montecarlo/output/U-238/output_efficiency_uranium.root";
    c.energy_bins_eff = {10, 100, 500, 1000, 2000};
    c.atoms = 5.04e17;
    return c;
}