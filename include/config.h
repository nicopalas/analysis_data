#pragma once
#include <vector>
#include <string>

enum class Sample {uranium, gold, cerium};

struct AnalysisConfig {
    Sample sample;
    std::string tree_name;
    std::string input_file;
    std::string output_tag;

    std::vector<double> energy_bins;

    std::string efficiency_file;
    std::vector<double> energy_bins_eff;
};

static AnalysisConfig makeUraniumConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::Uranium;
    c.tree_name       = "events_uranium";
    c.input_file      = "/Users/nico/Desktop/Tese/Macros/Macros/n_tof_cerium/ROOT_files/cathode_candidates/coincidences.root";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "output_efficiency_uranium.root";
    c.energy_bins_eff = {1, 10, 100, 500, 1000};
    return c;
}

static AnalysisConfig makeGoldConfig(
    const std::vector<double>& energy_bins,
    const std::string& tag = "nominal")
{
    AnalysisConfig c;
    c.sample          = Sample::Gold;
    c.tree_name       = "events_gold";
    c.input_file      = "/Users/nico/Desktop/Tese/Macros/Macros/n_tof_cerium/ROOT_files/cathode_candidates/coincidences.root";
    c.output_tag      = tag;
    c.energy_bins     = energy_bins;
    c.efficiency_file = "output_efficiency_gold.root";
    c.energy_bins_eff = {120, 300, 600, 1000};
    return c;
}