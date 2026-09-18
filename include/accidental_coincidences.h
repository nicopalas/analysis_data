// accidental_shape.h
#pragma once
#include "types.h"
#include "config.h"
#include "cuts.h"
#include "constants.h"
#include "TTree.h"
#include <cmath>

// Fill the angular shape of accidental coincidences using the SAME dt window
// and the SAME cuts as the ROI selection. Offsets must be the ones derived
// from the physics file -- they encode detector alignment, not sample stats.
static void fillAccidentalShape(
    TTree* tree_acc,
    const AnalysisConfig& cfg,
    Vec3D& counts_acc,
    double mean_x0, double mean_x1,
    double mean_y0, double mean_y1,
    double emin, double emax)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    double tof0, tof1, amp0, amp1, x0, x1, y0, y1, neutron_energy;
    int RunNumber;
    double cos_theta, cos_theta_det, phi, phi_det;
    double sumX0, sumX1, sumY0, sumY1;
    tree_acc->SetBranchAddress("tof0",  &tof0);
    tree_acc->SetBranchAddress("tof1",  &tof1);
    tree_acc->SetBranchAddress("amp0",  &amp0);
    tree_acc->SetBranchAddress("amp1",  &amp1);
    tree_acc->SetBranchAddress("cos_theta", &cos_theta);
    tree_acc->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tree_acc->SetBranchAddress("phi", &phi);
    tree_acc->SetBranchAddress("phi_det", &phi_det);
    tree_acc->SetBranchAddress("neutron_energy", &neutron_energy);
    tree_acc->SetBranchAddress("RunNumber", &RunNumber);
    tree_acc->SetBranchAddress("sumX0", &sumX0);
    tree_acc->SetBranchAddress("sumX1", &sumX1);
    tree_acc->SetBranchAddress("sumY1", &sumY1);
    tree_acc->SetBranchAddress("sumY0", &sumY0);

    Long64_t nentries = tree_acc->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree_acc->GetEntry(i);
        if (RunNumber == 118771 || RunNumber == 118789 || RunNumber == 118668 || RunNumber == 118587 || RunNumber == 118751 ||
    RunNumber == 118558 || RunNumber == 118565 || RunNumber == 118569 || RunNumber == 118570 || RunNumber == 118581 ||
    RunNumber == 118588 || RunNumber == 118601 || RunNumber == 118602 || RunNumber == 118603 || RunNumber == 118605 ||
    RunNumber == 118606 || RunNumber == 118607 || RunNumber == 118608 || RunNumber == 118609 || RunNumber == 118610 ||
    RunNumber == 118612 || RunNumber == 118614 || RunNumber == 118630 || RunNumber == 118645 || RunNumber == 118646 ||
    RunNumber == 118647 || RunNumber == 118649 || RunNumber == 118650 || RunNumber == 118651 || RunNumber == 118652 ||
    RunNumber == 118653 || RunNumber == 118654 || RunNumber == 118681 || RunNumber == 118683 || RunNumber == 118684 ||
    RunNumber == 118691 || RunNumber == 118695 || RunNumber == 118696 || RunNumber == 118743 || RunNumber == 118744 ||
    RunNumber == 118791 || RunNumber == 118792 || RunNumber == 118793 || RunNumber == 118794 || RunNumber == 118795) continue;
        if (neutron_energy < emin || neutron_energy > emax) continue;
        if (neutron_energy > 2000) continue;

        int e_bin = findBin(cfg.energy_bins, neutron_energy);
        if (e_bin < 0 || e_bin >= nbins) continue;
        if (cfg.sample == Sample::gold){
        if (sumY1<103 || sumX0<103 || sumX1>107 || sumY0>102) continue;
        }

        EventCuts c = getCuts(cfg.sample, neutron_energy);
        if (!passAmplitudeCut(amp0, amp1, c)) continue;   // same as ROI

        double dt = tof1 - tof0;
        if (dt < c.roi_min || dt > c.roi_max) continue;   // THE point

        if (cos_theta_det < 0.0) continue;
        if (std::fabs(cos_theta_det) > 1.0 || std::fabs(cos_theta) > 1.0) continue;

        int j  = int(std::fabs(cos_theta) / dcos_beam);
        int ii = int(cos_theta_det / dcos_det);
        if (j >= nbins_beam || ii >= nbins_det) continue;

        counts_acc[e_bin][j][ii] += 1.0;
    }
}