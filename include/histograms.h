#pragma once
#include "types.h" 
#include "config.h"
#include "cuts.h"
#include "constants.h"
#include "TTree.h"
#include "TH1D.h"
#include <vector>
#include <cmath>


// ---------------------- histograms.h --------------------------------
// 1. read the coincidences tree and apply cuts
// 2. store signals in roi or background vectors
// 3. fill dt histograms for each region to subtract background



static void fillHistograms(
    TTree* tree,
    const AnalysisConfig& cfg,
    std::vector<TH1D*>& hists_tof,
    Vec3D& counts_roi,
    Vec3D& counts_bkg,
    double emin,
    double emax )
    {
        int nbins = (int) cfg.energy_bins.size()-1;

        double cos_theta_det, cos_theta;
        double tof1,tof0, neutron_energy;
        float amp1, amp0;

        tree->SetBranchAddress("tof1", &tof1);
        tree->SetBranchAddress("tof0", &tof0);
        tree->SetBranchAddress("amp0", &amp0);
        tree->SetBranchAddress("amp1", &amp1);
        tree->SetBranchAddress("cos_theta",      &cos_theta);
        tree->SetBranchAddress("cos_theta_det",  &cos_theta_det);
        tree->SetBranchAddress("neutron_energy", &neutron_energy);

        Long64_t nentries = tree->GetEntries();
        
        for (Long64_t i = 0; i<nentries; i++){
            tree->GetEntry(i);
            if (neutron_energy<emin || neutron_energy>emax) continue;
            int e_bin = findBin (cfg.energy_bins, neutron_energy);
            if(e_bin<0 || e_bin>=nbins) continue;
            if (neutron_energy>1000) continue;
            EventCuts c = getCuts (cfg.sample, neutron_energy);
            if (!passAmplitudeCut(amp0, amp1, c)) continue;
            if (cos_theta_det<0.0) continue;
            if (cos_theta_det>1 || std::abs(cos_theta)>1) continue;

            double dt = tof1-tof0;
            hists_tof[e_bin]->Fill(dt);

            int j = int (std::abs(cos_theta)/dcos_beam);
            int ii = int (cos_theta_det/dcos_det);
            if (j>=nbins_beam) continue;
            if (ii>=nbins_det) continue;
            

            if (dt<c.roi_min || dt>c.roi_max)
                counts_bkg[e_bin][j][ii]++;
            else
                counts_roi[e_bin][j][ii]++;
        }
    }