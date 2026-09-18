#pragma once
#include "types.h"
#include "config.h"
#include "cuts.h"
#include "constants.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "TTree.h"
#include "TH1D.h"
#include <map>
#include <vector>
#include <cmath>


TFile *cut = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cut_gold_angle.root");
TCutG *cut0 = (TCutG*) cut->Get("cut0");

TFile *cut_u = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cut_uranium_angle.root");
TCutG *cut1 = (TCutG*) cut_u->Get("cut_uranium");
std::vector<int> bad_Runs = {118589, 118618, 118619, 118620, 118622, 118623, 118624, 118625, 118627, 118685, 118694, 118701, 118705, 118706, 118710, 118711, 118712, 118713, 118719, 118726, 118727, 118728, 118729, 118730, 118733, 118734, 118736, 118739, 118740, 118748, 118769, 118773, 118775, 118776, 118777, 118778, 118779, 118780, 118781, 118782, 118784, 118785, 118786, 118787, 118788, 118790};
using DatasetMap = std::map<int, RooDataSet*>;
static void compute_angles(
    double x0, double y0, double x1, double y1,
    double offsetx0, double offsetx1,
    double offsety0, double offsety1,
    double& cos_theta_det,
    double& phi_det,
    double& phi,
    double& cos_theta)
{
    double dx = (2*x1 +2.5) - (2*x0 - 2.5);
    double dy =  2*y1 - 2*y0 - (2*offsety1-2*offsety0);
    double dz =  5.0;

    double nd = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (nd > 0.0) {
        double sth = std::sqrt(dx*dx + dy*dy) / nd;
        cos_theta_det = dz / nd;
        phi_det       = std::atan2(dy, dx);

        double nx = (-sth * std::cos(phi_det) + cos_theta_det) / std::sqrt(2.0);
        double ny =   sth * std::sin(phi_det);
        double nz = ( sth * std::cos(phi_det) + cos_theta_det) / std::sqrt(2.0);
        double nb = std::sqrt(nx*nx + ny*ny + nz*nz);

        cos_theta = nz/nb;
        phi = std::atan2(ny,nx);
    } else {
        cos_theta_det = -999.;
        phi_det       = -999.;
        cos_theta     = -999.;
        phi = -999.;
    }
}

static void fillHistograms(
    TTree* tree,
    const AnalysisConfig& cfg,
    std::vector<TH1D*>& hists_tof,
    Vec3D& counts_roi,
    Vec3D& counts_bkg,
    Vec3D& counts_upeak,
    double emin,
    double emax,
    double& out_mean_x0, double& out_mean_x1,
    double& out_mean_y0, double& out_mean_y1)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    double cos_theta_det, cos_theta;
    double tof1, tof0, neutron_energy;
    double  amp0, amp1;
    double x0,x1,y0,y1;
    double phi, phi_det;
    int RunNumber;
    float PulseIntensity;
    double sumX0, sumY1, sumX1, sumY0;
    int full_position;
    tree->SetBranchAddress("tof1",          &tof1);
    tree->SetBranchAddress("tof0",          &tof0);
    tree->SetBranchAddress("amp0",          &amp0);
    tree->SetBranchAddress("amp1",          &amp1);
    tree->SetBranchAddress("x1",           &x1);
    tree->SetBranchAddress("x0",           &x0);
    tree->SetBranchAddress("y0",           &y0);
    tree->SetBranchAddress("y1",           &y1);
    tree->SetBranchAddress("sumX0", &sumX0);
    tree->SetBranchAddress("sumX1", &sumX1);
    tree->SetBranchAddress("sumY1", &sumY1);
    tree->SetBranchAddress("sumY0", &sumY0);
    tree->SetBranchAddress("neutron_energy",&neutron_energy);
    tree->SetBranchAddress("cos_theta", &cos_theta);
    tree->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("phi_det", &phi_det);
    tree->SetBranchAddress("RunNumber", &RunNumber);
    tree->SetBranchAddress("PulseIntensity", &PulseIntensity);
    tree->SetBranchAddress("full_position", &full_position);

    Long64_t nentries = tree->GetEntries();
    Long64_t n_mean = 0;
    double mean_x0 = 0.0, mean_x1 = 0.0, mean_y0 = 0.0, mean_y1 = 0.0;
    for (Long64_t i = 0; i < nentries; ++i){
        tree->GetEntry(i);
        if(neutron_energy < emin || neutron_energy > emax) continue;
        if (RunNumber == 118771 || RunNumber == 118789 || !full_position) continue;
        double ratio = (amp1 - amp0) / (amp0 + amp1);
        double dt    = tof1 - tof0; 
        EventCuts c = getCuts(cfg.sample, neutron_energy);
        if (!passAmplitudeCut(amp0, amp1, c)) continue;
        mean_x0 += x0;
        mean_x1 += x1;
        mean_y0 += y0;
        mean_y1 += y1;
        ++n_mean;
    }
    mean_x0 /= n_mean;
    mean_x1 /= n_mean;
    mean_y0 /= n_mean;
    mean_y1 /= n_mean;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (RunNumber == 118771 || RunNumber == 118789) continue;
        if(neutron_energy < emin || neutron_energy > emax) continue;
        int e_bin = findBin(cfg.energy_bins, neutron_energy);
        if(e_bin < 0 || e_bin >= nbins) continue;
        if(neutron_energy > 1000) continue;
        double ratio = (amp1 - amp0) / (amp0 + amp1);
        double dt    = tof1 - tof0;   // (or compute later, but same thing)
        EventCuts c = getCuts(cfg.sample, neutron_energy);
        if (!passAmplitudeCut(amp0, amp1, c)) continue;
        if (!full_position) continue;
        double cos_theta_corrected, cos_theta_det_corrected, phi_det, phi;
        compute_angles(x0, y0, x1, y1,
                       mean_x0, mean_x1, mean_y0, mean_y1,
                       cos_theta_det_corrected, phi_det, phi, cos_theta_corrected);
        if(cos_theta_det_corrected < 0.0) continue;
        if(std::fabs(cos_theta_det_corrected) > 1 || std::fabs(cos_theta_corrected) > 1) continue;
        double ampsum = amp0+amp1;
        hists_tof[e_bin]->Fill(dt);

        int j  = int(std::fabs(cos_theta_corrected) / dcos_beam);
        int ii = int(cos_theta_det_corrected / dcos_det);
        if(j  >= nbins_beam) continue;
        if(ii >= nbins_det)  continue;

        // In fillHistograms:
        if(dt >= c.roi_min && dt <= c.roi_max){
            counts_roi[e_bin][j][ii]++;
            } 
        else if(inUraniumPeak(dt, c) && cos_theta_det>0.6){
            counts_upeak[e_bin][j][ii]++;
        }
        else{
            counts_bkg[e_bin][j][ii]++;
        }
    }
}

