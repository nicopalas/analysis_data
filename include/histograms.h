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

using DatasetMap = std::map<int, RooDataSet*>;

static void fillHistograms(
    TTree* tree,
    const AnalysisConfig& cfg,
    std::vector<TH1D*>& hists_tof,
    Vec3D& counts_roi,
    Vec3D& counts_bkg,
    Vec3D& counts_upeak,
    double emin,
    double emax)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    double cos_theta_det, cos_theta;
    double tof1, tof0, neutron_energy;
    float  amp0, amp1;
    double x0,x1,y0,y1;

    tree->SetBranchAddress("tof1",          &tof1);
    tree->SetBranchAddress("tof0",          &tof0);
    tree->SetBranchAddress("amp0",          &amp0);
    tree->SetBranchAddress("amp1",          &amp1);
    tree->SetBranchAddress("x1",           &x1);
    tree->SetBranchAddress("x0",           &x0);
    tree->SetBranchAddress("y0",           &y0);
    tree->SetBranchAddress("y1",           &y1);
    tree->SetBranchAddress("neutron_energy",&neutron_energy);

    Long64_t nentries = tree->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        double cos_theta_det = 5.0 /(sqrt((x0-x1-5.0)*(x0-x1-5.0) + (y1-y0)*(y1-y0)+25));
        double phi_det = TMath::ATan2(y1-y0,x0-x1-5.0);
        double sin_theta_det = std::sqrt(1-cos_theta_det*cos_theta_det);
        double cos_theta = 1/sqrt(2) * (-sin_theta_det*std::cos(phi_det)+cos_theta_det);
        if(neutron_energy < emin || neutron_energy > emax) continue;
        int e_bin = findBin(cfg.energy_bins, neutron_energy);
        if(e_bin < 0 || e_bin >= nbins) continue;
        if(neutron_energy > 1000) continue;
        EventCuts c = getCuts(cfg.sample, neutron_energy);
        if(!passAmplitudeCut(amp0, amp1, c)) continue;
        if(cos_theta_det < 0.0) continue;
        if(std::fabs(cos_theta_det) > 1 || std::fabs(cos_theta) > 1) continue;

        double dt = tof1 - tof0;
        hists_tof[e_bin]->Fill(dt);

        int j  = int(std::fabs(cos_theta) / dcos_beam);
        int ii = int(cos_theta_det / dcos_det);
        if(j  >= nbins_beam) continue;
        if(ii >= nbins_det)  continue;

        if(dt >= c.roi_min && dt <= c.roi_max){
            counts_roi[e_bin][j][ii]++;
        } else if(inUraniumPeak(dt, c)){
            counts_upeak[e_bin][j][ii]++;
        } else {
            counts_bkg[e_bin][j][ii]++;
        }
    }
}

static void fillRooDatasets(
    TTree*                tree,
    const AnalysisConfig& cfg,
    DatasetMap&           datasets,
    RooRealVar&           dt_var,
    RooRealVar&           amp_sum_var,  // era amp0_var+amp1_var por separado — unificado
    double                emin,
    double                emax)
{
    int nbins = (int)cfg.energy_bins.size() - 1;

    double cos_theta, cos_theta_det, neutron_energy, tof1, tof0;
    float  amp0, amp1;

    tree->SetBranchAddress("cos_theta",     &cos_theta);
    tree->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tree->SetBranchAddress("neutron_energy",&neutron_energy);
    tree->SetBranchAddress("tof1",          &tof1);
    tree->SetBranchAddress("tof0",          &tof0);
    tree->SetBranchAddress("amp0",          &amp0);
    tree->SetBranchAddress("amp1",          &amp1);

    // observable set — amp_sum_var ya declarado como parámetro
    RooArgSet obs(dt_var, amp_sum_var);

    // crea un dataset vacío por cada bin de energía
    for(int i = 0; i < nbins; i++){
        datasets[i] = new RooDataSet(
            Form("data_ebin%d", i),
            Form("Data energy bin %d", i),
            obs);
    }

    Long64_t nentries = tree->GetEntries();

    for(Long64_t i = 0; i < nentries; i++){
        tree->GetEntry(i);

        // cortes de calidad — SIN passAmplitudeCut
        if(std::abs(cos_theta) > 1.0 || cos_theta_det > 1.0
                                      || cos_theta_det < 0.0) continue;
        if(neutron_energy > emax || neutron_energy < emin)    continue;
        if(neutron_energy > 1500) continue;

        int ebin = findBin(cfg.energy_bins, neutron_energy);
        if(ebin < 0 || ebin >= nbins) continue;  // nombre consistente: ebin

        double dt      = tof1 - tof0;
        double amp_sum = amp0 + amp1;

        // corte mínimo de amplitud — ruido electrónico
        if(amp_sum < 10e3) continue;  // <- punto y coma que faltaba

        // verifica rango de los observables
        if(dt      < dt_var.getMin()      || dt      > dt_var.getMax())      continue;
        if(amp_sum < amp_sum_var.getMin() || amp_sum > amp_sum_var.getMax()) continue;

        dt_var.setVal(dt);
        amp_sum_var.setVal(amp_sum);
        datasets[ebin]->add(obs);  // minúscula — RooDataSet::add()
    }

    for(int i = 0; i < nbins; i++)
        std::cout << Form("ebin %2d: %6d eventos\n",
                          i, datasets[i]->numEntries());
}