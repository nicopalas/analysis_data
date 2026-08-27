// ===========================================================================
//  systematics_cuts.C
//
//  Event-selection systematics for U-238(n,f) and Au-197(n,f).
//
//  Method
//  ------
//  For every cut variation the ENTIRE chain is re-run:
//      fill -> background fit -> signal -> efficiency -> anisotropy -> sigma
//  so the (large) part of the selection effect that cancels between the
//  numerator and the efficiency actually cancels. Varying only the signal pass
//  while keeping the nominal efficiency overestimates the systematic badly.
//
//  Each shift is tested against Barlow's criterion for correlated samples,
//  sigma(Delta)^2 = |sigma_var^2 - sigma_nom^2|, so statistical noise from the
//  events added/removed is not promoted to a systematic.
//
//  Shifts are accumulated into a covariance matrix (outer product of the shift
//  vectors), because a threshold change moves whole energy regions coherently
//  and must not be added bin-by-bin in quadrature.
//
//  Variations studied
//  ------------------
//    U-238 : amp_min, amp_max, ratio_max, roi_width, roi_shift,
//            e_threshold, cut_vs_E_smooth
//    Au-197: upeak_width, upeak_shift, amp_min, amp_max, ratio_max
//
//  Run
//  ---
//    root -l -b -q 'systematics_cuts.C+("U238","all")'
//    root -l -b -q 'systematics_cuts.C+("Au197","all")'
//    root -l -b -q 'systematics_cuts.C+("Au197","scan_upeak")'
//    root -l -b -q 'systematics_cuts.C+("U238","scan_ampmin")'
//    root -l -b -q 'systematics_cuts.C+("both","all")'
// ===========================================================================

#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/cut_variations.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/histograms.h"
#include "../include/background_subtraction.h"
#include "../include/signal_selection.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"
#include "../include/cross_section.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ===========================================================================
//  SAMPLE CONFIGURATION
// ===========================================================================
typedef AnalysisConfig (*ConfigMaker)(const std::vector<double>&,
                                      const std::string&);

struct SampleCfg {
    std::string  name;              // "U238" / "Au197"
    std::string  outdir;
    std::string  acceptance_file;
    ConfigMaker  makeConfig = nullptr;

    // event-energy window handed to fillHistograms
    double emin = 1.0, emax = 1000.0;

    // coarse bins for the efficiency pass
    std::vector<double> ebins_eff;

    // fine logarithmic bins for anisotropy / cross section
    int    nbins_aniso = 50;
    double aniso_lo    = 1.5;
    double aniso_hi    = 1000.0;

    // TOF histogram ranges
    double tof_lo_eff = -20, tof_hi_eff = 20;
    double tof_lo_ani = -30, tof_hi_ani = 30;

    // normalisation window for normalised_xs()
    double norm_lo = 8.0, norm_hi = 10.0, norm_ref = 2.006;

    // representative fine bins used by the plateau scans
    int probe[4] = {0,1,2,3};
};

static SampleCfg uraniumCfg(){
    SampleCfg c;
    c.name            = "U238";
    c.outdir          = "/Users/nico/Desktop/Tese/Analysis/cross_section/"
                        "output/U-238/systematics/";
    c.acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/"
                        "data/acceptance_coincidence.csv";
    c.makeConfig      = makeUraniumConfig;

    c.emin = 1.0;   c.emax = 1000.0;
    c.ebins_eff  = {1, 10, 100, 300, 600, 750, 1000};
    c.nbins_aniso = 50;
    c.aniso_lo    = 1.5;
    c.aniso_hi    = 1000.0;

    c.tof_lo_eff = -20; c.tof_hi_eff = 20;
    c.tof_lo_ani = -30; c.tof_hi_ani = 30;

    c.norm_lo = 8.0; c.norm_hi = 10.0; c.norm_ref = 2.006;

    c.probe[0] = 5; c.probe[1] = 18; c.probe[2] = 30; c.probe[3] = 42;
    return c;
}

static SampleCfg goldCfg(){
    SampleCfg c;
    c.name            = "Au197";
    c.outdir          = "/Users/nico/Desktop/Tese/Analysis/cross_section/"
                        "output/Au-197/systematics/";
    c.acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/"
                        "data/acceptance_coincidence.csv";
    c.makeConfig      = makeGoldConfig;

    // Au-197 fission has a high threshold: do not try to analyse below ~30 MeV.
    c.emin = 40.0;  c.emax = 1000.0;
    c.ebins_eff   = {40, 100, 300, 600, 1000};

    c.nbins_aniso = 5;          // fewer bins: gold has much lower statistics
    c.aniso_lo    = 40.0;
    c.aniso_hi    = 1000.0;

    // gold ROI is much narrower than uranium's; the uranium-contamination
    // window sits near -15..-2.5 ns, so the histogram must cover it.
    c.tof_lo_eff = -25; c.tof_hi_eff = 15;
    c.tof_lo_ani = -25; c.tof_hi_ani = 15;

    // adjust to the window where you have a reference Au cross section
    c.norm_lo = 80.0; c.norm_hi = 100.0; c.norm_ref = 1.0;

    c.probe[0] = 1; c.probe[1] = 2; c.probe[2] = 3; c.probe[3] = 4;
    return c;
}

// global knobs shared by both samples
namespace syscfg {
    const double barlow_z      = 2.0;   // significance threshold
    const bool   floor_on_fail = true;  // quote sigma(Delta) when Barlow fails
}

// ===========================================================================
//  RESULT CONTAINER
// ===========================================================================
struct ChainResult {
    bool        ok = false;
    std::string name;
    int nE = 0, nBeam = 0, nDet = 0;

    std::vector<double> sigma,      u_sigma;        // raw absolute XS
    std::vector<double> sigma_norm, u_sigma_norm;   // normalised XS
    std::vector<double> wfb,        u_wfb;          // W(0)/W(90)
    Vec2D               w,          u_w;            // full W(theta)
    std::vector<double> ecen;                       // bin centres
};

// ---------------------------------------------------------------------------
// acceptance does not depend on the cuts -> load once per path
// ---------------------------------------------------------------------------
static bool getAcceptance(const std::string& file, Vec2D& acceptance){
    static std::string cached_path;
    static Vec2D       cached;
    static bool        loaded = false;

    if(!loaded || cached_path != file){
        Vec2D fine;
        if(!loadAcceptanceCSV(file, fine)){
            std::cerr << "[syst] failed to load acceptance CSV: " << file << "\n";
            return false;
        }
        cached      = rebin(fine);
        cached_path = file;
        loaded      = true;
    }
    acceptance = cached;
    return true;
}

// ===========================================================================
//  FULL ANALYSIS CHAIN UNDER ONE CUT VARIATION
// ===========================================================================
static ChainResult runChain(const SampleCfg& S, const CutVariation& var, int uid){

    cutVar() = var;                  // <-- activates the variation globally

    ChainResult R;
    R.name  = var.name;
    R.nBeam = nbins_beam;
    R.nDet  = nbins_det;

    Vec2D acceptance;
    if(!getAcceptance(S.acceptance_file, acceptance)) return R;

    // ------------------- efficiency pass (coarse bins) ---------------------
    const std::vector<double>& eb = S.ebins_eff;
    const int nE_eff = (int)eb.size() - 1;
    if(nE_eff < 1){ std::cerr << "[syst] bad ebins_eff\n"; return R; }

    AnalysisConfig cfg_eff = S.makeConfig(eb, Form("eff_v%d", uid));

    TFile* fin = TFile::Open(cfg_eff.input_file.c_str());
    if(!fin || fin->IsZombie()){
        std::cerr << "[syst] cannot open " << cfg_eff.input_file << "\n"; return R;
    }
    TTree* tree = (TTree*)fin->Get(cfg_eff.tree_name.c_str());
    if(!tree){
        std::cerr << "[syst] tree " << cfg_eff.tree_name << " not found\n";
        fin->Close(); return R;
    }

    std::vector<TH1D*> h_eff(nE_eff, nullptr);
    for(int i = 0; i < nE_eff; ++i){
        h_eff[i] = new TH1D(Form("htof_eff_%s_v%d_%d", S.name.c_str(), uid, i),
                            "", 200, S.tof_lo_eff, S.tof_hi_eff);
        h_eff[i]->SetDirectory(0);
    }

    Vec3D roi_e(nE_eff, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D bkg_e(nE_eff, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D up_e (nE_eff, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree, cfg_eff, h_eff, roi_e, bkg_e, up_e, S.emin, S.emax);
    fin->Close();

    std::vector<BackgroundFit> bf_e(nE_eff);
    std::vector<double> cs_e(nE_eff, 0.0),     u_cs_e(nE_eff, 0.0);
    std::vector<double> up_cs_e(nE_eff, 0.0),  u_up_cs_e(nE_eff, 0.0);

    for(int i = 0; i < nE_eff; ++i){
        double Ec   = std::sqrt(eb[i] * eb[i+1]);
        EventCuts c = getCuts(cfg_eff.sample, Ec);
        bf_e[i] = fitBackground(cfg_eff, h_eff[i], c.roi_min, c.roi_max, i);
        if(!bf_e[i].func || !bf_e[i].hist_subtracted){
            std::cerr << "[syst] " << S.name << " / " << var.name
                      << " : efficiency fit " << i << " FAILED\n";
            for(auto* h : h_eff) delete h;
            return R;
        }
        cs_e[i]   = bf_e[i].counts_subtract_bkg;
        u_cs_e[i] = bf_e[i].u_counts_subtract_bkg;
    }

    Vec3D sig_e  (nE_eff, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_sig_e(nE_eff, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    computeSignal(cfg_eff, roi_e, bkg_e, up_e,
                  cs_e, u_cs_e, up_cs_e, u_up_cs_e, sig_e, u_sig_e);

    std::vector<EfficiencyResult> eff(nE_eff);
    for(int e = 0; e < nE_eff; ++e)
        eff[e] = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                   sig_e, u_sig_e, acceptance, e);

    for(auto* h : h_eff) delete h;
    for(auto& b : bf_e){ delete b.func; delete b.hist_subtracted; }

    // ------------- anisotropy / cross-section pass (fine bins) -------------
    const int nE = S.nbins_aniso;
    std::vector<double> ea = buildLogBins(nE, S.aniso_lo, S.aniso_hi);
    AnalysisConfig cfg = S.makeConfig(ea, Form("aniso_v%d", uid));

    TFile* fin2 = TFile::Open(cfg.input_file.c_str());
    if(!fin2 || fin2->IsZombie()) return R;
    TTree* tree2 = (TTree*)fin2->Get(cfg.tree_name.c_str());
    if(!tree2){ fin2->Close(); return R; }

    std::vector<TH1D*> h_a(nE, nullptr);
    for(int i = 0; i < nE; ++i){
        h_a[i] = new TH1D(Form("htof_ani_%s_v%d_%d", S.name.c_str(), uid, i),
                          "", 200, S.tof_lo_ani, S.tof_hi_ani);
        h_a[i]->SetDirectory(0);
    }

    Vec3D roi_a(nE, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D bkg_a(nE, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D up_a (nE, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    fillHistograms(tree2, cfg, h_a, roi_a, bkg_a, up_a, S.emin, S.emax);
    fin2->Close();

    std::vector<BackgroundFit> bf_a(nE);
    std::vector<double> cs_a(nE, 0.0),    u_cs_a(nE, 0.0);
    std::vector<double> up_cs_a(nE, 0.0), u_up_cs_a(nE, 0.0);

    for(int i = 0; i < nE; ++i){
        double Ec   = std::sqrt(ea[i] * ea[i+1]);
        EventCuts c = getCuts(cfg.sample, Ec);
        bf_a[i] = fitBackground(cfg, h_a[i], c.roi_min, c.roi_max, i);
        cs_a[i]   = bf_a[i].func ? bf_a[i].counts_subtract_bkg   : 0.0;
        u_cs_a[i] = bf_a[i].func ? bf_a[i].u_counts_subtract_bkg : 0.0;
    }

    Vec3D sig_a  (nE, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_sig_a(nE, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    computeSignal(cfg, roi_a, bkg_a, up_a,
                  cs_a, u_cs_a, up_cs_a, u_up_cs_a, sig_a, u_sig_a);

    R.nE = nE;
    R.w  .assign(nE, std::vector<double>(nbins_beam, 0.0));
    R.u_w.assign(nE, std::vector<double>(nbins_beam, 0.0));
    R.wfb  .assign(nE, 0.0);  R.u_wfb  .assign(nE, 0.0);
    R.sigma.assign(nE, 0.0);  R.u_sigma.assign(nE, 0.0);
    R.ecen .assign(nE, 0.0);

    std::vector<double> Elo(nE), Ehi(nE);
    for(int e = 0; e < nE; ++e){ Elo[e] = ea[e]; Ehi[e] = ea[e+1]; }

    for(int e = 0; e < nE; ++e){
        double Ec = std::sqrt(ea[e] * ea[e+1]);
        R.ecen[e] = Ec;

        int ie = findBin(S.ebins_eff, Ec);
        if(ie < 0)       ie = 0;
        if(ie >= nE_eff) ie = nE_eff - 1;

        AnisotropyResult A = anisotropy(nbins_beam, nbins_det, sig_a, u_sig_a,
                                        acceptance, e,
                                        eff[ie].eps, eff[ie].u_eps, cfg);
        for(int j = 0; j < nbins_beam; ++j){
            R.w[e][j]   = A.w[j];
            R.u_w[e][j] = A.u_w[j];
        }
        R.wfb  [e] = A.w  [nbins_beam - 1];
        R.u_wfb[e] = A.u_w[nbins_beam - 1];

        CrossSection X = cross_section_absolute(
            cfg.flux_file, cfg.flux_hist,
            nbins_beam, nbins_det,
            sig_a, u_sig_a,
            eff[ie].eps, eff[ie].u_eps,
            acceptance, e,
            Elo, Ehi, cfg.atoms, cfg);

        R.sigma  [e] = X.sigma;
        R.u_sigma[e] = X.u_sigma;
    }
    R.sigma_norm   = R.sigma;      // use the absolute cross section
    R.u_sigma_norm = R.u_sigma;

    
    for(auto* h : h_a) delete h;
    for(auto& b : bf_a){ delete b.func; delete b.hist_subtracted; }

    R.ok = true;
    return R;
}

// ===========================================================================
//  BARLOW TEST FOR CORRELATED (NESTED) SAMPLES
// ===========================================================================
struct BarlowPoint {
    double delta = 0.0;   // var - nom
    double sig_d = 0.0;   // sqrt(|u_var^2 - u_nom^2|)
    double z     = 0.0;   // delta / sig_d
    double syst  = 0.0;   // adopted systematic
    bool   sign  = false; // passed significance
};

static BarlowPoint barlow(double nom, double u_nom, double var, double u_var){
    BarlowPoint b;
    b.delta = var - nom;
    b.sig_d = std::sqrt(std::fabs(u_var * u_var - u_nom * u_nom));
    b.z     = (b.sig_d > 0.0) ? b.delta / b.sig_d : 0.0;
    b.sign  = std::fabs(b.z) > syscfg::barlow_z;
    if(b.sign)                     b.syst = std::fabs(b.delta);
    else if(syscfg::floor_on_fail) b.syst = b.sig_d;
    else                           b.syst = 0.0;
    return b;
}

// ===========================================================================
//  SYSTEMATIC SOURCE BOOKKEEPING
// ===========================================================================
struct SystSource {
    std::string name;
    std::vector<double> shift_sigma;
    std::vector<double> shift_wfb;
    std::vector<double> z_sigma;
    std::vector<double> z_wfb;
    bool correlated_in_E     = true;
    bool correlated_in_theta = true;
};

static std::vector<double> symmetrise(const std::vector<double>& nom,
                                      const std::vector<double>& up,
                                      const std::vector<double>& dn,
                                      bool has_dn){
    std::vector<double> d(nom.size(), 0.0);
    for(size_t i = 0; i < nom.size(); ++i){
        if(has_dn && i < dn.size() && i < up.size()) d[i] = 0.5 * (up[i] - dn[i]);
        else if(i < up.size())                       d[i] = up[i] - nom[i];
    }
    return d;
}

static TH2D* buildCovariance(const std::vector<SystSource>& src,
                             const std::vector<double>& ecen,
                             bool use_sigma,
                             const char* hname){
    const int n = (int)ecen.size();
    TH2D* cov = new TH2D(hname, ";E_{n} bin;E_{n} bin", n, 0, n, n, 0, n);
    cov->SetDirectory(0);
    for(const auto& s : src){
        const std::vector<double>& d = use_sigma ? s.shift_sigma : s.shift_wfb;
        if((int)d.size() < n) continue;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j){
                double v = s.correlated_in_E ? d[i] * d[j]
                                             : (i == j ? d[i] * d[i] : 0.0);
                cov->SetBinContent(i+1, j+1, cov->GetBinContent(i+1, j+1) + v);
            }
    }
    return cov;
}

// ===========================================================================
//  VARIATION LISTS
// ===========================================================================
struct VarPair {
    std::string  name;
    CutVariation up, dn;
    bool has_dn    = true;
    bool corrE     = true;
    bool corrTheta = true;
};

// ---- U-238: full selection study --------------------------------------------
static std::vector<VarPair> variationsUranium(){
    std::vector<VarPair> v;

    // amplitude threshold. Fragment energy loss in the deposit scales ~1/cos(t),
    // so this bias is angle dependent and does NOT cancel in W(theta)/W(90).
    v.push_back({"amp_min",   varAmpMin(1.05), varAmpMin(0.95), true, true, false});
    v.push_back({"amp_max",   varAmpMax(1.05), varAmpMax(0.95), true, true, true});

    // pulse-height asymmetry. The one-sided form is intentional (unequal
    // detector gains), so only the VALUE is varied, never the sign convention.
    v.push_back({"ratio_max", varRatio(+0.05), varRatio(-0.05), true, true, false});

    // ROI. Correlated with the background fit range through fitBackground():
    // treat width+shift as ONE source when quoting, not two independent ones.
    v.push_back({"roi_width", varRoiWidth(1.10), varRoiWidth(0.90), true, true, true});
    v.push_back({"roi_shift", varRoiShift(+0.5), varRoiShift(-0.5), true, true, true});

    // artefacts of the step function in getCuts vs E
    v.push_back({"e_threshold", varEThreshold(1.20), varEThreshold(0.83),
                 true, true, true});
    {   VarPair p; p.name = "cut_vs_E_smooth"; p.up = varSmoothE();
        p.has_dn = false; p.corrE = false; v.push_back(p); }

    return v;
}

// ---- Au-197: uranium-peak subtraction + amplitude cuts ----------------------
static std::vector<VarPair> variationsGold(){
    std::vector<VarPair> v;

    // ---- uranium contamination peak -------------------------------------
    // The Au deposit carries a U impurity whose fissions arrive in a separate
    // TOF window (nominally about -15 .. -2.5 ns). Two independent things can
    // go wrong: the window is too wide/narrow, or it is mispositioned.
    //
    // width : widen/narrow both edges symmetrically about the window centre
    v.push_back({"upeak_width", varUPeak(-1.5, +1.5), varUPeak(+1.5, -1.5),
                 true, true, true});
    // shift : slide the whole window without changing its width
    v.push_back({"upeak_shift", varUPeak(+1.0, +1.0), varUPeak(-1.0, -1.0),
                 true, true, true});
    // a large widening probes how much true Au signal the window eats
    {   VarPair p; p.name = "upeak_wide"; p.up = varUPeak(-3.0, +3.0);
        p.has_dn = false; v.push_back(p); }

    // ---- amplitude cuts ---------------------------------------------------
    v.push_back({"amp_min",   varAmpMin(1.05), varAmpMin(0.95), true, true, false});
    v.push_back({"amp_max",   varAmpMax(1.05), varAmpMax(0.95), true, true, true});

    // ---- amplitude ratio --------------------------------------------------
    // Au uses a much tighter ratio_max (0.10-0.15 at low E) than U (0.30), so
    // a smaller step is appropriate here.
    v.push_back({"ratio_max", varRatio(+0.03), varRatio(-0.03), true, true, false});

    return v;
}

static std::vector<VarPair> buildVariations(const SampleCfg& S){
    return (S.name == "Au197") ? variationsGold() : variationsUranium();
}

// ===========================================================================
//  PLATEAU SCANS
// ===========================================================================
static void plateauScan(const SampleCfg& S,
                        const ChainResult& nom,
                        const std::vector<double>& xvals,
                        const std::vector<CutVariation>& vars,
                        const char* label,
                        const char* xtitle){

    std::vector<ChainResult> R;
    for(size_t k = 0; k < vars.size(); ++k){
        std::cout << "\n=== SCAN " << label << "  point " << k
                  << " (" << xvals[k] << ") ===\n";
        R.push_back(runChain(S, vars[k], 100 + (int)k));
    }

    TFile* fo = TFile::Open(
        (S.outdir + "scan_" + label + "_" + S.name + ".root").c_str(), "RECREATE");

    TCanvas* c = new TCanvas(Form("c_scan_%s", label), label, 1100, 800);
    c->Divide(2, 2);

    for(int p = 0; p < 4; ++p){
        int e = S.probe[p];
        if(e >= nom.nE) continue;

        std::vector<double> x, ysig, eysig, yw, eyw;
        for(size_t k = 0; k < R.size(); ++k){
            if(!R[k].ok) continue;
            x.push_back(xvals[k]);

            double dn_w = (nom.wfb[e] != 0.0) ? std::fabs(nom.wfb[e]) : 1.0;
            yw.push_back(R[k].wfb[e] / dn_w);
            BarlowPoint bw = barlow(nom.wfb[e], nom.u_wfb[e],
                                    R[k].wfb[e], R[k].u_wfb[e]);
            eyw.push_back(bw.sig_d / dn_w);

            double dn_s = (nom.sigma_norm[e] != 0.0)
                        ? std::fabs(nom.sigma_norm[e]) : 1.0;
            ysig.push_back(R[k].sigma_norm[e] / dn_s);
            BarlowPoint bs = barlow(nom.sigma_norm[e], nom.u_sigma_norm[e],
                                    R[k].sigma_norm[e], R[k].u_sigma_norm[e]);
            eysig.push_back(bs.sig_d / dn_s);
        }
        if(x.empty()) continue;

        TGraphErrors* gw = new TGraphErrors((int)x.size(), x.data(), yw.data(),
                                            nullptr, eyw.data());
        gw->SetName(Form("scan_%s_W_ebin%d", label, e));
        gw->SetTitle(Form("E = %.1f MeV;%s;ratio to nominal", nom.ecen[e], xtitle));
        gw->SetMarkerStyle(20); gw->SetMarkerColor(kBlue+1);
        gw->SetLineColor(kBlue+1);

        TGraphErrors* gs = new TGraphErrors((int)x.size(), x.data(), ysig.data(),
                                            nullptr, eysig.data());
        gs->SetName(Form("scan_%s_sigma_ebin%d", label, e));
        gs->SetMarkerStyle(24); gs->SetMarkerColor(kRed+1);
        gs->SetLineColor(kRed+1);

        c->cd(p+1);
        gw->Draw("AP");
        gs->Draw("P SAME");
        if(p == 0){
            TLegend* l = new TLegend(0.15, 0.75, 0.55, 0.88);
            l->SetBorderSize(0); l->SetFillStyle(0);
            l->AddEntry(gw, "W(0)/W(90)",     "lp");
            l->AddEntry(gs, "#sigma_{norm}",  "lp");
            l->Draw();
        }
        gw->Write(); gs->Write();
    }

    c->SaveAs((S.outdir + "scan_" + label + "_" + S.name + ".pdf").c_str());
    c->Write();
    fo->Close();

    std::cout << "\n[syst] plateau scan '" << label << "' written for "
              << S.name << ".\n"
              << "       Quote the systematic as the spread over the FLAT "
                 "region, not the full scan range.\n";
    cutVar() = varNominal();
}

// ===========================================================================
//  RUN ONE SAMPLE
// ===========================================================================
static void runSample(const SampleCfg& S, const std::string& mode){

    gSystem->mkdir(S.outdir.c_str(), kTRUE);

    std::cout << "\n\n##################################################\n"
              << "###  SAMPLE " << S.name << "   mode = " << mode << "\n"
              << "##################################################\n";

    std::cout << "\n=== NOMINAL ===\n";
    ChainResult nom = runChain(S, varNominal(), 0);
    if(!nom.ok){
        std::cerr << "[syst] nominal chain FAILED for " << S.name << "\n";
        return;
    }
    const int nE = nom.nE;

    // ------------------------- scan modes ----------------------------------
    if(mode == "scan_ampmin"){
        std::vector<double>       xv = { 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20};
        std::vector<CutVariation> vv;
        for(double f : xv) vv.push_back(varAmpMin(f));
        plateauScan(S, nom, xv, vv, "ampmin", "amp_min / nominal");
        return;
    }
    if(mode == "scan_ampmax"){
        std::vector<double>       xv = {0.90, 0.95, 1.00, 1.05, 1.10};
        std::vector<CutVariation> vv;
        for(double f : xv) vv.push_back(varAmpMax(f));
        plateauScan(S, nom, xv, vv, "ampmax", "amp_max / nominal");
        return;
    }
    if(mode == "scan_ratio"){
        std::vector<double>       xv = {-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09};
        std::vector<CutVariation> vv;
        for(double d : xv) vv.push_back(varRatio(d));
        plateauScan(S, nom, xv, vv, "ratio", "#Delta ratio_max");
        return;
    }
    if(mode == "scan_upeak"){
        // widen/narrow the uranium-contamination window symmetrically
        std::vector<double>       xv = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
        std::vector<CutVariation> vv;
        for(double d : xv) vv.push_back(varUPeak(-d, +d));
        plateauScan(S, nom, xv, vv, "upeak", "window half-width change (ns)");
        return;
    }

    // ------------------------- full variation loop -------------------------
    std::vector<VarPair>    vars = buildVariations(S);
    std::vector<SystSource> sources;

    std::ofstream tab((S.outdir + "systematics_table_" + S.name + ".csv").c_str());
    tab << "sample,source,ebin,E_MeV,"
        << "sigma_nom,sigma_var,delta_sigma,sigd_sigma,z_sigma,syst_sigma,rel_sigma,"
        << "w_nom,w_var,delta_w,sigd_w,z_w,syst_w,rel_w\n";

    int uid = 1;
    for(auto& vp : vars){

        std::cout << "\n=== " << S.name << " : " << vp.name << " (up) ===\n";
        ChainResult up = runChain(S, vp.up, uid++);
        if(!up.ok){
            std::cerr << "[syst] " << vp.name << " up FAILED -> skipped\n";
            continue;
        }

        ChainResult dn;
        bool has_dn = vp.has_dn;
        if(has_dn){
            std::cout << "\n=== " << S.name << " : " << vp.name << " (down) ===\n";
            dn = runChain(S, vp.dn, uid++);
            if(!dn.ok){
                std::cerr << "[syst] " << vp.name
                          << " down FAILED -> treated as one-sided\n";
                has_dn = false;
            }
        }

        SystSource s;
        s.name                = vp.name;
        s.correlated_in_E     = vp.corrE;
        s.correlated_in_theta = vp.corrTheta;
        s.shift_sigma = symmetrise(nom.sigma_norm, up.sigma_norm,
                                   has_dn ? dn.sigma_norm : up.sigma_norm, has_dn);
        s.shift_wfb   = symmetrise(nom.wfb, up.wfb,
                                   has_dn ? dn.wfb : up.wfb, has_dn);
        s.z_sigma.assign(nE, 0.0);
        s.z_wfb  .assign(nE, 0.0);

        for(int e = 0; e < nE; ++e){
            BarlowPoint bs = barlow(nom.sigma_norm[e], nom.u_sigma_norm[e],
                                    up.sigma_norm[e],  up.u_sigma_norm[e]);
            BarlowPoint bw = barlow(nom.wfb[e], nom.u_wfb[e],
                                    up.wfb[e],  up.u_wfb[e]);
            s.z_sigma[e] = bs.z;
            s.z_wfb  [e] = bw.z;

            if(!bs.sign && !syscfg::floor_on_fail) s.shift_sigma[e] = 0.0;
            if(!bw.sign && !syscfg::floor_on_fail) s.shift_wfb  [e] = 0.0;

            double rel_s = (nom.sigma_norm[e] != 0.0)
                         ? std::fabs(s.shift_sigma[e] / nom.sigma_norm[e]) : 0.0;
            double rel_w = (nom.wfb[e] != 0.0)
                         ? std::fabs(s.shift_wfb[e]   / nom.wfb[e])        : 0.0;

            tab << S.name << "," << vp.name << "," << e << "," << nom.ecen[e] << ","
                << nom.sigma_norm[e] << "," << up.sigma_norm[e] << ","
                << bs.delta << "," << bs.sig_d << "," << bs.z << ","
                << std::fabs(s.shift_sigma[e]) << "," << rel_s << ","
                << nom.wfb[e] << "," << up.wfb[e] << ","
                << bw.delta << "," << bw.sig_d << "," << bw.z << ","
                << std::fabs(s.shift_wfb[e]) << "," << rel_w << "\n";
        }
        sources.push_back(s);

        double m_s = 0.0, m_w = 0.0;
        int    n_sig_s = 0, n_sig_w = 0;
        for(int e = 0; e < nE; ++e){
            if(nom.sigma_norm[e] != 0.0)
                m_s = std::max(m_s, std::fabs(s.shift_sigma[e] / nom.sigma_norm[e]));
            if(nom.wfb[e] != 0.0)
                m_w = std::max(m_w, std::fabs(s.shift_wfb[e] / nom.wfb[e]));
            if(std::fabs(s.z_sigma[e]) > syscfg::barlow_z) ++n_sig_s;
            if(std::fabs(s.z_wfb  [e]) > syscfg::barlow_z) ++n_sig_w;
        }
        std::cout << "  --> " << std::left << std::setw(14) << vp.name
                  << "  max|dsigma|/sigma = " << std::fixed << std::setprecision(3)
                  << std::setw(7) << 100.0 * m_s << " %"
                  << "  (" << n_sig_s << "/" << nE << " bins pass Barlow)"
                  << "   max|dW|/W = " << std::setw(7) << 100.0 * m_w << " %"
                  << "  (" << n_sig_w << "/" << nE << ")\n";
    }
    tab.close();

    // ------------------------- covariance + plots --------------------------
    TFile* fo = TFile::Open(
        (S.outdir + "systematics_cuts_" + S.name + ".root").c_str(), "RECREATE");

    TH2D* cov_s = buildCovariance(sources, nom.ecen, true,
                                  Form("cov_sigma_%s", S.name.c_str()));
    TH2D* cov_w = buildCovariance(sources, nom.ecen, false,
                                  Form("cov_w_%s", S.name.c_str()));

    TH2D* cor_s = (TH2D*)cov_s->Clone(Form("cor_sigma_%s", S.name.c_str()));
    for(int i = 1; i <= nE; ++i)
        for(int j = 1; j <= nE; ++j){
            double di = std::sqrt(cov_s->GetBinContent(i,i));
            double dj = std::sqrt(cov_s->GetBinContent(j,j));
            cor_s->SetBinContent(i, j, (di > 0 && dj > 0)
                                 ? cov_s->GetBinContent(i,j)/(di*dj) : 0.0);
        }

    std::vector<double> x, y, ex, ey_sys, ey_tot;
    for(int e = 0; e < nE; ++e){
        if(nom.sigma_norm[e] <= 0.0) continue;
        double sd = std::sqrt(std::max(0.0, cov_s->GetBinContent(e+1, e+1)));
        x.push_back(nom.ecen[e]);
        y.push_back(nom.sigma_norm[e]);
        ex.push_back(0.0);
        ey_sys.push_back(sd);
        ey_tot.push_back(std::sqrt(sd*sd +
                         nom.u_sigma_norm[e] * nom.u_sigma_norm[e]));
    }

    if(!x.empty()){
        TGraphErrors* g_sys = new TGraphErrors((int)x.size(), x.data(), y.data(),
                                               ex.data(), ey_sys.data());
        g_sys->SetName(Form("g_sigma_syst_%s", S.name.c_str()));
        g_sys->SetFillColorAlpha(kOrange+1, 0.45);

        TGraphErrors* g_tot = new TGraphErrors((int)x.size(), x.data(), y.data(),
                                               ex.data(), ey_tot.data());
        g_tot->SetName(Form("g_sigma_tot_%s", S.name.c_str()));
        g_tot->SetMarkerStyle(20);
        g_tot->SetTitle(Form("#sigma_{norm}(%s) with selection systematics;"
                             "E_{n} (MeV);#sigma (barn)", S.name.c_str()));

        TCanvas* cc = new TCanvas(Form("c_syst_%s", S.name.c_str()),
                                  "cut systematics", 1000, 700);
        cc->SetLogx();
        g_tot->Draw("AP");
        g_sys->Draw("E2 SAME");
        g_tot->Draw("P SAME");
        cc->SaveAs((S.outdir + "sigma_with_cut_systematics_" + S.name + ".pdf").c_str());
        g_sys->Write(); g_tot->Write(); cc->Write();
    }

    // per-source overlay: sigma and W on separate canvases
    for(int which = 0; which < 2; ++which){
        const bool useSigma = (which == 0);
        TCanvas* c2 = new TCanvas(
            Form("c_src_%s_%s", useSigma ? "sigma" : "W", S.name.c_str()),
            "per-source shifts", 1000, 700);
        c2->SetLogx();
        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.14, 0.60, 0.50, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);

        int col = 1;
        for(const auto& s : sources){
            std::vector<double> xx, yy;
            for(int e = 0; e < nE; ++e){
                double den = useSigma ? nom.sigma_norm[e] : nom.wfb[e];
                if(den <= 0.0) continue;
                double sh  = useSigma ? s.shift_sigma[e] : s.shift_wfb[e];
                xx.push_back(nom.ecen[e]);
                yy.push_back(100.0 * sh / den);
            }
            if(xx.empty()) continue;
            TGraph* g = new TGraph((int)xx.size(), xx.data(), yy.data());
            g->SetLineColor(col); g->SetLineWidth(2); g->SetMarkerColor(col);
            g->SetName(Form("shift_%s_%s_%s",
                            useSigma ? "sigma" : "W",
                            s.name.c_str(), S.name.c_str()));
            mg->Add(g, "L");
            leg->AddEntry(g, s.name.c_str(), "l");
            g->Write();
            if(++col == 5 || col == 10) ++col;
        }
        mg->SetTitle(Form("Relative shift of %s per variation (%s);"
                          "E_{n} (MeV);shift (%%)",
                          useSigma ? "#sigma_{norm}" : "W(0)/W(90)",
                          S.name.c_str()));
        mg->Draw("A");
        leg->Draw();
        c2->SaveAs((S.outdir + std::string("per_source_shifts_") +
                    (useSigma ? "sigma_" : "W_") + S.name + ".pdf").c_str());
        c2->Write();
    }

    cov_s->Write(); cov_w->Write(); cor_s->Write();
    fo->Close();

    cutVar() = varNominal();

    std::cout << "\n[syst] " << S.name << " done.\n"
              << "  table : " << S.outdir << "systematics_table_" << S.name << ".csv\n"
              << "  root  : " << S.outdir << "systematics_cuts_" << S.name << ".root\n"
              << "  NOTE  : cov_sigma_* is the CORRELATED covariance. Do not\n"
              << "          collapse it to a diagonal band before fitting a2/a4\n"
              << "          or computing EXFOR pulls.\n";
}

// ===========================================================================
//  ENTRY POINT
// ===========================================================================
void systematics_cuts(const char* sample = "U238", const char* mode = "all"){

    gStyle->SetOptStat(0);

    const std::string s(sample), m(mode);

    if(s == "U238" || s == "uranium" || s == "U"){
        runSample(uraniumCfg(), m);
    }
    else if(s == "Au197" || s == "gold" || s == "Au"){
        runSample(goldCfg(), m);
    }
    else if(s == "both" || s == "all"){
        runSample(uraniumCfg(), m);
        runSample(goldCfg(),    m);
    }
    else{
        std::cerr << "[syst] unknown sample '" << s << "'. "
                     "Use U238, Au197 or both.\n";
    }
}