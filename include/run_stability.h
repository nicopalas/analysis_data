#pragma once
#include "types.h"
#include "config.h"
#include "cuts.h"
#include "constants.h"
#include "histograms.h"          // compute_angles, findBin
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TMath.h"
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

// ===========================================================================
//  Run selection
// ===========================================================================
struct RunFilter {
    int           run_min = -1;          // -1 = no lower bound
    int           run_max = -1;          // -1 = no upper bound
    std::set<int> excluded;              // runs vetoed for known reasons
    std::set<int> allowed;               // if non-empty, ONLY these runs

    bool accept(int r) const {
        if (excluded.count(r))                  return false;
        if (!allowed.empty() && !allowed.count(r)) return false;
        if (run_min >= 0 && r < run_min)        return false;
        if (run_max >= 0 && r > run_max)        return false;
        return true;
    }
};

// list of distinct run numbers actually present in the tree, sorted
static std::vector<int> getRunList(TTree* tree, const RunFilter& base)
{
    int RunNumber;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("RunNumber", 1);
    tree->SetBranchAddress("RunNumber", &RunNumber);

    std::set<int> runs;
    const Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (base.accept(RunNumber)) runs.insert(RunNumber);
    }

    tree->SetBranchStatus("*", 1);
    tree->ResetBranchAddresses();

    return std::vector<int>(runs.begin(), runs.end());
}

// split the run list into consecutive blocks of `group_size` RUNS
// (blocks of runs actually present, so gaps in the numbering do not create
//  empty or unbalanced groups)
static std::vector<RunFilter> makeRunGroups(const std::vector<int>& runs,
                                            const RunFilter& base,
                                            int group_size)
{
    std::vector<RunFilter> groups;
    for (size_t i = 0; i < runs.size(); i += group_size) {
        RunFilter g   = base;
        g.run_min     = -1;
        g.run_max     = -1;
        g.allowed.clear();
        const size_t j_end = std::min(i + (size_t)group_size, runs.size());
        for (size_t j = i; j < j_end; ++j) g.allowed.insert(runs[j]);
        groups.push_back(g);
    }
    return groups;
}

// ===========================================================================
//  Per-group bookkeeping
// ===========================================================================
struct GroupDiag {
    int    run_lo   = 0;
    int    run_hi   = 0;
    int    n_runs   = 0;
    double n_events = 0.0;      // events passing all cuts
    double sum_pulse= 0.0;      // sum of PulseIntensity over accepted events
    double mean_x0  = 0.0;
    double mean_x1  = 0.0;
    double mean_y0  = 0.0;
    double mean_y1  = 0.0;
};

// ---------------------------------------------------------------------------
//  Beam-spot offsets computed ONCE over the whole dataset.
//
//  This matters: the original fillHistograms recomputes the mean positions
//  from whatever events it is given. If you let each group recompute its own
//  offsets you subtract away exactly the beam-position drift you are trying
//  to detect. So the offsets are frozen globally and the per-group means are
//  reported separately as a diagnostic.
// ---------------------------------------------------------------------------
static void computeGlobalOffsets(TTree* tree, const RunFilter& rf,
                                 double& off_x0, double& off_x1,
                                 double& off_y0, double& off_y1)
{
    double x0, x1, y0, y1;
    int RunNumber;
    tree->SetBranchAddress("x0", &x0);
    tree->SetBranchAddress("x1", &x1);
    tree->SetBranchAddress("y0", &y0);
    tree->SetBranchAddress("y1", &y1);
    tree->SetBranchAddress("RunNumber", &RunNumber);

    double sx0 = 0, sx1 = 0, sy0 = 0, sy1 = 0;
    Long64_t n = 0;
    const Long64_t nent = tree->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        tree->GetEntry(i);
        if (!rf.accept(RunNumber)) continue;
        sx0 += x0; sx1 += x1; sy0 += y0; sy1 += y1;
        ++n;
    }
    if (n == 0) { off_x0 = off_x1 = off_y0 = off_y1 = 0.0; return; }
    off_x0 = sx0 / n; off_x1 = sx1 / n;
    off_y0 = sy0 / n; off_y1 = sy1 / n;

    tree->ResetBranchAddresses();
}

// ---------------------------------------------------------------------------
//  Same selection as fillHistograms, but with an explicit run filter and
//  externally supplied offsets.
//
//  NOTE ON A DISCREPANCY IN THE ORIGINAL CODE:
//  fillHistograms computes cos_theta_corrected / cos_theta_det_corrected via
//  compute_angles and then never uses them: the cuts and the binning use the
//  branch values cos_theta / cos_theta_det instead, and the local phi/phi_det
//  shadow the branch variables. To make the stability test reflect the real
//  analysis, the default here reproduces that behaviour. Flip
//  use_corrected_angles to true only if the shadowing was a bug you want fixed
//  — but then fix it in fillHistograms too, or the two will not be comparable.
// ---------------------------------------------------------------------------
static void fillHistogramsRuns(
    TTree* tree,
    const AnalysisConfig& cfg,
    const RunFilter& rf,
    std::vector<TH1D*>& hists_tof,
    Vec3D& counts_roi,
    Vec3D& counts_bkg,
    Vec3D& counts_upeak,
    double emin, double emax,
    double off_x0, double off_x1, double off_y0, double off_y1,
    GroupDiag& diag,
    bool use_corrected_angles = false)
{
    const int nbins = (int)cfg.energy_bins.size() - 1;

    double cos_theta_det, cos_theta;
    double tof1, tof0, neutron_energy;
    double amp0, amp1;
    double x0, x1, y0, y1;
    double phi_br, phi_det_br;
    int    RunNumber;
    float  PSpulse;

    tree->SetBranchAddress("tof1",           &tof1);
    tree->SetBranchAddress("tof0",           &tof0);
    tree->SetBranchAddress("amp0",           &amp0);
    tree->SetBranchAddress("amp1",           &amp1);
    tree->SetBranchAddress("x0",             &x0);
    tree->SetBranchAddress("x1",             &x1);
    tree->SetBranchAddress("y0",             &y0);
    tree->SetBranchAddress("y1",             &y1);
    tree->SetBranchAddress("neutron_energy", &neutron_energy);
    tree->SetBranchAddress("cos_theta",      &cos_theta);
    tree->SetBranchAddress("cos_theta_det",  &cos_theta_det);
    tree->SetBranchAddress("phi",            &phi_br);
    tree->SetBranchAddress("phi_det",        &phi_det_br);
    tree->SetBranchAddress("RunNumber",      &RunNumber);
    tree->SetBranchAddress("PulseIntensity", &PSpulse);

    double sx0 = 0, sx1 = 0, sy0 = 0, sy1 = 0;
    std::set<int> runs_seen;

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (!rf.accept(RunNumber))                        continue;
        if (neutron_energy < emin || neutron_energy > emax) continue;
        if (neutron_energy > 2000)                        continue;

        const int e_bin = findBin(cfg.energy_bins, neutron_energy);
        if (e_bin < 0 || e_bin >= nbins) continue;

        const EventCuts c = getCuts(cfg.sample, neutron_energy);
        if (!passAmplitudeCut(amp0, amp1, c)) continue;

        double ct_corr, ctd_corr, phi_corr, phid_corr;
        compute_angles(x0, y0, x1, y1,
                       off_x0, off_x1, off_y0, off_y1,
                       ctd_corr, phid_corr, phi_corr, ct_corr);

        const double ct  = use_corrected_angles ? ct_corr  : cos_theta;
        const double ctd = use_corrected_angles ? ctd_corr : cos_theta_det;

        if (ctd < 0.0)                                     continue;
        if (std::fabs(ctd) > 1 || std::fabs(ct) > 1)       continue;

        const double dt = tof1 - tof0;
        hists_tof[e_bin]->Fill(dt);

        const int j  = int(std::fabs(ct) / dcos_beam);
        const int ii = int(ctd / dcos_det);
        if (j  >= nbins_beam) continue;
        if (ii >= nbins_det)  continue;

        if (dt >= c.roi_min && dt <= c.roi_max) counts_roi  [e_bin][j][ii]++;
        else if (inUraniumPeak(dt, c))          counts_upeak[e_bin][j][ii]++;
        else                                    counts_bkg  [e_bin][j][ii]++;

        diag.n_events  += 1.0;
        diag.sum_pulse += PSpulse;
        sx0 += x0; sx1 += x1; sy0 += y0; sy1 += y1;
        runs_seen.insert(RunNumber);
    }

    if (diag.n_events > 0) {
        diag.mean_x0 = sx0 / diag.n_events;
        diag.mean_x1 = sx1 / diag.n_events;
        diag.mean_y0 = sy0 / diag.n_events;
        diag.mean_y1 = sy1 / diag.n_events;
    }
    diag.n_runs = (int)runs_seen.size();
    if (!runs_seen.empty()) {
        diag.run_lo = *runs_seen.begin();
        diag.run_hi = *runs_seen.rbegin();
    }

    tree->ResetBranchAddresses();
}

// ===========================================================================
//  Comparison across groups
// ===========================================================================
struct StabilityTest {
    std::string name;
    double      weighted_mean = 0.0;
    double      u_mean        = 0.0;
    double      chi2          = 0.0;
    int         ndf           = 0;
    double      pvalue        = 1.0;
    std::vector<double> pull;      // (y_i - mean) / sigma_i
    std::vector<int>    outliers;  // indices with |pull| > 3
};

// Constant-hypothesis test: are all the groups compatible with one value?
static StabilityTest testConstant(const std::string& name,
                                  const std::vector<double>& y,
                                  const std::vector<double>& ey,
                                  double pull_threshold = 3.0)
{
    StabilityTest t;
    t.name = name;

    double sw = 0.0, swy = 0.0;
    int    n  = 0;
    for (size_t i = 0; i < y.size(); ++i) {
        if (!(ey[i] > 0.0) || !std::isfinite(y[i])) continue;
        const double w = 1.0 / (ey[i] * ey[i]);
        sw  += w;
        swy += w * y[i];
        ++n;
    }
    if (sw <= 0.0 || n < 2) return t;

    t.weighted_mean = swy / sw;
    t.u_mean        = std::sqrt(1.0 / sw);

    t.pull.assign(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) {
        if (!(ey[i] > 0.0) || !std::isfinite(y[i])) continue;
        const double p = (y[i] - t.weighted_mean) / ey[i];
        t.pull[i] = p;
        t.chi2   += p * p;
        if (std::fabs(p) > pull_threshold) t.outliers.push_back((int)i);
    }
    t.ndf    = n - 1;
    t.pvalue = (t.ndf > 0) ? TMath::Prob(t.chi2, t.ndf) : 1.0;
    return t;
}

static void printStabilityTest(const StabilityTest& t,
                               const std::vector<GroupDiag>& diags)
{
    std::cout << "\n--- " << t.name << " ---\n"
              << "  weighted mean = " << t.weighted_mean
              << " +/- " << t.u_mean << "\n"
              << "  chi2/ndf = " << t.chi2 << "/" << t.ndf;
    if (t.ndf > 0) std::cout << " = " << t.chi2 / t.ndf;
    std::cout << "   p = " << t.pvalue << "\n";

    for (size_t i = 0; i < t.pull.size(); ++i) {
        std::cout << "   g" << std::setw(2) << i
                  << " [" << diags[i].run_lo << "-" << diags[i].run_hi << "]"
                  << "  pull = " << std::setw(8) << std::fixed
                  << std::setprecision(2) << t.pull[i];
        if (std::fabs(t.pull[i]) > 3.0) std::cout << "   <== OUTLIER";
        std::cout << "\n";
    }
    std::cout.unsetf(std::ios::fixed);

    // An overall chi2/ndf well above 1 with no single outlier means a drift or
    // an underestimated per-group error, not a bad run. Look at the pull
    // sequence: monotonic = drift, scattered = extra unaccounted variance.
    if (t.ndf > 0 && t.chi2 / t.ndf > 2.0 && t.outliers.empty())
        std::cout << "   note: spread exceeds the quoted errors with no single"
                     " outlier -> drift or underestimated uncertainties\n";
}

// scatter plot of one observable vs group, with the weighted mean band
static void plotStability(const StabilityTest& t,
                          const std::vector<double>& y,
                          const std::vector<double>& ey,
                          const std::vector<GroupDiag>& diags,
                          const std::string& outfile)
{
    const int n = (int)y.size();
    std::vector<double> x(n), ex(n, 0.0);
    for (int i = 0; i < n; ++i) x[i] = i;

    TCanvas* c = new TCanvas(Form("c_stab_%s", t.name.c_str()), "", 900, 500);
    TGraphErrors* g = new TGraphErrors(n, x.data(),
                                       const_cast<double*>(y.data()),
                                       ex.data(),
                                       const_cast<double*>(ey.data()));
    g->SetTitle(Form("%s;run group;%s", t.name.c_str(), t.name.c_str()));
    g->SetMarkerStyle(20);
    g->Draw("AP");

    const double lo = g->GetXaxis()->GetXmin();
    const double hi = g->GetXaxis()->GetXmax();

    TBox* band = new TBox(lo, t.weighted_mean - t.u_mean,
                          hi, t.weighted_mean + t.u_mean);
    band->SetFillColorAlpha(kAzure + 1, 0.25);
    band->Draw("same");

    TLine* line = new TLine(lo, t.weighted_mean, hi, t.weighted_mean);
    line->SetLineColor(kAzure + 2);
    line->SetLineWidth(2);
    line->Draw("same");

    g->Draw("P same");
    c->Print(outfile.c_str());
    delete c;
}