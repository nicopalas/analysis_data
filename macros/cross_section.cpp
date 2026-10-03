// =============================================================================
//  fission_xs.C
//
//  Fission cross sections of Au-197 and U-238 from PPAC coincidence data.
//
//  Usage (from ROOT):
//      root -l fission_xs.C                       // runs both samples
//      root -l -e '.L fission_xs.C' -e 'gold_xs()'
//      root -l -e '.L fission_xs.C' -e 'uranium_xs()'
//
// -----------------------------------------------------------------------------
//  METHOD (per neutron-energy bin e)
//
//      N_corr(e) = sum_{b,d}  N(e,b,d) / ( A(b,d) * eps(e,d) )
//      sigma(e)  = N_corr(e) / ( Phi(e) * n_areal )          [cm^2 -> barn]
//
//      N(e,b,d)  counts in the ToF ROI, in beam-angle bin b (|cos theta|)
//                and detector-angle bin d (cos theta_det)
//      A(b,d)    acceptance of cell (b,d), from acceptance_coincidence.csv
//      eps(e,d)  detection efficiency (heff_ebin<k> histograms)
//      Phi(e)    number of neutrons in bin e (integrated evaluated flux)
//      n_areal   areal density of the target [atoms/cm^2]
//
//  Cells with A <= kAccMin or eps <= kEpsMin are dropped. The fraction of
//  counts kept is saved as h_kept_frac_<tag> so you can check that it does
//  not vary strongly with energy (if it does, the dropped angular region
//  biases the SHAPE of sigma(E), not only its scale).
//
//  The raw cross section is then scaled to reference data with a weighted
//  mean of the ratios sigma_ref / sigma_raw over all bins whose centre lies
//  INSIDE the reference energy range (no extrapolation). If there is no
//  overlap, no normalisation is applied and the absolute result is plotted,
//  with a clear warning.
//
// -----------------------------------------------------------------------------
//  ASSUMPTIONS TO VERIFY  (search for "VERIFY" in the code)
//
//   1. hEval_Abs: x axis in eV, content per unit lethargy (E dPhi/dE).
//      If it is dPhi/dE instead, set kFluxPerLethargy = false.
//   2. The flux histogram corresponds to the SAME runs used here (the bad
//      runs are excluded from the counts). If the evaluated flux is per
//      pulse / per nominal proton bunch, set kFluxScale accordingly.
//      Irrelevant after normalisation, essential for absolute results.
//   3. A(b,d) is a dimensionless detection probability. If it is a solid
//      angle in sr, N/A is dN/dOmega and must be multiplied by the solid
//      angle of the cell before summing.
//   4. n_areal is an areal density [atoms/cm^2], consistent with Phi being
//      a neutron count (not a fluence per cm^2).
//   5. heff_ebin<k> has nbins_det bins matching dcos_det, and k indexes
//      the edges given in SampleSetup::eff_bins (checked at run time for
//      the number of bins only).
//
// -----------------------------------------------------------------------------
//  CHANGES WITH RESPECT TO THE PREVIOUS VERSION
//
//   - Normalisation: weighted mean of ratios over the overlap region, no
//     TGraph::Eval extrapolation; guards against empty bins; the
//     normalisation uncertainty and chi2/ndf are reported and saved.
//   - U-238 reference values are already in barn (~1.6 b at 40-50 MeV);
//     the wrong mb -> barn factor was removed.
//   - Efficiency uncertainty: contributions of the same eps(d) are summed
//     linearly over beam bins before squaring (they are fully correlated).
//   - Flux: edge bins are weighted by their fractional overlap (in ln E),
//     so no flux bin is counted twice; coverage is checked.
//   - Geometric bin centres used consistently; graphs use asymmetric
//     x errors equal to the true bin edges.
//   - Removed dead code: unused beam-spot pre-pass, unused TCutG loading,
//     unused branches (PulseIntensity float/double mismatch risk gone).
//   - A single bad-run list for everything.
//   - Events at |cos| exactly equal to the upper grid edge are kept.
//   - Gold and uranium share one implementation; only SampleSetup differs.
//   - Efficiency-bin mapping warns if a fine bin straddles two eff bins.
// =============================================================================

#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/types.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/constants.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/acceptance.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/config.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/histograms.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cuts.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TParameter.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace xs {

// =============================================================================
//  Global settings
// =============================================================================
const std::string kBaseDir  = "/Users/nico/Desktop/Tese/Analysis/cross_section/";
const std::string kDataFile = kBaseDir + "data/coincidences_final.root";
const std::string kFluxFile = kBaseDir + "data/flux_data/evalFlux_prelim.root";
const std::string kFluxHist = "hEval_Abs";
const std::string kAccFile  = kBaseDir + "acceptance_coincidence.csv";

// Neutron-energy binning of the cross section (log bins)
constexpr int    kNbins = 20;
constexpr double kEmin  = 60.0;    // MeV
constexpr double kEmax  = 1000.0;  // MeV

constexpr double kMeVtoeV = 1.0e6;
constexpr double kBarn    = 1.0e-24;   // cm^2

// Flux interpretation                                             VERIFY (1,2)
constexpr bool   kFluxPerLethargy = true;  // content = E dPhi/dE
constexpr double kFluxScale       = 1.0;   // e.g. N_protons(used runs)/N_ref

// Cell-selection thresholds for the acceptance/efficiency correction
constexpr double kAccMin = 1.0e-6;   // minimum acceptance of a cell
constexpr double kEpsMin = 0.6;      // minimum reliable efficiency

// Runs excluded from the analysis (one list, used everywhere)
const std::set<int> kBadRuns = {118771, 118789, 118668, 118587, 118751};

// =============================================================================
//  Data structures
// =============================================================================

// Reference cross-section points. ALWAYS in barn, E strictly increasing.
struct ReferenceData {
    std::string         label;
    std::vector<double> E;        // MeV
    std::vector<double> sigma;    // barn
    std::vector<double> u_sigma;  // barn
};

// Everything that differs between the gold and uranium analyses.
struct SampleSetup {
    std::string         name;          // for plots, e.g. "Au-197"
    std::string         tag;           // for object names, e.g. "gold"
    Sample              sample;        // passed to getCuts()
    std::string         tree_name;     // tree in kDataFile
    std::string         eff_path;      // file with heff_ebin<k>
    std::vector<double> eff_bins;      // MeV edges of the efficiency bins
    double              n_areal;       // atoms/cm^2                  VERIFY (4)
    double              sumdiff0_min;  // keep if sumX0 - sumY0 >= this
    double              sumdiff1_max;  // keep if sumX1 - sumY1 <= this
    ReferenceData       ref;
    std::string         out_path;
};

// Energy binning with geometric bin centres (natural for log bins).
struct EnergyBinning {
    std::vector<double> edges;  // MeV, size n()+1
    int    n()         const { return (int)edges.size() - 1; }
    double lo(int e)   const { return edges[e]; }
    double hi(int e)   const { return edges[e + 1]; }
    double center(int e) const { return std::sqrt(edges[e] * edges[e + 1]); }
};

struct CorrectedCounts {
    std::vector<double> N, u_N;   // corrected counts and uncertainty
    std::vector<double> raw;      // counts in the ROI, all cells
    std::vector<double> kept;     // counts in cells used for the correction
    std::vector<bool>   has_eff;  // efficiency available for this bin
};

struct Normalisation {
    bool   ok      = false;
    double scale   = 1.0;
    double u_scale = 0.0;
    double chi2    = 0.0;
    int    n_used  = 0;
};

EnergyBinning makeBinning()
{
    EnergyBinning eb;
    eb.edges = buildLogBins(kNbins, kEmin, kEmax, 0.8);
    return eb;
}

// =============================================================================
//  Flux: number of neutrons per energy bin
// -----------------------------------------------------------------------------
//  Each flux-histogram bin [lo,hi] contributes only over its overlap
//  [max(lo,a), min(hi,b)] with the cross-section bin [a,b]:
//    lethargy histogram:  dPhi = content * ln(ohi/olo)   (exact for a
//                         histogram that is flat per unit lethargy)
//    dPhi/dE histogram:   dPhi = content * (ohi - olo)
//  Errors are added in quadrature, i.e. flux bins are treated as
//  uncorrelated. Evaluated-flux errors are usually correlated, so u_flux is
//  a lower bound; after normalisation the common part cancels anyway.
// =============================================================================
void integrateFlux(const TH1D &h, const EnergyBinning &eb,
                   std::vector<double> &flux, std::vector<double> &u_flux)
{
    flux.assign(eb.n(), 0.0);
    u_flux.assign(eb.n(), 0.0);
    const TAxis *ax = h.GetXaxis();

    for (int e = 0; e < eb.n(); ++e) {
        const double a = eb.lo(e) * kMeVtoeV;
        const double b = eb.hi(e) * kMeVtoeV;
        if (a < ax->GetXmin() || b > ax->GetXmax())
            std::cerr << "[flux] WARNING: bin " << e << " [" << eb.lo(e) << ", "
                      << eb.hi(e) << "] MeV not fully covered by " << kFluxHist << "\n";

        double sum = 0.0, var = 0.0;
        for (int k = 1; k <= h.GetNbinsX(); ++k) {
            const double olo = std::max(ax->GetBinLowEdge(k), a);
            const double ohi = std::min(ax->GetBinUpEdge(k),  b);
            if (ohi <= olo || olo <= 0.0) continue;

            const double width = kFluxPerLethargy ? std::log(ohi / olo) : (ohi - olo);
            sum += h.GetBinContent(k) * width;
            var += std::pow(h.GetBinError(k) * width, 2);
        }
        flux[e]   = sum * kFluxScale;
        u_flux[e] = std::sqrt(var) * kFluxScale;
    }
}

// =============================================================================
//  Event loop: fill N(e, b, d)
// =============================================================================
Vec3D countEvents(TTree *t, const SampleSetup &s, const EnergyBinning &eb)
{
    double tof0, tof1, E, amp0, amp1, cth, cthd, sx0, sx1, sy0, sy1;
    int run;

    // Read only the branches that are used (faster, and avoids type
    // mismatches on branches we do not need).
    const std::vector<std::string> used = {
        "tof0", "tof1", "neutron_energy", "amp0", "amp1", "cos_theta",
        "cos_theta_det", "sumX0", "sumX1", "sumY0", "sumY1", "RunNumber"};
    t->SetBranchStatus("*", 0);
    for (const auto &br : used) t->SetBranchStatus(br.c_str(), 1);

    t->SetBranchAddress("tof0",           &tof0);
    t->SetBranchAddress("tof1",           &tof1);
    t->SetBranchAddress("neutron_energy", &E);
    t->SetBranchAddress("amp0",           &amp0);
    t->SetBranchAddress("amp1",           &amp1);
    t->SetBranchAddress("cos_theta",      &cth);
    t->SetBranchAddress("cos_theta_det",  &cthd);
    t->SetBranchAddress("sumX0",          &sx0);
    t->SetBranchAddress("sumX1",          &sx1);
    t->SetBranchAddress("sumY0",          &sy0);
    t->SetBranchAddress("sumY1",          &sy1);
    t->SetBranchAddress("RunNumber",      &run);

    Vec3D counts(eb.n(), Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    const double cos_beam_max = nbins_beam * dcos_beam;  // upper edge of grid
    const double cos_det_max  = nbins_det  * dcos_det;

    const Long64_t nentries = t->GetEntries();
    Long64_t n_pass = 0;
    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        if (kBadRuns.count(run)) continue;
        if (E < eb.edges.front() || E >= eb.edges.back()) continue;
        const int ie = findBin(eb.edges, E);
        if (ie < 0 || ie >= eb.n()) continue;

        // Amplitude and ToF-ROI cuts (energy dependent)
        const EventCuts c = getCuts(s.sample, E);
        if (!passAmplitudeCut(amp0, amp1, c)) continue;
        const double dt = tof1 - tof0;
        if (dt < c.roi_min || dt > c.roi_max) continue;


        // Angular binning. The "!(x <= y)" form also rejects NaN.
        const double acth = std::fabs(cth);
        if (!(acth <= 1.0) || !(cthd >= 0.0 && cthd <= 1.0)) continue;

        int ib = int(acth / dcos_beam);
        int id = int(cthd / dcos_det);
        // Values exactly on the upper grid edge belong to the last bin
        if (ib == nbins_beam && std::fabs(acth - cos_beam_max) < 1e-12) --ib;
        if (id == nbins_det  && std::fabs(cthd - cos_det_max)  < 1e-12) --id;
        if (ib < 0 || ib >= nbins_beam || id < 0 || id >= nbins_det) continue;

        counts[ie][ib][id] += 1.0;
        ++n_pass;
    }
    std::cout << "[events] " << s.tree_name << ": " << n_pass << " / " << nentries
              << " entries pass all cuts\n";
    return counts;
}

// =============================================================================
//  Efficiency: eps(e, d) for every cross-section bin e
// -----------------------------------------------------------------------------
//  Each fine bin is mapped to the efficiency bin containing its geometric
//  centre. A warning is printed if the fine bin straddles two efficiency
//  bins. eps[e] is left empty when no efficiency is available; that bin is
//  then not corrected and gets no cross-section point.
// =============================================================================
void loadEfficiency(TFile *f, const SampleSetup &s, const EnergyBinning &eb,
                    std::vector<std::vector<double>> &eps,
                    std::vector<std::vector<double>> &u_eps)
{
    eps.assign(eb.n(), {});
    u_eps.assign(eb.n(), {});
    const int n_eff = (int)s.eff_bins.size() - 1;
    std::map<int, std::pair<std::vector<double>, std::vector<double>>> cache;
    std::set<int> missing;

    for (int e = 0; e < eb.n(); ++e) {
        const int k    = findBin(s.eff_bins, eb.center(e));
        const int k_lo = findBin(s.eff_bins, eb.lo(e));
        const int k_hi = findBin(s.eff_bins, eb.hi(e) * (1.0 - 1e-12));

        if (k < 0 || k >= n_eff) {
            std::cerr << "[eff] WARNING: E = " << eb.center(e) << " MeV outside the "
                      << "efficiency range, bin " << e << " not corrected\n";
            continue;
        }
        if (k_lo != k_hi)
            std::cerr << "[eff] note: bin " << e << " [" << eb.lo(e) << ", " << eb.hi(e)
                      << "] MeV straddles efficiency bins " << k_lo << "/" << k_hi
                      << ", using bin " << k << "\n";

        if (!cache.count(k)) {
            if (missing.count(k)) continue;
            TH1D *h = dynamic_cast<TH1D *>(f->Get(Form("heff_ebin%d", k)));
            if (!h) {
                std::cerr << "[eff] WARNING: heff_ebin" << k << " not found\n";
                missing.insert(k);
                continue;
            }
            if (h->GetNbinsX() != nbins_det) {                        // VERIFY (5)
                std::cerr << "[eff] WARNING: heff_ebin" << k << " has " << h->GetNbinsX()
                          << " bins, expected nbins_det = " << nbins_det << "\n";
                missing.insert(k);
                continue;
            }
            std::vector<double> v(nbins_det), u(nbins_det);
            for (int d = 0; d < nbins_det; ++d) {
                v[d] = h->GetBinContent(d + 1);
                u[d] = h->GetBinError(d + 1);
            }
            cache[k] = {v, u};
        }
        eps[e]   = cache[k].first;
        u_eps[e] = cache[k].second;
    }
}

// =============================================================================
//  Acceptance / efficiency correction
// -----------------------------------------------------------------------------
//  N = sum_{b,d} n_bd / (A_bd * eps_d)
//
//  Uncertainty:
//   - Poisson:    var += n_bd / (A_bd * eps_d)^2              (per cell)
//   - efficiency: eps_d is the SAME number for all beam bins b, so
//                   dN/deps_d = -sum_b n_bd / (A_bd * eps_d^2) = -S_d
//                 and var += (S_d * u_eps_d)^2                (per d)
//  Correlations between different d (eps(d) sharing one reference bin in
//  the efficiency computation) are still neglected.
// =============================================================================
CorrectedCounts correctCounts(const Vec3D &counts, const Vec2D &acc,
                              const std::vector<std::vector<double>> &eps,
                              const std::vector<std::vector<double>> &u_eps,
                              const EnergyBinning &eb)
{
    CorrectedCounts r;
    const int n = eb.n();
    r.N.assign(n, 0.0);   r.u_N.assign(n, 0.0);
    r.raw.assign(n, 0.0); r.kept.assign(n, 0.0);
    r.has_eff.assign(n, false);

    for (int e = 0; e < n; ++e) {
        for (int b = 0; b < nbins_beam; ++b)
            for (int d = 0; d < nbins_det; ++d) r.raw[e] += counts[e][b][d];

        if (eps[e].empty()) continue;
        r.has_eff[e] = true;

        double N = 0.0, var_stat = 0.0, var_eff = 0.0;
        for (int d = 0; d < nbins_det; ++d) {
            const double ep = eps[e][d];
            if (ep <= kEpsMin) continue;

            double S_d = 0.0;
            for (int b = 0; b < nbins_beam; ++b) {
                const double nc = counts[e][b][d];
                const double A  = acc[b][d];
                if (A <= kAccMin || nc <= 0.0) continue;

                const double w = A * ep;
                r.kept[e] += nc;
                N         += nc / w;
                var_stat  += nc / (w * w);
                S_d       += nc / (A * ep * ep);
            }
            var_eff += std::pow(S_d * u_eps[e][d], 2);
        }
        r.N[e]   = N;
        r.u_N[e] = std::sqrt(var_stat + var_eff);
    }
    return r;
}

// =============================================================================
//  Raw cross section  sigma = N_corr / (Phi * n_areal)  [barn]
// =============================================================================
TH1D *computeCrossSection(const CorrectedCounts &cc,
                          const std::vector<double> &flux,
                          const std::vector<double> &u_flux,
                          const EnergyBinning &eb, double n_areal,
                          const std::string &name)
{
    TH1D *h = new TH1D(name.c_str(), ";E_{n} (MeV);#sigma (barn)",
                       eb.n(), eb.edges.data());
    h->SetDirectory(nullptr);

    for (int e = 0; e < eb.n(); ++e) {
        if (flux[e] <= 0.0 || cc.N[e] <= 0.0) continue;
        const double sigma = cc.N[e] / (flux[e] * n_areal) / kBarn;
        const double rel   = std::sqrt(std::pow(cc.u_N[e] / cc.N[e], 2) +
                                       std::pow(u_flux[e] / flux[e], 2));
        h->SetBinContent(e + 1, sigma);
        h->SetBinError  (e + 1, sigma * rel);
    }
    return h;
}

// =============================================================================
//  Reference interpolation (linear, ONLY inside the reference range)
//  The uncertainty is interpolated linearly as well (conservative).
// =============================================================================
bool evalReference(const ReferenceData &r, double E, double &val, double &err)
{
    if (r.E.size() < 2 || E < r.E.front() || E > r.E.back()) return false;
    std::size_t j = std::upper_bound(r.E.begin(), r.E.end(), E) - r.E.begin();
    j = std::min(j, r.E.size() - 1);
    const std::size_t i = j - 1;
    const double t = (E - r.E[i]) / (r.E[j] - r.E[i]);
    val = (1.0 - t) * r.sigma[i]   + t * r.sigma[j];
    err = (1.0 - t) * r.u_sigma[i] + t * r.u_sigma[j];
    return true;
}

// =============================================================================
//  Normalisation to the reference
// -----------------------------------------------------------------------------
//  For every bin with data whose geometric centre lies inside the reference
//  range:  R_i = sigma_ref(E_i) / sigma_raw_i, with relative errors added in
//  quadrature. The scale is the weighted mean of the R_i, weights 1/u_R^2.
//
//  u_scale is a FULLY CORRELATED normalisation uncertainty common to all
//  points. It is printed and saved separately, not added to the point
//  errors. chi2/ndf of the ratios tests shape agreement in the overlap
//  region; a large value means the raw shape disagrees with the reference.
// =============================================================================
Normalisation computeScale(const TH1D *h, const EnergyBinning &eb,
                           const ReferenceData &ref)
{
    Normalisation res;
    std::vector<std::pair<double, double>> ratios;  // (R, weight)
    double sw = 0.0, swR = 0.0;

    for (int e = 0; e < eb.n(); ++e) {
        const double x  = h->GetBinContent(e + 1);
        const double ux = h->GetBinError(e + 1);
        if (x <= 0.0 || ux <= 0.0) continue;

        double r, ur;
        if (!evalReference(ref, eb.center(e), r, ur) || r <= 0.0) continue;

        const double R  = r / x;
        const double uR = R * std::sqrt(std::pow(ur / r, 2) + std::pow(ux / x, 2));
        const double w  = 1.0 / (uR * uR);
        ratios.push_back({R, w});
        sw  += w;
        swR += w * R;
        printf("[norm] E = %7.2f MeV  raw = %.4e  ref = %.4e  ratio = %.4e +/- %.4e\n",
               eb.center(e), x, r, R, uR);
    }
    if (ratios.empty()) return res;

    res.ok      = true;
    res.n_used  = (int)ratios.size();
    res.scale   = swR / sw;
    res.u_scale = 1.0 / std::sqrt(sw);
    for (const auto &p : ratios) res.chi2 += p.second * std::pow(p.first - res.scale, 2);
    return res;
}

// =============================================================================
//  Graph helpers
// =============================================================================
TGraphAsymmErrors *makeGraph(const TH1D *h, const EnergyBinning &eb,
                             const std::string &name)
{
    std::vector<double> x, y, exl, exh, ey;
    for (int e = 0; e < eb.n(); ++e) {
        const double v = h->GetBinContent(e + 1);
        if (v <= 0.0) continue;
        x.push_back(eb.center(e));
        y.push_back(v);
        exl.push_back(eb.center(e) - eb.lo(e));
        exh.push_back(eb.hi(e) - eb.center(e));
        ey.push_back(h->GetBinError(e + 1));
    }
    auto *g = new TGraphAsymmErrors((int)x.size(), x.data(), y.data(),
                                    exl.data(), exh.data(), ey.data(), ey.data());
    g->SetName(name.c_str());
    return g;
}

TGraphErrors *makeReferenceGraph(const ReferenceData &r, const std::string &name)
{
    std::vector<double> ex(r.E.size(), 0.0);
    auto *g = new TGraphErrors((int)r.E.size(), r.E.data(), r.sigma.data(),
                               ex.data(), r.u_sigma.data());
    g->SetName(name.c_str());
    g->SetTitle(r.label.c_str());
    return g;
}

// =============================================================================
//  Full analysis for one sample
// =============================================================================
void runCrossSection(const SampleSetup &s, const EnergyBinning &eb)
{
    std::cout << "\n================ " << s.name << " ================\n";

    // ── 1. Flux ────────────────────────────────────────────────────────────
    std::vector<double> flux, u_flux;
    {
        TFile *ff = TFile::Open(kFluxFile.c_str(), "READ");
        if (!ff || ff->IsZombie()) { std::cerr << "Cannot open flux file\n"; return; }
        TH1D *hf = dynamic_cast<TH1D *>(ff->Get(kFluxHist.c_str()));
        if (!hf) { std::cerr << "Flux histogram not found\n"; ff->Close(); return; }
        integrateFlux(*hf, eb, flux, u_flux);
        ff->Close();
        delete ff;
    }

    // ── 2. Acceptance ──────────────────────────────────────────────────────
    Vec2D acc_fine;
    if (!loadAcceptanceCSV(kAccFile, acc_fine)) {
        std::cerr << "Failed to load acceptance CSV\n"; return;
    }
    const Vec2D acceptance = rebin(acc_fine);                        // VERIFY (3)
    if ((int)acceptance.size() < nbins_beam || (int)acceptance[0].size() < nbins_det) {
        std::cerr << "Acceptance grid smaller than nbins_beam x nbins_det\n"; return;
    }

    // ── 3. Counts ──────────────────────────────────────────────────────────
    Vec3D counts;
    {
        TFile *fin = TFile::Open(kDataFile.c_str(), "READ");
        if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
        TTree *tin = dynamic_cast<TTree *>(fin->Get(s.tree_name.c_str()));
        if (!tin) { std::cerr << "Tree " << s.tree_name << " not found\n"; fin->Close(); return; }
        counts = countEvents(tin, s, eb);
        fin->Close();
        delete fin;
    }

    // ── 4. Efficiency ──────────────────────────────────────────────────────
    std::vector<std::vector<double>> eps, u_eps;
    {
        TFile *fe = TFile::Open(s.eff_path.c_str(), "READ");
        if (!fe || fe->IsZombie()) { std::cerr << "Cannot open " << s.eff_path << "\n"; return; }
        loadEfficiency(fe, s, eb, eps, u_eps);
        fe->Close();
        delete fe;
    }

    // ── 5. Corrected counts and raw cross section ─────────────────────────
    const CorrectedCounts cc = correctCounts(counts, acceptance, eps, u_eps, eb);
    TH1D *h_raw = computeCrossSection(cc, flux, u_flux, eb, s.n_areal,
                                      "h_cs_raw_" + s.tag);

    TH1D *h_kept = new TH1D(("h_kept_frac_" + s.tag).c_str(),
                            ";E_{n} (MeV);fraction of ROI counts used",
                            eb.n(), eb.edges.data());
    h_kept->SetDirectory(nullptr);
    TH1D *h_flux = new TH1D(("h_flux_" + s.tag).c_str(),
                            ";E_{n} (MeV);neutrons per bin",
                            eb.n(), eb.edges.data());
    h_flux->SetDirectory(nullptr);

    printf("\n ebin   E(MeV)    raw_cts  kept   N_corr +/- err          Phi            sigma_raw (b)\n");
    for (int e = 0; e < eb.n(); ++e) {
        const double kf = cc.raw[e] > 0 ? cc.kept[e] / cc.raw[e] : 0.0;
        h_kept->SetBinContent(e + 1, kf);
        h_flux->SetBinContent(e + 1, flux[e]);
        h_flux->SetBinError  (e + 1, u_flux[e]);
        printf(" %3d  %8.2f  %8.0f  %4.0f%%  %10.2f +/- %-9.2f  %.4e  %.4e +/- %.4e%s\n",
               e, eb.center(e), cc.raw[e], 100.0 * kf, cc.N[e], cc.u_N[e], flux[e],
               h_raw->GetBinContent(e + 1), h_raw->GetBinError(e + 1),
               cc.has_eff[e] ? "" : "   (no efficiency)");
    }

    // ── 6. Normalisation ───────────────────────────────────────────────────
    const Normalisation norm = computeScale(h_raw, eb, s.ref);
    TH1D *h_norm = nullptr;
    if (norm.ok) {
        h_norm = (TH1D *)h_raw->Clone(("h_cs_norm_" + s.tag).c_str());
        h_norm->SetDirectory(nullptr);
        h_norm->Scale(norm.scale);   // scales contents and errors
        const int ndf = norm.n_used - 1;
        printf("[norm] scale = %.4e +/- %.4e (%.1f%%, correlated)  from %d bin(s)",
               norm.scale, norm.u_scale, 100.0 * norm.u_scale / norm.scale, norm.n_used);
        if (ndf > 0) printf("  chi2/ndf = %.2f/%d", norm.chi2, ndf);
        printf("\n");
    } else {
        std::cerr << "[norm] WARNING: no bin of this measurement ["
                  << eb.edges.front() << ", " << eb.edges.back()
                  << "] MeV lies inside the reference range [";
        if (!s.ref.E.empty()) std::cerr << s.ref.E.front() << ", " << s.ref.E.back();
        std::cerr << "] MeV. NO normalisation applied; the ABSOLUTE cross "
                     "section is plotted. Add reference points in the "
                     "measured range to normalise.\n";
    }
    const TH1D *h_plot = norm.ok ? h_norm : h_raw;

    // ── 7. Graphs ──────────────────────────────────────────────────────────
    TGraphAsymmErrors *gr_exp = makeGraph(h_plot, eb, "g_cs_" + s.tag +
                                          (norm.ok ? "_norm" : "_raw"));
    gr_exp->SetTitle((s.name + " this work" +
                      (norm.ok ? " (normalised)" : " (absolute)")).c_str());
    gr_exp->SetMarkerStyle(20);
    gr_exp->SetMarkerColor(kBlue + 1);
    gr_exp->SetLineColor(kBlue + 1);
    gr_exp->SetLineWidth(2);

    TGraphErrors *gr_ref = makeReferenceGraph(s.ref, "g_" + s.tag + "_ref");
    gr_ref->SetMarkerStyle(21);
    gr_ref->SetMarkerColor(kRed + 1);
    gr_ref->SetLineColor(kRed + 1);
    gr_ref->SetLineWidth(2);

    // ── 8. Plot ────────────────────────────────────────────────────────────
    TCanvas *c = nullptr;
    if (gr_exp->GetN() > 0) {
        double ymin = 1e300, ymax = 0.0;
        auto update = [&](double y) { if (y > 0) { ymin = std::min(ymin, y); ymax = std::max(ymax, y); } };
        for (int i = 0; i < gr_exp->GetN(); ++i) update(gr_exp->GetY()[i]);
        for (int i = 0; i < gr_ref->GetN(); ++i) update(gr_ref->GetY()[i]);
        const double xmin = 0.8 * std::min(eb.edges.front(),
                                           s.ref.E.empty() ? eb.edges.front() : s.ref.E.front());
        const double xmax = 1.2 * eb.edges.back();

        c = new TCanvas(("c_" + s.tag).c_str(), (s.name + " cross section").c_str(), 900, 600);
        c->SetLogx();
        c->SetLogy();
        TH1F *frame = c->DrawFrame(xmin, 0.5 * ymin, xmax, 2.0 * ymax);
        frame->SetTitle((s.name + " fission cross section;E_{n} (MeV);#sigma (barn)").c_str());
        gr_ref->Draw("P same");
        gr_exp->Draw("P same");

        TLegend *leg = new TLegend(0.15, 0.72, 0.50, 0.88);
        leg->AddEntry(gr_exp, norm.ok ? "This work (normalised)" : "This work (absolute)", "lp");
        leg->AddEntry(gr_ref, "Reference", "lp");
        leg->Draw();
        c->Update();
    } else {
        std::cerr << "[plot] no non-empty cross-section bins, nothing to draw\n";
    }

    // ── 9. Save ────────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(s.out_path.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { std::cerr << "Cannot create " << s.out_path << "\n"; return; }
    fout->cd();
    if (c) c->Write();
    gr_ref->Write();
    gr_exp->Write();
    h_raw->Write();
    if (h_norm) h_norm->Write();
    h_kept->Write();
    h_flux->Write();
    TParameter<double>("norm_scale",   norm.scale).Write();
    TParameter<double>("norm_u_scale", norm.u_scale).Write();
    TParameter<int>   ("norm_ok",      norm.ok ? 1 : 0).Write();
    fout->Close();
    delete fout;

    std::cout << "[DONE] " << s.out_path << "\n";
}

}  // namespace xs

// =============================================================================
//  Entry points
// =============================================================================

void gold_xs()
{
    const xs::EnergyBinning eb = xs::makeBinning();

    xs::SampleSetup s;
    s.name      = "Au-197";
    s.tag       = "gold";
    s.sample    = Sample::gold;
    s.tree_name = "events_gold";
    s.eff_path  = xs::kBaseDir + "output/Au-197/output_efficiency_gold.root";
    s.eff_bins  = makeGoldConfig(eb.edges).energy_bins_eff;
    s.n_areal   = 9.2e17;          // atoms/cm^2                        VERIFY (4)
    s.out_path  = xs::kBaseDir + "output/cs_toy_gold.root";

    // Reference Au-197(n,f), given in mb and converted to barn here.
    s.ref.label = "Au-197 reference;E_{n} (MeV);#sigma (barn)";
    s.ref.E     = {46.3, 66.6, 73.9, 94.1, 132.9, 144.6, 173.3};
    const std::vector<double> sig_mb = {0.103, 0.81, 1.20, 2.81, 6.1, 8.1, 10.3};
    const std::vector<double> err_mb = {0.019, 0.12, 0.17, 0.39, 0.9, 1.2, 1.6};
    for (double v : sig_mb) s.ref.sigma.push_back(v * 1e-3);
    for (double v : err_mb) s.ref.u_sigma.push_back(v * 1e-3);

    xs::runCrossSection(s, eb);
}

void uranium_xs()
{
    const xs::EnergyBinning eb = xs::makeBinning();

    xs::SampleSetup s;
    s.name      = "U-238";
    s.tag       = "uranium";
    s.sample    = Sample::uranium;
    s.tree_name = "events_uranium";
    s.eff_path  = xs::kBaseDir + "output/U-238/output_efficiency_uranium.root";
    // Must match the binning used to produce heff_ebin<k> in the file.
    // No makeUraniumConfig() available, so hard-coded here.     VERIFY (5)
    s.eff_bins  = {1.0, 10.0, 100.0, 500.0, 1000.0};
    s.n_areal   = 6.67e17;         // atoms/cm^2                        VERIFY (4)
    s.out_path  = xs::kBaseDir + "output/cs_toy_uranium.root";

    // Reference U-238(n,f). These values are ALREADY in barn
    // (sigma_f ~ 1.6 b at 40-50 MeV); no unit conversion.
    // NOTE: they cover 42-50 MeV, below the 60 MeV threshold of this
    // analysis, so no normalisation will be applied until points in the
    // measured range (>= 60 MeV) are added here.
    s.ref.label   = "U-238 reference;E_{n} (MeV);#sigma (barn)";
    s.ref.E       = {0.50,0.52,0.54,0.57,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.94,0.96,0.98,
        1.00,1.10,1.25,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.60,4.00,
        4.50,4.70,5.00,5.30,5.50,5.80,6.00,6.20,6.50,7.00,7.50,7.75,8.00,8.50,
        9.00,10.00,11.00,11.50,12.00,13.00,14.00,14.50,15.00,16.00,17.00,18.00,
        19.00,20.00,21.00,22.00,23.00,24.00,25.00,26.00,27.00,28.00,29.00,30.00,
        32.00,34.00,36.00,38.00,40.00,42.00,44.00,46.00,48.00,50.00,52.00,54.00,
        56.00,58.00,60.00,64.00,68.00,72.00,76.00,80.00,84.00,88.00,92.00,96.00,
        100.00,104.00,108.00,112.00,116.00,120.00,128.00,136.00,144.00,152.00,
        160.00,168.00,176.00,184.00,192.00,200.00,300.00,400.00,500.00,600.00,
        700.00,800.00,900.00,1000.00};
    s.ref.sigma   = {0.00026980,0.00067119,0.00059598,0.00065116,0.00116531,0.00129202,0.00184360,
        0.00264695,0.00450784,0.00675016,0.01378015,0.01690077,0.01550547,0.01588961,
        0.01410576,0.02895309,0.03320903,0.18619023,0.41835630,0.48206710,0.53553997,
        0.54749355,0.54624349,0.54151562,0.53708920,0.52481337,0.54752327,0.55451845,
        0.55921914,0.55980956,0.54833662,0.55240729,0.54820060,0.56779697,0.61392195,
        0.68694589,0.82159525,0.95004294,0.99851630,0.99569271,1.01920402,1.01472289,
        1.01452311,1.01117324,1.00882950,1.00505332,0.98752794,1.03156514,1.14994058,
        1.19263970,1.24300610,1.32225423,1.32703520,1.32162220,1.35796397,1.40707482,
        1.51901925,1.55512900,1.60213378,1.54300822,1.57505132,1.58626925,1.56289817,
        1.62440482,1.61076480,1.65532032,1.70277087,1.69981759,1.65676440,1.63970023,
        1.66526066,1.65841277,1.65557833,1.68089904,1.64272747,1.61694393,1.64917876,
        1.61998653,1.63154717,1.62447590,1.58763664,1.57865896,1.52995539,1.50445837,
        1.51984211,1.49515595,1.49600313,1.44827525,1.40759496,1.39989479,1.42732347,
        1.38255921,1.42826710,1.35680338,1.36360367,1.33475089,1.32909427,1.29432732,
        1.28907594,1.32032323,1.31550749,1.29939456,1.33658872,1.34154384,1.31700244,
        1.32466100,1.44100000,1.49700000,1.46100000,1.45800000,1.47500000,1.48200000,
        1.46800000,1.46700001};
    s.ref.u_sigma = {0.00002049,0.00010438,0.00002249,0.00002441,0.00003694,0.00002922,0.00004329,
        0.00004879,0.00008038,0.00010401,0.00017266,0.00021365,0.00023027,0.00021990,
        0.00015131,0.00029119,0.00032341,0.00134738,0.00262918,0.00302871,0.00290432,
        0.00301741,0.00311783,0.00300113,0.00344929,0.00290484,0.00298915,0.00307938,
        0.00360237,0.00358805,0.00362927,0.00400621,0.00365471,0.00415490,0.00463600,
        0.00555936,0.00612914,0.00688432,0.00719061,0.00859903,0.00773614,0.00729852,
        0.00729594,0.00815660,0.00806416,0.00962545,0.00864000,0.00748159,0.00674434,
        0.00570974,0.00834980,0.01067722,0.01146856,0.01263919,0.01082061,0.01487632,
        0.01428978,0.02590652,0.01836317,0.02226853,0.01699061,0.02641331,0.02127290,
        0.02271004,0.02358701,0.01967756,0.02992725,0.02426689,0.02619886,0.03118292,
        0.02826156,0.02867969,0.02690437,0.03005004,0.02926416,0.03718286,0.03321005,
        0.03569087,0.04041503,0.03035458,0.03109424,0.02486835,0.02395078,0.02856810,
        0.02707488,0.03238298,0.03572723,0.03841870,0.03488576,0.03670170,0.03425680,
        0.03842244,0.04389951,0.04586037,0.04611334,0.04018981,0.03825088,0.03790471,
        0.03312034,0.03876146,0.03685174,0.04049074,0.03526862,0.05257003,0.03563506,
        0.03218347,0.07185272,0.10454333,0.08123140,0.08067812,0.08530600,0.06620452,
        0.06558076,0.06740687};

    xs::runCrossSection(s, eb);
}

void fission_xs()
{
    gold_xs();
    uranium_xs();
}