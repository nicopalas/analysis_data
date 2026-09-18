#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/histograms.h"
#include "../include/run_stability.h"
#include "../include/background_subtraction.h"
#include "../include/signal_selection.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"

// ===========================================================================
//  Run-by-run stability of the Au analysis.
//
//  Differences with respect to the uranium version
//  -----------------------------------------------
//  * The background model has TWO components: a flat continuum and the
//    uranium contamination peak (Crystal Ball). Both must be subtracted, and
//    both are passed to computeSignal. Dropping the U peak would count the
//    contamination as signal.
//  * The ratio (U peak / flat) is an extra stability observable, and a good
//    one: the contamination level should not depend on the run. If it does,
//    something changed in the target or in the selection.
//
//  What to trust
//  -------------
//  Every observable here is a ratio, so none of them depends on the
//  integrated beam intensity. The raw yield per group is NOT tested, because
//  the beam normalisation is not available in the event tree (summing
//  PulseIntensity over events is weighted by the event rate, which is what
//  you are trying to measure). With a pulse-level tree or a beam monitor the
//  yield becomes a real test too, and probably the most sensitive one.
// ===========================================================================
void uranium(int group_size = 30)
{
    const double emin = 100;
    const double emax = 2000;

    const std::string outdir =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/Au-197/";
    // NOTE: ROOT does not create directories. mkdir -p this path first or
    // every Print() below fails with "Cannot open file".

    const std::vector<double> ebins = {100, 300, 600, 900, 2000};
    const int nbins_e = (int)ebins.size() - 1;

    AnalysisConfig cfg0 = makeGoldConfig(ebins, "stab");

    // --- base run selection --------------------------------------------------
    RunFilter base;
    base.run_min  = 118558;
    base.run_max  = 118795;
    base.excluded = {118771, 118789};

    // --- open data -----------------------------------------------------------
    TFile* fin = TFile::Open(cfg0.input_file.c_str());
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening " << cfg0.input_file << "\n";
        return;
    }
    TTree* tree = (TTree*)fin->Get(cfg0.tree_name.c_str());
    if (!tree) { std::cerr << "Tree not found\n"; fin->Close(); return; }

    // --- run list and grouping -----------------------------------------------
    const std::vector<int> runs = getRunList(tree, base);
    std::cout << "Runs available: " << runs.size() << "  ["
              << (runs.empty() ? 0 : runs.front()) << " - "
              << (runs.empty() ? 0 : runs.back())  << "]\n";
    if ((int)runs.size() < 2 * group_size) {
        std::cerr << "Not enough runs for a meaningful comparison\n";
        fin->Close();
        return;
    }

    std::vector<RunFilter> groups = makeRunGroups(runs, base, group_size);
    std::cout << "Groups of " << group_size << " runs: " << groups.size() << "\n";

    // a short trailing group carries a huge error bar and drags the weighted
    // mean around; merge it into the previous one
    if (groups.size() > 1) {
        const int last = (int)runs.size() % group_size;
        if (last != 0 && last < group_size) {
            const size_t prev = groups.size() - 2;
            for (int r : groups.back().allowed) groups[prev].allowed.insert(r);
            groups.pop_back();
            std::cout << "Merged the short trailing group -> "
                      << groups.size() << " groups\n";
        }
    }
    const int ng = (int)groups.size();

    // --- global beam-spot offsets --------------------------------------------
    double off_x0, off_x1, off_y0, off_y1;
    computeGlobalOffsets(tree, base, off_x0, off_x1, off_y0, off_y1);
    std::cout << "Global offsets: x0=" << off_x0 << " x1=" << off_x1
              << " y0=" << off_y0 << " y1=" << off_y1 << "\n";

    // --- acceptance ----------------------------------------------------------
    Vec2D dOmega_fine, acceptance;
    if (!loadAcceptanceCSV(
            "/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv",
            dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV\n";
        fin->Close();
        return;
    }
    acceptance = rebin(dOmega_fine);

    // =======================================================================
    //  Pass 0: full dataset -> efficiency, applied identically to every group
    // =======================================================================
    std::vector<TH1D*> h_all(nbins_e, nullptr);
    for (int e = 0; e < nbins_e; ++e) {
        h_all[e] = new TH1D(Form("htof_all_%d", e), "", 100, -15, 15);
        h_all[e]->SetDirectory(0);
    }

    Vec3D roi_all  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D bkg_all  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D upeak_all(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D acc_all  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    GroupDiag diag_all;
    fillHistogramsRuns(tree, cfg0, base, h_all,
                       roi_all, bkg_all, upeak_all,
                       emin, emax, off_x0, off_x1, off_y0, off_y1, diag_all);

    TFile* f_acc = TFile::Open(cfg0.acc_file.c_str());
    TTree* t_acc = f_acc ? (TTree*)f_acc->Get(cfg0.acc_tree_name.c_str()) : nullptr;
    if (!t_acc) {
        std::cerr << "Accidentals tree not found\n";
        fin->Close();
        return;
    }
    // NOTE: fillAccidentalShape still needs a RunFilter argument, exactly as
    // fillHistogramsRuns has. Until then each group gets the global accidental
    // shape rescaled by its exposure (see below), which hides any run
    // dependence of the accidental rate.
    fillAccidentalShape(t_acc, cfg0, acc_all,
                        off_x0, off_x1, off_y0, off_y1, emin, emax);

    std::vector<BackgroundFit> bf_all(nbins_e);
    std::vector<double> cs_all   (nbins_e, 0.0), u_cs_all   (nbins_e, 0.0);
    std::vector<double> cs_up_all(nbins_e, 0.0), u_cs_up_all(nbins_e, 0.0);

    for (int e = 0; e < nbins_e; ++e) {
        const double    Ec = std::sqrt(ebins[e] * ebins[e + 1]);
        const EventCuts c  = getCuts(cfg0.sample, Ec);

        bf_all[e] = fitBackground(cfg0, h_all[e], c.roi_min, c.roi_max, e);

        // gold: flat continuum AND uranium contamination peak
        cs_all[e]      = bf_all[e].counts_subtract_bkg;
        u_cs_all[e]    = bf_all[e].u_counts_subtract_bkg;
        cs_up_all[e]   = bf_all[e].counts_subtract_upeak;
        u_cs_up_all[e] = bf_all[e].u_counts_subtract_upeak;

        std::cout << "full  e" << e
                  << "  chi2/ndf=" << bf_all[e].chi2ndf
                  << "  flat="     << cs_all[e]
                  << "  Upeak="    << cs_up_all[e] << "\n";
    }

    Vec3D sig_all  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_sig_all(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    computeSignal(cfg0, roi_all, acc_all, upeak_all,
                  cs_all, u_cs_all, cs_up_all, u_cs_up_all,
                  sig_all, u_sig_all);

    std::vector<EfficiencyResult> eff(nbins_e);
    for (int e = 0; e < nbins_e; ++e)
        eff[e] = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                   sig_all, u_sig_all, acceptance, e);

    // =======================================================================
    //  Pass 1: one full chain per group
    // =======================================================================
    std::vector<GroupDiag> diags(ng);

    // observables, indexed [energy bin][group]
    auto mat = [&]{ return std::vector<std::vector<double>>(
                        nbins_e, std::vector<double>(ng, 0.0)); };

    auto aniso_v = mat(), aniso_e = mat();   // W(0)/W(90)
    auto sb_v    = mat(), sb_e    = mat();   // signal / total background
    auto up_v    = mat(), up_e    = mat();   // U peak / flat continuum
    auto mu_v    = mat(), mu_e    = mat();   // ToF mean in the ROI
    auto sig_w_v = mat(), sig_w_e = mat();   // ToF RMS in the ROI
    auto chi2_v  = mat();                    // fit quality

    std::vector<std::vector<TH1D*>> h_group(ng, std::vector<TH1D*>(nbins_e, nullptr));

    for (int g = 0; g < ng; ++g) {

        AnalysisConfig cfg = makeGoldConfig(ebins, Form("stab_g%d", g));

        for (int e = 0; e < nbins_e; ++e) {
            h_group[g][e] = new TH1D(Form("htof_g%d_e%d", g, e), "", 100, -15, 15);
            h_group[g][e]->SetDirectory(0);
        }

        Vec3D roi  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D bkg  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D upeak(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D acc  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        fillHistogramsRuns(tree, cfg, groups[g], h_group[g],
                           roi, bkg, upeak,
                           emin, emax, off_x0, off_x1, off_y0, off_y1, diags[g]);

        fillAccidentalShape(t_acc, cfg, acc,
                            off_x0, off_x1, off_y0, off_y1, emin, emax);
        if (diag_all.n_events > 0) {
            const double frac = diags[g].n_events / diag_all.n_events;
            for (int e = 0; e < nbins_e; ++e)
                for (int j = 0; j < nbins_beam; ++j)
                    for (int ii = 0; ii < nbins_det; ++ii)
                        acc[e][j][ii] *= frac;
        }

        std::vector<double> cs   (nbins_e, 0.0), u_cs   (nbins_e, 0.0);
        std::vector<double> cs_up(nbins_e, 0.0), u_cs_up(nbins_e, 0.0);

        for (int e = 0; e < nbins_e; ++e) {
            const double    Ec = std::sqrt(ebins[e] * ebins[e + 1]);
            const EventCuts c  = getCuts(cfg.sample, Ec);

            BackgroundFit bf = fitBackground(cfg, h_group[g][e],
                                             c.roi_min, c.roi_max, e);

            cs[e]        = bf.counts_subtract_bkg;
            u_cs[e]      = bf.u_counts_subtract_bkg;
            cs_up[e]     = bf.counts_subtract_upeak;
            u_cs_up[e]   = bf.u_counts_subtract_upeak;
            chi2_v[e][g] = bf.chi2ndf;

            // --- ToF shape inside the ROI: detector-timing diagnostic -------
            const int b1 = h_group[g][e]->FindBin(c.roi_min);
            const int b2 = h_group[g][e]->FindBin(c.roi_max);
            h_group[g][e]->GetXaxis()->SetRange(b1, b2);
            const double n_roi = h_group[g][e]->Integral(b1, b2);
            mu_v[e][g]    = h_group[g][e]->GetMean();
            mu_e[e][g]    = h_group[g][e]->GetMeanError();
            sig_w_v[e][g] = h_group[g][e]->GetRMS();
            sig_w_e[e][g] = h_group[g][e]->GetRMSError();
            h_group[g][e]->GetXaxis()->SetRange(0, 0);

            // --- signal to TOTAL background (flat + U peak) -----------------
            const double bkg_tot   = cs[e] + cs_up[e];
            const double u_bkg_tot = std::sqrt(u_cs[e]   * u_cs[e] +
                                               u_cs_up[e] * u_cs_up[e]);
            const double s = n_roi - bkg_tot;
            if (bkg_tot > 0.0 && n_roi > 0.0 && std::fabs(s) > 0.0) {
                sb_v[e][g] = s / bkg_tot;
                sb_e[e][g] = std::fabs(sb_v[e][g]) *
                             std::sqrt(n_roi / (s * s) +
                                       (u_bkg_tot * u_bkg_tot) /
                                       (bkg_tot * bkg_tot));
            }

            // --- uranium contamination relative to the continuum ------------
            // A ratio, so it does not depend on the exposure of the group.
            if (cs[e] > 0.0 && cs_up[e] > 0.0) {
                up_v[e][g] = cs_up[e] / cs[e];
                up_e[e][g] = up_v[e][g] *
                             std::sqrt((u_cs_up[e] * u_cs_up[e]) /
                                       (cs_up[e] * cs_up[e]) +
                                       (u_cs[e] * u_cs[e]) / (cs[e] * cs[e]));
            }
        }

        Vec3D sig  (nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        Vec3D u_sig(nbins_e, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
        computeSignal(cfg, roi, acc, upeak,
                      cs, u_cs, cs_up, u_cs_up,
                      sig, u_sig);

        for (int e = 0; e < nbins_e; ++e) {
            AnisotropyResult a = anisotropy(nbins_beam, nbins_det,
                                            sig, u_sig, acceptance, e,
                                            eff[e].eps, eff[e].u_eps, cfg);
            aniso_v[e][g] = a.w  [nbins_beam - 1];
            aniso_e[e][g] = a.u_w[nbins_beam - 1];
        }

        std::cout << "group " << g
                  << "  runs " << diags[g].run_lo << "-" << diags[g].run_hi
                  << "  (" << diags[g].n_runs << " runs)"
                  << "  events=" << diags[g].n_events << "\n";
    }
    fin->Close();
    if (f_acc) f_acc->Close();

    // =======================================================================
    //  Comparison across groups
    // =======================================================================
    std::cout << "\n==================== STABILITY TESTS ====================\n";

    std::vector<StabilityTest> all_tests;

    for (int e = 0; e < nbins_e; ++e) {
        const std::string tag = Form("E %.0f-%.0f MeV", ebins[e], ebins[e + 1]);

        StabilityTest t1 = testConstant("W(0)/W(90)   " + tag, aniso_v[e], aniso_e[e]);
        StabilityTest t2 = testConstant("S/B          " + tag, sb_v[e],    sb_e[e]);
        StabilityTest t3 = testConstant("Upeak/flat   " + tag, up_v[e],    up_e[e]);
        StabilityTest t4 = testConstant("ToF mean     " + tag, mu_v[e],    mu_e[e]);
        StabilityTest t5 = testConstant("ToF RMS      " + tag, sig_w_v[e], sig_w_e[e]);

        printStabilityTest(t1, diags);
        printStabilityTest(t2, diags);
        printStabilityTest(t3, diags);
        printStabilityTest(t4, diags);
        printStabilityTest(t5, diags);

        plotStability(t1, aniso_v[e], aniso_e[e], diags,
                      outdir + Form("stab_aniso_gold_e%d.pdf", e));
        plotStability(t3, up_v[e], up_e[e], diags,
                      outdir + Form("stab_upeak_gold_e%d.pdf", e));
        plotStability(t4, mu_v[e], mu_e[e], diags,
                      outdir + Form("stab_tofmean_gold_e%d.pdf", e));

        all_tests.push_back(t1);
        all_tests.push_back(t2);
        all_tests.push_back(t3);
        all_tests.push_back(t4);
        all_tests.push_back(t5);
    }

    // --- background fit quality per group ------------------------------------
    std::cout << "\n--- Background fit chi2/ndf ---\n";
    for (int e = 0; e < nbins_e; ++e) {
        std::cout << "  e" << e << ":";
        for (int g = 0; g < ng; ++g)
            std::cout << "  " << std::setprecision(3) << chi2_v[e][g];
        std::cout << "\n";
    }

    // --- ToF shape comparison, each group against the full dataset -----------
    std::cout << "\n--- ToF shape vs full dataset (Chi2Test, UU NORM) ---\n";
    for (int e = 0; e < nbins_e; ++e) {
        for (int g = 0; g < ng; ++g) {
            if (h_group[g][e]->Integral() < 50) continue;
            const double p = h_group[g][e]->Chi2Test(h_all[e], "UU NORM");
            std::cout << "  e" << e << " g" << g << "  p = " << p;
            if (p < 0.01) std::cout << "   <== shape differs";
            std::cout << "\n";
        }
    }

    // --- beam spot drift ------------------------------------------------------
    std::cout << "\n--- Beam spot means per group (diagnostic, no errors) ---\n";
    for (int g = 0; g < ng; ++g)
        std::cout << "  g" << g
                  << "  x0=" << diags[g].mean_x0 << "  x1=" << diags[g].mean_x1
                  << "  y0=" << diags[g].mean_y0 << "  y1=" << diags[g].mean_y1
                  << "\n";

    // --- summary --------------------------------------------------------------
    std::cout << "\n==================== SUMMARY ====================\n";
    std::map<int, int> flags;
    for (const auto& t : all_tests)
        for (int idx : t.outliers) flags[idx]++;

    if (flags.empty()) {
        std::cout << "No group flagged at |pull| > 3.\n";
    } else {
        for (const auto& kv : flags)
            std::cout << "  group " << kv.first
                      << " [" << diags[kv.first].run_lo << "-"
                      << diags[kv.first].run_hi << "] fails "
                      << kv.second << " / " << all_tests.size() << " tests\n";
        std::cout << "\nOne failure out of many tests is expected by chance."
                     " Veto runs only when a group fails several unrelated"
                     " observables.\n";
    }

    // --- error inflation factors ---------------------------------------------
    // If a group is not to blame, the spread itself is the systematic: scale
    // the final uncertainty of each observable by sqrt(chi2/ndf) when > 1.
    std::cout << "\n--- Suggested error scale factors sqrt(chi2/ndf) ---\n";
    for (const auto& t : all_tests) {
        if (t.ndf <= 0) continue;
        const double r = t.chi2 / t.ndf;
        if (r > 1.0)
            std::cout << "  " << t.name << "  x" << std::sqrt(r) << "\n";
    }

    // --- cleanup --------------------------------------------------------------
    for (int e = 0; e < nbins_e; ++e) delete h_all[e];
    for (int g = 0; g < ng; ++g)
        for (int e = 0; e < nbins_e; ++e) delete h_group[g][e];
}