#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/histograms.h"
#include "../include/background_subtraction.h"
#include "../include/signal_selection.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"
#include "../include/fit_anisotropy.h"
#include "../include/cross_section.h"



void gold_analysis()
{
    const std::string outdir = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/Au-197/";

    // ================================================================
    // EFFICIENCY — coarse binning
    // ================================================================
    const double emin = 40.;
    const double emax = 1000.;
    const std::vector<double> energy_bins_eff = {40, 300, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    AnalysisConfig cfg_eff = makeGoldConfig(energy_bins_eff, "eff");

    // --- open data ---
    TFile* fin = TFile::Open(cfg_eff.input_file.c_str());
    if(!fin || fin->IsZombie()){
        std::cerr << "Error opening " << cfg_eff.input_file << "\n";
        return;
    }
    TTree* tree = (TTree*)fin->Get(cfg_eff.tree_name.c_str());
    if(!tree){
        std::cerr << "Tree not found\n";
        fin->Close();
        return;
    }
    std::cout << "Gold entries: " << tree->GetEntries() << "\n";

    // --- histograms ---
    std::vector<TH1D*> hists_tof_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        hists_tof_eff[i] = new TH1D(
            Form("htof_gold_eff_%d", i), "", 100, -15, 15);
        hists_tof_eff[i]->SetDirectory(0);
    }

    Vec3D counts_roi_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_eff(nbins_eff,                                  // NEW
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_acc_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        double mx0, mx1, my0, my1;
    fillHistograms(tree, cfg_eff, hists_tof_eff,
                   counts_roi_eff, counts_bkg_eff, counts_upeak_eff,
                   emin, emax, mx0, mx1, my0, my1);
    fin->Close();

    TFile* f_acc = TFile::Open(cfg_eff.acc_file.c_str());
    if(!f_acc || f_acc->IsZombie()){
        std::cerr << "Error opening " << cfg_eff.acc_file << "\n";
        return;
    }
    TTree* t_acc = (TTree*)f_acc->Get(cfg_eff.acc_tree_name.c_str());
    if(!t_acc){
        std::cerr << "Accidentals tree " << cfg_eff.acc_tree_name << " not found\n";
        f_acc->Close();
        return;
    }
    fillAccidentalShape(t_acc, cfg_eff, counts_acc_eff,
                        mx0, mx1, my0, my1, emin, emax);
    f_acc->Close();

    for(int i = 0; i < nbins_eff; ++i)
        std::cout << "Gold eff ebin " << i
                  << "  entries=" << hists_tof_eff[i]->GetEntries() << "\n";

    // --- fit background ---
    std::vector<BackgroundFit> bfs_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        double Ec = std::sqrt(energy_bins_eff[i] * energy_bins_eff[i+1]);
        EventCuts c = getCuts(cfg_eff.sample, Ec);
        bfs_eff[i] = fitBackground(cfg_eff, hists_tof_eff[i],
                                   c.roi_min, c.roi_max, i);
        std::cout << "Gold eff ebin " << i
                  << "  chi2/ndf=" << bfs_eff[i].chi2ndf << "\n";
        if(!bfs_eff[i].func || !bfs_eff[i].hist_subtracted){
            std::cerr << "Gold eff ebin " << i << " fit failed\n";
            return;
        }
    }

    // declara aquí, antes de computeSignal
    Vec3D counts_signal_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    std::vector<double> cs_bkg_eff(nbins_eff),    u_cs_bkg_eff(nbins_eff);
    std::vector<double> cs_upeak_eff(nbins_eff),   u_cs_upeak_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        cs_bkg_eff[i]     = bfs_eff[i].counts_subtract_bkg;
        u_cs_bkg_eff[i]   = bfs_eff[i].u_counts_subtract_bkg;
        cs_upeak_eff[i]   = bfs_eff[i].counts_subtract_upeak;
        u_cs_upeak_eff[i] = bfs_eff[i].u_counts_subtract_upeak;
    }

    computeSignal(cfg_eff,
                  counts_roi_eff, counts_acc_eff, counts_upeak_eff,
                  cs_bkg_eff,   u_cs_bkg_eff,
                  cs_upeak_eff, u_cs_upeak_eff,
                  counts_signal_eff, u_counts_signal_eff);

    // --- acceptance ---
    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);

    // --- efficiency ---
    std::vector<EfficiencyResult> eff(nbins_eff);
    for(int e = 0; e < nbins_eff; ++e)
        eff[e] = computeEfficiency(
            nbins_det -1,
            nbins_beam,
            nbins_det,
            counts_signal_eff,
            u_counts_signal_eff,
            acceptance,
            e);

    // --- save efficiency ---
    TFile* fout_eff = TFile::Open(
        (outdir + "output_efficiency_gold.root").c_str(), "RECREATE");
    if(!fout_eff || fout_eff->IsZombie()){
        std::cerr << "Error creating gold efficiency file\n";
        return;
    }
    for(int e = 0; e < nbins_eff; ++e){
        TH1D* heff = new TH1D(Form("heff_ebin%d", e), "", nbins_det, 0, 1);
        for(int i = 0; i < nbins_det; ++i){
            heff->SetBinContent(i+1, eff[e].eps[i]);
            heff->SetBinError(i+1,   eff[e].u_eps[i]);
        }
        heff->Write();
        hists_tof_eff[e]->Write();
        bfs_eff[e].hist_subtracted->Write();
    }
    fout_eff->Close();

    // --- plot efficiency ---
    std::vector<TF1*>  fits_eff(nbins_eff, nullptr);
    std::vector<TH1D*> hists_sub_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        fits_eff[i]      = bfs_eff[i].func;
        hists_sub_eff[i] = bfs_eff[i].hist_subtracted;
    }
    plotBackgroundFits(hists_tof_eff, hists_sub_eff, fits_eff, nbins_eff,
                       outdir + "background_subtraction_gold_eff.pdf");
    plotEfficiency(eff, nbins_eff, nbins_det, energy_bins_eff,
               outdir + "efficiency_gold.pdf");
    plotEfficiencyResolution(eff, nbins_eff, nbins_det, energy_bins_eff,
                         outdir + "efficiency_resolution_gold.pdf");

    // ================================================================
    // ANISOTROPY — 5 bins log entre 40 y 300 MeV
    // ================================================================
    const int nbins_aniso = 9;
    std::vector<double> energy_bins_aniso = buildLogBins(nbins_aniso, 40, 1000, 0.8);

    AnalysisConfig cfg_aniso = makeGoldConfig(energy_bins_aniso, "aniso");

    // --- open data ---
    TFile* fin2 = TFile::Open(cfg_aniso.input_file.c_str());
    if(!fin2 || fin2->IsZombie()){
        std::cerr << "Error opening file for gold anisotropy\n";
        return;
    }
    TTree* tree2 = (TTree*)fin2->Get(cfg_aniso.tree_name.c_str());
    if(!tree2){ fin2->Close(); return; }

    // --- histograms ---
    std::vector<TH1D*> hists_tof_aniso(nbins_aniso, nullptr);
    for(int i = 0; i < nbins_aniso; ++i){
        hists_tof_aniso[i] = new TH1D(
            Form("htof_gold_aniso_%d", i), "", 100, -20, 20);
        hists_tof_aniso[i]->SetDirectory(0);
    }

    Vec3D counts_roi_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_aniso(nbins_aniso,                              // NEW
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_acc_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

        double mx0_a, mx1_a, my0_a, my1_a;
    fillHistograms(tree2, cfg_aniso, hists_tof_aniso,
                   counts_roi_aniso, counts_bkg_aniso, counts_upeak_aniso,
                   emin, emax, mx0_a, mx1_a, my0_a, my1_a);
    fin2->Close();

    TFile* f_acc2 = TFile::Open(cfg_aniso.acc_file.c_str());
    if(!f_acc2 || f_acc2->IsZombie()){
        std::cerr << "Error opening accidentals for anisotropy\n";
        return;
    }
    TTree* t_acc2 = (TTree*)f_acc2->Get(cfg_aniso.acc_tree_name.c_str());
    if(!t_acc2){
        std::cerr << "Accidentals tree not found (aniso)\n";
        f_acc2->Close();
        return;
    }
    fillAccidentalShape(t_acc2, cfg_aniso, counts_acc_aniso,
                        mx0_a, mx1_a, my0_a, my1_a, emin, emax);
    f_acc2->Close();

    // --- fit background ---
    std::vector<BackgroundFit> bfs_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        double Ec = std::sqrt(energy_bins_aniso[i] * energy_bins_aniso[i+1]);
        EventCuts c = getCuts(cfg_aniso.sample, Ec);
        bfs_aniso[i] = fitBackground(cfg_aniso, hists_tof_aniso[i],
                                     c.roi_min, c.roi_max, i);
        std::cout << "Gold aniso ebin " << i
                  << "  chi2/ndf=" << bfs_aniso[i].chi2ndf << "\n";
        if(!bfs_aniso[i].func || !bfs_aniso[i].hist_subtracted){
            std::cerr << "Gold aniso ebin " << i << " fit failed\n";
            return;
        }
    }

    // --- signal ---
    Vec3D counts_signal_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    // --- signal aniso ---
    std::vector<double> cs_bkg_aniso(nbins_aniso),   u_cs_bkg_aniso(nbins_aniso);
    std::vector<double> cs_upeak_aniso(nbins_aniso),  u_cs_upeak_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        cs_bkg_aniso[i]    = bfs_aniso[i].counts_subtract_bkg;
        u_cs_bkg_aniso[i]  = bfs_aniso[i].u_counts_subtract_bkg;
        cs_upeak_aniso[i]  = bfs_aniso[i].counts_subtract_upeak;
        u_cs_upeak_aniso[i]= bfs_aniso[i].u_counts_subtract_upeak;
    }
        computeSignal(cfg_aniso,
                  counts_roi_aniso, counts_acc_aniso, counts_upeak_aniso,
                  cs_bkg_aniso,     u_cs_bkg_aniso,
                  cs_upeak_aniso,   u_cs_upeak_aniso,
                  counts_signal_aniso, u_counts_signal_aniso);


       // --- anisotropy: puntos + ajuste de Legendre ---
    const bool fit_a4 = false;

    std::vector<AnisotropyResult> aniso(nbins_aniso);
    std::vector<LegendreResult>   leg(nbins_aniso);

    for(int e = 0; e < nbins_aniso; ++e){
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);

        int e_eff = findBin(energy_bins_eff, Ec);
        if(e_eff < 0)          e_eff = 0;
        if(e_eff >= nbins_eff) e_eff = nbins_eff - 1;

        aniso[e] = anisotropy(
            nbins_beam, nbins_det,
            counts_signal_aniso, u_counts_signal_aniso,
            acceptance, e,
            eff[e_eff].eps, eff[e_eff].u_eps,
            cfg_aniso);

        leg[e] = legendre_fit(
            nbins_beam, nbins_det,
            counts_signal_aniso,
            acceptance, e,
            eff[e_eff].eps, eff[e_eff].u_eps,
            cfg_aniso, fit_a4);

        std::cout << "Gold aniso ebin " << e << "  E=" << Ec << " MeV";
        if(leg[e].valid)
            std::cout << "  W0/W90(fit)=" << leg[e].anisotropy
                      << " +/- " << leg[e].u_anisotropy
                      << "  a2=" << leg[e].a2 << " +/- " << leg[e].u_a2
                      << "  chi2/ndf=" << leg[e].chi2ndf;
        else
            std::cout << "  fit failed";
        std::cout << "\n";
    }

    // --- save anisotropy ---
    TFile* fout_aniso = TFile::Open(
        (outdir + "output_anisotropy_gold.root").c_str(), "RECREATE");
    if(!fout_aniso || fout_aniso->IsZombie()){
        std::cerr << "Error creating gold anisotropy file\n";
        return;
    }

    // W(cos)/W(90) del ajuste
    auto W_over_W90 = [](double* x, double* p){
        double c = x[0];
        double N = 1. + p[0]*legP2(c)  + p[1]*legP4(c);
        double D = 1. + p[0]*legP2(0.) + p[1]*legP4(0.);
        return N / D;
    };

    for(int e = 0; e < nbins_aniso; ++e){
        const LegendreResult& L = leg[e];
        hists_tof_aniso[e]->Write();
        if(!L.valid) continue;

        int n = (int)L.w.size();
        std::vector<double> y(n), ey(n), ex(n, 0.0);
        for(int k = 0; k < n; ++k){
            y[k]  = L.w[k]   / L.W90;
            ey[k] = L.u_w[k] / L.W90;
        }

        TGraphErrors* g = new TGraphErrors(
            n, L.cos_theta.data(), y.data(), ex.data(), ey.data());
        g->SetName(Form("anisotropy_gold_ebin%d", e));
        g->SetTitle(Form("W(#theta)/W(90) %.1f-%.1f MeV  (#chi^{2}/ndf=%.2f);"
                         "cos(#theta_{beam});W(#theta)/W(90)",
                         energy_bins_aniso[e], energy_bins_aniso[e+1], L.chi2ndf));
        g->SetMarkerStyle(20);

        TF1* fn = new TF1(Form("fit_gold_ebin%d", e), W_over_W90, 0., 1., 2);
        fn->SetParameters(L.a2, L.a4);
        fn->SetLineColor(kRed+1);
        g->GetListOfFunctions()->Add(fn);

        g->Write();
    }

    // --- R = W(0)/W(90) y a2 frente a energía ---
    std::vector<double> xR, exR, yR, eyR, ya2, eya2;
    for(int e = 0; e < nbins_aniso; ++e){
        if(!leg[e].valid) continue;
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);
        xR.push_back(Ec);
        exR.push_back(0.0);
        yR.push_back(leg[e].anisotropy);
        eyR.push_back(leg[e].u_anisotropy);
        ya2.push_back(leg[e].a2);
        eya2.push_back(leg[e].u_a2);
    }

    TGraphErrors* g_R = new TGraphErrors(
        (int)xR.size(), xR.data(), yR.data(), exR.data(), eyR.data());
    g_R->SetName("anisotropy_ratio_fit_gold");
    g_R->SetTitle("Au-197: W(0^{#circ})/W(90^{#circ}) from Legendre fit;"
                  "E_{n} (MeV);W(0)/W(90)");
    g_R->SetMarkerStyle(20);
    g_R->SetMarkerColor(kBlue+1);
    g_R->SetLineColor(kBlue+1);
    g_R->Write();

    TGraphErrors* g_a2 = new TGraphErrors(
        (int)xR.size(), xR.data(), ya2.data(), exR.data(), eya2.data());
    g_a2->SetName("a2_fit_gold");
    g_a2->SetTitle("Au-197: a_{2};E_{n} (MeV);a_{2}");
    g_a2->SetMarkerStyle(20);
    g_a2->Write();

    fout_aniso->Close();

    // --- plot anisotropy ---
    plotAnisotropy(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_gold.pdf");

    plotAnisotropyFit(leg, energy_bins_aniso,
                      outdir + "anisotropy_fit_gold.pdf");          // los 4 bins

    plotAnisotropyRatioFit(leg, energy_bins_aniso,
                      outdir + "anisotropy_ratio_fit_gold", "^{197}Au(n,f)");

const int nbins_cs = nbins_aniso;
const std::vector<double>& energy_bins_cs = energy_bins_aniso;
// E_low y E_high en MeV para readFlux
std::vector<double> E_low_cs(nbins_cs), E_high_cs(nbins_cs);
for (int e = 0; e < nbins_cs; ++e) {
    E_low_cs[e]  = energy_bins_cs[e];
    E_high_cs[e] = energy_bins_cs[e+1];
}
    // reutilizamos counts_signal_aniso, u_counts_signal_aniso y eff ya calculados
std::vector<CrossSection> cs_abs(nbins_cs);
for (int e = 0; e < nbins_cs; ++e) {

    double Ec   = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    int e_eff   = findBin(energy_bins_eff, Ec);
    if (e_eff < 0)          e_eff = 0;
    if (e_eff >= nbins_eff) e_eff = nbins_eff - 1;

    cs_abs[e] = cross_section_absolute(
        cfg_aniso.flux_file,     
        cfg_aniso.flux_hist,  
        nbins_beam,
        nbins_det,
        counts_signal_aniso,
        u_counts_signal_aniso,
        eff[e_eff].eps,
        eff[e_eff].u_eps,
        acceptance,
        e,
        E_low_cs,
        E_high_cs,
        cfg_aniso.atoms,                
        cfg_aniso);

    std::cout << "CS_abs ebin " << e
              << "  E=" << Ec << " MeV"
              << "  sigma=" << cs_abs[e].sigma
              << " +/- "    << cs_abs[e].u_sigma << " barn\n";
}

// ── normalización a la referencia de Au-197 ──────────────────────────────
const std::vector<RefPoint> ref_gold = {   // mb → barn
    { 46.3,  0.103e-3, 0.019e-3},
    { 66.6,  0.81e-3,  0.12e-3 },
    { 73.9,  1.20e-3,  0.17e-3 },
    { 94.1,  2.81e-3,  0.39e-3 },
    {132.9,  6.1e-3,   0.9e-3  },
    {144.6,  8.1e-3,   1.2e-3  },
    {173.3, 10.3e-3,   1.6e-3  }
};

NormalisedCS cs_norm = normalise_to_reference(
    cs_abs, energy_bins_cs, ref_gold, 173.2);

if (!cs_norm.ok)
    std::cerr << "[WARN] normalización fallida, sigo con la sección eficaz absoluta\n";

// --- save cross section ---
TFile* fout_cs = TFile::Open(
    (outdir + "output_cross_section_gold.root").c_str(), "RECREATE");
if (!fout_cs || fout_cs->IsZombie()) {
    std::cerr << "Error creating cross section output file\n";
    return;
}

// TGraphErrors con sigma vs energia
std::vector<double> x_cs, y_cs, ex_cs, ey_cs;
for (int e = 0; e < nbins_cs; ++e) {
    if (cs_abs[e].sigma <= 0.0) continue;
    double Ec = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    x_cs.push_back(Ec);
    y_cs.push_back(cs_abs[e].sigma);
    ex_cs.push_back((E_high_cs[e] - E_low_cs[e]) / 2.0);
    ey_cs.push_back(cs_abs[e].u_sigma);
}

TGraphErrors* g_cs = new TGraphErrors(
    (int)x_cs.size(),
    x_cs.data(), y_cs.data(),
    ex_cs.data(), ey_cs.data());
g_cs->SetName("g_sigma_abs_Au197");
g_cs->SetTitle("#sigma_{abs}(Au-197);E_{n} (MeV);#sigma (barn)");
g_cs->SetMarkerStyle(20);
g_cs->SetMarkerColor(kRed+1);
g_cs->SetLineColor(kRed+1);
g_cs->SetLineWidth(2);
g_cs->Write();

if (cs_norm.ok) {
    std::vector<double> x_n, y_n, ex_n, ey_n, ey_n_stat;
    for (int e = 0; e < nbins_cs; ++e) {
        if (cs_norm.sigma[e] <= 0.0) continue;
        double Ec = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
        x_n .push_back(Ec);
        y_n .push_back(cs_norm.sigma[e]);
        ex_n.push_back((E_high_cs[e] - E_low_cs[e]) / 2.0);
        ey_n.push_back(cs_norm.u_sigma[e]);
        ey_n_stat.push_back(cs_norm.u_sigma_stat[e]);

        printf("CS_norm ebin %d  E=%.1f MeV  sigma=%.4e +/- %.4e (tot) "
               "+/- %.4e (stat) barn\n",
               e, Ec, cs_norm.sigma[e], cs_norm.u_sigma[e], cs_norm.u_sigma_stat[e]);
    }

    TGraphErrors* g_cs_norm = new TGraphErrors(
        (int)x_n.size(), x_n.data(), y_n.data(), ex_n.data(), ey_n.data());
    g_cs_norm->SetName("g_sigma_norm_Au197");
    g_cs_norm->SetTitle("#sigma_{norm}(Au-197);E_{n} (MeV);#sigma (barn)");
    g_cs_norm->SetMarkerStyle(20);
    g_cs_norm->SetMarkerColor(kBlue+1);
    g_cs_norm->SetLineColor(kBlue+1);
    g_cs_norm->SetLineWidth(2);
    g_cs_norm->Write();

    // banda con error sólo estadístico, por si quieres separar sistemático
    TGraphErrors* g_cs_norm_stat = new TGraphErrors(
        (int)x_n.size(), x_n.data(), y_n.data(), ex_n.data(), ey_n_stat.data());
    g_cs_norm_stat->SetName("g_sigma_norm_Au197_stat");
    g_cs_norm_stat->SetTitle("#sigma_{norm}(Au-197) stat only;E_{n} (MeV);#sigma (barn)");
    g_cs_norm_stat->Write();

    // referencia
    std::vector<double> xr, yr, exr, eyr;
    for (const auto& p : ref_gold) {
        xr.push_back(p.E); yr.push_back(p.sigma);
        exr.push_back(0.0); eyr.push_back(p.u_sigma);
    }
    TGraphErrors* g_ref = new TGraphErrors(
        (int)xr.size(), xr.data(), yr.data(), exr.data(), eyr.data());
    g_ref->SetName("g_gold_ref");
    g_ref->SetTitle("Au-197 reference;E_{n} (MeV);#sigma (barn)");
    g_ref->SetMarkerStyle(21);
    g_ref->SetMarkerColor(kRed+1);
    g_ref->SetLineColor(kRed+1);
    g_ref->Write();
}

fout_cs->Close();


    for(int i = 0; i < nbins_eff;   ++i) delete hists_tof_eff[i];
    for(int i = 0; i < nbins_aniso; ++i) delete hists_tof_aniso[i];
}