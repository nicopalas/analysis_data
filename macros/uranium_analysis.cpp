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
#include "../include/cross_section.h"
#include "../include/fit_anisotropy.h"
#include "TF1.h"
#include "TH2D.h"

void uranium_analysis(){

    const double emin = 1.;
    const double emax = 1000;

    const std::string outdir = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/";

    // ================================================================
    // EFFICIENCY — coarse binning
    // ================================================================
    const std::vector<double> energy_bins_eff = {1, 10, 100, 500, 1000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    AnalysisConfig cfg_eff = makeUraniumConfig(energy_bins_eff, "eff");

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
    std::cout << "Entries: " << tree->GetEntries() << "\n";

    // --- histograms ---
    std::vector<TH1D*> hists_tof_eff(nbins_eff, nullptr);
    for(int i = 0; i < nbins_eff; ++i){
        hists_tof_eff[i] = new TH1D(Form("htof_eff_%d", i), "", 100, -15, 15);
        hists_tof_eff[i]->SetDirectory(0);
    }

    Vec3D counts_roi_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_eff(nbins_eff,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_eff(nbins_eff,
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
        for(int e = 0; e < nbins_eff; ++e){
        double tot = 0.0;
        for(int j = 0; j < nbins_beam; ++j)
            for(int ii = 0; ii < nbins_det; ++ii)
                tot += counts_acc_eff[e][j][ii];
        std::cout << "ACC eff ebin " << e << "  entries=" << tot << "\n";
    }
    f_acc->Close();

    for(int i = 0; i < nbins_eff; ++i)
        std::cout << "Eff ebin " << i
                  << "  entries=" << hists_tof_eff[i]->GetEntries() << "\n";

    // --- fit background ---
    std::vector<BackgroundFit> bfs_eff(nbins_eff);
    for(int i = 0; i < nbins_eff; ++i){
        double Ec = std::sqrt(energy_bins_eff[i] * energy_bins_eff[i+1]);
        EventCuts c = getCuts(cfg_eff.sample, Ec);
        bfs_eff[i] = fitBackground(cfg_eff, hists_tof_eff[i], c.roi_min, c.roi_max, i);
        std::cout << "Eff ebin " << i << "  chi2/ndf=" << bfs_eff[i].chi2ndf << "\n";
        if(!bfs_eff[i].func || !bfs_eff[i].hist_subtracted){
            std::cerr << "Eff ebin " << i << " fit failed\n";
            return;
        }
    }

    std::vector<double> cs_eff(nbins_eff),            u_cs_eff(nbins_eff);
    std::vector<double> cs_upeak_eff(nbins_eff, 0.0), u_cs_upeak_eff(nbins_eff, 0.0);
    for(int i = 0; i < nbins_eff; ++i){
        cs_eff[i]   = bfs_eff[i].counts_subtract_bkg;
        u_cs_eff[i] = bfs_eff[i].u_counts_subtract_bkg;
    }

    

    Vec3D counts_signal_eff(nbins_eff,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_eff(nbins_eff,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    computeSignal(cfg_eff,
                  counts_roi_eff, counts_acc_eff, counts_upeak_eff,
                  cs_eff,         u_cs_eff,
                  cs_upeak_eff,   u_cs_upeak_eff,
                  counts_signal_eff, u_counts_signal_eff);

    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);
    // --- efficiency ---
    std::vector<EfficiencyResult> eff(nbins_eff);
    for(int e = 0; e < nbins_eff; ++e){
        eff[e] = computeEfficiency(
            nbins_det-1,
            nbins_beam,
            nbins_det,
            counts_signal_eff,
            u_counts_signal_eff,
            acceptance,
            e);
        }

    // --- save efficiency ---
    TFile* fout_eff = TFile::Open(
        (outdir + "output_efficiency_uranium.root").c_str(), "RECREATE");
    if(!fout_eff || fout_eff->IsZombie()){
        std::cerr << "Error creating efficiency output file\n";
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
                       outdir + "background_subtraction_uranium_eff.pdf");
    plotEfficiency(eff, nbins_eff, nbins_det, energy_bins_eff,
                   outdir + "efficiency_uranium.pdf");
    plotEfficiencyResolution(eff, nbins_eff, nbins_det, energy_bins_eff,
                    outdir + "efficiency_resolution_uranium.pdf");

    // ================================================================
    // ANISOTROPY — fine logarithmic binning
    // ================================================================
    const int nbins_aniso = 50;
    std::vector<double> energy_bins_aniso = buildLogBins(nbins_aniso, 1.2, 1000.0, 1.5);

    AnalysisConfig cfg_aniso = makeUraniumConfig(energy_bins_aniso, "aniso");

    TFile* fin2 = TFile::Open(cfg_aniso.input_file.c_str());
    if(!fin2 || fin2->IsZombie()){
        std::cerr << "Error opening file for anisotropy\n";
        return;
    }
    TTree* tree2 = (TTree*)fin2->Get(cfg_aniso.tree_name.c_str());

    std::vector<TH1D*> hists_tof_aniso(nbins_aniso, nullptr);
    for(int i = 0; i < nbins_aniso; ++i){
        hists_tof_aniso[i] = new TH1D(Form("htof_aniso_%d", i), "", 100, -15, 15);
        hists_tof_aniso[i]->SetDirectory(0);
    }

    Vec3D counts_roi_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_bkg_aniso(nbins_aniso,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D counts_upeak_aniso(nbins_aniso,
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
        for(int e = 0; e < nbins_aniso; ++e){
        double tot = 0.0;
        for(int j = 0; j < nbins_beam; ++j)
            for(int ii = 0; ii < nbins_det; ++ii)
                tot += counts_acc_aniso[e][j][ii];
        std::cout << "ACC aniso ebin " << e << "  entries=" << tot << "\n";
    }
    f_acc2->Close();

    // --- fit background ---
    std::vector<BackgroundFit> bfs_aniso(nbins_aniso);
    for(int i = 0; i < nbins_aniso; ++i){
        double Ec = std::sqrt(energy_bins_aniso[i] * energy_bins_aniso[i+1]);
        EventCuts c = getCuts(cfg_aniso.sample, Ec);
        bfs_aniso[i] = fitBackground(cfg_aniso, hists_tof_aniso[i],
                                      c.roi_min, c.roi_max, i);
        std::cout << "Aniso ebin " << i << "  chi2/ndf=" << bfs_aniso[i].chi2ndf << "\n";
    }




    std::vector<double> cs_aniso(nbins_aniso),            u_cs_aniso(nbins_aniso);
    std::vector<double> cs_upeak_aniso(nbins_aniso, 0.0), u_cs_upeak_aniso(nbins_aniso, 0.0);
    for(int i = 0; i < nbins_aniso; ++i){
        cs_aniso[i]   = bfs_aniso[i].counts_subtract_bkg;
        u_cs_aniso[i] = bfs_aniso[i].u_counts_subtract_bkg;
    }
    Vec3D counts_signal_aniso(nbins_aniso,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts_signal_aniso(nbins_aniso,
    Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    computeSignal(cfg_aniso,
                  counts_roi_aniso, counts_acc_aniso, counts_upeak_aniso,
                  cs_aniso,         u_cs_aniso,
                  cs_upeak_aniso,   u_cs_upeak_aniso,
                  counts_signal_aniso, u_counts_signal_aniso);
        // ================================================================
    // ANISOTROPY — puntos + ajuste de Legendre
    // ================================================================
    const bool fit_a4 = false;

    std::vector<AnisotropyResult> aniso(nbins_aniso);
    std::vector<LegendreResult>   leg(nbins_aniso);

    for(int e = 0; e < nbins_aniso; ++e){
        double Ec = std::sqrt(energy_bins_aniso[e] * energy_bins_aniso[e+1]);

        int e_eff = findBin(energy_bins_eff, Ec);
        if(e_eff < 0)          e_eff = 0;
        if(e_eff >= nbins_eff) e_eff = nbins_eff - 1;

        // distribución angular por puntos (se mantiene para los plots antiguos)
        aniso[e] = anisotropy(
            nbins_beam, nbins_det,
            counts_signal_aniso, u_counts_signal_aniso,
            acceptance, e,
            eff[e_eff].eps, eff[e_eff].u_eps,
            cfg_aniso);

        // ajuste W(cos) = A0 (1 + a2 P2 + a4 P4) y R = W(0)/W(90) del ajuste
        leg[e] = legendre_fit(
            nbins_beam, nbins_det,
            counts_signal_aniso,
            acceptance, e,
            eff[e_eff].eps, eff[e_eff].u_eps,
            cfg_aniso, fit_a4);

        std::cout << "Aniso ebin " << e << "  E=" << Ec << " MeV";
        if(leg[e].valid)
            std::cout << "  W0/W90(fit)=" << leg[e].anisotropy
                      << " +/- " << leg[e].u_anisotropy
                      << "  a2=" << leg[e].a2 << " +/- " << leg[e].u_a2;
        else
            std::cout << "  fit failed";
        std::cout << "\n";
    }

    // --- save anisotropy ---
    TFile* fout_aniso = TFile::Open(
        (outdir + "output_anisotropy_uranium_complete.root").c_str(), "RECREATE");
    if(!fout_aniso || fout_aniso->IsZombie()){
        std::cerr << "Error creating anisotropy output file\n";
        return;
    }

    // Curva normalizada a W(90) del propio ajuste: (1 + a2 P2 + a4 P4) / (1 + a2 P2(0) + a4 P4(0))
    auto W_over_W90 = [](double* x, double* p){
        double c = x[0];
        double N = 1. + p[0]*legP2(c)  + p[1]*legP4(c);
        double D = 1. + p[0]*legP2(0.) + p[1]*legP4(0.);
        return N / D;
    };

    for(int e = 0; e < nbins_aniso; ++e){
        const LegendreResult& L = leg[e];

        if(e % 10 == 0 && L.valid){
            int n = (int)L.w.size();
            std::vector<double> y(n), ey(n), ex(n, 0.0);
            for(int k = 0; k < n; ++k){
                y[k]  = L.w[k]   / L.W90;   // puntos normalizados con el W(90) ajustado
                ey[k] = L.u_w[k] / L.W90;
            }

            TGraphErrors* gpt = new TGraphErrors(
                n, L.cos_theta.data(), y.data(), ex.data(), ey.data());
            gpt->SetName(Form("anisotropy_ebin%d", e));
            gpt->SetTitle(Form(
                "W(#theta)/W(90) %.2f-%.2f MeV  (#chi^{2}/ndf=%.2f);"
                "cos(#theta_{beam});W(#theta)/W(90)",
                energy_bins_aniso[e], energy_bins_aniso[e+1], L.chi2ndf));
            gpt->SetMinimum(0.4);
            gpt->SetMaximum(3.5);
            gpt->SetMarkerStyle(20);

            TF1* fn = new TF1(Form("fit_ebin%d", e), W_over_W90, 0., 1., 2);
            fn->SetParameters(L.a2, L.a4);
            fn->SetLineColor(kRed+1);
            gpt->GetListOfFunctions()->Add(fn);   // se dibuja y guarda con el grafo

            gpt->Write();
        }
        hists_tof_aniso[e]->Write();
    }

    // --- R = W(0)/W(90) y a2 frente a energía, del ajuste ---
    std::vector<double> xR, yR, exR, eyR, ya2, eya2;
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

    TGraphErrors* g_a2 = new TGraphErrors(
        (int)xR.size(), xR.data(), ya2.data(), exR.data(), eya2.data());
    g_a2->SetName("a2_fit");
    g_a2->SetTitle("a_{2};E_{n} (MeV);a_{2}");
    g_a2->SetMarkerStyle(20);
    g_a2->Write();

    // --- plots ---

    plotAnisotropy(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                   outdir + "anisotropy_uranium_complete.pdf");

    std::vector<ExforSource> sources = {
        {"/Users/nico/Downloads/13709003.csv"},
        {"/Users/nico/Downloads/14660003 (1).csv"},
        {"/Users/nico/Downloads/41756002 (2).csv"}
    };
    plotAnisotropyFit(leg, energy_bins_aniso,
                      outdir + "anisotropy_fit_uranium.pdf", 5);   // 1 de cada 10 bins

    TGraphErrors* g_R = plotAnisotropyRatioFit(leg, energy_bins_aniso,
                      outdir + "anisotropy_ratio_fit_uranium", "^{238}U(n,f)");

    if(g_R){
        plotAnisoVsExfor(g_R, sources, outdir + "aniso_vs_exfor.pdf");
        plotPullsVsExfor(g_R, sources, outdir + "pulls_vs_exfor_grid.pdf");
        plotPullsOverlay(g_R, sources, outdir + "pulls_vs_exfor_overlay.pdf");
        plotAnisoVsExforIndividual(g_R, sources, outdir + "aniso_vs_exfor_individual.pdf");
    }

    fout_aniso->cd();
    fout_aniso->Close();


    // ================================================================
// CROSS SECTION ABSOLUTE — fine logarithmic binning
// ================================================================
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
        nbins_beam, nbins_det,
        counts_signal_aniso,      u_counts_signal_aniso,
        eff[e_eff].eps,           eff[e_eff].u_eps,
        acceptance,
        e,
        E_low_cs, E_high_cs,
        cfg_aniso.atoms,
        cfg_aniso);

}

// ── signal counts diagnostic ──────────────────────────────────────────
std::cout << "\n=== Signal counts per energy bin ===\n";
for (int e = 0; e < nbins_cs; ++e) {
    double Ec = std::sqrt(energy_bins_cs[e] * energy_bins_cs[e+1]);
    double total_signal = 0.;
    for (int j = 0; j < nbins_beam; ++j)
        for (int ii = 0; ii < nbins_det; ++ii)
            total_signal += counts_signal_aniso[e][j][ii];
}
// --- save cross section ---
TFile* fout_cs = TFile::Open(
    (outdir + "output_cross_section_uranium.root").c_str(), "RECREATE");
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
g_cs->SetName("g_sigma_abs_U238");
g_cs->SetTitle("#sigma_{abs}(U-238);E_{n} (MeV);#sigma (barn)");
g_cs->SetMarkerStyle(20);
g_cs->SetMarkerColor(kRed+1);
g_cs->SetLineColor(kRed+1);
g_cs->SetLineWidth(2);
g_cs->Write();

// ── normalised cross section ──────────────────────────────────────────
TH1D* h_cs_raw = new TH1D("h_cs_raw", "", nbins_cs, energy_bins_cs.data());
for (int e = 0; e < nbins_cs; ++e){
    h_cs_raw->SetBinContent(e+1, cs_abs[e].sigma);
    h_cs_raw->SetBinError  (e+1, cs_abs[e].u_sigma);
}
TH1D* h_cs_norm = (TH1D*)h_cs_raw->Clone("h_cs_norm");
h_cs_norm->Reset();
normalised_xs(h_cs_raw, 8.0, 10.0, 2.006, h_cs_norm);

TCanvas* c_cs = new TCanvas("c_cs", "Cross section U-238", 900, 600);
c_cs->SetLogx();
h_cs_norm->SetLineColor(kBlue+1);
h_cs_norm->SetMarkerColor(kBlue+1);
h_cs_norm->SetMarkerStyle(20);
h_cs_norm->GetXaxis()->SetTitle("E_{n} (MeV)");
h_cs_norm->GetYaxis()->SetTitle("#sigma (barn)");
h_cs_norm->SetTitle("#sigma_{norm}(U-238)");
h_cs_norm->Draw("E");

fout_cs->cd();
g_cs->Write();
h_cs_norm->Write();
c_cs->Write();
fout_cs->Close();

for(int i = 0; i < nbins_eff;   ++i) delete hists_tof_eff[i];
for(int i = 0; i < nbins_aniso; ++i) delete hists_tof_aniso[i];
}