#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/constants.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/types.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/acceptance.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/config.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/anisotropy.h"
#include <map>
#include <utility>

// =====================================================================
//  Estilo "paper de física de partículas" (mismo que en efficiency_toy.C,
//  duplicado aquí porque cada macro de ROOT se ejecuta como una unidad
//  independiente).
// =====================================================================
void SetPhysicsStyleAniso(TStyle *st = gStyle){
    st->SetOptStat(0);
    st->SetOptTitle(0);
    st->SetCanvasColor(kWhite);
    st->SetPadColor(kWhite);
    st->SetFrameFillColor(kWhite);
    st->SetPadTopMargin(0.07);
    st->SetPadBottomMargin(0.14);
    st->SetPadLeftMargin(0.14);
    st->SetPadRightMargin(0.06);
    st->SetTitleFont(42, "XYZ");
    st->SetLabelFont(42, "XYZ");
    st->SetTextFont(42);
    st->SetTitleSize(0.050, "XYZ");
    st->SetTitleOffset(1.10, "X");
    st->SetTitleOffset(1.35, "Y");
    st->SetLabelSize(0.042, "XYZ");
    st->SetLabelOffset(0.010, "XYZ");
    st->SetNdivisions(508, "XYZ");
    st->SetFrameBorderMode(0);
    st->SetFrameLineWidth(2);
    st->SetCanvasBorderMode(0);
    st->SetPadBorderMode(0);
    st->SetPadTickX(1);
    st->SetPadTickY(1);
    st->SetLegendBorderSize(0);
    st->SetLegendFillColor(0);
    st->SetLegendFont(42);
}

void DrawTagAniso(double x, double y, const char* label, double size = 0.045){
    TLatex *lx = new TLatex();
    lx->SetNDC();
    lx->SetTextFont(42);
    lx->SetTextSize(size);
    lx->DrawLatex(x, y, Form("#font[62]{Gold} #font[52]{Preliminary}  %s", label));
}

void anisotropy(){

    SetPhysicsStyleAniso();

    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*)fin->Get("events_gold");
    if (!tin) { std::cerr << "Tree not found\n"; return; }

    TFile *fcut1 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/gold.root", "READ");
    if (!fcut1 || fcut1->IsZombie()) { std::cerr << "Cannot open cut1 file\n"; return; }
    TCutG *cut1 = (TCutG*)fcut1->Get("cut1");
    TCutG *cut2 = (TCutG*)fcut1->Get("cut2");
    if (!cut1 || !cut2) { std::cerr << "TCutG cut1/cut2 not found\n"; return; }

    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);

    const int nbins = 3;
    std::vector<double> energy_bins = {40, 240, 450, 1000};
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // =====================================================================
    // CONFIGURACIÓN DEL ANÁLISIS (gold)
    // -----------------------------------------------------------------------
    // makeGoldConfig() fija cfg.efficiency_file y, sobre todo,
    // cfg.energy_bins_eff = {30, 300, 600, 1000} -- SOLO 3 bins de energía
    // para la eficiencia, distintos de los 10 bins finos usados aquí para
    // la anisotropía. Por eso NO se puede usar el mismo índice "ebin" para
    // leer "eff_ebin_%d"; hay que mapear cada bin fino al bin grueso de
    // eficiencia que le corresponde (se hace más abajo con findBin()).
    //
    // cos_det_cut y n_toys no aparecen en makeGoldConfig(); si no tienen
    // valor por defecto en la definición de AnalysisConfig (que no tengo
    // a la vista), hay que fijarlos aquí explícitamente -- ajusta estos dos
    // valores a los reales de tu análisis.
    // =====================================================================
    AnalysisConfig cfg = makeUraniumConfig(energy_bins);
    cfg.cos_det_cut = 0.6;   // TODO: corte real en cos(theta_det)
    cfg.n_toys       = 2000; // TODO: nº de toys del bootstrap

    // TODO: confirmar que efficiency_file vive en este directorio; el
    // string guardado en config.h es solo el nombre del fichero, no la ruta
    // completa.
    std::string eff_path = cfg.efficiency_file;
    TFile *eff_file = TFile::Open(eff_path.c_str(), "READ");
    if (!eff_file || eff_file->IsZombie()) { std::cerr << "Cannot open efficiency file: " << eff_path << "\n"; return; }

    // ===== EVENT LOOP =====
    double tof1, tof0, neutron_energy;
    float amp0, amp1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("tof1", &tof1);
    tin->SetBranchAddress("tof0", &tof0);
    tin->SetBranchAddress("amp0", &amp0);
    tin->SetBranchAddress("amp1", &amp1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);
    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);

    Vec3D counts(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Vec3D u_counts(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));

    Long64_t nentries = tin->GetEntries();
    for (int i = 0; i < nentries; i++) {
        tin->GetEntry(i);
        double dt = tof1 - tof0;
        if (std::fabs(cos_theta) > 1 || std::fabs(cos_theta_det) > 1) continue;
        if (neutron_energy < 40.0 || neutron_energy > 1000.0) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;
        int bin_beam = int(std::fabs(cos_theta) / dcos_beam);
        int bin_det = int(std::fabs(cos_theta_det) / dcos_det);

        if (cut1->IsInside(amp1 + amp0, dt) && cut2->IsInside((amp1 - amp0) / (amp0 + amp1), dt)) {
            counts[bin][bin_beam][bin_det]++;
        }
    }

    for (int e = 0; e < nbins; e++)
        for (int b = 0; b < nbins_beam; b++)
            for (int d = 0; d < nbins_det; d++)
                u_counts[e][b][d] = std::sqrt(counts[e][b][d]);

    fin->Close();
    fcut1->Close();

    // =====================================================================
    // CALCULAR ANISOTROPÍA POR BIN DE ENERGÍA
    // -----------------------------------------------------------------------
    // Para cada bin FINO de energía (los 10 de "energy_bins"), se busca a
    // qué bin GRUESO de eficiencia pertenece según cfg.energy_bins_eff
    // (3 bins: [30,300), [300,600), [600,1000]) y se carga ese histograma
    // "eff_ebin_%d" -- varios bins finos comparten la misma curva de
    // eficiencia si caen dentro del mismo bin grueso. Los histogramas se
    // cachean para no releerlos del fichero más de una vez por bin grueso.
    // =====================================================================
    std::vector<AnisotropyResult> all_anisotropies(nbins);
    std::map<int, std::pair<std::vector<double>, std::vector<double>>> eff_cache;

    for (int ebin = 0; ebin < nbins; ebin++) {
        double E_center = std::sqrt(E_low[ebin] * E_high[ebin]);
        int eff_bin = findBin(cfg.energy_bins_eff, E_center);
        if (eff_bin < 0 || eff_bin >= (int)cfg.energy_bins_eff.size() - 1) {
            std::cerr << "Warning: E=" << E_center
                      << " MeV fuera del rango de energy_bins_eff, se salta bin " << ebin << "\n";
            continue;
        }

        if (eff_cache.find(eff_bin) == eff_cache.end()) {
            TH1D *h_eff = (TH1D*)eff_file->Get(Form("eff_ebin_%d", eff_bin));
            if (!h_eff) {
                std::cerr << "Warning: eff_ebin_" << eff_bin << " not found in efficiency file, skipping bin " << ebin << "\n";
                continue;
            }
            std::vector<double> eps(nbins_det, 0.0), u_eps(nbins_det, 0.0);
            for (int d = 0; d < nbins_det; d++) {
                eps[d]   = h_eff->GetBinContent(d + 1);
                u_eps[d] = h_eff->GetBinError(d + 1);
            }
            eff_cache[eff_bin] = {eps, u_eps};
        }

        std::vector<double> eps   = eff_cache[eff_bin].first;
        std::vector<double> u_eps = eff_cache[eff_bin].second;

        all_anisotropies[ebin] = anisotropy(
            nbins_beam, nbins_det,
            counts, u_counts,
            acceptance,
            ebin,
            eps, u_eps,
            cfg
        );
    }
    eff_file->Close();

    // =====================================================================
    // GRÁFICO: ratio W(0 deg) / W(90 deg) vs energía del neutrón
    // -----------------------------------------------------------------------
    // result.w ya está normalizado al bin de referencia (bin_90 = primer bin
    // de cos_theta_beam, ~90 grados). Por tanto el último bin de result.w
    // (el más cercano a cos_theta_beam = 1, es decir ~0 grados, dirección
    // del haz) ES DIRECTAMENTE el ratio W(0)/W(90) para ese bin de energía
    // -- no hace falta ningún cálculo adicional, solo leerlo.
    // =====================================================================
    std::vector<double> xE, yRatio, eyRatio;
    for (int e = 0; e < nbins; e++) {
        if (all_anisotropies[e].w.empty()) continue;
        double val = all_anisotropies[e].w.back();
        double err = all_anisotropies[e].u_w.back();
        if (!std::isfinite(val)) continue;
        xE.push_back(std::sqrt(E_low[e] * E_high[e]));
        yRatio.push_back(val);
        eyRatio.push_back(std::isfinite(err) ? err : 0.0);
    }

    TCanvas *c_ratio = new TCanvas("c_anisotropy_ratio", "Anisotropy ratio W(0)/W(90) vs Energy", 1100, 850);
    c_ratio->SetLogx();

    TGraphErrors *g_ratio = new TGraphErrors((int)xE.size(), xE.data(), yRatio.data(), nullptr, eyRatio.data());
    g_ratio->SetName("anisotropy_ratio_0_90_vs_E");
    g_ratio->SetTitle(";E_{n} (MeV);W(0^{o}) / W(90^{o})");
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerColor(kAzure + 2);
    g_ratio->SetLineColor(kAzure + 2);
    g_ratio->SetMarkerSize(1.3);
    g_ratio->SetLineWidth(1);
    g_ratio->Draw("AP");
    g_ratio->GetXaxis()->SetTitleSize(0.05);
    g_ratio->GetYaxis()->SetTitleSize(0.05);
    g_ratio->GetXaxis()->SetLabelSize(0.045);
    g_ratio->GetYaxis()->SetLabelSize(0.045);

    TLine *iso_line = new TLine(g_ratio->GetXaxis()->GetXmin(), 1.0, g_ratio->GetXaxis()->GetXmax(), 1.0);
    iso_line->SetLineStyle(9);
    iso_line->SetLineColor(kGray + 2);
    iso_line->SetLineWidth(1);
    iso_line->Draw("same");

    DrawTagAniso(0.16, 0.90, "");
    c_ratio->Update();
    c_ratio->SaveAs("anisotropy_ratio_0_90_vs_energy.png");
    c_ratio->SaveAs("anisotropy_ratio_0_90_vs_energy.pdf");

    // ===== GUARDAR =====
    TFile *fout = new TFile("/Users/nico/Desktop/Tese/Analysis/cross_section/anisotropy_au.root", "RECREATE");
    g_ratio->Write();
    c_ratio->Write();

    // Guardar también la distribución angular completa por bin de energía,
    // por si hace falta revisarla más adelante (no se dibuja aquí).
    for (int e = 0; e < nbins; e++) {
        if (all_anisotropies[e].w.empty()) continue;
        std::vector<double> cosb(nbins_beam);
        for (int b = 0; b < nbins_beam; b++) cosb[b] = (b + 0.5) * dcos_beam;
        TGraphErrors *g_ang = new TGraphErrors(nbins_beam, cosb.data(),
                                                all_anisotropies[e].w.data(), nullptr,
                                                all_anisotropies[e].u_w.data());
        g_ang->SetName(Form("W_costheta_beam_ebin_%d", e));
        g_ang->SetTitle(";cos(#theta_{beam});W(#theta) / W(90^{o})");
        g_ang->Write();
    }

    fout->Close();

    std::cout << "\n=== DONE ===\n";
    std::cout << "Plots saved:\n";
    std::cout << "  - anisotropy_ratio_0_90_vs_energy.png / .pdf\n";
    std::cout << "  - anisotropy_au.root\n";
}