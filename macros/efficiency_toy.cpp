#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/constants.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/types.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/acceptance.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/efficiency.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/config.h"

// =====================================================================
//  Estilo "paper de física de partículas" (LHCb/ATLAS-like)
// =====================================================================
void SetPhysicsStyle(TStyle *st = gStyle){
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
    st->SetPalette(kBird);
    st->SetNumberContours(255);
}

// Marca tipo "ATLAS-style" en la esquina del pad
void DrawTag(double x, double y, const char* label, double size = 0.045){
    TLatex *lx = new TLatex();
    lx->SetNDC();
    lx->SetTextFont(42);
    lx->SetTextSize(size);
    lx->DrawLatex(x, y, Form("#font[62]{U^{238}} %s", label));
}

void efficiency_toy(){

    SetPhysicsStyle();

    // Paleta cualitativa fija (evita el rojo/verde puro, más "paper-friendly")
    const int nColors = 10;
    int colors[] = {
        kBlack, kAzure + 2, kRed + 1, kSpring - 1, kOrange + 7,
        kViolet - 3, kTeal + 3, kMagenta + 2, kGray + 2, kYellow + 2
    };
    int markers[] = {20, 21, 22, 23, 33, 34, 24, 25, 26, 27};

    // ===== CARGAR DATOS =====
    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*)fin->Get("events_gold");
    if (!tin) { std::cerr << "Tree not found\n"; return; }
    TFile *fcut1 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/gold.root", "READ");
    if (!fcut1 || fcut1->IsZombie()) { std::cerr << "Cannot open cut1 file\n"; return; }
    TCutG *cut1 = (TCutG*)fcut1->Get("cut1");
    TCutG *cut2 = (TCutG*)fcut1->Get("cut2");
    //TCutG *cut3 = (TCutG*)fcut1->Get("cut3");
    if (!cut1 || !cut2) { std::cerr << "TCutG cut1/cut2 not found\n"; return; }

    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);

    const int nbins = 4;
    std::vector<double> energy_bins = {40, 200, 400, 700, 1000};
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // ===== EVENT LOOP =====
    double tof1, tof0, neutron_energy;
    double amp0, amp1;
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

        if (bool passCut = amp0>7e3 && (amp0 + amp1 >= 17000.0)) {
            counts[bin][bin_beam][bin_det]++;
        }
    }

    for (int e = 0; e < nbins; e++)
        for (int b = 0; b < nbins_beam; b++)
            for (int d = 0; d < nbins_det; d++)
                u_counts[e][b][d] = std::sqrt(counts[e][b][d]);

    fin->Close();
    //fcut1->Close();

    // ===== CALCULAR EFICIENCIAS =====
    std::vector<EfficiencyResult> all_efficiencies(nbins);
    std::vector<int> ref_bins(nbins);

    for (int ebin = 0; ebin < nbins; ebin++) {
        ref_bins[ebin] = -1;
        for (int d = nbins_det - 1; d >= 0; d--) {
            bool has_counts = false;
            for (int b = 0; b < nbins_beam; b++) {
                if (counts[ebin][b][d] > 0 && acceptance[b][d] > 0) { has_counts = true; break; }
            }
            if (has_counts) { ref_bins[ebin] = d; break; } // first that has counts
        }

        if (ref_bins[ebin] == -1) {
            all_efficiencies[ebin].eps.resize(nbins_det, 0.0);
            all_efficiencies[ebin].u_eps.resize(nbins_det, 0.0);
            continue;
        }

        all_efficiencies[ebin] = computeEfficiency(
            ref_bins[ebin], nbins_beam, nbins_det, counts, u_counts, acceptance, ebin);
    }

    // =====================================================================
    // EFICIENCIA VS ÁNGULO DE HAZ, DERIVADA DE LA EFICIENCIA EN SISTEMA DETECTOR
    // -----------------------------------------------------------------------
    // La eficiencia intrínseca del detector es función de cos_theta_det (el
    // ángulo con el que el neutrón llega al detector), no del ángulo de haz.
    // No tiene sentido re-determinarla de forma independiente en el eje de
    // haz con el mismo método de "bin de referencia" -- eso daría dos
    // determinaciones de la misma cantidad física por caminos distintos,
    // sin que una sea más correcta que la otra.
    //
    // En vez de eso, aquí se PROYECTA eps_det(E, cos_theta_det), ya
    // calculada, sobre cada bin de haz: para un bin de haz b, se promedian
    // los eps_det(E,d) de los bins de detector d que ese bin de haz alcanza,
    // ponderando por la DISTRIBUCIÓN DE CUENTAS REAL counts[e][b][d] (en vez
    // de por la aceptancia geométrica pura). Como N(e,b,d) ~ Phi(e,b) *
    // acceptance(b,d) * eps_det(e,d), este peso ya lleva embebidas tanto la
    // aceptancia como la eficiencia (y la forma real de la física dentro del
    // bin, no una suposición de reparto uniforme sobre el ángulo sólido).
    //
    //   eps_beam(E,b)   = Σ_d N(E,b,d) * eps_det(E,d)  /  Σ_d N(E,b,d)
    //   u_eps_beam(E,b) = sqrt( Σ_d (N(E,b,d)/ΣN)^2 * u_eps_det(E,d)^2 )
    // =====================================================================
    std::vector<EfficiencyResult> all_eff_beam(nbins);

    for (int ebin = 0; ebin < nbins; ebin++) {
        all_eff_beam[ebin].eps.resize(nbins_beam, 0.0);
        all_eff_beam[ebin].u_eps.resize(nbins_beam, 0.0);

        for (int b = 0; b < nbins_beam; b++) {
            double wsum = 0.0, wsum_eps = 0.0, var = 0.0;
            for (int d = 0; d < nbins_det; d++) {
                double w = counts[ebin][b][d];
                if (w <= 0) continue;
                wsum     += w;
                wsum_eps += w * all_efficiencies[ebin].eps[d];
                var      += w * w * all_efficiencies[ebin].u_eps[d] * all_efficiencies[ebin].u_eps[d];
            }
            if (wsum > 0) {
                all_eff_beam[ebin].eps[b]   = wsum_eps / wsum;
                all_eff_beam[ebin].u_eps[b] = std::sqrt(var) / wsum;
            }
        }
    }

    // ===== HISTOGRAMAS 1D (eficiencia vs cos_theta_det, por bin de energía) =====
    std::vector<TH1D*> efficiencies(nbins, nullptr);
    for (int i = 0; i < nbins; i++) {
        efficiencies[i] = new TH1D(Form("eff_ebin_%d", i), "", nbins_det, 0, 1);
        efficiencies[i]->SetDirectory(0);
        efficiencies[i]->SetLineWidth(2);
        efficiencies[i]->SetMarkerSize(1.1);
        for (int d = 0; d < nbins_det; d++) {
            efficiencies[i]->SetBinContent(d + 1, all_efficiencies[i].eps[d]);
            efficiencies[i]->SetBinError(d + 1, all_efficiencies[i].u_eps[d]);
        }
    }

    std::vector<TH1D*> efficiencies_beam(nbins, nullptr);
    for (int i = 0; i < nbins; i++) {
        efficiencies_beam[i] = new TH1D(Form("eff_beam_ebin_%d", i), "", nbins_beam, 0, 1);
        efficiencies_beam[i]->SetDirectory(0);
        efficiencies_beam[i]->SetLineWidth(2);
        efficiencies_beam[i]->SetMarkerSize(1.1);
        for (int b = 0; b < nbins_beam; b++) {
            efficiencies_beam[i]->SetBinContent(b + 1, all_eff_beam[i].eps[b]);
            efficiencies_beam[i]->SetBinError(b + 1, all_eff_beam[i].u_eps[b]);
        }
    }

    // =====================================================================
    // EFICIENCIA PROMEDIO POR BIN DE ENERGÍA (colapsando beam y det bins)
    // -----------------------------------------------------------------------
    // Igual que eps_beam, pero ahora se pesa por TODAS las cuentas del bin
    // de energía (suma sobre b y d), no solo sobre d dentro de un b fijo.
    //
    //   eps_avg(E)   = Σ_{b,d} N(E,b,d) * eps_det(E,d)  /  Σ_{b,d} N(E,b,d)
    //   u_eps_avg(E) = sqrt( Σ_{b,d} (N(E,b,d)/ΣN)^2 * u_eps_det(E,d)^2 )
    // =====================================================================
    std::vector<double> eps_avg(nbins, 0.0), u_eps_avg(nbins, 0.0);

    for (int e = 0; e < nbins; e++) {
        double wsum = 0.0, wsum_eps = 0.0, var = 0.0;
        for (int b = 0; b < nbins_beam; b++) {
            for (int d = 0; d < nbins_det; d++) {
                double w = counts[e][b][d];
                if (w <= 0) continue;
                wsum     += w;
                wsum_eps += w * all_efficiencies[e].eps[d];
                var      += w * w * all_efficiencies[e].u_eps[d] * all_efficiencies[e].u_eps[d];
            }
        }
        if (wsum > 0) {
            eps_avg[e]   = wsum_eps / wsum;
            u_eps_avg[e] = std::sqrt(var) / wsum;
        }
        std::cout << "Bin " << e << " (" << E_low[e] << "-" << E_high[e]
                   << " MeV): <eff> = " << eps_avg[e] << " +/- " << u_eps_avg[e] << "\n";
    }

    // =====================================================================
    // EFICIENCIA PROMEDIO POR BIN DE ENERGÍA, PONDERADA POR ACEPTANCIA
    // -----------------------------------------------------------------------
    // A diferencia de la versión anterior (ponderada por counts), aquí el
    // peso es la aceptancia geométrica pura dOmega(b,d), independiente de
    // cuántos eventos cayeron realmente en cada bin. Esto evita el sesgo
    // circular de sobre-pesar los bins donde eps ya es alta (porque ahí se
    // detectan más cuentas). Es la versión correcta para usar como factor
    // de corrección al extraer flujo / sección eficaz.
    //
    //   eps_avg(E)   = Σ_{b,d} dOmega(b,d) * eps_det(E,d)  /  Σ_{b,d} dOmega(b,d)
    //   u_eps_avg(E) = sqrt( Σ_{b,d} (dOmega(b,d)/ΣdOmega)^2 * u_eps_det(E,d)^2 )
    //
    // Solo se suman los bins de detector donde eps_det fue efectivamente
    // determinada (d <= ref_bins[e]); más allá del bin de referencia no hay
    // información y no debe contribuir con eps=0 artificial.
    // =====================================================================
    std::vector<double> eps_avg_acc(nbins, 0.0), u_eps_avg_acc(nbins, 0.0);

    for (int e = 0; e < nbins; e++) {
        double wsum = 0.0, wsum_eps = 0.0, var = 0.0;

        int dmax = (ref_bins[e] >= 0) ? ref_bins[e] : -1;
        for (int b = 0; b < nbins_beam; b++) {
            for (int d = 0; d <= dmax; d++) {
                double w = acceptance[b][d];
                if (w <= 0) continue;
                wsum     += w;
                wsum_eps += w * all_efficiencies[e].eps[d];
                var      += w * w * all_efficiencies[e].u_eps[d] * all_efficiencies[e].u_eps[d];
            }
        }
        if (wsum > 0) {
            eps_avg_acc[e]   = wsum_eps / wsum;
            u_eps_avg_acc[e] = std::sqrt(var) / wsum;
        }
        std::cout << "Bin " << e << " (" << E_low[e] << "-" << E_high[e]
                   << " MeV): <eff>_acc = " << eps_avg_acc[e] << " +/- " << u_eps_avg_acc[e] << "\n";
    }

    // ===== GRÁFICO 1: Eficiencia vs cos(theta_det), overlay de todos los ebins =====
    TCanvas *c1 = new TCanvas("c_efficiencies", "Efficiency vs cos(#theta_{det})", 1100, 850);

    TLegend *leg1 = new TLegend(0.62, 0.14, 0.94, 0.50);
    leg1->SetNColumns(2);
    leg1->SetTextSize(0.026);
    leg1->SetHeader("Energy bins", "C");

    double max_eff1 = 1.1;
    for (int e = 0; e < nbins; e++)
        for (int d = 1; d <= nbins_det; d++)
            max_eff1 = std::max(max_eff1, efficiencies[e]->GetBinContent(d));
    max_eff1 = std::max(1.1, max_eff1 * 1.1);

    for (int e = 0; e < nbins; e++) {
        efficiencies[e]->SetLineColor(colors[e % nColors]);
        efficiencies[e]->SetMarkerColor(colors[e % nColors]);
        efficiencies[e]->SetMarkerStyle(markers[e % nColors]);
        efficiencies[e]->SetMinimum(0.0);
        efficiencies[e]->SetMaximum(max_eff1);

        if (e == 0) {
            efficiencies[e]->Draw("E1");
            efficiencies[e]->GetXaxis()->SetTitle("cos(#theta_{det})");
            efficiencies[e]->GetYaxis()->SetTitle("Efficiency");
        } else {
            efficiencies[e]->Draw("E1 same");
        }
        leg1->AddEntry(efficiencies[e], Form("%.1f-%.1f MeV", E_low[e], E_high[e]), "lp");
    }

    TLine *line1 = new TLine(0, 1, 1, 1);
    line1->SetLineStyle(9);
    line1->SetLineColor(kGray + 2);
    line1->SetLineWidth(1);
    line1->Draw("same");

    leg1->Draw();
    DrawTag(0.14, 0.93, "");
    c1->Update();
    c1->SaveAs("efficiency_vs_costheta.png");
    c1->SaveAs("efficiency_vs_costheta.pdf");

    // ===== GRÁFICO 2: Eficiencia vs Energía, overlay de todos los det-bins =====
    TCanvas *c2 = new TCanvas("c_eff_vs_energy", "Efficiency vs Energy", 1100, 850);
    c2->SetLogx();

    std::vector<TGraphErrors*> graphs_eff_vs_E(nbins_det, nullptr);
    std::vector<double> x_vals(nbins), y_vals(nbins), ey_vals(nbins);

    double max_eff2 = 1.1;
    for (int d = 0; d < nbins_det; d++)
        for (int e = 0; e < nbins; e++)
            max_eff2 = std::max(max_eff2, all_efficiencies[e].eps[d]);
    max_eff2 = std::max(1.15, max_eff2 * 1.15);

    TLegend *leg2 = new TLegend(0.62, 0.10, 0.94, 0.52);
    leg2->SetNColumns(2);
    leg2->SetTextSize(0.022);
    leg2->SetHeader("cos(#theta_{det})", "C");

    for (int d = 0; d < nbins_det; d++) {
        for (int e = 0; e < nbins; e++) {
            x_vals[e]  = std::sqrt(E_low[e] * E_high[e]);
            y_vals[e]  = all_efficiencies[e].eps[d];
            ey_vals[e] = all_efficiencies[e].u_eps[d];
        }

        TGraphErrors *g = new TGraphErrors(nbins, x_vals.data(), y_vals.data(), nullptr, ey_vals.data());
        g->SetName(Form("eff_vs_E_detbin_%d", d));
        g->SetTitle(";E_{n} (MeV);Efficiency");
        g->SetMarkerStyle(markers[d % nColors]);
        g->SetMarkerColor(colors[d % nColors]);
        g->SetLineColor(colors[d % nColors]);
        g->SetMarkerSize(1.1);
        g->SetLineWidth(1);
        graphs_eff_vs_E[d] = g;

        if (d == 0) {
            g->Draw("AP");
            g->GetYaxis()->SetRangeUser(0.0, max_eff2);
            g->GetXaxis()->SetTitle("E_{n} (MeV)");
            g->GetYaxis()->SetTitle("Efficiency");
        } else {
            g->Draw("P same");
        }

        double cos_lo = d * dcos_det, cos_hi = (d + 1) * dcos_det;
        leg2->AddEntry(g, Form("%.2f-%.2f", cos_lo, cos_hi), "p");
    }

    leg2->Draw();
    DrawTag(0.16, 0.93, "");
    c2->Update();
    c2->SaveAs("efficiency_vs_energy.png");
    c2->SaveAs("efficiency_vs_energy.pdf");

    // ===== GRÁFICO 7: Eficiencia PROMEDIO vs Energía =====
    TCanvas *c7 = new TCanvas("c_eff_avg_vs_energy", "Average efficiency vs Energy", 1100, 850);
    c7->SetLogx();

    TGraphErrors *g_avg = new TGraphErrors(nbins, x_vals.data(), eps_avg.data(), nullptr, u_eps_avg.data());
    g_avg->SetName("eff_avg_vs_E");
    g_avg->SetTitle(";E_{n} (MeV);#LTEfficiency#GT");
    g_avg->SetMarkerStyle(20);
    g_avg->SetMarkerColor(kAzure + 2);
    g_avg->SetLineColor(kAzure + 2);
    g_avg->SetMarkerSize(1.3);
    g_avg->SetLineWidth(2);
    g_avg->Draw("AP");
    g_avg->GetYaxis()->SetRangeUser(0.0, max_eff2);
    g_avg->GetXaxis()->SetTitle("E_{n} (MeV)");
    g_avg->GetYaxis()->SetTitle("#LTEfficiency#GT");

    TLine *line7 = new TLine(g_avg->GetXaxis()->GetXmin(), 1.0, g_avg->GetXaxis()->GetXmax(), 1.0);
    line7->SetLineStyle(9);
    line7->SetLineColor(kGray + 2);
    line7->Draw("same");

    DrawTag(0.16, 0.93, "");
    c7->Update();
    c7->SaveAs("efficiency_avg_vs_energy.png");
    c7->SaveAs("efficiency_avg_vs_energy.pdf");

    // =====================================================================
    // GRÁFICO 3: Eficiencia vs cos(theta_beam), overlay de todos los ebins
    // -----------------------------------------------------------------------
    // Proyección de eps_det(E, cos_theta_det) ponderada por la distribución
    // real de cuentas (ver bloque de cálculo más arriba). Es la eficiencia
    // media que tuvieron realmente los eventos detectados en cada bin de
    // haz -- un diagnóstico de la medida, NO un factor de corrección para
    // extraer flujo/sección eficaz (para eso hace falta la versión ponderada
    // por aceptancia geométrica, ver comentario arriba).
    // =====================================================================
    TCanvas *c3 = new TCanvas("c_eff_vs_costheta_beam", "Efficiency vs cos(#theta_{beam})", 1100, 850);

    TLegend *leg3 = new TLegend(0.62, 0.14, 0.94, 0.50);
    leg3->SetNColumns(2);
    leg3->SetTextSize(0.026);
    leg3->SetHeader("Energy bins", "C");

    double max_eff3 = 1.1;
    for (int e = 0; e < nbins; e++)
        for (int b = 1; b <= nbins_beam; b++)
            max_eff3 = std::max(max_eff3, efficiencies_beam[e]->GetBinContent(b));
    max_eff3 = std::max(1.1, max_eff3 * 1.1);

    for (int e = 0; e < nbins; e++) {
        efficiencies_beam[e]->SetLineColor(colors[e % nColors]);
        efficiencies_beam[e]->SetMarkerColor(colors[e % nColors]);
        efficiencies_beam[e]->SetMarkerStyle(markers[e % nColors]);
        efficiencies_beam[e]->SetMinimum(0.0);
        efficiencies_beam[e]->SetMaximum(max_eff3);

        if (e == 0) {
            efficiencies_beam[e]->Draw("E1");
            efficiencies_beam[e]->GetXaxis()->SetTitle("cos(#theta_{beam})");
            efficiencies_beam[e]->GetYaxis()->SetTitle("Efficiency");
        } else {
            efficiencies_beam[e]->Draw("E1 same");
        }
        leg3->AddEntry(efficiencies_beam[e], Form("%.1f-%.1f MeV", E_low[e], E_high[e]), "lp");
    }

    TLine *line3 = new TLine(0, 1, 1, 1);
    line3->SetLineStyle(9);
    line3->SetLineColor(kGray + 2);
    line3->SetLineWidth(1);
    line3->Draw("same");

    leg3->Draw();
    DrawTag(0.14, 0.93, "");
    c3->Update();
    c3->SaveAs("efficiency_vs_costheta_beam.png");
    c3->SaveAs("efficiency_vs_costheta_beam.pdf");

    // =====================================================================
    // GRÁFICO 4: Comparativa "small multiples" — un pad por bin de cos_theta_det,
    // eficiencia vs energía, para comparar la forma de la curva entre bins
    // =====================================================================
    int nCols = (int)std::ceil(std::sqrt((double)nbins_det));
    int nRows = (int)std::ceil((double)nbins_det / nCols);

    TCanvas *c4 = new TCanvas("c_eff_vs_E_grid", "Efficiency vs Energy - per cos(theta_det) bin",
                               320 * nCols, 280 * nRows);
    c4->Divide(nCols, nRows, 0.001, 0.001);

    for (int d = 0; d < nbins_det; d++) {
        TPad *pad = (TPad*)c4->cd(d + 1);
        pad->SetLogx();
        pad->SetTopMargin(0.10);
        pad->SetBottomMargin(0.18);
        pad->SetLeftMargin(0.18);
        pad->SetRightMargin(0.04);

        TGraphErrors *g = graphs_eff_vs_E[d];
        g->SetMarkerColor(kAzure + 2);
        g->SetLineColor(kAzure + 2);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.9);
        g->Draw("AP");
        g->GetYaxis()->SetRangeUser(0.0, max_eff2);
        g->GetXaxis()->SetTitleSize(0.075);
        g->GetYaxis()->SetTitleSize(0.075);
        g->GetXaxis()->SetLabelSize(0.065);
        g->GetYaxis()->SetLabelSize(0.065);
        g->GetXaxis()->SetTitle("E_{n} (MeV)");
        g->GetYaxis()->SetTitle("Efficiency");
        g->GetXaxis()->SetTitleOffset(1.05);
        g->GetYaxis()->SetTitleOffset(1.1);

        double cos_lo = d * dcos_det, cos_hi = (d + 1) * dcos_det;
        TLatex *lp = new TLatex();
        lp->SetNDC();
        lp->SetTextFont(62);
        lp->SetTextSize(0.09);
        lp->DrawLatex(0.22, 0.88, Form("cos#theta_{det} #in [%.2f, %.2f]", cos_lo, cos_hi));

        TLine *ref = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
        ref->SetLineStyle(9);
        ref->SetLineColor(kGray + 1);
        ref->Draw("same");
    }
    c4->cd(0);
    c4->Update();
    c4->SaveAs("efficiency_vs_energy_per_costheta_bin.png");
    c4->SaveAs("efficiency_vs_energy_per_costheta_bin.pdf");

    // =====================================================================
    // GRÁFICO 5: Comparativa "small multiples" — un pad por bin de cos_theta_beam,
    // eficiencia vs energía. Análogo a c4 pero para el eje de haz; sirve para
    // comparar directamente forma/nivel de la eficiencia entre ambos ejes.
    // =====================================================================
    std::vector<TGraphErrors*> graphs_eff_vs_E_beam(nbins_beam, nullptr);
    std::vector<double> yb_vals(nbins), eyb_vals(nbins);

    for (int b = 0; b < nbins_beam; b++) {
        for (int e = 0; e < nbins; e++) {
            yb_vals[e]  = all_eff_beam[e].eps[b];
            eyb_vals[e] = all_eff_beam[e].u_eps[b];
        }
        TGraphErrors *g = new TGraphErrors(nbins, x_vals.data(), yb_vals.data(), nullptr, eyb_vals.data());
        g->SetName(Form("eff_vs_E_beambin_%d", b));
        g->SetTitle(";E_{n} (MeV);Efficiency");
        graphs_eff_vs_E_beam[b] = g;
    }

    int nColsB = (int)std::ceil(std::sqrt((double)nbins_beam));
    int nRowsB = (int)std::ceil((double)nbins_beam / nColsB);

    TCanvas *c5 = new TCanvas("c_eff_vs_E_grid_beam", "Efficiency vs Energy - per cos(theta_beam) bin",
                               320 * nColsB, 280 * nRowsB);
    c5->Divide(nColsB, nRowsB, 0.001, 0.001);

    for (int b = 0; b < nbins_beam; b++) {
        TPad *pad = (TPad*)c5->cd(b + 1);
        pad->SetLogx();
        pad->SetTopMargin(0.10);
        pad->SetBottomMargin(0.18);
        pad->SetLeftMargin(0.18);
        pad->SetRightMargin(0.04);

        TGraphErrors *g = graphs_eff_vs_E_beam[b];
        g->SetMarkerColor(kOrange + 7);
        g->SetLineColor(kOrange + 7);
        g->SetMarkerStyle(21);
        g->SetMarkerSize(0.9);
        g->Draw("AP");
        g->GetYaxis()->SetRangeUser(0.0, max_eff2);
        g->GetXaxis()->SetTitleSize(0.075);
        g->GetYaxis()->SetTitleSize(0.075);
        g->GetXaxis()->SetLabelSize(0.065);
        g->GetYaxis()->SetLabelSize(0.065);
        g->GetXaxis()->SetTitle("E_{n} (MeV)");
        g->GetYaxis()->SetTitle("Efficiency");
        g->GetXaxis()->SetTitleOffset(1.05);
        g->GetYaxis()->SetTitleOffset(1.1);

        double cos_lo = b * dcos_beam, cos_hi = (b + 1) * dcos_beam;
        TLatex *lp = new TLatex();
        lp->SetNDC();
        lp->SetTextFont(62);
        lp->SetTextSize(0.09);
        lp->DrawLatex(0.20, 0.88, Form("cos#theta_{beam} #in [%.2f, %.2f]", cos_lo, cos_hi));

        TLine *ref = new TLine(g->GetXaxis()->GetXmin(), 1.0, g->GetXaxis()->GetXmax(), 1.0);
        ref->SetLineStyle(9);
        ref->SetLineColor(kGray + 1);
        ref->Draw("same");
    }
    c5->cd(0);
    c5->Update();
    c5->SaveAs("efficiency_vs_energy_per_costheta_beam_bin.png");
    c5->SaveAs("efficiency_vs_energy_per_costheta_beam_bin.pdf");

    // =====================================================================
    // GRÁFICO 6: Incertidumbre relativa (u_eps/eps) vs cos(theta_det), por
    // bin de energía. Esta es la información que un heatmap de eficiencia
    // NO te da: dónde tus resultados son estadísticamente robustos y dónde
    // son ruido. Útil para decidir qué bins reportar / cuáles fusionar.
    // =====================================================================
    TCanvas *c6 = new TCanvas("c_rel_unc", "Relative uncertainty vs cos(#theta_{det})", 1100, 850);
    c6->SetLogy();

    TLegend *leg6 = new TLegend(0.62, 0.55, 0.94, 0.91);
    leg6->SetNColumns(2);
    leg6->SetTextSize(0.026);
    leg6->SetHeader("Energy bins", "C");

    std::vector<TH1D*> rel_unc(nbins, nullptr);
    double max_runc = 0.05;
    for (int e = 0; e < nbins; e++) {
        rel_unc[e] = new TH1D(Form("relunc_ebin_%d", e), "", nbins_det, 0, 1);
        rel_unc[e]->SetDirectory(0);
        rel_unc[e]->SetLineWidth(2);
        rel_unc[e]->SetMarkerSize(1.1);
        for (int d = 0; d < nbins_det; d++) {
            double eps = all_efficiencies[e].eps[d];
            double relerr = (eps > 0) ? all_efficiencies[e].u_eps[d] / eps : 0.0;
            rel_unc[e]->SetBinContent(d + 1, relerr);
            max_runc = std::max(max_runc, relerr);
        }
    }
    max_runc *= 1.3;

    for (int e = 0; e < nbins; e++) {
        rel_unc[e]->SetLineColor(colors[e % nColors]);
        rel_unc[e]->SetMarkerColor(colors[e % nColors]);
        rel_unc[e]->SetMarkerStyle(markers[e % nColors]);
        rel_unc[e]->SetMinimum(1e-3);
        rel_unc[e]->SetMaximum(max_runc);

        if (e == 0) {
            rel_unc[e]->Draw("PL");
            rel_unc[e]->GetXaxis()->SetTitle("cos(#theta_{det})");
            rel_unc[e]->GetYaxis()->SetTitle("#sigma_{#varepsilon} / #varepsilon");
        } else {
            rel_unc[e]->Draw("PL same");
        }
        leg6->AddEntry(rel_unc[e], Form("%.1f-%.1f MeV", E_low[e], E_high[e]), "lp");
    }

    leg6->Draw();
    DrawTag(0.14, 0.93, "");
    c6->Update();
    c6->SaveAs("efficiency_relative_uncertainty.png");
    c6->SaveAs("efficiency_relative_uncertainty.pdf");

    // ===== GRÁFICO 8: Eficiencia PROMEDIO (ponderada por aceptancia) vs Energía =====
    TCanvas *c8 = new TCanvas("c_eff_avg_acc_vs_energy", "Average efficiency (acceptance-weighted) vs Energy", 1100, 850);
    c8->SetLogx();

    TGraphErrors *g_avg_acc = new TGraphErrors(nbins, x_vals.data(), eps_avg_acc.data(), nullptr, u_eps_avg_acc.data());
    g_avg_acc->SetName("eff_avg_acc_vs_E");
    g_avg_acc->SetTitle(";E_{n} (MeV);#LTEfficiency#GT_{acc}");
    g_avg_acc->SetMarkerStyle(20);
    g_avg_acc->SetMarkerColor(kOrange + 7);
    g_avg_acc->SetLineColor(kOrange + 7);
    g_avg_acc->SetMarkerSize(1.3);
    g_avg_acc->SetLineWidth(2);
    g_avg_acc->Draw("AP");
    g_avg_acc->GetYaxis()->SetRangeUser(0.0, max_eff2);
    g_avg_acc->GetXaxis()->SetTitle("E_{n} (MeV)");
    g_avg_acc->GetYaxis()->SetTitle("#LTEfficiency#GT (acceptance-weighted)");

    TLine *line8 = new TLine(g_avg_acc->GetXaxis()->GetXmin(), 1.0, g_avg_acc->GetXaxis()->GetXmax(), 1.0);
    line8->SetLineStyle(9);
    line8->SetLineColor(kGray + 2);
    line8->Draw("same");

    DrawTag(0.16, 0.93, "");
    c8->Update();
    c8->SaveAs("efficiency_avg_acceptance_weighted_vs_energy.png");
    c8->SaveAs("efficiency_avg_acceptance_weighted_vs_energy.pdf");

    // ===== GUARDAR =====
    TFile *fout = new TFile("/Users/nico/Desktop/Tese/Analysis/cross_section/efficiencies_u_toy.root", "RECREATE");

    for (int i = 0; i < nbins; i++) efficiencies[i]->Write();
    for (int i = 0; i < nbins; i++) efficiencies_beam[i]->Write();
    for (int i = 0; i < nbins; i++) rel_unc[i]->Write();
    for (int d = 0; d < nbins_det; d++) if (graphs_eff_vs_E[d]) graphs_eff_vs_E[d]->Write();
    for (int b = 0; b < nbins_beam; b++) if (graphs_eff_vs_E_beam[b]) graphs_eff_vs_E_beam[b]->Write();

    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    c6->Write();
    g_avg->Write();
    c7->Write();
    g_avg_acc->Write();
    c8->Write();

    fout->Close();

    std::cout << "\n=== DONE ===\n";
    std::cout << "Plots saved:\n";
    std::cout << "  - efficiency_vs_costheta.png / .pdf              (eff vs cos_theta_det)\n";
    std::cout << "  - efficiency_vs_energy.png / .pdf                (eff vs E, overlay per det bin)\n";
    std::cout << "  - efficiency_vs_costheta_beam.png / .pdf         (eff vs cos_theta_beam, consistency check)\n";
    std::cout << "  - efficiency_vs_energy_per_costheta_bin.png/.pdf (small multiples, det bins)\n";
    std::cout << "  - efficiency_vs_energy_per_costheta_beam_bin.*   (small multiples, beam bins)\n";
    std::cout << "  - efficiency_relative_uncertainty.png / .pdf     (statistical reliability map)\n";
    std::cout << "  - efficiencies_u_toy.root\n";
}