#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"

// ============================================================
//  Ajuste gaussiano en ventana FWTM (full width at tenth max)
//  Devuelve {mean, sigma, mean_err, sigma_err} vía referencias
// ============================================================
bool FitFWTM(TH1D* h, double& mean, double& sigma, double& mean_err, double& sigma_err)
{
    if (!h || h->GetEntries() < 20) return false;

    int nbins_h   = h->GetNbinsX();
    int maxbin    = h->GetMaximumBin();
    double peak   = h->GetBinCenter(maxbin);
    double maxval = h->GetBinContent(maxbin);
    if (maxval <= 0) return false;

    double tenth = maxval / 10.0;

    double left  = h->GetBinCenter(1);
    for (int b = maxbin; b >= 1; b--) {
        if (h->GetBinContent(b) > tenth) continue;
        left = h->GetBinCenter(b);
        break;
    }

    double right = h->GetBinCenter(nbins_h);
    for (int b = maxbin; b <= nbins_h; b++) {
        if (h->GetBinContent(b) > tenth) continue;
        right = h->GetBinCenter(b);
        break;
    }

    if (right <= left) return false;

    TF1 f("fwtm_gaus", "gaus", left, right);
    f.SetParameter(1, peak);
    TFitResultPtr r = h->Fit(&f, "QRSN", "", left, right);
    if (!r.Get() || r->Status() != 0) {
        // seguimos adelante igualmente, pero avisamos
        std::cout << "[WARN] fit no convergio bien en " << h->GetName() << std::endl;
    }

    mean      = f.GetParameter(1);
    sigma     = f.GetParameter(2);
    mean_err  = f.GetParError(1);
    sigma_err = f.GetParError(2);
    return true;
}

// ============================================================
//  Overlay de un conjunto de histogramas (uno por ebin) en un
//  único canvas, y lo guarda + escribe
// ============================================================
void DrawOverlay(std::vector<TH1D*>& hists, int ebins,
                  std::vector<double>& E_low, std::vector<double>& E_high,
                  const char* title, const char* fname, TFile* fout)
{
    const int colors[] = {kBlack, kAzure+2, kRed+1, kSpring-1, kOrange+7,
                           kViolet-3, kTeal+3, kMagenta+2, kGray+2, kYellow+2};
    const int markers[] = {20,21,22,23,33,34,24,25,26,27};

    TCanvas* c = new TCanvas(Form("c_%s", fname), title, 1000, 800);
    TLegend* leg = new TLegend(0.65, 0.55, 0.94, 0.90);
    leg->SetNColumns(1);
    leg->SetTextSize(0.025);
    leg->SetBorderSize(0);

    // normalizar cada histograma a area unitaria (comparar forma, no cuentas absolutas)
    std::vector<TH1D*> hnorm(ebins, nullptr);
    for (int i = 0; i < ebins; i++) {
        if (!hists[i] || hists[i]->GetEntries() == 0) continue;
        hnorm[i] = (TH1D*) hists[i]->Clone(Form("%s_norm", hists[i]->GetName()));
        hnorm[i]->SetDirectory(0);
        double integral = hnorm[i]->Integral();
        if (integral > 0) hnorm[i]->Scale(1.0 / integral);
    }

    double maxY = 0;
    for (int i = 0; i < ebins; i++)
        if (hnorm[i]) maxY = std::max(maxY, hnorm[i]->GetMaximum());

    bool first = true;
    for (int i = 0; i < ebins; i++) {
        if (!hnorm[i]) continue;
        hnorm[i]->SetLineColor(colors[i % 10]);
        hnorm[i]->SetMarkerColor(colors[i % 10]);
        hnorm[i]->SetMarkerStyle(markers[i % 10]);
        hnorm[i]->SetMaximum(maxY * 1.2);
        hnorm[i]->SetTitle(title);
        hnorm[i]->GetYaxis()->SetTitle("Normalized counts");
        if (first) { hnorm[i]->Draw("E1"); first = false; }
        else        hnorm[i]->Draw("E1 same");
        leg->AddEntry(hnorm[i], Form("%.1f-%.1f MeV", E_low[i], E_high[i]), "lp");
    }
    leg->Draw();
    c->SaveAs(Form("%s.png", fname));
    fout->cd();
    c->Write();
    for (auto* h : hnorm) if (h) h->Write();
}

void sums_study(){

    const int    ebins = 2;
    const int    bins  = 100;
    const double emin  = 40.;
    const double emax  = 1000.;
    const double sum_lo = 80., sum_hi = 120.;   // rango de los sumX/sumY

    TFile *fin = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root");
    if (!fin || fin->IsZombie()) { std::cerr << "No se pudo abrir coincidences.root\n"; return; }

    TTree *tin_u = (TTree*) fin->Get("events_gold");
    if (!tin_u) { std::cerr << "No se encontro events_uranium\n"; return; }

    Long64_t nentries_u = tin_u->GetEntries();

    double sumX0_u, sumX1_u, sumY1_u, sumY0_u;
    double neutron_energy_u;
    float  amp0_u, amp1_u;

    tin_u->SetBranchAddress("sumX0", &sumX0_u);
    tin_u->SetBranchAddress("sumY0", &sumY0_u);
    tin_u->SetBranchAddress("sumX1", &sumX1_u);
    tin_u->SetBranchAddress("sumY1", &sumY1_u);
    tin_u->SetBranchAddress("neutron_energy", &neutron_energy_u);
    tin_u->SetBranchAddress("amp0", &amp0_u);
    tin_u->SetBranchAddress("amp1", &amp1_u);

    std::vector<double> neutron_energy_bins = {40, 300, 1000};
    std::vector<double> E_low(ebins), E_high(ebins);
    for (int e = 0; e < ebins; ++e) {
        E_low[e]  = neutron_energy_bins[e];
        E_high[e] = neutron_energy_bins[e+1];
    }

    std::vector<TH1D*> hists_sumx0(ebins, nullptr), hists_sumx1(ebins, nullptr);
    std::vector<TH1D*> hists_sumy0(ebins, nullptr), hists_sumy1(ebins, nullptr);
    std::vector<TH1D*> hists_sumx0_b(ebins, nullptr), hists_sumx1_b(ebins, nullptr);
    std::vector<TH1D*> hists_sumy0_b(ebins, nullptr), hists_sumy1_b(ebins, nullptr);

    for (int i = 0; i < ebins; i++){
        hists_sumx0[i]   = new TH1D(Form("hists_x0_ebin%d",   i), "", bins, sum_lo, sum_hi);
        hists_sumx1[i]   = new TH1D(Form("hists_x1_ebin%d",   i), "", bins, sum_lo, sum_hi);
        hists_sumy0[i]   = new TH1D(Form("hists_y0_ebin%d",   i), "", bins, sum_lo, sum_hi);
        hists_sumy1[i]   = new TH1D(Form("hists_y1_ebin%d",   i), "", bins, sum_lo, sum_hi);
        hists_sumx0_b[i] = new TH1D(Form("hists_x0_b_ebin%d", i), "", bins, sum_lo, sum_hi);
        hists_sumx1_b[i] = new TH1D(Form("hists_x1_b_ebin%d", i), "", bins, sum_lo, sum_hi);
        hists_sumy0_b[i] = new TH1D(Form("hists_y0_b_ebin%d", i), "", bins, sum_lo, sum_hi);
        hists_sumy1_b[i] = new TH1D(Form("hists_y1_b_ebin%d", i), "", bins, sum_lo, sum_hi);
        for (auto* h : {hists_sumx0[i], hists_sumx1[i], hists_sumy0[i], hists_sumy1[i],
                        hists_sumx0_b[i], hists_sumx1_b[i], hists_sumy0_b[i], hists_sumy1_b[i]})
            h->SetDirectory(0);
    }

    for (Long64_t i = 0; i < nentries_u; i++){
        tin_u->GetEntry(i);
        int ebin = findBin(neutron_energy_bins, neutron_energy_u);
        if (ebin < 0 || ebin >= ebins) continue;

        bool is_signal = (amp0_u + amp1_u > 15e3) && (amp0_u+amp1_u)<35e3 &&
                          ((amp1_u - amp0_u) / (amp0_u + amp1_u) < 0.2);

        if (is_signal){
            hists_sumx0[ebin]->Fill(sumX0_u);
            hists_sumy0[ebin]->Fill(sumY0_u);
            hists_sumx1[ebin]->Fill(sumX1_u);
            hists_sumy1[ebin]->Fill(sumY1_u);
        } else {
            hists_sumx0_b[ebin]->Fill(sumX0_u);
            hists_sumy0_b[ebin]->Fill(sumY0_u);
            hists_sumx1_b[ebin]->Fill(sumX1_u);
            hists_sumy1_b[ebin]->Fill(sumY1_u);
        }
    }
    fin->Close();

    TFile* fout = new TFile("sums_study_gold.root", "RECREATE");

    // ---- overlays por variable, señal y fondo ----
    DrawOverlay(hists_sumx0,   ebins, E_low, E_high, "sumX0 signal;sumX0;counts",     "sumx0_signal",     fout);
    DrawOverlay(hists_sumx1,   ebins, E_low, E_high, "sumX1 signal;sumX1;counts",     "sumx1_signal",     fout);
    DrawOverlay(hists_sumy0,   ebins, E_low, E_high, "sumY0 signal;sumY0;counts",     "sumy0_signal",     fout);
    DrawOverlay(hists_sumy1,   ebins, E_low, E_high, "sumY1 signal;sumY1;counts",     "sumy1_signal",     fout);
    DrawOverlay(hists_sumx0_b, ebins, E_low, E_high, "sumX0 background;sumX0;counts", "sumx0_background", fout);
    DrawOverlay(hists_sumx1_b, ebins, E_low, E_high, "sumX1 background;sumX1;counts", "sumx1_background", fout);
    DrawOverlay(hists_sumy0_b, ebins, E_low, E_high, "sumY0 background;sumY0;counts", "sumy0_background", fout);
    DrawOverlay(hists_sumy1_b, ebins, E_low, E_high, "sumY1 background;sumY1;counts", "sumy1_background", fout);

    // ---- fits FWTM y evolucion centroide/sigma vs energia ----
    // agrupamos señal y fondo de la MISMA variable para superponerlos
    struct VarPair { std::vector<TH1D*>* sig; std::vector<TH1D*>* bkg; const char* tag; };
    std::vector<VarPair> varpairs = {
        {&hists_sumx0, &hists_sumx0_b, "sumX0"},
        {&hists_sumx1, &hists_sumx1_b, "sumX1"},
        {&hists_sumy0, &hists_sumy0_b, "sumY0"},
        {&hists_sumy1, &hists_sumy1_b, "sumY1"},
    };

    std::vector<double> Ecenter(ebins);
    for (int e = 0; e < ebins; e++) Ecenter[e] = std::sqrt(E_low[e] * E_high[e]);

    auto fitAll = [&](std::vector<TH1D*>* hists,
                       std::vector<double>& mean, std::vector<double>& sigma,
                       std::vector<double>& mean_err, std::vector<double>& sigma_err){
        for (int e = 0; e < ebins; e++) {
            double m=0, s=0, me=0, se=0;
            if (FitFWTM((*hists)[e], m, s, me, se)) {
                mean[e]=m; sigma[e]=s; mean_err[e]=me; sigma_err[e]=se;
            }
        }
    };

    for (auto& vp : varpairs) {
        std::vector<double> mean_s(ebins,0), sigma_s(ebins,0), mean_err_s(ebins,0), sigma_err_s(ebins,0);
        std::vector<double> mean_b(ebins,0), sigma_b(ebins,0), mean_err_b(ebins,0), sigma_err_b(ebins,0);
        std::vector<double> ex(ebins,0);

        fitAll(vp.sig, mean_s, sigma_s, mean_err_s, sigma_err_s);
        fitAll(vp.bkg, mean_b, sigma_b, mean_err_b, sigma_err_b);

        TGraphErrors* g_mean_s  = new TGraphErrors(ebins, Ecenter.data(), mean_s.data(),  ex.data(), mean_err_s.data());
        TGraphErrors* g_sigma_s = new TGraphErrors(ebins, Ecenter.data(), sigma_s.data(), ex.data(), sigma_err_s.data());
        TGraphErrors* g_mean_b  = new TGraphErrors(ebins, Ecenter.data(), mean_b.data(),  ex.data(), mean_err_b.data());
        TGraphErrors* g_sigma_b = new TGraphErrors(ebins, Ecenter.data(), sigma_b.data(), ex.data(), sigma_err_b.data());

        g_mean_s->SetName(Form("centroid_vs_E_%s_signal", vp.tag));
        g_sigma_s->SetName(Form("sigma_vs_E_%s_signal", vp.tag));
        g_mean_b->SetName(Form("centroid_vs_E_%s_background", vp.tag));
        g_sigma_b->SetName(Form("sigma_vs_E_%s_background", vp.tag));

        // señal: azul, circulo lleno | fondo: rojo, cuadrado hueco
        g_mean_s->SetMarkerStyle(20);  g_mean_s->SetMarkerColor(kAzure+2);  g_mean_s->SetLineColor(kAzure+2);
        g_mean_b->SetMarkerStyle(25);  g_mean_b->SetMarkerColor(kRed+1);   g_mean_b->SetLineColor(kRed+1);
        g_sigma_s->SetMarkerStyle(20); g_sigma_s->SetMarkerColor(kAzure+2); g_sigma_s->SetLineColor(kAzure+2);
        g_sigma_b->SetMarkerStyle(25); g_sigma_b->SetMarkerColor(kRed+1);  g_sigma_b->SetLineColor(kRed+1);

        // --- canvas centroide: señal + fondo superpuestos ---
        TCanvas* cm = new TCanvas(Form("c_centroid_%s", vp.tag), vp.tag, 900, 700);
        cm->SetLogx();
        double ymin_m = std::min(TMath::MinElement(ebins, mean_s.data()), TMath::MinElement(ebins, mean_b.data()));
        double ymax_m = std::max(TMath::MaxElement(ebins, mean_s.data()), TMath::MaxElement(ebins, mean_b.data()));
        double pad_m  = 0.1 * (ymax_m - ymin_m + 1e-9);
        g_mean_s->SetTitle(Form(";E_{n} (MeV);Centroid (%s)", vp.tag));
        g_mean_s->GetYaxis()->SetRangeUser(ymin_m - pad_m, ymax_m + pad_m);
        g_mean_s->Draw("AP");
        g_mean_b->Draw("P same");
        TLegend* leg_m = new TLegend(0.65, 0.75, 0.92, 0.90);
        leg_m->SetBorderSize(0);
        leg_m->AddEntry(g_mean_s, "signal", "lp");
        leg_m->AddEntry(g_mean_b, "background", "lp");
        leg_m->Draw();
        cm->SaveAs(Form("centroid_vs_E_%s.png", vp.tag));

        // --- canvas sigma: señal + fondo superpuestos ---
        TCanvas* cs = new TCanvas(Form("c_sigma_%s", vp.tag), vp.tag, 900, 700);
        cs->SetLogx();
        double ymin_s = std::min(TMath::MinElement(ebins, sigma_s.data()), TMath::MinElement(ebins, sigma_b.data()));
        double ymax_s = std::max(TMath::MaxElement(ebins, sigma_s.data()), TMath::MaxElement(ebins, sigma_b.data()));
        double pad_s  = 0.1 * (ymax_s - ymin_s + 1e-9);
        g_sigma_s->SetTitle(Form(";E_{n} (MeV);Sigma (%s)", vp.tag));
        g_sigma_s->GetYaxis()->SetRangeUser(std::max(0.0, ymin_s - pad_s), ymax_s + pad_s);
        g_sigma_s->Draw("AP");
        g_sigma_b->Draw("P same");
        TLegend* leg_s = new TLegend(0.65, 0.75, 0.92, 0.90);
        leg_s->SetBorderSize(0);
        leg_s->AddEntry(g_sigma_s, "signal", "lp");
        leg_s->AddEntry(g_sigma_b, "background", "lp");
        leg_s->Draw();
        cs->SaveAs(Form("sigma_vs_E_%s.png", vp.tag));

        fout->cd();
        g_mean_s->Write();  g_mean_b->Write();
        g_sigma_s->Write(); g_sigma_b->Write();
        cm->Write();
        cs->Write();
    }

    fout->Close();
    std::cout << "sums_study terminado. Output: sums_study_uranium.root" << std::endl;
}