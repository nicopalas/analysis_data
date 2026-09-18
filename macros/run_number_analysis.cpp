#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/histograms.h"

// ---------------------------------------------------------------------------
// Robust trend: median (and 25/75 quantiles as a band) of each run-slice.
// FF amplitude spectra are broad and asymmetric, so a gaussian centroid is
// the wrong summary here — quantiles need no fit and cannot fail to converge.
// ---------------------------------------------------------------------------
static void quantile_trend(TH2D* h2, const char* base, const char* ytitle,
                           TGraphErrors*& g_med, TGraph*& g_q25, TGraph*& g_q75,
                           Long64_t min_counts = 200)
{
    std::vector<double> x, med, emed, q25, q75, ex;
    const double probs[3] = {0.25, 0.50, 0.75};
    double q[3];

    for (int bx = 1; bx <= h2->GetNbinsX(); bx++) {
        TH1D* s = h2->ProjectionY(Form("_p_%s_%d", base, bx), bx, bx);
        s->SetDirectory(0);
        if (s->GetEntries() >= min_counts) {
            s->GetQuantiles(3, q, (double*)probs);
            double n = s->GetEntries();
            // approximate error on the median for a roughly normal core
            double sigma_eff = (q[2] - q[0]) / 1.349;
            x   .push_back(h2->GetXaxis()->GetBinCenter(bx));
            ex  .push_back(0.0);
            med .push_back(q[1]);
            emed.push_back(1.2533 * sigma_eff / std::sqrt(n));
            q25 .push_back(q[0]);
            q75 .push_back(q[2]);
        }
        delete s;
    }

    g_med = new TGraphErrors(x.size(), x.data(), med.data(), ex.data(), emed.data());
    g_med->SetName(Form("g_%s_median", base));
    g_med->SetTitle(Form(";run number;%s", ytitle));
    g_med->SetMarkerStyle(20); g_med->SetMarkerSize(0.6);

    g_q25 = new TGraph(x.size(), x.data(), q25.data());
    g_q25->SetName(Form("g_%s_q25", base));
    g_q25->SetLineStyle(2); g_q25->SetMarkerStyle(1);

    g_q75 = new TGraph(x.size(), x.data(), q75.data());
    g_q75->SetName(Form("g_%s_q75", base));
    g_q75->SetLineStyle(2); g_q75->SetMarkerStyle(1);
}

// ---------------------------------------------------------------------------
void run_number_analysis(const char* tree_name = "events_gold",
                         double amp_lo = 0., double amp_hi = 50e3, int nampbins = 250,
                         double cos_min = 0.0)      // optional angular restriction
{
    TFile *fin = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection_test.root");
    if (!fin || fin->IsZombie()) { std::cerr << "[ERROR] cannot open input\n"; return; }
    TTree *tin = (TTree*) fin->Get(tree_name);
    if (!tin) { std::cerr << "[ERROR] tree " << tree_name << " not found\n"; fin->Close(); return; }
    if (!tin->GetBranch("RunNumber")) { std::cerr << "[ERROR] no RunNumber branch\n"; fin->Close(); return; }

    double  cos_theta, cos_theta_det;
    double  tof0, tof1;
    double  amp0, amp1;                  // /D in events_selection
    float   amp0_0, amp1_0;              // amp_0_<s>, booked /F — the reference
    int     RunNumber, det0, det1;

    tin->SetBranchAddress("amp0",          &amp0);
    tin->SetBranchAddress("amp1",          &amp1);
    tin->SetBranchAddress("amp_0_0",       &amp0_0);
    tin->SetBranchAddress("amp_0_1",       &amp1_0);
    tin->SetBranchAddress("tof0",          &tof0);
    tin->SetBranchAddress("tof1",          &tof1);
    tin->SetBranchAddress("cos_theta",     &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tin->SetBranchAddress("RunNumber",     &RunNumber);
    tin->SetBranchAddress("det0",          &det0);
    tin->SetBranchAddress("det1",          &det1);

    const int run_min = (int) tin->GetMinimum("RunNumber");
    const int run_max = (int) tin->GetMaximum("RunNumber");
    const int nruns   = run_max - run_min + 1;
    std::cout << "[INFO] " << tree_name << ": " << tin->GetEntries()
              << " entries, runs " << run_min << "-" << run_max << "\n";

    // one histogram per detector actually present in this tree
    std::map<int, TH2D*> h2_amp, h2_raw;
    auto book = [&](int d) {
        if (h2_amp.count(d)) return;
        h2_amp[d] = new TH2D(Form("h2_amp_det%d", d), Form("det %d calibrated;run number;amplitude", d),
                             nruns, run_min-0.5, run_max+0.5, nampbins, amp_lo, amp_hi);
        h2_raw[d] = new TH2D(Form("h2_raw_det%d", d), Form("det %d raw (amp_0);run number;amplitude", d),
                             nruns, run_min-0.5, run_max+0.5, nampbins, amp_lo, amp_hi);
        h2_amp[d]->SetDirectory(0); h2_raw[d]->SetDirectory(0);
    };

    // sum of the two fragments — the quantity your other cuts use
    TH2D* h2_sum = new TH2D("h2_ampsum", ";run number;amp0 + amp1",
                            nruns, run_min-0.5, run_max+0.5, nampbins, amp_lo, 2*amp_hi);
    h2_sum->SetDirectory(0);

    Long64_t n = tin->GetEntries(), nkept = 0;
    for (Long64_t i = 0; i < n; i++) {
        tin->GetEntry(i);
        if (i % 200000 == 0) std::cout << "[INFO] " << i << "/" << n << "\r" << std::flush;

        if (amp0 < -9000. || amp1 < -9000.) continue;
        if (cos_min > 0. && std::fabs(cos_theta_det) < cos_min) continue;

        const int    det[2] = {det0, det1};
        const double a  [2] = {amp0, amp1};
        const double a0 [2] = {(double)amp0_0, (double)amp1_0};

        for (int s = 0; s < 2; s++) {
            int d = det[s];
            if (d < 0 || d >= 10) continue;
            book(d);
            h2_amp[d]->Fill(RunNumber, a[s]);
            if (a0[s] > -9000.) h2_raw[d]->Fill(RunNumber, a0[s]);
        }
        h2_sum->Fill(RunNumber, amp0 + amp1);
        nkept++;
    }
    std::cout << "\n[INFO] kept " << nkept << " entries\n";
    fin->Close();

    TFile* fout = new TFile(Form("run_number_analysis_%s.root", tree_name), "RECREATE");

    for (auto& kv : h2_amp) {
        int d = kv.first;
        TH2D* hA = kv.second;
        TH2D* hR = h2_raw[d];
        if (hA->GetEntries() < 500) continue;

        hA->Write(); hR->Write();
        hA->ProjectionY(Form("h_amp_det%d_all", d))->Write();

        TGraphErrors *gA; TGraph *gA25, *gA75;
        TGraphErrors *gR; TGraph *gR25, *gR75;
        quantile_trend(hA, Form("amp_det%d", d), "amplitude (calibrated)", gA, gA25, gA75);
        quantile_trend(hR, Form("raw_det%d", d), "amplitude (raw)",        gR, gR25, gR75);
        gA->Write(); gA25->Write(); gA75->Write();
        gR->Write(); gR25->Write(); gR75->Write();

        // relative drift: median normalised to its own average over all runs
        auto relative = [&](TGraphErrors* g, const char* nm, int col) -> TGraph* {
            int N = g->GetN(); if (N == 0) return nullptr;
            double ref = TMath::Mean(N, g->GetY());
            if (ref == 0) return nullptr;
            std::vector<double> y(N);
            for (int k = 0; k < N; k++) y[k] = g->GetY()[k] / ref;
            TGraph* gr = new TGraph(N, g->GetX(), y.data());
            gr->SetName(nm); gr->SetTitle(Form("det %d;run number;median / <median>", d));
            gr->SetMarkerStyle(20); gr->SetMarkerSize(0.6);
            gr->SetMarkerColor(col); gr->SetLineColor(col);
            return gr;
        };
        TGraph* grA = relative(gA, Form("g_relamp_det%d", d), kAzure+2);
        TGraph* grR = relative(gR, Form("g_relraw_det%d", d), kRed+1);
        if (grA) grA->Write();
        if (grR) grR->Write();

        TCanvas* c = new TCanvas(Form("c_det%d", d), Form("det %d", d), 1400, 900);
        c->Divide(2, 2);
        c->cd(1); hA->Draw("COLZ");
        c->cd(2); hR->Draw("COLZ");
        c->cd(3);
            gA->SetMarkerColor(kAzure+2); gA->SetLineColor(kAzure+2);
            gA->Draw("AP");
            gA25->Draw("L same"); gA75->Draw("L same");
        c->cd(4);
            if (grA && grR) {
                gPad->DrawFrame(run_min-0.5, 0.85, run_max+0.5, 1.15)
                    ->SetTitle(Form("det %d;run number;median / <median>", d));
                grR->Draw("P same");   // raw = reference
                grA->Draw("P same");   // calibrated
                TLegend* l = new TLegend(0.62, 0.78, 0.92, 0.92);
                l->SetBorderSize(0);
                l->AddEntry(grA, "calibrated", "p");
                l->AddEntry(grR, "raw (amp_0)", "p");
                l->Draw();
            }
        c->SaveAs(Form("run_number_analysis_%s_det%d.png", tree_name, d));
        c->Write();

        std::cout << "det " << d << ": " << (Long64_t)hA->GetEntries() << " hits, "
                  << gA->GetN() << " runs with a median\n";
    }

    h2_sum->Write();
    TGraphErrors* gS; TGraph* gS25, *gS75;
    quantile_trend(h2_sum, "ampsum", "amp0 + amp1", gS, gS25, gS75);
    gS->Write(); gS25->Write(); gS75->Write();

    fout->Close();
    std::cout << "[DONE] run_number_analysis_" << tree_name << ".root\n";
}