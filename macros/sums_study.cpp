// sums_vs_run.cpp — sum distributions per detector vs run number
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include <unordered_set>

static const int NDET = 10;

bool FitFWTM(TH1D* h, double& mean, double& sigma, double& mean_err, double& sigma_err)
{
    if (!h || h->GetEntries() < 50) return false;
    int nbins_h = h->GetNbinsX();
    int maxbin  = h->GetMaximumBin();
    double peak = h->GetBinCenter(maxbin);
    double maxval = h->GetBinContent(maxbin);
    if (maxval <= 0) return false;
    double tenth = maxval / 10.0;

    double left = h->GetBinCenter(1);
    for (int b = maxbin; b >= 1; b--) { if (h->GetBinContent(b) > tenth) continue; left = h->GetBinCenter(b); break; }
    double right = h->GetBinCenter(nbins_h);
    for (int b = maxbin; b <= nbins_h; b++) { if (h->GetBinContent(b) > tenth) continue; right = h->GetBinCenter(b); break; }
    if (right <= left) return false;

    TF1 f("fwtm_gaus", "gaus", left, right);
    f.SetParameter(1, peak);
    h->Fit(&f, "QRSN", "", left, right);
    mean = f.GetParameter(1);   sigma = f.GetParameter(2);
    mean_err = f.GetParError(1); sigma_err = f.GetParError(2);
    return (sigma > 0 && sigma < (right - left));
}

// slice a (run, sum) TH2 run by run, fit each slice, return centroid or sigma vs run
static TGraphErrors* trend(TH2D* h2, const char* name, const char* title, bool want_sigma)
{
    std::vector<double> x, y, ex, ey;
    for (int bx = 1; bx <= h2->GetNbinsX(); bx++) {
        TH1D* slice = h2->ProjectionY(Form("_px_%s_%d", name, bx), bx, bx);
        slice->SetDirectory(0);
        double m, s, me, se;
        if (FitFWTM(slice, m, s, me, se)) {
            x .push_back(h2->GetXaxis()->GetBinCenter(bx));
            ex.push_back(0.0);
            y .push_back(want_sigma ? s  : m);
            ey.push_back(want_sigma ? se : me);
        }
        delete slice;
    }
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(name); g->SetTitle(title);
    g->SetMarkerStyle(20); g->SetMarkerSize(0.6);
    return g;
}

void sums_study(const char* tree_name = "events_uranium",
                 double sum_lo = 80., double sum_hi = 150., int nsumbins = 140,
                 bool dedup = true)
{
    TFile* fin = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection_test.root");
    if (!fin || fin->IsZombie()) { std::cerr << "[ERROR] cannot open input file\n"; return; }
    TTree* tin = (TTree*) fin->Get(tree_name);
    if (!tin) { std::cerr << "[ERROR] tree " << tree_name << " not found\n"; fin->Close(); return; }
    if (!tin->GetBranch("RunNumber")) { std::cerr << "[ERROR] no RunNumber branch\n"; fin->Close(); return; }

    const int run_min = (int) tin->GetMinimum("RunNumber");
    const int run_max = (int) tin->GetMaximum("RunNumber");
    const int nruns   = run_max - run_min + 1;
    std::cout << "[INFO] " << tree_name << ": " << tin->GetEntries()
              << " entries, runs " << run_min << " to " << run_max << "\n";

    Int_t    RunNumber, coincID, det0, det1;
    Double_t sumX0, sumY0, sumX1, sumY1, neutron_energy;
    Float_t  amp0, amp1;   // amp<s> is booked as /D in events_selection — see note below

    tin->SetBranchAddress("RunNumber",      &RunNumber);
    tin->SetBranchAddress("coincID",        &coincID);
    tin->SetBranchAddress("det0",           &det0);
    tin->SetBranchAddress("det1",           &det1);
    tin->SetBranchAddress("sumX0",          &sumX0);
    tin->SetBranchAddress("sumY0",          &sumY0);
    tin->SetBranchAddress("sumX1",          &sumX1);
    tin->SetBranchAddress("sumY1",          &sumY1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);

    std::vector<TH2D*> h2X(NDET, nullptr), h2Y(NDET, nullptr);
    for (int d = 0; d < NDET; d++) {
        h2X[d] = new TH2D(Form("h2_sumX_det%d", d), Form("det %d;run number;sumX [ns]", d),
                          nruns, run_min - 0.5, run_max + 0.5, nsumbins, sum_lo, sum_hi);
        h2Y[d] = new TH2D(Form("h2_sumY_det%d", d), Form("det %d;run number;sumY [ns]", d),
                          nruns, run_min - 0.5, run_max + 0.5, nsumbins, sum_lo, sum_hi);
        h2X[d]->SetDirectory(0); h2Y[d]->SetDirectory(0);
    }

    // the same anode hit appears in the pairs (d-1,d) and (d,d+1); count it once
    std::unordered_set<unsigned long long> seen;
    auto key = [](int run, int cid, int det) {
        return ((unsigned long long)(unsigned)run << 36) ^
               ((unsigned long long)(unsigned)cid << 4)  ^ (unsigned long long)det;
    };

    Long64_t n = tin->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        tin->GetEntry(i);
        if (i % 200000 == 0) std::cout << "[INFO] " << i << "/" << n << "\r" << std::flush;

        const int  det[2]  = {det0, det1};
        const double sX[2] = {sumX0, sumX1};
        const double sY[2] = {sumY0, sumY1};

        for (int s = 0; s < 2; s++) {
            int d = det[s];
            if (d < 0 || d >= NDET) continue;
            if (sX[s] < -9000. || sY[s] < -9000.) continue;              // DUMMY
            if (dedup && !seen.insert(key(RunNumber, coincID, d)).second) continue;
            h2X[d]->Fill(RunNumber, sX[s]);
            h2Y[d]->Fill(RunNumber, sY[s]);
        }
    }
    std::cout << "\n";
    fin->Close();

    TFile* fout = new TFile(Form("sums_vs_run_%s.root", tree_name), "RECREATE");
    for (int d = 0; d < NDET; d++) {
        if (h2X[d]->GetEntries() == 0) continue;

        h2X[d]->Write(); h2Y[d]->Write();

        TH1D* pX = h2X[d]->ProjectionY(Form("h_sumX_det%d_all", d)); pX->Write();
        TH1D* pY = h2Y[d]->ProjectionY(Form("h_sumY_det%d_all", d)); pY->Write();

        TGraphErrors* gpx = trend(h2X[d], Form("g_peakX_vs_run_det%d",  d), Form("det %d;run number;sumX centroid [ns]", d), false);
        TGraphErrors* gsx = trend(h2X[d], Form("g_sigmaX_vs_run_det%d", d), Form("det %d;run number;sumX sigma [ns]",    d), true);
        TGraphErrors* gpy = trend(h2Y[d], Form("g_peakY_vs_run_det%d",  d), Form("det %d;run number;sumY centroid [ns]", d), false);
        TGraphErrors* gsy = trend(h2Y[d], Form("g_sigmaY_vs_run_det%d", d), Form("det %d;run number;sumY sigma [ns]",    d), true);
        gpx->Write(); gsx->Write(); gpy->Write(); gsy->Write();

        TCanvas* c = new TCanvas(Form("c_det%d", d), Form("det %d", d), 1400, 900);
        c->Divide(2, 2);
        c->cd(1); h2X[d]->Draw("COLZ");
        c->cd(2); h2Y[d]->Draw("COLZ");
        c->cd(3); gpx->SetMarkerColor(kAzure+2); gpx->SetLineColor(kAzure+2); gpx->Draw("AP");
                    gpx->GetYaxis()->SetRangeUser(100., 110.);
                  gpy->SetMarkerColor(kRed+1);   gpy->SetLineColor(kRed+1);
                  gpy->SetMarkerStyle(25);       gpy->Draw("P same");
        c->cd(4); gsx->SetMarkerColor(kAzure+2); gsx->SetLineColor(kAzure+2); gsx->Draw("AP");
                  gsy->SetMarkerColor(kRed+1);   gsy->SetLineColor(kRed+1);
                  gsy->SetMarkerStyle(25);       gsy->Draw("P same");
        c->SaveAs(Form("sums_vs_run_%s_det%d.png", tree_name, d));
        c->Write();

        std::cout << "det " << d << ": " << (Long64_t)h2X[d]->GetEntries()
                  << " hits, X centroid over " << gpx->GetN() << " runs\n";
    }
    fout->Close();
    std::cout << "[DONE] sums_vs_run_" << tree_name << ".root\n";
}