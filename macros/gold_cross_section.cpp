#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"

void gold_cross_section(){
    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*)fin->Get("events_gold");
    TTree *tin_u = (TTree*) fin->Get("events_uranium");
    if (!tin) { std::cerr << "Tree not found\n"; return; }

    TFile *flux_file = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/flux_data/evalFlux_prelim.root",
        "READ");
    if (!flux_file || flux_file->IsZombie()) {
        std::cerr << "Cannot open flux file\n"; return;
    }
    TH1D *hist_flux = (TH1D*)flux_file->Get("hEval_Abs");
    if (!hist_flux) { std::cerr << "Flux histogram not found\n"; return; }
    hist_flux = (TH1D*)hist_flux->Clone("h_flux_clone");
    hist_flux->SetDirectory(nullptr);
    flux_file->Close();

    // ── cuts ─────────────────────────────────────────────────────────────────
    TFile *fcut0 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/sum_amps_gold.root", "READ");
    if (!fcut0 || fcut0->IsZombie()) { std::cerr << "Cannot open cut0 file\n"; return; }
    TCutG *cut0 = (TCutG*)fcut0->Get("sum_amps");
    if (!cut0) { std::cerr << "TCutG gold0 not found\n"; return; }

    TFile *fcut1 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/amps_gold.root", "READ");
    if (!fcut1 || fcut1->IsZombie()) { std::cerr << "Cannot open cut1 file\n"; return; }
    TCutG *cut1 = (TCutG*)fcut1->Get("amps");
    if (!cut1) { std::cerr << "TCutG gold1 not found\n"; return; }



    // ── energy binning ────────────────────────────────────────────────────────
    const int nbins = 10;
    std::vector<double> energy_bins = buildLogBins(nbins, 50.0, 1000.0);
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // ── flux ─────────────────────────────────────────────────────────────────
    std::vector<double> flux(nbins, 0.0), u_flux(nbins, 0.0);
    for (int e = 0; e < nbins; ++e) {
        Int_t bin1 = hist_flux->FindBin(E_low[e]  * 1e6 + 1e-9);
        Int_t bin2 = hist_flux->FindBin(E_high[e] * 1e6 - 1e-9);
        bin1 = std::max(bin1, 1);
        bin2 = std::min(bin2, hist_flux->GetNbinsX());
        if (bin1 > bin2) continue;
        double integral = 0.0, err2 = 0.0;
        for (int b = bin1; b <= bin2; ++b) {
            double Ec  = hist_flux->GetBinCenter(b);
            double w   = hist_flux->GetBinWidth(b);
            double val = hist_flux->GetBinContent(b);
            double err = hist_flux->GetBinError(b);
            if (Ec <= 0.0) continue;
            integral += (val / Ec) * w;
            err2     += TMath::Power((err / Ec) * w, 2);
        }
        flux[e]   = integral;
        u_flux[e] = std::sqrt(err2);
    }
    delete hist_flux;

    // ── event loop ───────────────────────────────────────────────────────────
    double tof1, tof0, neutron_energy;
    double  amp0, amp1;
    tin->SetBranchAddress("tof1",           &tof1);
    tin->SetBranchAddress("tof0",           &tof0);
    tin->SetBranchAddress("amp0",           &amp0);
    tin->SetBranchAddress("amp1",           &amp1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);

    std::vector<double> counts(nbins, 0.0);
    Long64_t nentries = tin->GetEntries();

    for (Long64_t i = 0; i < nentries; i++){
        tin->GetEntry(i);
        if (neutron_energy < 10.0 || neutron_energy > 1000.0) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;

        double dt = tof1 - tof0;
        if (cut0->IsInside(amp0+amp1, dt) && cut1->IsInside(amp1, amp0))
            counts[bin]++;
    }
    fin->Close();

    // ── compute raw cross section ─────────────────────────────────────────────
    const double barn    = 1.0e-24;
    const double N_atoms = 9.2e17;

    TH1D *h_cs_raw = new TH1D("h_cs_raw", "", nbins, energy_bins.data());
    for (int e = 0; e < nbins; ++e) {
        double Phi   = flux[e];
        double u_Phi = u_flux[e];
        if (Phi <= 0.0 || counts[e] <= 0.0) continue;

        double c       = counts[e];
        double uc      = std::sqrt(c);
        double sigma   = c / (Phi * N_atoms * barn);
        double u_sigma = sigma * std::sqrt((uc*uc)/(c*c) + (u_Phi*u_Phi)/(Phi*Phi));

        h_cs_raw->SetBinContent(e+1, sigma);
        h_cs_raw->SetBinError  (e+1, u_sigma);

        double Ec = std::sqrt(E_low[e] * E_high[e]);
        printf("ebin=%2d  E=%.3f MeV  counts=%.0f  sigma=%.4e +/- %.4e barn\n",
               e, Ec, c, sigma, u_sigma);
    }

    // ── reference gold data ───────────────────────────────────────────────────
    std::vector<double> E_ref  = {46.3, 66.6, 73.9, 94.1, 132.9, 144.6, 173.3};
    std::vector<double> sig_ref= {0.103, 0.81, 1.20, 2.81, 6.1,   8.1,   10.3};
    std::vector<double> u_ref  = {0.019, 0.12, 0.17, 0.39, 0.9,   1.2,   1.6};
    std::vector<double> ex_ref(E_ref.size(), 0.0);
    for (auto &s : sig_ref) s *= 1e-3;   // mb → barn
    for (auto &u : u_ref)   u *= 1e-3;

    TGraphErrors *gr_ref = new TGraphErrors(
        (int)E_ref.size(),
        E_ref.data(), sig_ref.data(),
        ex_ref.data(), u_ref.data());
    gr_ref->SetName("g_gold_ref");
    gr_ref->SetTitle("Au-197 reference;E_{n} (MeV);#sigma (barn)");
    gr_ref->SetMarkerStyle(21);
    gr_ref->SetMarkerColor(kRed+1);
    gr_ref->SetLineColor(kRed+1);
    gr_ref->SetLineWidth(2);

    // ── normalise experimental to reference at first bin ─────────────────────
    // find first non-empty bin and scale exp to match gr_ref->Eval at that energy
    TH1D *h_cs_norm = (TH1D*)h_cs_raw->Clone("h_cs_norm");
    h_cs_norm->SetDirectory(nullptr);

    double scale = 1.0;
    for (int e = 0; e < nbins; ++e) {
        double exp_val = h_cs_raw->GetBinContent(e+1);
        if (exp_val <= 0.0) continue;
        double Ec      = h_cs_raw->GetBinCenter(e+1);
        double ref_val = gr_ref->Eval(Ec);
        if (ref_val <= 0.0) continue;
        scale = ref_val / exp_val;
        printf("Normalisation at first bin: Ec=%.1f MeV  exp=%.4e  ref=%.4e  scale=%.4e\n",
               Ec, exp_val, ref_val, scale);
        break;
    }

    for (int e = 1; e <= h_cs_norm->GetNbinsX(); ++e) {
        h_cs_norm->SetBinContent(e, h_cs_norm->GetBinContent(e) * scale);
        h_cs_norm->SetBinError  (e, h_cs_norm->GetBinError(e)   * scale);
    }

    // ── TGraphErrors from normalised histogram ────────────────────────────────
    std::vector<double> x_exp, y_exp, ex_exp, ey_exp;
    for (int e = 1; e <= h_cs_norm->GetNbinsX(); ++e) {
        if (h_cs_norm->GetBinContent(e) <= 0.0) continue;
        x_exp.push_back(h_cs_norm->GetBinCenter(e));
        y_exp.push_back(h_cs_norm->GetBinContent(e));
        ex_exp.push_back(h_cs_norm->GetBinWidth(e) / 2.0);
        ey_exp.push_back(h_cs_norm->GetBinError(e));
    }
    TGraphErrors *gr_exp = new TGraphErrors(
        (int)x_exp.size(),
        x_exp.data(), y_exp.data(),
        ex_exp.data(), ey_exp.data());
    gr_exp->SetName("g_cs_gold_norm");
    gr_exp->SetTitle("Au-197 this work (normalised);E_{n} (MeV);#sigma (barn)");
    gr_exp->SetMarkerStyle(20);
    gr_exp->SetMarkerColor(kBlue+1);
    gr_exp->SetLineColor(kBlue+1);
    gr_exp->SetLineWidth(2);

    // ── save and plot ─────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_gold.root",
        "RECREATE");

    TCanvas *c = new TCanvas("c_gold", "Au-197 cross section", 900, 600);
    c->SetLogy();
    gr_exp->Draw("AP");
    gr_exp->GetXaxis()->SetTitle("E_{n} (MeV)");
    gr_exp->GetYaxis()->SetTitle("#sigma (barn)");
    gr_exp->SetTitle("Au-197 cross section");

    gr_ref->Draw("P same");
    TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.88);
    leg->AddEntry(gr_exp, "This work", "lp");
    leg->AddEntry(gr_ref, "Reference", "lp");
    leg->Draw();

    c->Write();
    gr_ref->Write();
    gr_exp->Write();
    h_cs_raw->Write();
    h_cs_norm->Write();

    fout->Close();
    std::cout << "[DONE] cs_toy_gold.root\n";
}