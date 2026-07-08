#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"

void cross_section_toy(){

    // ── input data ────────────────────────────────────────────────────────────
    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*)fin->Get("events_uranium");
    if (!tin) { std::cerr << "Tree not found\n"; return; }

    // ── fission cut ───────────────────────────────────────────────────────────
    TFile *fcut = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/fission.root", "READ");
    if (!fcut || fcut->IsZombie()) { std::cerr << "Cannot open cut file\n"; return; }
    TCutG *cut = (TCutG*)fcut->Get("fission");
    if (!cut) { std::cerr << "TCutG 'fission' not found\n"; return; }

    // ── flux ──────────────────────────────────────────────────────────────────
    // Read flux ONCE here, not inside the cross section loop
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

    // ── energy binning ────────────────────────────────────────────────────────
    const int nbins = 30;
    std::vector<double> energy_bins = buildLogBins(nbins, 1.0, 1000.0);
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // ── read flux into vectors ────────────────────────────────────────────────
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

    // ── branch addresses ──────────────────────────────────────────────────────
    double tof1, tof0, neutron_energy;
    double  amp0, amp1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("tof1",           &tof1);
    tin->SetBranchAddress("tof0",           &tof0);
    tin->SetBranchAddress("amp0",           &amp0);
    tin->SetBranchAddress("amp1",           &amp1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);
    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);

    // ── count events inside TCutG per energy bin ──────────────────────────────
    std::vector<double> counts(nbins, 0.0);

    Long64_t nentries = tin->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tin->GetEntry(i);
        if (neutron_energy < 1.0 || neutron_energy > 1000.0) continue;
        if (std::fabs(cos_theta)>1 || std::fabs(cos_theta_det)>1) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;
        double sum_amp = amp1;

        if (amp0>8000)
            counts[bin]++;
    }

    // ── cross section ─────────────────────────────────────────────────────────
    // sigma = counts / (Phi * N_atoms * barn)
    // Poisson uncertainty on counts: u_counts = sqrt(counts)
    const double barn    = 1.0e-24;
    const double N_atoms = 5.05e17;   // adjust to your target

    TH1D *counts_hist = new TH1D("counts_cut", "", nbins, 1.0, 1000.);
    for (int i = 0; i<nbins; i++){
        counts_hist->Fill(counts[i]);
    }
    std::cout << "Cheguei 7" << std::endl;
    std::vector<double> counts_cut(nbins, 0.0);
    for (int i = 0; i<nentries; i++){
        tin->GetEntry(i);
        double dt = tof1-tof0;
        double amp_sum = amp0+amp1;
        if (neutron_energy < 1.0 || neutron_energy > 1000.0) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;
        if (amp0<8000) continue;
        counts_cut[bin]++;
    }
    TH1D *counts_hist_cut = new TH1D("counts_cut", "", nbins, 1.0, 1000.);
    for (int i = 0; i<nbins; i++){
        counts_hist_cut->Fill(counts_cut[i]);
    }
    fin->Close();



    std::vector<double> x_cs, y_cs, ex_cs, ey_cs;
    std::vector<double> x_cs_cut, y_cs_cut, ex_cs_cut, ey_cs_cut;

    std::cout << "Cheguei 6" << std::endl;
    std::cout << "\n=== Cross section ===\n";
    for (int e = 0; e < nbins; ++e) {
        double Phi   = flux[e];
        double u_Phi = u_flux[e];
        if (Phi <= 0.0 || counts[e] <= 0.0) continue;

        double c    = counts[e];
        double c_cut = counts_cut[e];
        double uc   = std::sqrt(c);   // Poisson
        double uc_cut = std::sqrt(counts_cut[e]);

        double sigma   = c / (Phi * N_atoms * barn);
        double sigma_cut =  c_cut / (Phi * N_atoms * barn);
        double u_sigma = sigma * std::sqrt(
            (uc*uc)/(c*c) + (u_Phi*u_Phi)/(Phi*Phi));
        double u_sigma_cut = sigma_cut * std::sqrt(
            (uc_cut*uc_cut)/(c_cut*c_cut) + (u_Phi*u_Phi)/(Phi*Phi));

        double Ec = std::sqrt(E_low[e] * E_high[e]);
        printf("ebin=%2d  E=%.3f MeV  counts=%.0f  Phi=%.4e  sigma=%.4e +/- %.4e barn\n",
               e, Ec, c, Phi, sigma, u_sigma);

        x_cs.push_back(Ec);
        y_cs.push_back(sigma);
        ex_cs.push_back((E_high[e] - E_low[e]) / 2.0);
        ey_cs.push_back(u_sigma);
        x_cs_cut.push_back(Ec);
        y_cs_cut.push_back(sigma_cut);
        ex_cs_cut.push_back((E_high[e] - E_low[e]) / 2.0);
        ey_cs_cut.push_back(u_sigma_cut);
    }

std::cout << "Cheguei 5" << std::endl;

    // ── save ──────────────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy.root",
        "RECREATE");
    TGraphErrors *g = new TGraphErrors(
        (int)x_cs.size(),
        x_cs.data(), y_cs.data(),
        ex_cs.data(), ey_cs.data());
    g->SetName("g_cs_toy");
    g->SetTitle("U-238 cross section (toy);E_{n} (MeV);#sigma (barn)");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+1);
    g->SetLineColor(kBlue+1);
    g->Write();

    TGraphErrors *g_cut = new TGraphErrors(
        (int)x_cs.size(),
        x_cs_cut.data(), y_cs_cut.data(),
        ex_cs_cut.data(), ey_cs_cut.data());
    g_cut->SetName("g_cs_cut_toy");
    g_cut->SetTitle("U-238 cross section (toy);E_{n} (MeV);#sigma (barn)");
    g_cut->SetMarkerStyle(20);
    g_cut->SetMarkerColor(kBlue+1);
    g_cut->SetLineColor(kBlue+1);
    g_cut->Write();
std::cout << "Cheguei 4" << std::endl;
TH1D* h_cs_raw = new TH1D("h_cs_raw", "", nbins, energy_bins.data());
for (int e = 0; e < nbins; ++e){
    h_cs_raw->SetBinContent(e+1, y_cs[e]);
    h_cs_raw->SetBinError  (e+1, ey_cs[e]);
}
TH1D* h_cs_raw_cut = new TH1D("h_cs_raw_cut", "", nbins, energy_bins.data());
for (int e = 0; e < nbins; ++e){
    h_cs_raw_cut->SetBinContent(e+1, y_cs_cut[e]);
    h_cs_raw_cut->SetBinError  (e+1, ey_cs_cut[e]);
}
std::cout << "Cheguei 3" << std::endl;
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
std::cout << "Cheguei 2" << std::endl;
TH1D* h_cs_norm2 = (TH1D*)h_cs_raw_cut->Clone("h_cs_norm2");
h_cs_norm2->Reset();
normalised_xs(h_cs_raw_cut, 8.0, 10.0, 2.006, h_cs_norm2);
TCanvas* c_cs2 = new TCanvas("c_cs2", "Cross section U-238", 900, 600);
c_cs2->SetLogx();
h_cs_norm2->SetLineColor(kBlue+1);
h_cs_norm2->SetMarkerColor(kBlue+1);
h_cs_norm2->SetMarkerStyle(20);
h_cs_raw->Write();
h_cs_norm2->GetXaxis()->SetTitle("E_{n} (MeV)");
h_cs_norm2->GetYaxis()->SetTitle("#sigma (barn)");
h_cs_norm2->SetTitle("#sigma_{norm}(U-238)");
h_cs_norm2->Draw("E");
std::cout << "Cheguei 1" << std::endl;
fout->cd();
h_cs_norm->Write();
h_cs_norm2->Write();
c_cs->Write();
fout->Close();
std::cout << "[DONE] cs_toy.root\n";
}