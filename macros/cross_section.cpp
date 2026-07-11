#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/types.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/constants.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/acceptance.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/config.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"
#include <map>
#include <utility>

void gold_xs(){
    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root", "READ");
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
    TFile *fcut0 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/gold.root", "READ");
    if (!fcut0 || fcut0->IsZombie()) { std::cerr << "Cannot open cut0 file\n"; return; }
    TCutG *cut0 = (TCutG*)fcut0->Get("cut1");
    if (!cut0) { std::cerr << "TCutG gold0 not found\n"; return; }

    TCutG *cut1 = (TCutG*)fcut0->Get("cut2");
    if (!cut1) { std::cerr << "TCutG gold1 not found\n"; return; }

    // ── acceptance (solid angle per beam/det cell) ──────────────────────────
    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);

    // ── energy binning ────────────────────────────────────────────────────────
    const int nbins = 20;
    std::vector<double> energy_bins = buildLogBins(nbins, 40.0, 1000.0);
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // ── efficiency setup (gold) ─────────────────────────────────────────────
    // makeGoldConfig() fija cfg.efficiency_file y cfg.energy_bins_eff
    // ({30,300,600,1000}, 3 bins) -- distinto del binning fino (15 bins) de
    // esta sección eficaz. Igual que en anisotropy.C, cada uno de los 15
    // bins finos se mapea al bin de eficiencia que le corresponde.
    // NOTA: energy_bins_eff empieza en 30 MeV, pero esta sección eficaz
    // arranca en 10 MeV -- los bins de 10-30 MeV NO tienen eficiencia medida
    // y se dejan sin corregir (contribución nula a la sección eficaz), ver
    // aviso más abajo.
    AnalysisConfig cfg = makeGoldConfig(energy_bins);
    std::string eff_path = cfg.efficiency_file;
    TFile *eff_file = TFile::Open(eff_path.c_str(), "READ");
    if (!eff_file || eff_file->IsZombie()) { std::cerr << "Cannot open efficiency file: " << eff_path << "\n"; return; }

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

    // ── event loop (ahora también por ángulo de haz y de detector) ─────────
    double tof1, tof0, neutron_energy;
    float  amp0, amp1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("tof1",           &tof1);
    tin->SetBranchAddress("tof0",           &tof0);
    tin->SetBranchAddress("amp0",           &amp0);
    tin->SetBranchAddress("amp1",           &amp1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);
    tin->SetBranchAddress("cos_theta",      &cos_theta);
    tin->SetBranchAddress("cos_theta_det",  &cos_theta_det);

    Vec3D counts(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Long64_t nentries = tin->GetEntries();

    for (Long64_t i = 0; i < nentries; i++){
        tin->GetEntry(i);
        if (neutron_energy < 40.0 || neutron_energy > 1000.0) continue;
        if (std::fabs(cos_theta) > 1 || std::fabs(cos_theta_det) > 1) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;

        int bin_beam = int(std::fabs(cos_theta) / dcos_beam);
        int bin_det  = int(std::fabs(cos_theta_det) / dcos_det);
        if (bin_beam < 0 || bin_beam >= nbins_beam) continue;
        if (bin_det  < 0 || bin_det  >= nbins_det)  continue;

        double dt = tof1 - tof0;
        if (cut0->IsInside(amp0+amp1, dt) && cut1->IsInside((amp1-amp0)/(amp0+amp1), dt))
            counts[bin][bin_beam][bin_det]++;
    }
    fin->Close();

    // ── corregir cuentas por eficiencia y ángulo sólido ─────────────────────
    // N_corr(e) = Σ_{b,d} N(e,b,d) / ( Ω(b,d) · eps(e,d) )
    // -- estimador de máxima verosimilitud del nº de eventos físicos reales,
    // sumando sobre todas las celdas (b,d) con aceptancia y eficiencia
    // fiables. Propagación de error: contribución Poissoniana en N más
    // incertidumbre de la eficiencia, sumadas en cuadratura celda a celda
    // (las eps(e,d) de una misma fila de energía están correlacionadas por
    // venir del mismo bin de referencia en computeEfficiency -- esto asume
    // independencia, ver caveat de conversaciones previas).
    std::vector<double> counts_corr(nbins, 0.0), u_counts_corr(nbins, 0.0);
    std::map<int, std::pair<std::vector<double>, std::vector<double>>> eff_cache;
    const double eps_min = 0.6; // umbral mínimo de eficiencia fiable

    for (int e = 0; e < nbins; ++e) {
        double E_center = std::sqrt(E_low[e] * E_high[e]);
        int eff_bin = findBin(cfg.energy_bins_eff, E_center);
        if (eff_bin < 0 || eff_bin >= (int)cfg.energy_bins_eff.size() - 1) {
            std::cerr << "Warning: E=" << E_center
                      << " MeV fuera del rango de eficiencia medida (energy_bins_eff), "
                      << "bin " << e << " queda sin corregir (sigma=0)\n";
            continue;
        }

        if (eff_cache.find(eff_bin) == eff_cache.end()) {
            TH1D *h_eff = (TH1D*)eff_file->Get(Form("eff_ebin_%d", eff_bin));
            if (!h_eff) {
                std::cerr << "Warning: eff_ebin_" << eff_bin << " not found, bin " << e << " sin corregir\n";
                continue;
            }
            std::vector<double> eps(nbins_det, 0.0), u_eps(nbins_det, 0.0);
            for (int d = 0; d < nbins_det; d++) {
                eps[d]   = h_eff->GetBinContent(d + 1);
                u_eps[d] = h_eff->GetBinError(d + 1);
            }
            eff_cache[eff_bin] = {eps, u_eps};
        }
        const std::vector<double> &eps   = eff_cache[eff_bin].first;
        const std::vector<double> &u_eps = eff_cache[eff_bin].second;

        double N = 0.0, varN = 0.0;
        for (int b = 0; b < nbins_beam; ++b) {
            for (int d = 0; d < nbins_det; ++d) {
                double n  = counts[e][b][d];
                double om = acceptance[b][d];
                double ep = eps[d];
                if (om <= 1e-6 || ep <= eps_min || n <= 0.0) continue;

                double w = om * ep;
                N    += n / w;
                varN += n / (w * w);                                   // parte Poissoniana
                varN += std::pow(n * u_eps[d] / (om * ep * ep), 2);    // parte de incertidumbre de eficiencia
            }
        }
        counts_corr[e]   = N;
        u_counts_corr[e] = std::sqrt(varN);
    }
    eff_file->Close();

    // ── compute raw cross section (con cuentas corregidas) ─────────────────
    const double barn    = 1.0e-24;
    const double N_atoms = 9.2e17;

    TH1D *h_cs_raw = new TH1D("h_cs_raw", "", nbins, energy_bins.data());
    for (int e = 0; e < nbins; ++e) {
        double Phi   = flux[e];
        double u_Phi = u_flux[e];
        if (Phi <= 0.0 || counts_corr[e] <= 0.0) continue;

        double c       = counts_corr[e];
        double uc      = u_counts_corr[e];
        double sigma   = c / (Phi * N_atoms * barn);
        double u_sigma = sigma * std::sqrt((uc*uc)/(c*c) + (u_Phi*u_Phi)/(Phi*Phi));

        h_cs_raw->SetBinContent(e+1, sigma);
        h_cs_raw->SetBinError  (e+1, u_sigma);

        double Ec = std::sqrt(E_low[e] * E_high[e]);
        printf("ebin=%2d  E=%.3f MeV  N_corr=%.2f +/- %.2f  sigma=%.4e +/- %.4e barn\n",
               e, Ec, c, uc, sigma, u_sigma);
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
    TH1D *h_cs_norm = (TH1D*)h_cs_raw->Clone("h_cs_norm");
    h_cs_norm->SetDirectory(nullptr);

    double scale = 1.0;
    for (int e = 0; e < nbins; ++e) {
        double Ec      = h_cs_raw->GetBinCenter(e+1);
        double exp_val = h_cs_raw->GetBinContent(e+1);
        if (Ec<73.9) continue;
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



void uranium_xs(){
    TFile *fin = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root", "READ");
    if (!fin || fin->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }
    TTree *tin = (TTree*) fin->Get("events_uranium");
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
    TFile *fcut0 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/uranium.root", "READ");
    if (!fcut0 || fcut0->IsZombie()) { std::cerr << "Cannot open cut0 file\n"; return; }
    TCutG *cut0 = (TCutG*)fcut0->Get("cut1");
    if (!cut0) { std::cerr << "TCutG gold0 not found\n"; return; }

    TCutG *cut1 = (TCutG*)fcut0->Get("cut2");
    if (!cut1) { std::cerr << "TCutG gold1 not found\n"; return; }

    // ── acceptance (solid angle per beam/det cell) ──────────────────────────
    std::string acceptance_file = "/Users/nico/Desktop/Tese/Analysis/cross_section/data/acceptance_coincidence.csv";
    Vec2D acceptance, dOmega_fine;
    if (!loadAcceptanceCSV(acceptance_file, dOmega_fine)) {
        std::cerr << "Failed to load acceptance CSV" << std::endl;
        return;
    }
    acceptance = rebin(dOmega_fine);

    // ── energy binning ────────────────────────────────────────────────────────
    const int nbins = 20;
    std::vector<double> energy_bins = buildLogBins(nbins, 40.0, 1000.0);
    std::vector<double> E_low(nbins), E_high(nbins);
    for (int e = 0; e < nbins; ++e) {
        E_low[e]  = energy_bins[e];
        E_high[e] = energy_bins[e+1];
    }

    // ── efficiency setup (uranium) ──────────────────────────────────────────
    // NOTA / SUPUESTO: no tengo una makeUraniumConfig() a la vista en
    // config.h, así que reconstruyo aquí el binning de eficiencia tal como
    // se generó en el efficiency_toy.C original para uranio: 10 bins
    // logarítmicos entre 1 y 1000 MeV, fichero "efficiencies_u_toy.root".
    // Si tienes una makeUraniumConfig() real, dímelo y cambio esto por
    // cfg.energy_bins_eff / cfg.efficiency_file como en gold_xs().
    const int nbins_eff_u = 4;
    std::vector<double> energy_bins_eff_u = {1, 10, 100, 500, 1000};
    std::string eff_path_u = "/Users/nico/Desktop/Tese/Analysis/cross_section/efficiencies_u_toy.root";
    TFile *eff_file = TFile::Open(eff_path_u.c_str(), "READ");
    if (!eff_file || eff_file->IsZombie()) { std::cerr << "Cannot open efficiency file: " << eff_path_u << "\n"; return; }

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

    // ── event loop (ahora también por ángulo de haz y de detector) ─────────
    double tof1, tof0, neutron_energy;
    float  amp0, amp1;
    double cos_theta, cos_theta_det;
    tin->SetBranchAddress("tof1",           &tof1);
    tin->SetBranchAddress("tof0",           &tof0);
    tin->SetBranchAddress("amp0",           &amp0);
    tin->SetBranchAddress("amp1",           &amp1);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);
    tin->SetBranchAddress("cos_theta",      &cos_theta);
    tin->SetBranchAddress("cos_theta_det",  &cos_theta_det);

    Vec3D counts(nbins, Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    Long64_t nentries = tin->GetEntries();

    for (Long64_t i = 0; i < nentries; i++){
        tin->GetEntry(i);
        if (neutron_energy < 40.0 || neutron_energy > 1000.0) continue;
        if (std::fabs(cos_theta) > 1 || std::fabs(cos_theta_det) > 1) continue;

        int bin = findBin(energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;

        int bin_beam = int(std::fabs(cos_theta) / dcos_beam);
        int bin_det  = int(std::fabs(cos_theta_det) / dcos_det);
        if (bin_beam < 0 || bin_beam >= nbins_beam) continue;
        if (bin_det  < 0 || bin_det  >= nbins_det)  continue;

        double dt = tof1 - tof0;
        if (cut0->IsInside(amp0+amp1, dt) && cut1->IsInside((amp1-amp0)/(amp0+amp1), dt))
            counts[bin][bin_beam][bin_det]++;
    }
    fin->Close();

    // ── corregir cuentas por eficiencia y ángulo sólido ─────────────────────
    // Ver comentario detallado en gold_xs(); misma fórmula:
    // N_corr(e) = Σ_{b,d} N(e,b,d) / ( Ω(b,d) · eps(e,d) )
    std::vector<double> counts_corr(nbins, 0.0), u_counts_corr(nbins, 0.0);
    std::map<int, std::pair<std::vector<double>, std::vector<double>>> eff_cache;
    const double eps_min = 0.6;

    for (int e = 0; e < nbins; ++e) {
        double E_center = std::sqrt(E_low[e] * E_high[e]);
        int eff_bin = findBin(energy_bins_eff_u, E_center);
        if (eff_bin < 0 || eff_bin >= nbins_eff_u) {
            std::cerr << "Warning: E=" << E_center
                      << " MeV fuera del rango de eficiencia medida, bin " << e << " sin corregir\n";
            continue;
        }

        if (eff_cache.find(eff_bin) == eff_cache.end()) {
            TH1D *h_eff = (TH1D*)eff_file->Get(Form("eff_ebin_%d", eff_bin));
            if (!h_eff) {
                std::cerr << "Warning: eff_ebin_" << eff_bin << " not found, bin " << e << " sin corregir\n";
                continue;
            }
            std::vector<double> eps(nbins_det, 0.0), u_eps(nbins_det, 0.0);
            for (int d = 0; d < nbins_det; d++) {
                eps[d]   = h_eff->GetBinContent(d + 1);
                u_eps[d] = h_eff->GetBinError(d + 1);
            }
            eff_cache[eff_bin] = {eps, u_eps};
        }
        const std::vector<double> &eps   = eff_cache[eff_bin].first;
        const std::vector<double> &u_eps = eff_cache[eff_bin].second;

        double N = 0.0, varN = 0.0;
        for (int b = 0; b < nbins_beam; ++b) {
            for (int d = 0; d < nbins_det; ++d) {
                double n  = counts[e][b][d];
                double om = acceptance[b][d];
                double ep = eps[d];
                if (om <= 1e-6 || ep <= eps_min || n <= 0.0) continue;

                double w = om * ep;
                N    += n / w;
                varN += n / (w * w);
                varN += std::pow(n * u_eps[d] / (om * ep * ep), 2);
            }
        }
        counts_corr[e]   = N;
        u_counts_corr[e] = std::sqrt(varN);
    }
    eff_file->Close();

    // ── compute raw cross section (con cuentas corregidas) ─────────────────
    const double barn    = 1.0e-24;
    const double N_atoms = 6.67e17;

    TH1D *h_cs_raw = new TH1D("h_cs_raw_u", "", nbins, energy_bins.data());
    for (int e = 0; e < nbins; ++e) {
        double Phi   = flux[e];
        double u_Phi = u_flux[e];
        if (Phi <= 0.0 || counts_corr[e] <= 0.0) continue;

        double c       = counts_corr[e];
        double uc      = u_counts_corr[e];
        double sigma   = c / (Phi * N_atoms * barn);
        double u_sigma = sigma * std::sqrt((uc*uc)/(c*c) + (u_Phi*u_Phi)/(Phi*Phi));

        h_cs_raw->SetBinContent(e+1, sigma);
        h_cs_raw->SetBinError  (e+1, u_sigma);

        double Ec = std::sqrt(E_low[e] * E_high[e]);
        printf("ebin=%2d  E=%.3f MeV  N_corr=%.2f +/- %.2f  sigma=%.4e +/- %.4e barn\n",
               e, Ec, c, uc, sigma, u_sigma);
    }

    // ── reference uranium data ───────────────────────────────────────────────
    std::vector<double> E_ref  = {42, 44, 46, 48, 50};
    std::vector<double> sig_ref= {1.65841277, 1.65557833, 1.68089904, 1.64272747, 1.61694393};
    std::vector<double> u_ref  = {0.02867969, 0.02690437, 0.03005004, 0.02926416, 0.03718286};
    std::vector<double> ex_ref(E_ref.size(), 0.0);
    for (auto &s : sig_ref) s *= 1e-3;   // mb → barn
    for (auto &u : u_ref)   u *= 1e-3;

    TGraphErrors *gr_ref = new TGraphErrors(
        (int)E_ref.size(),
        E_ref.data(), sig_ref.data(),
        ex_ref.data(), u_ref.data());
    gr_ref->SetName("g_u_ref");
    gr_ref->SetTitle("U-238 reference;E_{n} (MeV);#sigma (barn)");
    gr_ref->SetMarkerStyle(21);
    gr_ref->SetMarkerColor(kRed+1);
    gr_ref->SetLineColor(kRed+1);
    gr_ref->SetLineWidth(2);

    // ── normalise experimental to reference at first bin ─────────────────────
    TH1D *h_cs_norm = (TH1D*)h_cs_raw->Clone("h_cs_norm_u");
    h_cs_norm->SetDirectory(nullptr);

    double scale = 1.0;
    for (int e = 0; e < nbins; ++e) {
        double Ec      = h_cs_raw->GetBinCenter(e+1);
        double exp_val = h_cs_raw->GetBinContent(e+1);
        if (Ec<46.0) continue;
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
    gr_exp->SetName("g_u_gold_norm");
    gr_exp->SetTitle("U-238 this work (normalised);E_{n} (MeV);#sigma (barn)");
    gr_exp->SetMarkerStyle(20);
    gr_exp->SetMarkerColor(kBlue+1);
    gr_exp->SetLineColor(kBlue+1);
    gr_exp->SetLineWidth(2);

    // ── save and plot ─────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_uranium.root",
        "RECREATE");

    // (corregido: nombre y título del canvas decían "gold"/"Au-197" por
    // copy-paste; ahora reflejan que esto es la sección eficaz del uranio)
    TCanvas *c = new TCanvas("c_uranium", "U-238 cross section", 900, 600);
    c->SetLogy();
    gr_exp->Draw("AP");
    gr_exp->GetXaxis()->SetTitle("E_{n} (MeV)");
    gr_exp->GetYaxis()->SetTitle("#sigma (barn)");
    gr_exp->SetTitle("U-238 cross section");

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
    std::cout << "[DONE] cs_toy_uranium.root\n";
}