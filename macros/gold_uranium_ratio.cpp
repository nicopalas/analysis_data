#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"
#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/cross_section.h"

void gold_uranium_ratio(){

    // ── open the two histograms produced by gold_xs() and uranium_xs() ────────
    TFile *f_gold = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_gold.root", "READ");
    if (!f_gold || f_gold->IsZombie()) { std::cerr << "Cannot open gold output file\n"; return; }
    TH1D *h_gold = (TH1D*)f_gold->Get("h_cs_raw");
    if (!h_gold) { std::cerr << "h_cs_raw (gold) not found\n"; return; }
    h_gold = (TH1D*)h_gold->Clone("h_cs_raw_gold_clone");
    h_gold->SetDirectory(nullptr);
    f_gold->Close();

    TFile *f_u = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_uranium.root", "READ");
    if (!f_u || f_u->IsZombie()) { std::cerr << "Cannot open uranium output file\n"; return; }
    TH1D *h_u = (TH1D*)f_u->Get("h_cs_raw_u");
    if (!h_u) { std::cerr << "h_cs_raw_u (uranium) not found\n"; return; }
    h_u = (TH1D*)h_u->Clone("h_cs_raw_u_clone");
    h_u->SetDirectory(nullptr);
    f_u->Close();

    if (h_gold->GetNbinsX() != h_u->GetNbinsX()) {
        std::cerr << "Gold and uranium histograms have different binning\n"; return;
    }
    const int nbins = h_gold->GetNbinsX();

    // ── raw bin-by-bin ratio Au/U (unscaled) ───────────────────────────────────
    TH1D *h_ratio_raw = (TH1D*)h_gold->Clone("h_ratio_raw");
    h_ratio_raw->SetDirectory(nullptr);
    h_ratio_raw->Reset();

    for (int e = 1; e <= nbins; ++e) {
        double g   = h_gold->GetBinContent(e);
        double ug  = h_gold->GetBinError(e);
        double u   = h_u->GetBinContent(e);
        double uu  = h_u->GetBinError(e);
        if (g <= 0.0 || u <= 0.0) continue;

        double r   = g / u;
        double ur  = r * std::sqrt((ug*ug)/(g*g) + (uu*uu)/(u*u));

        h_ratio_raw->SetBinContent(e, r);
        h_ratio_raw->SetBinError  (e, ur);
    }

    // ── reference gold data (Au-197), same as in gold_xs() ─────────────────────
    std::vector<double> E_ref_au   = {46.3, 66.6, 73.9, 94.1, 132.9, 144.6, 173.3};
    std::vector<double> sig_ref_au = {0.103, 0.81, 1.20, 2.81, 6.1,   8.1,   10.3}; // mb
    std::vector<double> u_ref_au   = {0.019, 0.12, 0.17, 0.39, 0.9,   1.2,   1.6};  // mb
    for (auto &s : sig_ref_au) s *= 1e-3;   // mb -> barn
    for (auto &u : u_ref_au)   u *= 1e-3;

    // ── reference uranium data (U-238(n,f)), full ENDF-style table ─────────────
    std::vector<double> E_ref_u = {
        0.50,0.52,0.54,0.57,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.94,0.96,0.98,
        1.00,1.10,1.25,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.60,4.00,
        4.50,4.70,5.00,5.30,5.50,5.80,6.00,6.20,6.50,7.00,7.50,7.75,8.00,8.50,
        9.00,10.00,11.00,11.50,12.00,13.00,14.00,14.50,15.00,16.00,17.00,18.00,
        19.00,20.00,21.00,22.00,23.00,24.00,25.00,26.00,27.00,28.00,29.00,30.00,
        32.00,34.00,36.00,38.00,40.00,42.00,44.00,46.00,48.00,50.00,52.00,54.00,
        56.00,58.00,60.00,64.00,68.00,72.00,76.00,80.00,84.00,88.00,92.00,96.00,
        100.00,104.00,108.00,112.00,116.00,120.00,128.00,136.00,144.00,152.00,
        160.00,168.00,176.00,184.00,192.00,200.00,300.00,400.00,500.00,600.00,
        700.00,800.00,900.00,1000.00
    };
    std::vector<double> sig_ref_u = {
        0.00026980,0.00067119,0.00059598,0.00065116,0.00116531,0.00129202,0.00184360,
        0.00264695,0.00450784,0.00675016,0.01378015,0.01690077,0.01550547,0.01588961,
        0.01410576,0.02895309,0.03320903,0.18619023,0.41835630,0.48206710,0.53553997,
        0.54749355,0.54624349,0.54151562,0.53708920,0.52481337,0.54752327,0.55451845,
        0.55921914,0.55980956,0.54833662,0.55240729,0.54820060,0.56779697,0.61392195,
        0.68694589,0.82159525,0.95004294,0.99851630,0.99569271,1.01920402,1.01472289,
        1.01452311,1.01117324,1.00882950,1.00505332,0.98752794,1.03156514,1.14994058,
        1.19263970,1.24300610,1.32225423,1.32703520,1.32162220,1.35796397,1.40707482,
        1.51901925,1.55512900,1.60213378,1.54300822,1.57505132,1.58626925,1.56289817,
        1.62440482,1.61076480,1.65532032,1.70277087,1.69981759,1.65676440,1.63970023,
        1.66526066,1.65841277,1.65557833,1.68089904,1.64272747,1.61694393,1.64917876,
        1.61998653,1.63154717,1.62447590,1.58763664,1.57865896,1.52995539,1.50445837,
        1.51984211,1.49515595,1.49600313,1.44827525,1.40759496,1.39989479,1.42732347,
        1.38255921,1.42826710,1.35680338,1.36360367,1.33475089,1.32909427,1.29432732,
        1.28907594,1.32032323,1.31550749,1.29939456,1.33658872,1.34154384,1.31700244,
        1.32466100,1.44100000,1.49700000,1.46100000,1.45800000,1.47500000,1.48200000,
        1.46800000,1.46700001
    };
    std::vector<double> u_ref_u = {
        0.00002049,0.00010438,0.00002249,0.00002441,0.00003694,0.00002922,0.00004329,
        0.00004879,0.00008038,0.00010401,0.00017266,0.00021365,0.00023027,0.00021990,
        0.00015131,0.00029119,0.00032341,0.00134738,0.00262918,0.00302871,0.00290432,
        0.00301741,0.00311783,0.00300113,0.00344929,0.00290484,0.00298915,0.00307938,
        0.00360237,0.00358805,0.00362927,0.00400621,0.00365471,0.00415490,0.00463600,
        0.00555936,0.00612914,0.00688432,0.00719061,0.00859903,0.00773614,0.00729852,
        0.00729594,0.00815660,0.00806416,0.00962545,0.00864000,0.00748159,0.00674434,
        0.00570974,0.00834980,0.01067722,0.01146856,0.01263919,0.01082061,0.01487632,
        0.01428978,0.02590652,0.01836317,0.02226853,0.01699061,0.02641331,0.02127290,
        0.02271004,0.02358701,0.01967756,0.02992725,0.02426689,0.02619886,0.03118292,
        0.02826156,0.02867969,0.02690437,0.03005004,0.02926416,0.03718286,0.03321005,
        0.03569087,0.04041503,0.03035458,0.03109424,0.02486835,0.02395078,0.02856810,
        0.02707488,0.03238298,0.03572723,0.03841870,0.03488576,0.03670170,0.03425680,
        0.03842244,0.04389951,0.04586037,0.04611334,0.04018981,0.03825088,0.03790471,
        0.03312034,0.03876146,0.03685174,0.04049074,0.03526862,0.05257003,0.03563506,
        0.03218347,0.07185272,0.10454333,0.08123140,0.08067812,0.08530600,0.06620452,
        0.06558076,0.06740687
    };
    // sig_ref_u, u_ref_u already in barn -- no conversion needed

    std::vector<double> ex_ref_au(E_ref_au.size(), 0.0);
    std::vector<double> ex_ref_u(E_ref_u.size(), 0.0);

    TGraphErrors *gr_ref_au = new TGraphErrors(
        (int)E_ref_au.size(), E_ref_au.data(), sig_ref_au.data(),
        ex_ref_au.data(), u_ref_au.data());
    gr_ref_au->SetName("g_ref_au");

    TGraphErrors *gr_ref_u = new TGraphErrors(
        (int)E_ref_u.size(), E_ref_u.data(), sig_ref_u.data(),
        ex_ref_u.data(), u_ref_u.data());
    gr_ref_u->SetName("g_ref_u");

    // ── reference ratio at the anchor point: Au(46.3 MeV) / U(46.0 MeV) ────────
    // use the tabulated points directly (both are close to 46 MeV and already
    // exist explicitly in each table, so no interpolation needed)
    double ref_au_46 = sig_ref_au[0];   // 46.3 MeV -> 6.1 mb entry
    double ref_u_46  = sig_ref_u[73];   // 46.00 MeV -> 1.68089904 barn entry
    double ref_ratio_46 = ref_au_46 / ref_u_46;

    printf("Reference anchor: sigma_Au(46.3 MeV)=%.4e barn, sigma_U(46.0 MeV)=%.4e barn, ref_ratio=%.4e\n",
           ref_au_46, ref_u_46, ref_ratio_46);

    // ── find the experimental bin closest to 46 MeV ────────────────────────────
    int best_bin = -1;
    double best_dE = 1e9;
    for (int e = 1; e <= nbins; ++e) {
        if (h_ratio_raw->GetBinContent(e) <= 0.0) continue;
        double Ec = h_ratio_raw->GetBinCenter(e);
        double dE = std::fabs(Ec - 46.0);
        if (dE < best_dE) { best_dE = dE; best_bin = e; }
    }
    if (best_bin < 0) { std::cerr << "No valid experimental ratio bin found near 46 MeV\n"; return; }

    double exp_ratio_46 = h_ratio_raw->GetBinContent(best_bin);
    double Ec_anchor    = h_ratio_raw->GetBinCenter(best_bin);

    // ── compute and print the scale ────────────────────────────────────────────
    double scale = ref_ratio_46 / exp_ratio_46;
    printf("Anchor bin: Ec=%.2f MeV  exp_ratio=%.4e  ref_ratio=%.4e  SCALE=%.6e\n",
           Ec_anchor, exp_ratio_46, ref_ratio_46, scale);

    // ── apply scale to the full ratio histogram ────────────────────────────────
    TH1D *h_ratio_norm = (TH1D*)h_ratio_raw->Clone("h_ratio_norm");
    h_ratio_norm->SetDirectory(nullptr);
    for (int e = 1; e <= nbins; ++e) {
        h_ratio_norm->SetBinContent(e, h_ratio_norm->GetBinContent(e) * scale);
        h_ratio_norm->SetBinError  (e, h_ratio_norm->GetBinError(e)   * scale);
    }

    // ── build TGraphErrors from the normalised ratio ────────────────────────────
    std::vector<double> x_r, y_r, ex_r, ey_r;
    for (int e = 1; e <= nbins; ++e) {
        if (h_ratio_norm->GetBinContent(e) <= 0.0) continue;
        x_r.push_back(h_ratio_norm->GetBinCenter(e));
        y_r.push_back(h_ratio_norm->GetBinContent(e));
        ex_r.push_back(h_ratio_norm->GetBinWidth(e) / 2.0);
        ey_r.push_back(h_ratio_norm->GetBinError(e));
    }
    TGraphErrors *gr_ratio = new TGraphErrors(
        (int)x_r.size(), x_r.data(), y_r.data(), ex_r.data(), ey_r.data());
    gr_ratio->SetName("g_ratio_au_u");
    gr_ratio->SetTitle("#sigma_{Au}/#sigma_{U} (normalised);E_{n} (MeV);Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetMarkerColor(kBlue+1);
    gr_ratio->SetLineColor(kBlue+1);
    gr_ratio->SetLineWidth(2);

    // ── reference ratio curve, evaluated at the uranium table energies
    //    within the gold reference's covered range ───────────────────────────────
    std::vector<double> x_ref_r, y_ref_r;
    double au_min = E_ref_au.front(), au_max = E_ref_au.back();
    for (size_t k = 0; k < E_ref_au.size(); ++k) {
        double E = E_ref_au[k];
        if (E < au_min || E > au_max) continue;
        double u_val = gr_ref_u->Eval(E);
        double au_val  = sig_ref_au[k];
        if (u_val <= 0.0) continue;
        x_ref_r.push_back(E);
        y_ref_r.push_back(au_val / u_val);
    }
    TGraphErrors *gr_ratio_ref = new TGraphErrors((int)x_ref_r.size());
    for (size_t k = 0; k < x_ref_r.size(); ++k)
        gr_ratio_ref->SetPoint(k, x_ref_r[k], y_ref_r[k]);
    gr_ratio_ref->SetName("g_ratio_ref");
    gr_ratio_ref->SetTitle("Reference ratio;E_{n} (MeV);Ratio");
    gr_ratio_ref->SetMarkerStyle(21);
    gr_ratio_ref->SetMarkerColor(kRed+1);
    gr_ratio_ref->SetLineColor(kRed+1);
    gr_ratio_ref->SetLineWidth(2);

    // ── save and plot ────────────────────────────────────────────────────────
    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_ratio.root",
        "RECREATE");

    TCanvas *c = new TCanvas("c_ratio", "Au/U cross section ratio", 900, 600);
    gr_ratio->Draw("AP");
    gr_ratio->GetXaxis()->SetTitle("E_{n} (MeV)");
    gr_ratio->GetYaxis()->SetTitle("#sigma_{Au}/#sigma_{U}");
    gr_ratio_ref->Draw("P same");

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(gr_ratio,     "This work", "lp");
    leg->AddEntry(gr_ratio_ref, "Reference", "lp");
    leg->Draw();

    c->Write();
    gr_ratio->Write();
    gr_ratio_ref->Write();
    h_ratio_raw->Write();
    h_ratio_norm->Write();

    fout->Close();
    std::cout << "[DONE] cs_toy_ratio.root\n";
}