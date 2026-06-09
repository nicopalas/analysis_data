#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers to load data
// ─────────────────────────────────────────────────────────────────────────────

// Load EXFOR cross-section CSV (cols[3]=DATA in barn*1e3→mb, cols[4]=err%,
// cols[11]=EN in eV→MeV)
TGraphErrors* load_exfor_cs(const char* csvfile)
{
    std::ifstream f(csvfile);
    std::string line;
    std::getline(f, line);   // skip header
    std::vector<double> x, y, ey;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')) cols.push_back(token);
        if (cols.size() < 12) continue;
        double en  = std::stod(cols[11]) / 1e6;          // eV → MeV
        double cs  = std::stod(cols[3])  * 1e3;          // barn → mb
        double err = cs * std::stod(cols[4]) / 100.0;
        x.push_back(en); y.push_back(cs); ey.push_back(err);
    }
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, ey.data());
    g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlue+1); g->SetLineColor(kBlue+1);
    return g;
}

// Load EXFOR anisotropy CSV (cols[3]=DATA NO-DIM, cols[4]=err%, cols[11]=EN eV→MeV)
TGraphErrors* load_exfor_aniso(const char* csvfile)
{
    std::ifstream f(csvfile);
    std::string line;
    std::getline(f, line);   // skip header
    std::vector<double> x, y, ey;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')) cols.push_back(token);
        if (cols.size() < 12) continue;
        double en   = std::stod(cols[7]) / 1e6;         // eV → MeV
        double val  = std::stod(cols[3]);                 // NO-DIM
        double err  = val * std::stod(cols[4]) / 100.0;
        x.push_back(en); y.push_back(val); ey.push_back(err);
    }
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, ey.data());
    g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlue+1); g->SetLineColor(kBlue+1);
    return g;
}

// Load cross section from ROOT file (h_cs_norm in barn → convert to mb)
TGraphErrors* load_cs_this_work()
{
    TFile* fmc = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/output_cross_section_uranium.root");
    if (!fmc || fmc->IsZombie()) { std::cerr << "[ERROR] Cannot open cs file\n"; return nullptr; }
    TH1D* h = (TH1D*)fmc->Get("h_cs_norm");
    if (!h) { std::cerr << "[ERROR] h_cs_norm not found\n"; fmc->Close(); return nullptr; }
    h->SetDirectory(0);
    fmc->Close();

    std::vector<double> x, y, ex, ey;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double sigma = h->GetBinContent(i) * 1e3;   // barn → mb
        if (sigma <= 0.) continue;
        x .push_back(h->GetBinCenter(i));
        y .push_back(sigma);
        ex.push_back(h->GetBinWidth(i) / 2.);
        ey.push_back(h->GetBinError(i) * 1e3);
    }
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerStyle(21); g->SetMarkerSize(0.7);
    g->SetMarkerColor(kRed+1); g->SetLineColor(kRed+1);
    delete h;
    return g;
}

TGraphErrors* load_aniso_this_work()
{
    TFile* fa = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/anisotropy_ratio_uranium_complete.root");
    if (!fa || fa->IsZombie()) { std::cerr << "[ERROR] Cannot open anisotropy file\n"; return nullptr; }
    TGraphErrors* g = (TGraphErrors*)fa->Get("anisotropy_ratio")->Clone();
    fa->Close();
    g->SetMarkerStyle(21); g->SetMarkerSize(0.7);
    g->SetMarkerColor(kRed+1); g->SetLineColor(kRed+1);
    return g;
}

// ─────────────────────────────────────────────────────────────────────────────
// Plot 1: cross section (left y-axis, mb) + anisotropy (right y-axis, NO-DIM)
// ─────────────────────────────────────────────────────────────────────────────
void plot_cs_with_anisotropy(const std::string& outdir)
{
    TGraphErrors* g_cs    = load_cs_this_work();
    TGraphErrors* g_aniso = load_aniso_this_work();
    if (!g_cs || !g_aniso) return;

    TCanvas* c = new TCanvas("c_cs_aniso", "Cross section + Anisotropy", 900, 600);
    c->SetLogx();

    // ── Left pad for cross section ────────────────────────────────────
    TPad* pad1 = new TPad("pad1", "", 0, 0, 1, 1);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->SetFillStyle(0);
    pad1->Draw();
    pad1->cd();

    TH1F* frame1 = pad1->DrawFrame(0.1, 100., 1000., 5000.);
    frame1->GetXaxis()->SetTitle("E_{n} (MeV)");
    frame1->GetYaxis()->SetTitle("#sigma (mb)");
    frame1->GetYaxis()->SetTitleColor(kRed+1);
    frame1->GetYaxis()->SetLabelColor(kRed+1);
    g_cs->Draw("P SAME");

    // ── Right pad for anisotropy ──────────────────────────────────────
    c->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0, 1, 1);
    pad2->SetLogx();
    pad2->SetFillStyle(0);
    pad2->SetFrameFillStyle(0);
    pad2->Draw();
    pad2->cd();

    // Find anisotropy range
    double aniso_min = 1e9, aniso_max = -1e9;
    for (int i = 0; i < g_aniso->GetN(); ++i) {
        double v = g_aniso->GetY()[i];
        double e = g_aniso->GetEY()[i];
        aniso_min = std::min(aniso_min, v - e);
        aniso_max = std::max(aniso_max, v + e);
    }
    double margin = (aniso_max - aniso_min) * 0.2;
    aniso_min -= margin; aniso_max += margin;

    TH1F* frame2 = pad2->DrawFrame(0.1, aniso_min, 1000., aniso_max);
    frame2->GetXaxis()->SetLabelSize(0);
    frame2->GetXaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetLabelSize(0);
    frame2->GetYaxis()->SetTickLength(0);

    // Draw anisotropy scaled to right axis
    g_aniso->Draw("P SAME");

    // Right axis
    TGaxis* axis_right = new TGaxis(1000., aniso_min, 1000., aniso_max,
                                    aniso_min, aniso_max, 510, "+L");
    axis_right->SetTitle("W(0^{#circ})/W(90^{#circ})");
    axis_right->SetTitleColor(kBlue+1);
    axis_right->SetLabelColor(kBlue+1);
    axis_right->SetLineColor(kBlue+1);
    axis_right->SetTitleOffset(1.2);
    axis_right->Draw();

    // ── Legend ────────────────────────────────────────────────────────
    TLegend* leg = new TLegend(0.55, 0.78, 0.88, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.032);
    leg->AddEntry(g_cs,    "^{238}U(n,f) #sigma — this work",         "lp");
    leg->AddEntry(g_aniso, "W(0^{#circ})/W(90^{#circ}) — this work",  "lp");
    leg->Draw();

    c->Update();
    c->SaveAs((outdir + "cs_with_anisotropy.pdf").c_str());

    TFile* fout = new TFile((outdir + "cs_with_anisotropy.root").c_str(), "RECREATE");
    g_cs->Write("g_cs");
    g_aniso->Write("g_aniso");
    fout->Close();
    std::cout << "[INFO] Plot 1 saved.\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Plot 2: cross section this work vs EXFOR
// ─────────────────────────────────────────────────────────────────────────────
void plot_cs_vs_exfor(const char* csvfile_cs, const std::string& outdir)
{
    TGraphErrors* g_exfor = load_exfor_cs(csvfile_cs);
    TGraphErrors* g_this  = load_cs_this_work();
    if (!g_exfor || !g_this) return;

    TCanvas* c = new TCanvas("c_cs_exfor", "Cross section vs EXFOR", 900, 600);
    c->SetLogx();
    c->SetLogy();

    TH1F* frame = c->DrawFrame(0.1, 100., 1000., 5000.);
    frame->GetXaxis()->SetTitle("E_{n} (MeV)");
    frame->GetYaxis()->SetTitle("#sigma (mb)");

    g_exfor->Draw("P SAME");
    g_this->Draw("P SAME");

    TLegend* leg = new TLegend(0.62, 0.78, 0.92, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.032);
    leg->AddEntry(g_exfor, "EXFOR 41756",          "lp");
    leg->AddEntry(g_this,  "This work",            "lp");
    leg->Draw();

    c->Update();
    c->SaveAs((outdir + "cs_vs_exfor.pdf").c_str());

    TFile* fout = new TFile((outdir + "cs_vs_exfor.root").c_str(), "RECREATE");
    g_exfor->Write("g_exfor_cs");
    g_this->Write("g_cs_this_work");
    fout->Close();
    std::cout << "[INFO] Plot 2 saved.\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Plot 3: anisotropy this work vs EXFOR
// ─────────────────────────────────────────────────────────────────────────────
void plot_aniso_vs_exfor(const char* csvfile_aniso, const std::string& outdir)
{
    TGraphErrors* g_exfor = load_exfor_aniso(csvfile_aniso);
    TGraphErrors* g_this  = load_aniso_this_work();
    if (!g_exfor || !g_this) return;

    // Find y range across both graphs
    double ymin = 1e9, ymax = -1e9;
    auto update_range = [&](TGraphErrors* g) {
        for (int i = 0; i < g->GetN(); ++i) {
            double v = g->GetY()[i], e = g->GetEY()[i];
            ymin = std::min(ymin, v - e);
            ymax = std::max(ymax, v + e);
        }
    };
    update_range(g_exfor);
    update_range(g_this);
    double margin = (ymax - ymin) * 0.15;
    ymin -= margin; ymax += margin;

    // Isotropic reference line
    TLine* line = new TLine(0.1, 1.0, 1000., 1.0);
    line->SetLineStyle(2); line->SetLineColor(kGray+1); line->SetLineWidth(1);

    TCanvas* c = new TCanvas("c_aniso_exfor", "Anisotropy vs EXFOR", 900, 600);
    c->SetLogx();

    TH1F* frame = c->DrawFrame(0.1, ymin, 1000., ymax);
    frame->GetXaxis()->SetTitle("E_{n} (MeV)");
    frame->GetYaxis()->SetTitle("W(0^{#circ})/W(90^{#circ})");

    line->Draw();
    g_exfor->Draw("P SAME");
    g_this->Draw("P SAME");

    TLegend* leg = new TLegend(0.55, 0.78, 0.92, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextFont(42);  leg->SetTextSize(0.032);
    leg->AddEntry(g_exfor, "EXFOR 41756",                              "lp");
    leg->AddEntry(g_this,  "This work",                               "lp");
    leg->AddEntry(line,    "Isotropic",                               "l");
    leg->Draw();

    c->Update();
    c->SaveAs((outdir + "aniso_vs_exfor.pdf").c_str());

    TFile* fout = new TFile((outdir + "aniso_vs_exfor.root").c_str(), "RECREATE");
    g_exfor->Write("g_exfor_aniso");
    g_this->Write("g_aniso_this_work");
    fout->Close();
    std::cout << "[INFO] Plot 3 saved.\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function — runs all three plots
// ─────────────────────────────────────────────────────────────────────────────
void plot_exfor(
    const char* csvfile_cs    = "/Users/nico/Downloads/417560032.csv",
    const char* csvfile_aniso = "/Users/nico/Downloads/41756002 (2).csv",
    const char* outdir        = "/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/")
{
    std::string od(outdir);
    plot_cs_with_anisotropy(od);
    plot_cs_vs_exfor(csvfile_cs, od);
    plot_aniso_vs_exfor(csvfile_aniso, od);
    std::cout << "[DONE] All three plots saved to " << od << "\n";
}