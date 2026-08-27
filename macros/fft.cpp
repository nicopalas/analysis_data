// ============================================================================
//  Delay-line propagation-velocity measurement via FFT of the diff-time
//  spectra, for both samples and both PPACs, X and Y cathodes each.
//
//  Refinements over the original:
//    - zero-padding (finer frequency grid, better peak interpolation)
//    - parabolic 3-point sub-bin peak interpolation (central value)
//    - local Gaussian fit around the peak with per-bin errors estimated
//      from the spectrum's own noise floor -> statistical sigma on f0
//    - frequency-axis bin centers made to coincide exactly with the FFT
//      output frequencies (removes a half-bin systematic in f0)
//
//  Physics:  dt spectrum is in ns  ->  FFT frequency axis is in GHz.
//            v = pitch * f0   [mm/ns], with pitch the strip pitch.
// ============================================================================
#include <TVirtualFFT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace Setup {
    const double V_NOM_MMNS   = 1.0;      // nominal delay-line velocity [mm/ns]
    const double PITCH_MM     = 2.0;      // strip pitch [mm]

    const int    N_BINS       = 2000;     // dt histogram bins
    const double T_MIN        = -100.0;   // [ns]
    const double T_MAX        =  100.0;   // [ns]

    const int    PAD_FACTOR   = 1;        // FFT zero-padding factor
    const double AMP_PERCENTILE = 0.8;    // keep only top 10% amplitudes

    const double F_SEARCH_MIN = 0.3;      // peak search band [GHz]
    const double F_SEARCH_MAX = 0.8;
    const double F_NOISE_MIN  = 1.0;      // noise-floor estimation band [GHz]
    const double F_NOISE_MAX  = 6.0;
    const double FIT_HALFWIN_GHZ = 0.10;  // half-window of the local Gaussian fit

    const char*  FILE_PATH    =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root";

    const std::string OUTDIR  = "./";     // where PDFs go
}

// ---------------------------------------------------------------------------
// FFT with zero-padding. Frequency-axis bin centers are placed exactly at the
// discrete FFT frequencies f_i = i / (nfft * dt), so GetBinCenter() returns
// true frequencies with no half-bin offset.
// ---------------------------------------------------------------------------
static TH1D* computeFFT(TH1D* hInput, const std::string& name, int padFactor)
{
    int nData = hInput->GetNbinsX();
    int nfft  = padFactor * nData;
    double dt = hInput->GetBinWidth(1);

     // --- resta del nivel DC / baseline antes de la transformada ---
    double mean = 0.0;
    for (int i = 0; i < nData; ++i) mean += hInput->GetBinContent(i + 1);
    mean /= nData;

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &nfft, "R2C ES K");
    for (int i = 0; i < nData; ++i)
        fft->SetPoint(i, hInput->GetBinContent(i + 1) - mean);
    for (int i = nData; i < nfft; ++i) fft->SetPoint(i, 0.0);   // zero-pad
    fft->Transform();

    TH1* hMag = nullptr;
    hMag = TH1::TransformHisto(fft, hMag, "MAG");

    int    nHalf = nfft / 2;
    double df    = 1.0 / (nfft * dt);            // frequency resolution [GHz]
    double xlo   = -0.5 * df;                    // so bin i center = i*df
    double xhi   = (nHalf - 0.5) * df;

    TH1D* hFreq = new TH1D(name.c_str(),
        (name + ";Frequency (GHz);Magnitude (arb. units)").c_str(),
        nHalf, xlo, xhi);

    for (int i = 0; i < nHalf; ++i)
        hFreq->SetBinContent(i + 1, hMag->GetBinContent(i + 1));

    hFreq->SetBinContent(1, 0.0);   // suppress DC

    delete hMag;
    delete fft;
    return hFreq;
}

// find the bin of the maximum inside the search band
static int findPeakBin(TH1D* h, double fmin, double fmax)
{
    int b0 = h->FindBin(fmin);
    int b1 = h->FindBin(fmax);
    int kmax = b0; double vmax = 0.0;
    for (int i = b0; i <= b1; ++i)
        if (h->GetBinContent(i) > vmax) { vmax = h->GetBinContent(i); kmax = i; }
    return kmax;
}

// parabolic 3-point interpolation -> sub-bin peak frequency
static double parabolicPeak(TH1D* h, int kmax)
{
    if (kmax <= 1 || kmax >= h->GetNbinsX()) return h->GetBinCenter(kmax);
    double y0 = h->GetBinContent(kmax - 1);
    double y1 = h->GetBinContent(kmax);
    double y2 = h->GetBinContent(kmax + 1);
    double denom = (y0 - 2*y1 + y2);
    if (std::fabs(denom) < 1e-12) return h->GetBinCenter(kmax);
    double delta = 0.5 * (y0 - y2) / denom;          // in bin-width units
    return h->GetBinCenter(kmax) + delta * h->GetBinWidth(kmax);
}

// noise floor: RMS of magnitude in a control band away from the peak
static double estimateNoiseRMS(TH1D* h, double fmin, double fmax)
{
    int b0 = h->FindBin(fmin);
    int b1 = h->FindBin(fmax);
    double sum = 0, sum2 = 0; int n = 0;
    for (int i = b0; i <= b1; ++i) {
        double v = h->GetBinContent(i);
        sum += v; sum2 += v*v; ++n;
    }
    if (n < 2) return 1.0;
    double mean = sum / n;
    double var  = sum2 / n - mean*mean;
    return (var > 0) ? std::sqrt(var) : 1.0;
}

struct PeakResult {
    double f0        = 0;   // GHz (Gaussian fit, falls back to parabolic)
    double sigma_f0  = 0;   // GHz (statistical, from fit)
    double f0_parab  = 0;   // GHz (parabolic interpolation, cross-check)
    double v         = 0;   // mm/ns
    double sigma_v   = 0;   // mm/ns
    double period    = 0;   // ns
    double chi2ndf   = 0;
    bool   fitOk     = false;
};

// local Gaussian fit around the peak; per-bin error = noise RMS
static PeakResult fitPeak(TH1D* h, double fmin, double fmax,
                          double fnoise_min, double fnoise_max,
                          double halfwin_ghz)
{
    PeakResult r;
    int kmax     = findPeakBin(h, fmin, fmax);
    r.f0_parab   = parabolicPeak(h, kmax);

    double noise = estimateNoiseRMS(h, fnoise_min, fnoise_max);
    for (int i = 1; i <= h->GetNbinsX(); ++i) h->SetBinError(i, noise);

    double xc = h->GetBinCenter(kmax);
    TF1 g("gpeak", "[0]*TMath::Gaus(x,[1],[2]) + [3]",
          xc - halfwin_ghz, xc + halfwin_ghz);
    g.SetParameters(h->GetBinContent(kmax), xc, h->GetBinWidth(kmax)*3, 0.0);
    g.SetParLimits(2, h->GetBinWidth(kmax)*0.5, halfwin_ghz);

    int status = h->Fit(&g, "QRN");   // Quiet, Range, No-draw
    r.fitOk = (status == 0);

    if (r.fitOk) {
        r.f0       = g.GetParameter(1);
        r.sigma_f0 = g.GetParError(1);
        r.chi2ndf  = (g.GetNDF() > 0) ? g.GetChisquare()/g.GetNDF() : 0.0;
    } else {                          // fall back to parabolic if fit fails
        r.f0       = r.f0_parab;
        r.sigma_f0 = h->GetBinWidth(kmax) / std::sqrt(12.0);  // bin-width proxy
    }

    r.period  = (r.f0 > 0) ? 1.0 / r.f0 : 0.0;
    r.v       = Setup::PITCH_MM * r.f0;
    r.sigma_v = Setup::PITCH_MM * r.sigma_f0;
    return r;
}

// draw one FFT spectrum with peak + fit annotations
static void drawFFT(TH1D* h, const PeakResult& p, const std::string& title,
                    const std::string& outFile)
{
    TCanvas* c = new TCanvas(("c_" + outFile).c_str(), title.c_str(), 1000, 550);
    c->SetLeftMargin(0.11);
    c->SetBottomMargin(0.13);

    h->GetXaxis()->SetRangeUser(0.2, 2.0);
    h->SetLineColor(kAzure + 2);
    h->SetLineWidth(2);
    h->SetTitle((title + ";Frequency (GHz);Magnitude (arb. units)").c_str());
    h->Draw("HIST");

    double yMax = h->GetMaximum();
    TLine* l = new TLine(p.f0, 0, p.f0, yMax * 1.05);
    l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->SetLineStyle(2);
    l->Draw();

    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.033);
    tex.SetTextColor(kRed + 1);
    tex.DrawLatex(0.55, 0.84, Form("f_{0} = %.4f #pm %.4f GHz", p.f0, p.sigma_f0));
    tex.DrawLatex(0.55, 0.77, Form("v = %.4f #pm %.4f mm/ns", p.v, p.sigma_v));
    tex.SetTextColor(kBlack);
    tex.DrawLatex(0.55, 0.70, Form("T = %.4f ns", p.period));
    tex.DrawLatex(0.55, 0.63, Form("parabolic f_{0} = %.4f GHz", p.f0_parab));
    tex.DrawLatex(0.55, 0.56, Form("#chi^{2}/ndf = %.2f", p.chi2ndf));
    tex.DrawLatex(0.55, 0.49, Form("(pitch = %.1f mm)", Setup::PITCH_MM));

    c->SaveAs((Setup::OUTDIR + outFile).c_str());
    std::cout << "  Saved: " << outFile << "\n";
}

// ---------------------------------------------------------------------------
// per-sample configuration
// ---------------------------------------------------------------------------
struct SampleConfig {
    std::string tree_name;
    std::string ppac_label[2];   // label for PPAC index 0 and 1
    std::string tag;             // short tag for output filenames
};

// process one sample: 2 PPACs x {X,Y} cathodes = 4 velocities
static void processSample(TFile* f, const SampleConfig& cfg)
{
    std::cout << "\n================================================\n";
    std::cout << "  SAMPLE: " << cfg.tree_name << "\n";
    std::cout << "================================================\n";

    TTree* tree = (TTree*) f->Get(cfg.tree_name.c_str());
    if (!tree) { std::cerr << "  ERROR: tree not found\n"; return; }

    double diffX0=0, diffX1=0, diffY0=0, diffY1=0;
    double sumX0=0, sumX1=0, sumY0=0, sumY1=0;
    double  amp0=0, amp1=0;
    double neutron_energy;
    tree->SetBranchAddress("diffX0", &diffX0);
    tree->SetBranchAddress("diffX1", &diffX1);
    tree->SetBranchAddress("diffY0", &diffY0);
    tree->SetBranchAddress("diffY1", &diffY1);
    tree->SetBranchAddress("sumX0", &sumX0);
    tree->SetBranchAddress("sumX1", &sumX1);
    tree->SetBranchAddress("sumY0", &sumY0);
    tree->SetBranchAddress("sumY1", &sumY1);
    tree->SetBranchAddress("amp0",   &amp0);
    tree->SetBranchAddress("amp1",   &amp1);
    tree->SetBranchAddress("neutron_energy", &neutron_energy);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "  Entries: " << nEntries << "\n";

    // --- pass 1: amplitude percentiles (per PPAC) ---
    std::vector<float> a0, a1;
    a0.reserve(nEntries); a1.reserve(nEntries);
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        float rat_amp = (amp1-amp0)/(amp0+amp1);
        if (sumX0<90 || sumX1<90 || sumY0<90 || sumY1<90) continue;
        if (neutron_energy>1000 || rat_amp>0.3 || amp0+amp1<20000) continue;
        if (amp0 > 0) a0.push_back(amp0);
        if (amp1 > 0) a1.push_back(amp1);
    }
    std::sort(a0.begin(), a0.end());
    std::sort(a1.begin(), a1.end());
    auto pct = [](std::vector<float>& v, double p) -> float {
        if (v.empty()) return 0;
        return v[std::max(0, (int)(p * v.size()) - 1)];
    };
    float cut0 = pct(a0, Setup::AMP_PERCENTILE);
    float cut1 = pct(a1, Setup::AMP_PERCENTILE);
    std::cout << "  amp0 cut = " << cut0 << " | amp1 cut = " << cut1 << "\n";

    // --- histograms: [ppac][cathode], cathode 0=X, 1=Y ---
    TH1D* h[2][2];
    const char* cath[2] = {"X", "Y"};
    for (int p = 0; p < 2; ++p)
        for (int cth = 0; cth < 2; ++cth)
            h[p][cth] = new TH1D(
                Form("h_%s_%s_%s", cfg.tag.c_str(), cfg.ppac_label[p].c_str(), cath[cth]),
                "", Setup::N_BINS, Setup::T_MIN, Setup::T_MAX);

    // --- pass 2: fill (X and Y of a PPAC share that PPAC's amplitude cut) ---
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        float rat_amp = (amp1-amp0)/(amp1+amp0);
        if (neutron_energy>1000 || amp0+amp1<18000 || rat_amp>0.0) continue;
        if (sumX0<90 || sumX1<90 || sumY0<90 || sumY1<90) continue;
        if (amp0 > cut0) { h[0][0]->Fill(diffX0); h[0][1]->Fill(diffY0); }
        if (amp1 > cut1) { h[1][0]->Fill(diffX1); h[1][1]->Fill(diffY1); }
    }

    // --- FFT + peak fit + velocity for each of the 4 spectra ---
    for (int p = 0; p < 2; ++p) {
        for (int cth = 0; cth < 2; ++cth) {
            std::string tag = cfg.ppac_label[p] + "_" + cath[cth];

            TH1D* hf = computeFFT(h[p][cth], "fft_" + tag, Setup::PAD_FACTOR);
            PeakResult r = fitPeak(hf,
                Setup::F_SEARCH_MIN, Setup::F_SEARCH_MAX,
                Setup::F_NOISE_MIN,  Setup::F_NOISE_MAX,
                Setup::FIT_HALFWIN_GHZ);

            std::cout << "\n  -- " << cfg.ppac_label[p]
                      << " cathode " << cath[cth] << " --\n";
            std::cout << "     f0      = " << r.f0 << " +/- " << r.sigma_f0
                      << " GHz  (parab " << r.f0_parab << ")\n";
            std::cout << "     v       = " << r.v << " +/- " << r.sigma_v
                      << " mm/ns\n";
            std::cout << "     T       = " << r.period << " ns\n";
            std::cout << "     chi2/ndf= " << r.chi2ndf
                      << (r.fitOk ? "" : "   [FIT FAILED -> parabolic fallback]")
                      << "\n";
            std::cout << "     discrep = "
                      << 100.0 * std::fabs(r.v - Setup::V_NOM_MMNS) / Setup::V_NOM_MMNS
                      << " % vs nominal\n";

            drawFFT(hf, r,
                    Form("FFT %s  cathode %s", cfg.ppac_label[p].c_str(), cath[cth]),
                    Form("fft_%s_%s_%s.pdf", cfg.tag.c_str(),
                         cfg.ppac_label[p].c_str(), cath[cth]));
        }
    }
}

int main()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    TFile* f = TFile::Open(Setup::FILE_PATH, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << Setup::FILE_PATH << "\n";
        return 1;
    }

    // gold: PPAC index 0 = PPAC7, index 1 = PPAC8
    SampleConfig gold {"events_gold", {"PPAC7", "PPAC8"}, "gold"};
    // uranium: PPAC index 0 = PPAC8, index 1 = PPAC9
    SampleConfig uranium {"events_uranium", {"PPAC8", "PPAC9"}, "uranium"};

    processSample(f, gold);
    processSample(f, uranium);

    f->Close();
    std::cout << "\nDone.\n";
    return 0;
}