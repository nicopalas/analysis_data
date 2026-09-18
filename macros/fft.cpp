// ============================================================================
//  Delay-line propagation-velocity measurement via FFT of the diff-time
//  spectra, run directly on out_cathodes.root (the hadd of every
//  anode_cathode_<run>.root produced by the preliminary matching).
//
//  WHY THAT INPUT
//    The strip pitch imprints a periodic modulation on the time-difference
//    spectrum, and the FFT reads its frequency.  That modulation is a property
//    of the delay line, not of the selection, so the RAW pairs from the
//    matching stage are the right sample: they carry far more statistics than
//    the selected ones, and no selection bias.  A constant offset in diff (the
//    electronic zero, `delta`) shifts the phase but NOT the frequency, so the
//    velocity is unaffected and no diff calibration is needed here.
//
//  WHAT CHANGED W.R.T. THE events_selection VERSION
//    [1] Tree is `cathode_coincidences`, not `events_gold`/`events_uranium`.
//    [2] x_diff<d> / y_diff<d> are std::vector<double> holding EVERY candidate
//        pair for that anode hit, not one scalar per event.  All of them are
//        histogrammed: the wrong pairings add a smooth background, which the
//        FFT ignores, while the periodic strip structure survives.
//    [3] The detectors are addressed by their own index (0..9) instead of by
//        sample.  There is no gold/uranium split at this stage.
//    [4] amp<d> is a per-detector scalar (the anode), so the amplitude
//        percentile cut is computed per detector over the hits where that
//        detector fired.
//    [5] The sum window is applied per CANDIDATE (each pair has its own sum),
//        not per event, and is a proper window [100,140] rather than sum>90.
//    [6] Noise band and peak search are clamped to the Nyquist frequency.  With
//        200 ns over 2000 bins, dt = 0.1 ns and Nyquist = 5 GHz, so the old
//        F_NOISE_MAX = 6 GHz ran into the overflow bin.
//
//  Physics:  the dt spectrum is in ns  ->  the FFT frequency axis is in GHz.
//            v = pitch * f0  [mm/ns], with pitch the strip pitch.
//            For a 20 cm active length spanning +-100 ns, v = 1 mm/ns and a
//            2 mm pitch gives f0 = 0.5 GHz.
// ============================================================================
#include <TVirtualFFT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace Setup {
    const int    NDET         = 10;
    // detectors to process; empty-looking entries are simply skipped if the
    // detector never fires
    const int    DET_FIRST    = 0;
    const int    DET_LAST     = 9;

    const double V_NOM_MMNS   = 1.0;      // nominal delay-line velocity [mm/ns]
    const double PITCH_MM     = 2.0;      // strip pitch [mm]

    const int    N_BINS       = 2000;     // dt histogram bins  -> dt = 0.1 ns
    const double T_MIN        = -100.0;   // [ns]
    const double T_MAX        =  100.0;   // [ns]

    const int    PAD_FACTOR   = 4;        // FFT zero-padding factor
    const double AMP_PERCENTILE = 0.8;    // keep the top 20% of anode amplitudes

    const double F_SEARCH_MIN = 0.3;      // peak search band [GHz]
    const double F_SEARCH_MAX = 0.8;
    const double F_NOISE_MIN  = 1.0;      // noise-floor estimation band [GHz]
    const double F_NOISE_MAX  = 4.5;      // below Nyquist (5 GHz for dt=0.1 ns)
    const double FIT_HALFWIN_GHZ = 0.10;  // half-window of the local gaussian fit

    // per-candidate sum window: the propagation time is L + 2*t_cable and lies
    // between 100 and 140 ns.  This removes the pairs that are not on the line.
    const double SUM_MIN      = 100.0;
    const double SUM_MAX      = 140.0;
    const double DIFF_ABS_MAX = 100.0;    // delay line length [ns]

    const char*  FILE_PATH    =
        "/nucl_lustre/n_tof_INTC_P_665/Analysis/anode_cathode_matching/out_cathodes.root";
    const char*  TREE_NAME    = "cathode_coincidences";

    const std::string OUTDIR  = "./";     // where the PDFs go
}

// ---------------------------------------------------------------------------
// FFT with zero-padding.  The frequency-axis bin centers are placed exactly at
// the discrete FFT frequencies f_i = i / (nfft * dt), so GetBinCenter() returns
// true frequencies with no half-bin offset.
// ---------------------------------------------------------------------------
static TH1D* computeFFT(TH1D* hInput, const std::string& name, int padFactor)
{
    int    nData = hInput->GetNbinsX();
    int    nfft  = padFactor * nData;
    double dt    = hInput->GetBinWidth(1);

    // subtract the baseline before transforming, otherwise the DC term leaks
    double mean = 0.0;
    for (int i = 0; i < nData; ++i) mean += hInput->GetBinContent(i + 1);
    mean /= nData;

    TVirtualFFT* fft = TVirtualFFT::FFT(1, &nfft, "R2C ES K");
    if (!fft) { std::cerr << "  ERROR: cannot create FFT\n"; return nullptr; }
    for (int i = 0; i < nData; ++i)
        fft->SetPoint(i, hInput->GetBinContent(i + 1) - mean);
    for (int i = nData; i < nfft; ++i) fft->SetPoint(i, 0.0);   // zero-pad
    fft->Transform();

    TH1* hMag = nullptr;
    hMag = TH1::TransformHisto(fft, hMag, "MAG");
    if (!hMag) { delete fft; return nullptr; }

    int    nHalf = nfft / 2;
    double df    = 1.0 / (nfft * dt);            // frequency resolution [GHz]
    double xlo   = -0.5 * df;                    // so that bin i center = i*df
    double xhi   = (nHalf - 0.5) * df;

    TH1D* hFreq = new TH1D(name.c_str(),
        (name + ";Frequency (GHz);Magnitude (arb. units)").c_str(),
        nHalf, xlo, xhi);
    hFreq->SetDirectory(nullptr);

    for (int i = 0; i < nHalf; ++i)
        hFreq->SetBinContent(i + 1, hMag->GetBinContent(i + 1));

    hFreq->SetBinContent(1, 0.0);   // suppress DC

    delete hMag;
    delete fft;
    return hFreq;
}

// bin of the maximum inside the search band (clamped to the histogram)
static int findPeakBin(TH1D* h, double fmin, double fmax)
{
    int nb = h->GetNbinsX();
    int b0 = std::max(1,  h->FindBin(fmin));
    int b1 = std::min(nb, h->FindBin(fmax));
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

// noise floor: RMS of the magnitude in a control band away from the peak
static double estimateNoiseRMS(TH1D* h, double fmin, double fmax)
{
    int nb = h->GetNbinsX();
    int b0 = std::max(1,  h->FindBin(fmin));
    int b1 = std::min(nb, h->FindBin(fmax));      // never reach the overflow bin
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
    double f0        = 0;   // GHz (gaussian fit, falls back to parabolic)
    double sigma_f0  = 0;   // GHz (statistical, from the fit)
    double f0_parab  = 0;   // GHz (parabolic interpolation, cross-check)
    double v         = 0;   // mm/ns
    double sigma_v   = 0;   // mm/ns
    double period    = 0;   // ns
    double chi2ndf   = 0;
    bool   fitOk     = false;
};

// local gaussian fit around the peak; per-bin error = the noise RMS
static PeakResult fitPeak(TH1D* h, double fmin, double fmax,
                          double fnoise_min, double fnoise_max,
                          double halfwin_ghz)
{
    PeakResult r;
    int kmax   = findPeakBin(h, fmin, fmax);
    r.f0_parab = parabolicPeak(h, kmax);

    double noise = estimateNoiseRMS(h, fnoise_min, fnoise_max);
    for (int i = 1; i <= h->GetNbinsX(); ++i) h->SetBinError(i, noise);

    double xc = h->GetBinCenter(kmax);
    TF1 g("gpeak", "[0]*TMath::Gaus(x,[1],[2]) + [3]",
          xc - halfwin_ghz, xc + halfwin_ghz);
    g.SetParameters(h->GetBinContent(kmax), xc, h->GetBinWidth(kmax) * 3, 0.0);
    g.SetParLimits(2, h->GetBinWidth(kmax) * 0.5, halfwin_ghz);

    int status = h->Fit(&g, "QRN");   // quiet, use the range, do not draw
    r.fitOk = (status == 0);

    if (r.fitOk) {
        r.f0       = g.GetParameter(1);
        r.sigma_f0 = g.GetParError(1);
        r.chi2ndf  = (g.GetNDF() > 0) ? g.GetChisquare() / g.GetNDF() : 0.0;
    } else {                          // fall back to the parabolic estimate
        r.f0       = r.f0_parab;
        r.sigma_f0 = h->GetBinWidth(kmax) / std::sqrt(12.0);
    }

    r.period  = (r.f0 > 0) ? 1.0 / r.f0 : 0.0;
    r.v       = Setup::PITCH_MM * r.f0;
    r.sigma_v = Setup::PITCH_MM * r.sigma_f0;
    return r;
}

// draw one FFT spectrum with the peak and fit annotations
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
    std::cout << "     saved: " << outFile << "\n";
}

// ---------------------------------------------------------------------------
static void processFile(const char* path, const char* tree_name)
{
    std::cout << "\n================================================\n"
              << "  FFT delay-line velocity from " << path << "\n"
              << "================================================\n";

    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) { std::cerr << "  ERROR: cannot open the file\n"; return; }
    TTree* tree = (TTree*) f->Get(tree_name);
    if (!tree) { std::cerr << "  ERROR: tree " << tree_name << " not found\n"; f->Close(); return; }

    const Long64_t nEntries = tree->GetEntries();
    std::cout << "  entries: " << nEntries << "\n";

    // ── branches: only what is needed, the merged file is large ─────────────
    Int_t   mult[Setup::NDET];
    Float_t amp[Setup::NDET];
    std::vector<double> *x_diff[Setup::NDET] = {}, *x_sum[Setup::NDET] = {};
    std::vector<double> *y_diff[Setup::NDET] = {}, *y_sum[Setup::NDET] = {};

    tree->SetBranchStatus("*", 0);
    for (int d = 0; d < Setup::NDET; ++d) {
        const char* names[] = {"mult%d", "amp%d", "x_diff%d", "x_sum%d",
                               "y_diff%d", "y_sum%d"};
        for (const char* n : names) tree->SetBranchStatus(Form(n, d), 1);

        tree->SetBranchAddress(Form("mult%d", d),   &mult[d]);
        tree->SetBranchAddress(Form("amp%d", d),    &amp[d]);
        tree->SetBranchAddress(Form("x_diff%d", d), &x_diff[d]);
        tree->SetBranchAddress(Form("x_sum%d", d),  &x_sum[d]);
        tree->SetBranchAddress(Form("y_diff%d", d), &y_diff[d]);
        tree->SetBranchAddress(Form("y_sum%d", d),  &y_sum[d]);
    }

    // ── pass 1: anode amplitude percentile, per detector ────────────────────
    // A high-amplitude cut sharpens the strip modulation: those are the hits
    // with the best timing, so the periodic structure is less washed out.
    std::vector<std::vector<float>> aBuf(Setup::NDET);
    for (int d = 0; d < Setup::NDET; ++d) aBuf[d].reserve(nEntries / 4);

    std::cout << "  pass 1: amplitude percentiles...\n";
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (i % 200000 == 0)
            std::cout << "    " << i << "/" << nEntries << "\r" << std::flush;
        for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d)
            if (mult[d] >= 1 && amp[d] > 0) aBuf[d].push_back(amp[d]);
    }

    std::vector<float> ampCut(Setup::NDET, 0.f);
    auto pct = [](std::vector<float>& v, double p) -> float {
        if (v.empty()) return 0.f;
        size_t k = (size_t) std::max(0.0, p * v.size() - 1);
        std::nth_element(v.begin(), v.begin() + k, v.end());
        return v[k];
    };
    std::cout << "\n  anode amplitude cuts (top "
              << (int) std::lround(100 * (1 - Setup::AMP_PERCENTILE)) << "%):\n";
    for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d) {
        ampCut[d] = pct(aBuf[d], Setup::AMP_PERCENTILE);
        std::cout << "    det " << d << "  n=" << aBuf[d].size()
                  << "  cut=" << ampCut[d] << "\n";
        aBuf[d].clear(); aBuf[d].shrink_to_fit();
    }

    // ── histograms: [det][cathode], cathode 0 = X, 1 = Y ────────────────────
    TH1D* h[Setup::NDET][2] = {};
    const char* cath[2] = {"X", "Y"};
    for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d)
        for (int c = 0; c < 2; ++c) {
            h[d][c] = new TH1D(Form("h_det%d_%s", d, cath[c]),
                               Form("det %d cathode %s;t_{1}-t_{2} [ns];pairs", d, cath[c]),
                               Setup::N_BINS, Setup::T_MIN, Setup::T_MAX);
            h[d][c]->SetDirectory(nullptr);
        }

    // ── pass 2: fill, per CANDIDATE PAIR ────────────────────────────────────
    // Every pair in the vector is histogrammed if it passes its own sum window.
    // The wrong pairings contribute a smooth background that the FFT ignores.
    std::cout << "\n  pass 2: filling the diff spectra...\n";
    long long nFilled = 0, nRejSum = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (i % 200000 == 0)
            std::cout << "    " << i << "/" << nEntries << "\r" << std::flush;

        for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d) {
            if (mult[d] < 1)          continue;
            if (amp[d] <= ampCut[d])  continue;

            if (x_diff[d] && x_sum[d]) {
                size_t n = std::min(x_diff[d]->size(), x_sum[d]->size());
                for (size_t j = 0; j < n; ++j) {
                    const double s = (*x_sum[d])[j], dd = (*x_diff[d])[j];
                    if (s < Setup::SUM_MIN || s > Setup::SUM_MAX) { nRejSum++; continue; }
                    if (std::fabs(dd) >= Setup::DIFF_ABS_MAX)     continue;
                    h[d][0]->Fill(dd); nFilled++;
                }
            }
            if (y_diff[d] && y_sum[d]) {
                size_t n = std::min(y_diff[d]->size(), y_sum[d]->size());
                for (size_t j = 0; j < n; ++j) {
                    const double s = (*y_sum[d])[j], dd = (*y_diff[d])[j];
                    if (s < Setup::SUM_MIN || s > Setup::SUM_MAX) { nRejSum++; continue; }
                    if (std::fabs(dd) >= Setup::DIFF_ABS_MAX)     continue;
                    h[d][1]->Fill(dd); nFilled++;
                }
            }
        }
    }
    std::cout << "\n  pairs histogrammed: " << nFilled
              << "   rejected by the sum window: " << nRejSum << "\n";

    // ── FFT + peak fit + velocity for every spectrum ────────────────────────
    std::cout << "\n  det pl        f0 [GHz]          v [mm/ns]      T [ns]"
                 "   chi2/ndf   discrep\n";
    for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d) {
        for (int c = 0; c < 2; ++c) {
            if (!h[d][c] || h[d][c]->GetEntries() < 1000) {
                std::cout << std::setw(5) << d << "  " << cath[c]
                          << "   too few entries, skipped\n";
                continue;
            }
            std::string tag = Form("det%d_%s", d, cath[c]);
            TH1D* hf = computeFFT(h[d][c], "fft_" + tag, Setup::PAD_FACTOR);
            if (!hf) continue;

            PeakResult r = fitPeak(hf,
                Setup::F_SEARCH_MIN, Setup::F_SEARCH_MAX,
                Setup::F_NOISE_MIN,  Setup::F_NOISE_MAX,
                Setup::FIT_HALFWIN_GHZ);

            std::cout << std::fixed << std::setprecision(4)
                      << std::setw(5) << d << "  " << cath[c]
                      << std::setw(9)  << r.f0 << " +/- " << std::setw(7) << r.sigma_f0
                      << std::setw(10) << r.v  << " +/- " << std::setw(7) << r.sigma_v
                      << std::setw(10) << r.period
                      << std::setw(10) << std::setprecision(2) << r.chi2ndf
                      << std::setw(9)  << std::setprecision(2)
                      << 100.0 * std::fabs(r.v - Setup::V_NOM_MMNS) / Setup::V_NOM_MMNS << " %"
                      << (r.fitOk ? "" : "   [FIT FAILED -> parabolic]") << "\n";

            drawFFT(hf, r, Form("FFT det %d cathode %s", d, cath[c]),
                    Form("fft_%s.pdf", tag.c_str()));
            delete hf;
        }
    }
    std::cout << std::defaultfloat;

    for (int d = Setup::DET_FIRST; d <= Setup::DET_LAST; ++d)
        for (int c = 0; c < 2; ++c) delete h[d][c];
    f->Close();
}

// ---------------------------------------------------------------------------
void delay_line_velocity(const char* path = Setup::FILE_PATH,
                         const char* tree_name = Setup::TREE_NAME)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    processFile(path, tree_name);
    std::cout << "\nDone.\n";
}

#ifndef __CLING__
void fft(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : Setup::FILE_PATH;
    delay_line_velocity(path, Setup::TREE_NAME);
    return 0;
}
#endif