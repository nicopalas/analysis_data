// acceptance_mc.C
//
// Geometric acceptance of a PPAC pair for back-to-back fission fragments.
//
// Works in the DETECTOR frame, where the 45-degree tilt is undone and the
// PPAC planes are perpendicular to the local z axis.
//
// FRAME CONSTRUCTION
// ------------------
// The PPAC volumes are rotated +45 deg about Y in the lab, so
//     v_lab = R_y(45) v_local   =>   v_local = R_y(-45) v_lab
// A pure lab-z displacement Delta (how the PPACs are stacked:
// z_ppac = z0 + i*d_PPACs, with x_ppac = 0) therefore maps to
//     x_det = -Delta/sqrt2 ,  z_det = +Delta/sqrt2        (so x_det = -z_det)
//
// With the target midway (d_target_PPAC = 2.5*sqrt2, d_PPACs = 5*sqrt2):
//     front PPAC centre:  x_det = -2.5 ,  z_det = +2.5
//     back  PPAC centre:  x_det = +2.5 ,  z_det = -2.5
//
// The CATHODES, however, are stacked along the PPAC NORMAL: CreatePPAC puts
// them at local (0, 0, +-(gas_gap + mylar_thickness)).  A displacement that is
// already purely local-z maps to detector-frame z with NO x component, so both
// cathodes of a PPAC share the same x_det; only z_det differs by +-0.32017 cm.
//
// COINCIDENCE CONDITION
// ---------------------
// A fragment must cross BOTH cathode planes of its PPAC inside the active
// area (a signal is needed on the X plane and on the Y plane), and the same
// for the back-to-back partner in the opposite PPAC: four intersections.
//
// Geometry from create_ntof_geo.C:
//   target  TGeoEltu semi-axes (7.8*sqrt2/2, 7.8/2) cm, U thickness 0.411e-4 cm
//   PPAC    TGeoBBox half-size 10 x 10 cm, gas_gap 0.32, mylar 1.7e-4
//
// TRUE-vs-RECONSTRUCTED ANGLE CHECK
// ----------------------------------
// dfx,dfy,dfz is the TRUE emission direction, already expressed in the
// detector frame (by construction dfz = cth). The reconstructed direction
// (dx,dy,dz -> cos_theta_det, cos_theta) instead comes from the cathode
// intersection points (xf,yf) and (xb,yb), i.e. it carries the effect of
// the finite gap between the two cathode planes of each PPAC. Comparing
// the two isolates that reconstruction bias from the true angular
// distribution, in both the detector frame and the beam frame.

#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

void acceptance_mc()
{
    // ── configuration ──────────────────────────────────────────────────────
    const Long64_t nevents = 50000000;      // 5e7: omega = counts/nevents

    const double target_a         = 7.8*TMath::Sqrt(2)/2.0;  // 5.5154 cm semi-axis x
    const double target_b         = 7.8/2.0;                 // 3.9    cm semi-axis y
    const double target_thickness = 0.41105294e-4;           // cm, U target

    const double ppac_half = 10.0;          // half active area [cm]

    // internal PPAC layout (CreatePPAC), along the PPAC normal
    const double gas_gap         = 0.32;
    const double mylar_thickness = 1.7e-4;
    const double cath_off        = gas_gap + mylar_thickness;   // 0.32017 cm
    // set to (gas_gap + mylar_thickness)/2 to use the GAS-VOLUME centres
    // instead of the cathode planes; the two differ by ~1.6 mm in z_det.

    // PPAC centres in the detector frame, origin at the target
    const double z_front = +2.5,  x_front = -2.5;
    const double z_back  = -2.5,  x_back  = +2.5;
    // (empirical values from the data alignment: -2.355 / +2.340)

    // Cathode planes.  Real internal layout per PPAC is X, anode, Y (in that
    // order along the PPAC normal). The forward fragment travels in +z_det
    // and meets the low-z layer of the front PPAC first (cathode X), then
    // cathode Y. The backward fragment travels in -z_det and meets the
    // HIGH-z layer of the back PPAC first, i.e. cathode Y, then cathode X --
    // the PPACs are all added with the same rotation, so the back one is
    // not flipped.
    const double zf_first  = z_front - cath_off;   // front cathode_x
    const double zf_second = z_front + cath_off;   // front cathode_y
    const double zb_first  = z_back  + cath_off;   // back  cathode_y
    const double zb_second = z_back  - cath_off;   // back  cathode_x

    const int    nbins_beam = 100;          // matches cos_theta_center in the CSV
    const int    nbins_det  = 20;           // matches cos_theta_det_center
    const double dcos_beam  = 1.0/nbins_beam;
    const double dcos_det   = 1.0/nbins_det;

    const double inv_sqrt2 = 1.0/TMath::Sqrt(2.0);

    // ── bookkeeping ────────────────────────────────────────────────────────
    std::vector<std::vector<double>> cell_counts(
        nbins_beam, std::vector<double>(nbins_det, 0.0));

    Long64_t counts_forward = 0, counts_backward = 0;
    Long64_t coincidence    = 0, missed          = 0;
    Long64_t coinc_midplane = 0;   // would pass a single-plane-per-PPAC test
    Long64_t lost_2cathode  = 0;   // passes mid-plane but fails both cathodes

    TRandom3 rng(0);

    TH2D* hist_thetas  = new TH2D("theta_det_beam", ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    TH2D* hist_emitted = new TH2D("emitted",        ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    TH2D* eff_beam     = new TH2D("eff_beam", ";cos#theta;#phi",
                                  100, 0, 1, 100, -TMath::Pi(), TMath::Pi());
    TH2D* eff_det      = new TH2D("eff_det",  ";cos#theta_{det};#phi_{det} (deg)",
                                  100, 0, 1, 100, -180, 180);

    TH1D* hist_cos_theta_det     = new TH1D("cos_theta_det",     "", nbins_det,  0, 1);
    TH1D* hist_cos_theta         = new TH1D("cos_theta",         "", nbins_beam, 0, 1);
    TH1D* hist_cos_theta_emitted = new TH1D("cos_theta_emitted", "", nbins_beam, 0, 1);
    TH1D* hist_x_front = new TH1D("x_front", ";x_{front} (cm)", 200, -20, 20);
    TH1D* hist_x_back  = new TH1D("x_back",  ";x_{back} (cm)",  200, -20, 20);
    TH1D* hist_y_back  = new TH1D("y_back",  ";y_{back} (cm)",  200, -20, 20);

    // ── NEW: true vs reconstructed angle check ────────────────────────────
    // Beam frame and detector frame, each with:
    //   - a 1D "true" and "reconstructed" angle distribution (deg), overlaid
    //   - a 2D (theta_true, theta_reco - theta_true) correlation
    const int    nbins_theta   = 1000;
    const double dtheta_range  = 15.0;   // deg; widen if the 2D map overflows

    TH1D* hist_theta_true_beam = new TH1D("theta_true_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 90);
    TH1D* hist_theta_reco_beam = new TH1D("theta_reco_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 90);
    TH1D* hist_theta_true_det  = new TH1D("theta_true_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);
    TH1D* hist_theta_reco_det  = new TH1D("theta_reco_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);

    TH2D* hist_dtheta_vs_theta_beam = new TH2D("dtheta_vs_theta_beam",
        ";#theta_{true,beam} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 90, 120, -dtheta_range, dtheta_range);
    TH2D* hist_dtheta_vs_theta_det = new TH2D("dtheta_vs_theta_det",
        ";#theta_{true,det} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 90, 120, -dtheta_range, dtheta_range);

    // same check, but in cos(theta) instead of theta (deg) -- same events,
    // same true/reco quantities, just the other common variable
    const int    nbins_costheta  = 100;
    const double dcostheta_range = 0.05;   // widen if the 2D map overflows

    TH1D* hist_costheta_true_beam = new TH1D("costheta_true_beam",
        ";cos#theta_{beam};counts", nbins_costheta, 0, 1);
    TH1D* hist_costheta_reco_beam = new TH1D("costheta_reco_beam",
        ";cos#theta_{beam};counts", nbins_costheta, 0, 1);
    TH1D* hist_costheta_true_det  = new TH1D("costheta_true_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);
    TH1D* hist_costheta_reco_det  = new TH1D("costheta_reco_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);

    TH2D* hist_dcostheta_vs_costheta_beam = new TH2D("dcostheta_vs_costheta_beam",
        ";cos#theta_{true,beam};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, 0, 1, 120, -dcostheta_range, dcostheta_range);
    TH2D* hist_dcostheta_vs_costheta_det = new TH2D("dcostheta_vs_costheta_det",
        ";cos#theta_{true,det};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, 0, 1, 120, -dcostheta_range, dcostheta_range);

    // ── event loop ─────────────────────────────────────────────────────────
    for (Long64_t i = 0; i < nevents; ++i) {

        if (i % 5000000 == 0) std::cout << "  event " << i << std::endl;

        // uniform over the elliptical target: uniform in the unit disk, then scale
        const double u_r   = TMath::Sqrt(rng.Uniform(0., 1.));
        const double u_phi = TMath::TwoPi() * rng.Uniform(0., 1.);
        const double ox = target_a * u_r * TMath::Cos(u_phi);
        const double oy = target_b * u_r * TMath::Sin(u_phi);
        const double oz = (0.5 - rng.Uniform(0., 1.)) * target_thickness;

        // isotropic into the forward hemisphere; the partner is exactly
        // back-to-back, so one event covers the full sphere as a pair
        const double cth = rng.Uniform(0., 1.);
        const double sth = TMath::Sqrt(1. - cth*cth);
        const double phd = TMath::TwoPi() * rng.Uniform(0., 1.);

        const double dfx = TMath::Cos(phd)*sth;
        const double dfy = TMath::Sin(phd)*sth;
        const double dfz = cth;

        if (dfz <= 0.) continue;            // parallel to the planes

        // TRUE emission direction, mapped to the beam frame with the same
        // R_y(45) rotation used below for the reconstructed direction.
        // (dfx,dfy,dfz) is already a unit vector in the detector frame, so
        // no renormalisation is needed here.
        const double cos_theta_true_beam = (-dfx + dfz) * inv_sqrt2;
        // cos_theta_true_det is just dfz == cth, already at hand.

        // intersection of a ray from (ox,oy,oz) with the plane z = zp
        auto cross = [&](double zp, double vx, double vy, double vz,
                         double& X, double& Y) -> bool {
            if (vz == 0.) return false;
            const double t = (zp - oz) / vz;
            if (t < 0.) return false;
            X = ox + t*vx;
            Y = oy + t*vy;
            return true;
        };
        auto inArea = [&](double X, double Y, double xoff) -> bool {
            return TMath::Abs(X - xoff) <= ppac_half
                && TMath::Abs(Y)        <= ppac_half;
        };

        // ---- forward fragment through the two front cathodes --------------
        double xf1, yf1, xf2, yf2;
        const bool okf1 = cross(zf_first,  dfx, dfy, dfz, xf1, yf1);
        const bool okf2 = cross(zf_second, dfx, dfy, dfz, xf2, yf2);
        const bool hit_front = okf1 && okf2
                            && inArea(xf1, yf1, x_front)
                            && inArea(xf2, yf2, x_front);

        // ---- backward fragment through the two back cathodes --------------
        double xb1, yb1, xb2, yb2;
        const bool okb1 = cross(zb_first,  -dfx, -dfy, -dfz, xb1, yb1);
        const bool okb2 = cross(zb_second, -dfx, -dfy, -dfz, xb2, yb2);
        const bool hit_back = okb1 && okb2
                           && inArea(xb1, yb1, x_back)
                           && inArea(xb2, yb2, x_back);

        if (!okf1 || !okf2 || !okb1 || !okb2) continue;

        // measured point per PPAC: each coordinate comes from its OWN cathode
        // plane, not an average of both. Real layout per PPAC is X, anode, Y,
        // so the forward-going fragment meets the X plane first, then Y.
        //   front (going +z): X @ zf_first -> xf1 ; Y @ zf_second -> yf2
        //   back  (going -z): Y @ zb_first -> yb1 ; X @ zb_second -> xb2
        const double xf = xf1, yf = yf2;
        const double xb = xb2, yb = yb1;

        // diagnostic: single-plane-per-PPAC test, for comparison
        double xfm, yfm, xbm, ybm;
        cross(z_front,  dfx,  dfy,  dfz, xfm, yfm);
        cross(z_back,  -dfx, -dfy, -dfz, xbm, ybm);
        const bool hit_mid = inArea(xfm, yfm, x_front) && inArea(xbm, ybm, x_back);

        // axis from the two reconstructed points.  Both fragments share a
        // vertex and are exactly back-to-back, so this reproduces the emission
        // direction exactly -- a useful internal check.
        const double dx = xf - xb;
        const double dy = yf - yb;
        const double dz = z_front - z_back;
        const double nn = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

        const double cos_theta_det = dz / nn;
        const double sin_theta_det = TMath::Sqrt(1. - cos_theta_det*cos_theta_det);
        const double phi_det       = TMath::ATan2(dy, dx);

        // back to the beam frame: v_lab = R_y(45) v_local
        //   => cos_theta_beam = (-vx + vz)/sqrt2
        const double vx = sin_theta_det * TMath::Cos(phi_det);
        const double vy = sin_theta_det * TMath::Sin(phi_det);
        const double vz = cos_theta_det;

        const double nx = ( vx + vz) * inv_sqrt2;
        const double ny =   vy;
        const double nz = (-vx + vz) * inv_sqrt2;
        const double nb = TMath::Sqrt(nx*nx + ny*ny + nz*nz);
        if (nb <= 0.) continue;

        const double cos_theta = nz / nb;
        const double phi       = TMath::ATan2(ny, nx);

        hist_emitted->Fill(TMath::Abs(cos_theta), cos_theta_det);
        hist_cos_theta_emitted->Fill(TMath::Abs(cos_theta));

        if (hit_front) ++counts_forward;
        if (hit_back)  ++counts_backward;
        if (!hit_front || !hit_back) ++missed;
        if (hit_mid)   ++coinc_midplane;
        if (hit_mid && !(hit_front && hit_back)) ++lost_2cathode;

        if (hit_front && hit_back) {
            ++coincidence;

            const int bin_beam = int(TMath::Abs(cos_theta)     / dcos_beam);
            const int bin_det  = int(TMath::Abs(cos_theta_det) / dcos_det);
            if (bin_beam >= nbins_beam || bin_det >= nbins_det) continue;

            cell_counts[bin_beam][bin_det] += 1.;

            hist_thetas->Fill(TMath::Abs(cos_theta), cos_theta_det);
            hist_cos_theta->Fill(TMath::Abs(cos_theta));
            hist_cos_theta_det->Fill(cos_theta_det);
            hist_x_front->Fill(xf);
            hist_x_back->Fill(xb);
            hist_y_back->Fill(yb);
            eff_beam->Fill(cos_theta, phi);
            eff_det->Fill(cos_theta_det, phi_det*TMath::RadToDeg());

            // ── true vs reconstructed angle, only where we have a full
            //    coincidence (i.e. an actual reconstructed point pair) ──────
            const double theta_true_beam_deg =
                TMath::ACos(TMath::Abs(cos_theta_true_beam)) * TMath::RadToDeg();
            const double theta_reco_beam_deg =
                TMath::ACos(TMath::Abs(cos_theta))            * TMath::RadToDeg();
            const double theta_true_det_deg  = TMath::ACos(cth)           * TMath::RadToDeg();
            const double theta_reco_det_deg  = TMath::ACos(cos_theta_det) * TMath::RadToDeg();

            hist_theta_true_beam->Fill(theta_true_beam_deg);
            hist_theta_reco_beam->Fill(theta_reco_beam_deg);
            hist_theta_true_det->Fill(theta_true_det_deg);
            hist_theta_reco_det->Fill(theta_reco_det_deg);

            hist_dtheta_vs_theta_beam->Fill(theta_true_beam_deg,
                                             theta_reco_beam_deg - theta_true_beam_deg);
            hist_dtheta_vs_theta_det->Fill(theta_true_det_deg,
                                            theta_reco_det_deg  - theta_true_det_deg);

            // ── same check in cos(theta) instead of theta (deg) ─────────────
            const double costheta_true_beam = TMath::Abs(cos_theta_true_beam);
            const double costheta_reco_beam = TMath::Abs(cos_theta);
            const double costheta_true_det  = cth;            // already in [0,1]
            const double costheta_reco_det  = cos_theta_det;   // already in [0,1]

            hist_costheta_true_beam->Fill(costheta_true_beam);
            hist_costheta_reco_beam->Fill(costheta_reco_beam);
            hist_costheta_true_det->Fill(costheta_true_det);
            hist_costheta_reco_det->Fill(costheta_reco_det);

            hist_dcostheta_vs_costheta_beam->Fill(costheta_true_beam,
                                                   costheta_reco_beam - costheta_true_beam);
            hist_dcostheta_vs_costheta_det->Fill(costheta_true_det,
                                                  costheta_reco_det  - costheta_true_det);
        }
    }

    // ── summary ────────────────────────────────────────────────────────────
    std::cout << "\n---------- RESULTS ----------\n"
              << "  events              = " << nevents << "\n"
              << "  hits front (2 cath) = " << 100.*counts_forward /nevents << " %\n"
              << "  hits back  (2 cath) = " << 100.*counts_backward/nevents << " %\n"
              << "  coincidences        = " << 100.*coincidence    /nevents << " %\n"
              << "  missed              = " << 100.*missed         /nevents << " %\n"
              << "\n  single-plane test   = " << 100.*coinc_midplane/nevents << " %\n"
              << "  lost by requiring both cathodes = "
              << 100.*lost_2cathode/nevents << " %  ("
              << (coinc_midplane > 0
                  ? 100.*lost_2cathode/double(coinc_midplane) : 0.)
              << " % of single-plane coincidences)\n"
              << "\n  <theta_reco - theta_true> (beam) = "
              << hist_dtheta_vs_theta_beam->GetMean(2) << " deg,  RMS = "
              << hist_dtheta_vs_theta_beam->GetRMS(2)  << " deg\n"
              << "  <theta_reco - theta_true> (det)  = "
              << hist_dtheta_vs_theta_det->GetMean(2)  << " deg,  RMS = "
              << hist_dtheta_vs_theta_det->GetRMS(2)   << " deg\n"
              << "\n  <costheta_reco - costheta_true> (beam) = "
              << hist_dcostheta_vs_costheta_beam->GetMean(2) << ",  RMS = "
              << hist_dcostheta_vs_costheta_beam->GetRMS(2)  << "\n"
              << "  <costheta_reco - costheta_true> (det)  = "
              << hist_dcostheta_vs_costheta_det->GetMean(2)  << ",  RMS = "
              << hist_dcostheta_vs_costheta_det->GetRMS(2)   << "\n"
              << std::endl;

    // ── CSV, same format as acceptance_coincidence.csv ─────────────────────
    std::ofstream csv("acceptance_coincidence.csv");
    csv << "cos_theta_center,cos_theta_det_center,counts,omega\n";
    csv << std::fixed;
    for (int j = 0; j < nbins_beam; ++j) {
        const double c_beam = (j + 0.5) * dcos_beam;
        for (int k = 0; k < nbins_det; ++k) {
            const double c_det  = (k + 0.5) * dcos_det;
            const double counts = cell_counts[j][k];
            csv << std::setprecision(10) << c_beam << ','
                << std::setprecision(10) << c_det  << ','
                << std::setprecision(0)  << counts << ','
                << std::setprecision(10) << counts/double(nevents) << '\n';
        }
    }
    csv.close();
    std::cout << "wrote acceptance_coincidence.csv ("
              << nbins_beam*nbins_det << " rows)\n";

    // ── ROOT output ────────────────────────────────────────────────────────
    TFile* fout = new TFile("mc_acceptance.root", "RECREATE");

    hist_emitted->Write();
    eff_beam->Write();
    eff_det->Write();
    hist_x_front->Write();
    hist_x_back->Write();
    hist_y_back->Write();
    hist_cos_theta_det->Write();

    TH1D* acceptance = new TH1D("acceptance", ";cos#theta_{det};acceptance",
                                nbins_det, 0, 1);
    for (int k = 0; k < nbins_det; ++k) {
        double emitted = 0., hits = 0.;
        for (int j = 0; j < nbins_beam; ++j) {
            emitted += hist_emitted->GetBinContent(j+1, k+1);
            hits    += cell_counts[j][k];
        }
        const double a = (emitted > 0.) ? hits/emitted : 0.;
        acceptance->SetBinContent(k+1, a);
        acceptance->SetBinError  (k+1, (emitted > 0.)
                                       ? TMath::Sqrt(a*(1.-a)/emitted) : 0.);
    }
    acceptance->Write();

    // clone before dividing, so the raw coincidence map survives
    TH2D* hist_ratio = (TH2D*)hist_thetas->Clone("acceptance_2d");
    hist_ratio->Divide(hist_emitted);
    hist_ratio->Write();
    hist_thetas->Write();

    TH1D* h_cos_ratio = (TH1D*)hist_cos_theta->Clone("cos_theta_acceptance");
    h_cos_ratio->Divide(hist_cos_theta_emitted);
    h_cos_ratio->Write();
    hist_cos_theta->Write();
    hist_cos_theta_emitted->Write();

    // ── NEW: write the true-vs-reconstructed angle histograms ──────────────
    hist_theta_true_beam->Write();
    hist_theta_reco_beam->Write();
    hist_theta_true_det->Write();
    hist_theta_reco_det->Write();
    hist_dtheta_vs_theta_beam->Write();
    hist_dtheta_vs_theta_det->Write();

    // ── NEW: canvas with correlation + overlaid distributions, both frames ─
    gStyle->SetOptStat(1111);
    gStyle->SetPalette(kBird);

    TCanvas* c_angle_check = new TCanvas("c_angle_check",
                                          "true vs reconstructed angle", 1200, 900);
    c_angle_check->Divide(2, 2);

    c_angle_check->cd(1);
    hist_dtheta_vs_theta_beam->SetTitle("beam frame: reco-true vs true");
    hist_dtheta_vs_theta_beam->Draw("COLZ");

    c_angle_check->cd(2);
    hist_theta_true_beam->SetLineColor(kBlue+1);
    hist_theta_true_beam->SetLineWidth(3);
    hist_theta_true_beam->SetLineStyle(1);
    hist_theta_reco_beam->SetLineColor(kRed+1);
    hist_theta_reco_beam->SetLineWidth(2);
    hist_theta_reco_beam->SetLineStyle(2);
    hist_theta_true_beam->SetTitle("beam frame: true vs reconstructed");
    hist_theta_true_beam->SetStats(0);
    hist_theta_reco_beam->SetStats(0);
    {
        const double ymax = 1.1 * TMath::Max(hist_theta_true_beam->GetMaximum(),
                                              hist_theta_reco_beam->GetMaximum());
        hist_theta_true_beam->SetMaximum(ymax);
    }
    hist_theta_true_beam->Draw("HIST");
    hist_theta_reco_beam->Draw("HIST SAME");
    TLegend* leg_beam = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg_beam->SetBorderSize(0);
    leg_beam->AddEntry(hist_theta_true_beam, "true", "l");
    leg_beam->AddEntry(hist_theta_reco_beam, "reconstructed", "l");
    leg_beam->Draw();

    c_angle_check->cd(3);
    hist_dtheta_vs_theta_det->SetTitle("detector frame: reco-true vs true");
    hist_dtheta_vs_theta_det->Draw("COLZ");

    c_angle_check->cd(4);
    hist_theta_true_det->SetLineColor(kBlue+1);
    hist_theta_true_det->SetLineWidth(3);
    hist_theta_true_det->SetLineStyle(1);
    hist_theta_reco_det->SetLineColor(kRed+1);
    hist_theta_reco_det->SetLineWidth(2);
    hist_theta_reco_det->SetLineStyle(2);
    hist_theta_true_det->SetTitle("detector frame: true vs reconstructed");
    hist_theta_true_det->SetStats(0);
    hist_theta_reco_det->SetStats(0);
    {
        const double ymax = 1.1 * TMath::Max(hist_theta_true_det->GetMaximum(),
                                              hist_theta_reco_det->GetMaximum());
        hist_theta_true_det->SetMaximum(ymax);
    }
    hist_theta_true_det->Draw("HIST");
    hist_theta_reco_det->Draw("HIST SAME");
    TLegend* leg_det = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg_det->SetBorderSize(0);
    leg_det->AddEntry(hist_theta_true_det, "true", "l");
    leg_det->AddEntry(hist_theta_reco_det, "reconstructed", "l");
    leg_det->Draw();

    c_angle_check->Write();
    c_angle_check->SaveAs("angle_true_vs_reco.png");
    c_angle_check->SaveAs("angle_true_vs_reco.pdf");

    // ── NEW: write the true-vs-reconstructed cos(theta) histograms ─────────
    hist_costheta_true_beam->Write();
    hist_costheta_reco_beam->Write();
    hist_costheta_true_det->Write();
    hist_costheta_reco_det->Write();
    hist_dcostheta_vs_costheta_beam->Write();
    hist_dcostheta_vs_costheta_det->Write();

    // ── NEW: same 2x2 canvas layout, but in cos(theta) ──────────────────────
    TCanvas* c_costheta_check = new TCanvas("c_costheta_check",
                                             "true vs reconstructed cos(theta)", 1200, 900);
    c_costheta_check->Divide(2, 2);

    c_costheta_check->cd(1);
    hist_dcostheta_vs_costheta_beam->SetTitle("beam frame: reco-true vs true");
    hist_dcostheta_vs_costheta_beam->Draw("COLZ");

    c_costheta_check->cd(2);
    hist_costheta_true_beam->SetLineColor(kBlue+1);
    hist_costheta_true_beam->SetLineWidth(3);
    hist_costheta_true_beam->SetLineStyle(1);
    hist_costheta_reco_beam->SetLineColor(kRed+1);
    hist_costheta_reco_beam->SetLineWidth(2);
    hist_costheta_reco_beam->SetLineStyle(2);
    hist_costheta_true_beam->SetTitle("beam frame: true vs reconstructed");
    hist_costheta_true_beam->SetStats(0);
    hist_costheta_reco_beam->SetStats(0);
    {
        const double ymax = 1.1 * TMath::Max(hist_costheta_true_beam->GetMaximum(),
                                              hist_costheta_reco_beam->GetMaximum());
        hist_costheta_true_beam->SetMaximum(ymax);
    }
    hist_costheta_true_beam->Draw("HIST");
    hist_costheta_reco_beam->Draw("HIST SAME");
    TLegend* leg_costheta_beam = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg_costheta_beam->SetBorderSize(0);
    leg_costheta_beam->AddEntry(hist_costheta_true_beam, "true", "l");
    leg_costheta_beam->AddEntry(hist_costheta_reco_beam, "reconstructed", "l");
    leg_costheta_beam->Draw();

    c_costheta_check->cd(3);
    hist_dcostheta_vs_costheta_det->SetTitle("detector frame: reco-true vs true");
    hist_dcostheta_vs_costheta_det->Draw("COLZ");

    c_costheta_check->cd(4);
    hist_costheta_true_det->SetLineColor(kBlue+1);
    hist_costheta_true_det->SetLineWidth(3);
    hist_costheta_true_det->SetLineStyle(1);
    hist_costheta_reco_det->SetLineColor(kRed+1);
    hist_costheta_reco_det->SetLineWidth(2);
    hist_costheta_reco_det->SetLineStyle(2);
    hist_costheta_true_det->SetTitle("detector frame: true vs reconstructed");
    hist_costheta_true_det->SetStats(0);
    hist_costheta_reco_det->SetStats(0);
    {
        const double ymax = 1.1 * TMath::Max(hist_costheta_true_det->GetMaximum(),
                                              hist_costheta_reco_det->GetMaximum());
        hist_costheta_true_det->SetMaximum(ymax);
    }
    hist_costheta_true_det->Draw("HIST");
    hist_costheta_reco_det->Draw("HIST SAME");
    TLegend* leg_costheta_det = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg_costheta_det->SetBorderSize(0);
    leg_costheta_det->AddEntry(hist_costheta_true_det, "true", "l");
    leg_costheta_det->AddEntry(hist_costheta_reco_det, "reconstructed", "l");
    leg_costheta_det->Draw();

    c_costheta_check->Write();
    c_costheta_check->SaveAs("costheta_true_vs_reco.png");
    c_costheta_check->SaveAs("costheta_true_vs_reco.pdf");

    fout->Close();
    std::cout << "wrote mc_acceptance.root\n";
    std::cout << "wrote angle_true_vs_reco.png / .pdf\n";
    std::cout << "wrote costheta_true_vs_reco.png / .pdf\n";
}