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

#include "TRandom3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
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

    // Cathode planes.  The forward fragment travels in +z_det and meets the
    // low-z layer of the front PPAC first (cathode_y), then cathode_x.
    // The backward fragment travels in -z_det and meets the HIGH-z layer of
    // the back PPAC first, i.e. cathode_x, then cathode_y -- the PPACs are all
    // added with the same rotation, so the back one is not flipped.
    const double zf_first  = z_front - cath_off;   // front cathode_y
    const double zf_second = z_front + cath_off;   // front cathode_x
    const double zb_first  = z_back  + cath_off;   // back  cathode_x
    const double zb_second = z_back  - cath_off;   // back  cathode_y

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

        // measured point per PPAC: midpoint of the crossing (both cathodes)
        const double xf = 0.5*(xf1 + xf2), yf = 0.5*(yf1 + yf2);
        const double xb = 0.5*(xb1 + xb2), yb = 0.5*(yb1 + yb2);

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

    fout->Close();
    std::cout << "wrote mc_acceptance.root\n";
}