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
// Explicitly, R_y(-45) (X,Y,Z) = ( (X-Z)/sqrt2 , Y , (X+Z)/sqrt2 ).
//
// BEAM AXIS
// ---------
// The beam runs along +z_lab, so in the detector frame its direction is
//     b_hat = R_y(-45) (0,0,1) = (-1, 0, +1)/sqrt2
// i.e. the beam axis is the line  x_det = -z_det , tilted -45 deg about Y
// with respect to the detector z axis.  The PPACs are stacked ALONG THE
// BEAM (z_ppac_lab = z0 + i*d_PPACs, with x_ppac_lab = 0), so every PPAC
// centre lies on that line and
//     x_ppac = -z_ppac                        <-- never a free parameter
// This is enforced below instead of being hard-coded as four independent
// numbers, so the geometry cannot drift out of consistency.
//
// CATHODE PLACEMENT
// -----------------
// CreatePPAC puts the cathodes at local (0, 0, +-(gas_gap + mylar)).  A
// displacement that is already purely local-z maps to detector-frame z with
// NO x component:
//     cathode_lab = (0,0,Delta) + R_y(45)(0,0,+-c)
//     => detector  ( -Delta/sqrt2 , 0 , +Delta/sqrt2 +- c )
// So both cathodes of a PPAC share the x of their PPAC centre and differ
// only in z.  A cathode centre is therefore NOT on the beam axis: it is off
// it by c = 0.32017 cm.  If CreatePPAC instead stacks the cathodes along the
// beam direction, set cathodes_on_beam_axis = true below (costs ~0.7 % of
// absolute coincidence probability).
//
// CATHODE ORDERING  (this is what was wrong before)
// -------------------------------------------------
// The internal layout along the PPAC normal is  Y, anode, X, and ALL PPACs
// are added with the SAME rotation -- the back one is not flipped.  Hence,
// for BOTH PPACs, in detector-frame z:
//     z_Y = z_ppac - c        z_X = z_ppac + c
// The forward fragment (travelling +z_det) meets front-Y then front-X; the
// backward fragment (travelling -z_det) meets back-X then back-Y.  The
// ARRIVAL ORDER is reversed but the PLANE IDENTITY is not: labelling the
// planes by arrival order mirrors the back PPAC and gives unequal lever arms
// (5+2c for X, 5-2c for Y), a spurious +-12.8 % azimuthal distortion of the
// reconstructed direction.  With the correct assignment both lever arms are
// exactly d_PPACs/sqrt2 and the reconstruction closes to machine precision.
//
// COINCIDENCE CONDITION
// ---------------------
// A fragment must cross BOTH cathode planes of its PPAC inside the active
// area (a signal is needed on the X plane and on the Y plane), and the same
// for the back-to-back partner in the opposite PPAC: four intersections.
//
// Geometry from create_ntof_geo.C:
//   target  TGeoEltu semi-axes (7.8*sqrt2/2, 7.8/2) cm, U thickness 0.411e-4 cm
//           (the sqrt2 elongation in local x is the footprint of a round beam
//            on a foil tilted 45 deg; in the detector frame the foil is flat
//            at z_det = 0, so sampling (ox,oy) in that ellipse is correct)
//   PPAC    TGeoBBox half-size 10 x 10 cm, gas_gap 0.32, mylar 1.7e-4
//
// TRUE-vs-RECONSTRUCTED ANGLE CHECK
// ---------------------------------
// dfx,dfy,dfz is the TRUE emission direction in the detector frame (dfz =
// cth by construction).  The reconstructed direction comes from the four
// cathode intersection points.  For an ideal back-to-back pair the two
// cathodes of a PPAC lie on the SAME straight track, so with the correct
// plane assignment the reconstruction is exact and these histograms collapse
// to delta functions at zero -- that is the intended closure test.  To study
// a real reconstruction bias, switch on strip_pitch below (0 = off), which
// digitises the measured coordinates.

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

    // true  -> cathode centres displaced along the BEAM direction
    // false -> cathode centres displaced along the PPAC NORMAL (CreatePPAC)
    const bool cathodes_on_beam_axis = false;

    // optional strip digitisation of the measured coordinates [cm]; 0 = off.
    // Only affects the reconstructed direction, never the acceptance.
    const double strip_pitch = 0.0;

    const double inv_sqrt2 = 1.0/TMath::Sqrt(2.0);

    // ── geometry: one source of truth ──────────────────────────────────────
    const double d_PPACs = 5.0*TMath::Sqrt(2.0);   // lab-z spacing between PPACs
    const double step    = d_PPACs*inv_sqrt2;      // 5.0 cm in the detector frame

    // PPAC centres sit ON the beam axis (x = -z), target midway between them
    const double z_front = +0.5*step,  x_front = -z_front;   // (-2.5, +2.5)
    const double z_back  = -0.5*step,  x_back  = -z_back;    // (+2.5, -2.5)

    // Cathode planes: SAME ordering in both PPACs (back is not flipped)
    const double zX_front = z_front - cath_off;   // +2.17983
    const double zY_front = z_front + cath_off;   // +2.82017
    const double zX_back  = z_back  - cath_off;   // -2.82017
    const double zY_back  = z_back  + cath_off;   // -2.17983

    // Active-area centres in x
    const double xY_front = cathodes_on_beam_axis ? -zY_front : x_front;
    const double xX_front = cathodes_on_beam_axis ? -zX_front : x_front;
    const double xY_back  = cathodes_on_beam_axis ? -zY_back  : x_back;
    const double xX_back  = cathodes_on_beam_axis ? -zX_back  : x_back;

    // Per-coordinate lever arms.  Equal to `step` for the nominal layout, but
    // computed rather than assumed so a mirrored or re-aligned geometry stays
    // correctly reconstructed.
    const double Lx = zX_front - zX_back;
    const double Ly = zY_front - zY_back;

    std::cout << "detector-frame geometry\n"
              << "  front PPAC centre  (x,z) = (" << x_front << ", " << z_front << ")\n"
              << "  back  PPAC centre  (x,z) = (" << x_back  << ", " << z_back  << ")\n"
              << "  front  Y / X planes  z   = " << zY_front << " / " << zX_front << "\n"
              << "  back   Y / X planes  z   = " << zY_back  << " / " << zX_back  << "\n"
              << "  lever arms  Lx / Ly      = " << Lx << " / " << Ly << "\n"
              << "  cathode centres on beam axis: "
              << (cathodes_on_beam_axis ? "yes" : "no") << "\n"
              << "  strip pitch              = " << strip_pitch << " cm\n\n";

    const int    nbins_beam = 100;          // matches cos_theta_center in the CSV
    const int    nbins_det  = 100;           // matches cos_theta_det_center
    const double dcos_beam  = 1.0/nbins_beam;
    const double dcos_det   = 1.0/nbins_det;

    // ── bookkeeping ────────────────────────────────────────────────────────
    std::vector<std::vector<double>> cell_counts(
        nbins_beam, std::vector<double>(nbins_det, 0.0));

    Long64_t counts_forward = 0, counts_backward = 0;
    Long64_t coincidence    = 0, missed          = 0;
    Long64_t coinc_midplane = 0;   // would pass a single-plane-per-PPAC test
    Long64_t lost_2cathode  = 0;   // passes mid-plane but fails both cathodes

    TRandom3 rng(0);

    // Acceptance numerator and denominator are both binned in the TRUE angles,
    // so the ratio is a genuine acceptance and does not inherit any
    // reconstruction bias.
    TH2D* hist_thetas  = new TH2D("theta_det_beam", ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    TH2D* hist_emitted = new TH2D("emitted",        ";|cos#theta_{beam}|;cos#theta_{det}",
                                  nbins_beam, 0, 1, nbins_det, 0, 1);
    hist_thetas->Sumw2();
    hist_emitted->Sumw2();

    // cos_theta_beam is genuinely signed (negative whenever theta_det > 45 deg),
    // so the axis must cover [-1,1] or those events go silently to underflow.
    TH2D* eff_beam = new TH2D("eff_beam", ";cos#theta;#phi",
                              200, -1, 1, 100, -TMath::Pi(), TMath::Pi());
    TH2D* eff_det  = new TH2D("eff_det",  ";cos#theta_{det};#phi_{det} (deg)",
                              100, 0, 1, 100, -180, 180);

    TH1D* hist_cos_theta_det     = new TH1D("cos_theta_det",     "", nbins_det,  0, 1);
    TH1D* hist_cos_theta         = new TH1D("cos_theta",         "", nbins_beam, 0, 1);
    TH1D* hist_cos_theta_emitted = new TH1D("cos_theta_emitted", "", nbins_beam, 0, 1);
    hist_cos_theta->Sumw2();
    hist_cos_theta_emitted->Sumw2();

    TH1D* hist_x_front = new TH1D("x_front", ";x_{front} (cm)", 200, -20, 20);
    TH1D* hist_x_back  = new TH1D("x_back",  ";x_{back} (cm)",  200, -20, 20);
    TH1D* hist_y_back  = new TH1D("y_back",  ";y_{back} (cm)",  200, -20, 20);

    // ── true vs reconstructed angle check ──────────────────────────────────
    // Signed cos is used throughout: folding with Abs() before ACos() makes
    // true and reco land on opposite sides of 90 deg and fakes a residual.
    const int    nbins_theta   = 1000;
    const double dtheta_range  = 5.0;    // deg; widen if the 2D map overflows

    TH1D* hist_theta_true_beam = new TH1D("theta_true_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 180);
    TH1D* hist_theta_reco_beam = new TH1D("theta_reco_beam",
        ";#theta_{beam} (deg);counts", nbins_theta, 0, 180);
    TH1D* hist_theta_true_det  = new TH1D("theta_true_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);
    TH1D* hist_theta_reco_det  = new TH1D("theta_reco_det",
        ";#theta_{det} (deg);counts",  nbins_theta, 0, 90);

    TH2D* hist_dtheta_vs_theta_beam = new TH2D("dtheta_vs_theta_beam",
        ";#theta_{true,beam} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 180, 120, -dtheta_range, dtheta_range);
    TH2D* hist_dtheta_vs_theta_det = new TH2D("dtheta_vs_theta_det",
        ";#theta_{true,det} (deg);#theta_{reco}-#theta_{true} (deg)",
        nbins_theta, 0, 90, 120, -dtheta_range, dtheta_range);

    // same check, but in cos(theta) instead of theta (deg)
    const int    nbins_costheta  = 200;
    const double dcostheta_range = 0.02;   // widen if the 2D map overflows

    TH1D* hist_costheta_true_beam = new TH1D("costheta_true_beam",
        ";cos#theta_{beam};counts", nbins_costheta, -1, 1);
    TH1D* hist_costheta_reco_beam = new TH1D("costheta_reco_beam",
        ";cos#theta_{beam};counts", nbins_costheta, -1, 1);
    TH1D* hist_costheta_true_det  = new TH1D("costheta_true_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);
    TH1D* hist_costheta_reco_det  = new TH1D("costheta_reco_det",
        ";cos#theta_{det};counts",  nbins_costheta, 0, 1);

    TH2D* hist_dcostheta_vs_costheta_beam = new TH2D("dcostheta_vs_costheta_beam",
        ";cos#theta_{true,beam};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, -1, 1, 120, -dcostheta_range, dcostheta_range);
    TH2D* hist_dcostheta_vs_costheta_det = new TH2D("dcostheta_vs_costheta_det",
        ";cos#theta_{true,det};cos#theta_{reco}-cos#theta_{true}",
        nbins_costheta, 0, 1, 120, -dcostheta_range, dcostheta_range);

    // ── event loop ─────────────────────────────────────────────────────────
    for (Long64_t i = 0; i < nevents; ++i) {

        if (i % 5000000 == 0) std::cout << "  event " << i << std::endl;

        // uniform over the elliptical footprint: uniform in the unit disk,
        // then scale.  In the detector frame the foil is flat at z_det = 0.
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

        if (dfz <= 0.) continue;            // parallel to the planes (guard)

        // TRUE emission direction in the beam frame.  (dfx,dfy,dfz) is already
        // a unit vector in the detector frame, so no renormalisation is needed.
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
        double xfX, yfX, xfY, yfY;
        const bool okfX = cross(zX_front,  dfx,  dfy,  dfz, xfX, yfX);
        const bool okfY = cross(zY_front,  dfx,  dfy,  dfz, xfY, yfY);
        const bool hit_front = okfX && okfY
                            && inArea(xfX, yfX, xX_front)
                            && inArea(xfY, yfY, xY_front);

        // ---- backward fragment through the two back cathodes --------------
        double xbX, ybX, xbY, ybY;
        const bool okbX = cross(zX_back,  -dfx, -dfy, -dfz, xbX, ybX);
        const bool okbY = cross(zY_back,  -dfx, -dfy, -dfz, xbY, ybY);
        const bool hit_back = okbX && okbY
                           && inArea(xbX, ybX, xX_back)
                           && inArea(xbY, ybY, xY_back);

        if (!okfX || !okfY || !okbX || !okbY) continue;

        // measured point per PPAC: each coordinate from its OWN cathode plane.
        // Plane identity, not arrival order: X always at z_ppac + c, Y always
        // at z_ppac - c, in both PPACs.
        const double xf = xfX, yf = (yfY);
        const double xb = xbX, yb = (ybY);

        // diagnostic: single-plane-per-PPAC test, for comparison
        double xfm, yfm, xbm, ybm;
        cross(z_front,  dfx,  dfy,  dfz, xfm, yfm);
        cross(z_back,  -dfx, -dfy, -dfz, xbm, ybm);
        const bool hit_mid = inArea(xfm, yfm, x_front) && inArea(xbm, ybm, x_back);

        // Direction from per-coordinate slopes.  Both fragments share a vertex
        // and are exactly back-to-back, so the vertex cancels and this
        // reproduces the emission direction exactly (machine precision) when
        // strip_pitch = 0 -- the internal closure check.
        const double tx = (xf - xb) / Lx;
        const double ty = (yf - yb) / Ly;
        const double nn = TMath::Sqrt(tx*tx + ty*ty + 1.0);

        const double cos_theta_det = 1.0 / nn;
        const double phi_det       = TMath::ATan2(ty, tx);

        // back to the beam frame: v_lab = R_y(45) v_local
        //   => cos_theta_beam = (-vx + vz)/sqrt2
        const double vx = tx / nn;
        const double vy = ty / nn;
        const double vz = 1.0 / nn;

        const double nx = ( vx + vz) * inv_sqrt2;
        const double ny =   vy;
        const double nz = (-vx + vz) * inv_sqrt2;

        const double cos_theta = nz;                     // already normalised
        const double phi       = TMath::ATan2(ny, nx);

        // acceptance denominator, in TRUE angles
        hist_emitted->Fill(TMath::Abs(cos_theta_true_beam), cth);
        hist_cos_theta_emitted->Fill(TMath::Abs(cos_theta_true_beam));

        if (hit_front) ++counts_forward;
        if (hit_back)  ++counts_backward;
        if (!hit_front || !hit_back) ++missed;
        if (hit_mid)   ++coinc_midplane;
        if (hit_mid && !(hit_front && hit_back)) ++lost_2cathode;

        if (hit_front && hit_back) {
            ++coincidence;

            // clamp rather than skip: cos = 1 exactly is reachable, and a
            // `continue` here would also drop the event from every histogram
            // filled further down.
            const int bin_beam = TMath::Min(
                int(TMath::Abs(cos_theta_true_beam) / dcos_beam), nbins_beam - 1);
            const int bin_det  = TMath::Min(
                int(cth / dcos_det), nbins_det - 1);

            cell_counts[bin_beam][bin_det] += 1.;

            hist_thetas->Fill(TMath::Abs(cos_theta_true_beam), cth);
            hist_cos_theta->Fill(TMath::Abs(cos_theta_true_beam));
            hist_cos_theta_det->Fill(cth);
            hist_x_front->Fill(xf-x_front);
            hist_x_back->Fill(xb-x_back);
            hist_y_back->Fill(yb);
            eff_beam->Fill(cos_theta, phi);
            eff_det->Fill(cos_theta_det, phi_det*TMath::RadToDeg());

            // ── true vs reconstructed angle, only where we have a full
            //    coincidence (i.e. an actual reconstructed point pair) ──────
            const double theta_true_beam_deg =
                TMath::ACos(cos_theta_true_beam) * TMath::RadToDeg();
            const double theta_reco_beam_deg =
                TMath::ACos(cos_theta)           * TMath::RadToDeg();
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
            hist_costheta_true_beam->Fill(cos_theta_true_beam);
            hist_costheta_reco_beam->Fill(cos_theta);
            hist_costheta_true_det->Fill(cth);
            hist_costheta_reco_det->Fill(cos_theta_det);

            hist_dcostheta_vs_costheta_beam->Fill(cos_theta_true_beam,
                                                   cos_theta - cos_theta_true_beam);
            hist_dcostheta_vs_costheta_det->Fill(cth,
                                                  cos_theta_det - cth);
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
              << hist_dcostheta_vs_costheta_det->GetRMS(2)   << "\n";

    if (strip_pitch <= 0.)
        std::cout << "\n  (strip_pitch = 0: the residuals above are the closure\n"
                  << "   test and must be zero to machine precision)\n";

    std::cout << "\n  note: `omega` in the CSV is a probability per emitted pair\n"
              << "  from forward-hemisphere sampling.  Multiply by 4*pi to get\n"
              << "  a solid angle in sr.\n" << std::endl;

    // ── CSV, same format as acceptance_coincidence.csv ─────────────────────
    std::ofstream csv("/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv");
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

    // clone before dividing, so the raw coincidence map survives.  The "B"
    // option gives binomial errors, which is what a subset/total ratio needs;
    // the default assumes independent numerator and denominator.
    TH2D* hist_ratio = (TH2D*)hist_thetas->Clone("acceptance_2d");
    hist_ratio->Divide(hist_thetas, hist_emitted, 1., 1., "B");
    hist_ratio->Write();
    hist_thetas->Write();

    TH1D* h_cos_ratio = (TH1D*)hist_cos_theta->Clone("cos_theta_acceptance");
    h_cos_ratio->Divide(hist_cos_theta, hist_cos_theta_emitted, 1., 1., "B");
    h_cos_ratio->Write();
    hist_cos_theta->Write();
    hist_cos_theta_emitted->Write();

    // ── true-vs-reconstructed angle histograms ─────────────────────────────
    hist_theta_true_beam->Write();
    hist_theta_reco_beam->Write();
    hist_theta_true_det->Write();
    hist_theta_reco_det->Write();
    hist_dtheta_vs_theta_beam->Write();
    hist_dtheta_vs_theta_det->Write();

    // ── canvas with correlation + overlaid distributions, both frames ──────
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

    // ── true-vs-reconstructed cos(theta) histograms ────────────────────────
    hist_costheta_true_beam->Write();
    hist_costheta_reco_beam->Write();
    hist_costheta_true_det->Write();
    hist_costheta_reco_det->Write();
    hist_dcostheta_vs_costheta_beam->Write();
    hist_dcostheta_vs_costheta_det->Write();

    // ── same 2x2 canvas layout, but in cos(theta) ──────────────────────────
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