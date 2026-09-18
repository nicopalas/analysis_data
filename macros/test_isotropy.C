// ================================================================
//  test_isotropy.C
//
//  Closure test del procedimiento angular:
//    - genera eventos isotropos en 4pi (plano en cos_theta)
//    - planos en tiempo entre -5 y 5 ns
//    - los pasa por computeEfficiency() + anisotropy()
//    - comprueba que W(theta)/W(90) == 1 dentro de errores
//
//  TEST A (toy_closure) : no toca fillHistograms(), llena directamente
//                         los Vec3D de cuentas. Testea la parte angular
//                         pura (acceptance -> eficiencia -> anisotropia).
//                         No depende de nombres de ramas ni de cortes.
//
//  TEST B (toy_tree)    : escribe un TTree sintetico con la misma
//                         estructura que tus datos para testear ademas
//                         fillHistograms() + fitBackground() +
//                         computeSignal(). Hay que adaptar los nombres
//                         de las ramas (bloque marcado mas abajo).
//
//  Uso:
//     root -l 'test_isotropy.C+("A")'
//     root -l 'test_isotropy.C+("B")'
// ================================================================

#include "../include/types.h"
#include "../include/constants.h"
#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/utils.h"
#include "../include/acceptance.h"
#include "../include/efficiency.h"
#include "../include/anisotropy.h"
#include "../include/plotting.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// ================================================================
//  PARAMETROS DEL TOY
// ================================================================
namespace toy {

    const std::string acc_csv =
        "/Users/nico/Desktop/Tese/Analysis/acceptance_coincidence.csv";
    const std::string outdir =
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/toy/";

    const UInt_t seed = 20240607;

    // --- anisotropia inyectada:  W(cb) = 1 + a2*cb^2 + a4*cb^4
    //     a2 = a4 = 0  ->  ISOTROPA  (el test que pediste)
    double a2 = 0.0;
    double a4 = 0.0;

    // --- eficiencia inyectada:  eps(cd) = eps0 * (1 - slope*(1-cd))
    //     slope = 0 -> eficiencia plana
    double eps0  = 1.0;
    double slope = 0.0;

    // --- estadistica
    const long N_per_ebin_aniso = 200000;
    const long N_per_ebin_eff   = 2000000;

    // --- TEST B: tiempo
    const double tof_lo   = -5.0,  tof_hi   =  5.0;   // senal PLANA
    const double tof_full = 15.0;                     // rango del histograma
    const double f_bkg    = 0.0;                      // fraccion de fondo plano
    const double f_acc    = 0.0;                      // fraccion de accidentales
    const long   N_tree   = 5000000;                  // eventos totales del arbol
    const double emin_gen = 1.6, emax_gen = 2000.0;
}

// distribuciones inyectadas ---------------------------------------
static double injW(double cb){
    return 1.0 + toy::a2*cb*cb + toy::a4*cb*cb*cb*cb;
}
static double injEps(double cd){
    double e = toy::eps0 * (1.0 - toy::slope*(1.0 - cd));
    return std::max(0.0, std::min(1.0, e));
}

// ================================================================
//  Muestreador de celdas (cos_beam, cos_det)
//
//  Clave del test: NO invento la geometria. Muestreo directamente
//  con probabilidad proporcional a dOmega del propio CSV de
//  aceptancia, multiplicada por W(cb) y eps(cd).
//  Asi, si el analisis divide por esa misma aceptancia, tiene que
//  devolver exactamente W(cb) -> plano si a2=a4=0.
//  (Se asume que ambos ejes del CSV cubren [0,1].)
// ================================================================
struct CellSampler {
    int nb = 0, nd = 0;
    std::vector<double> cum;
    double total = 0.0;
    bool   inject = true;

    CellSampler(const Vec2D& acc, bool inject_ = true) : inject(inject_) {
        nb = (int)acc.size();
        nd = nb ? (int)acc[0].size() : 0;
        cum.resize((size_t)nb*nd, 0.0);
        double s = 0.0;
        for(int j = 0; j < nb; ++j){
            for(int i = 0; i < nd; ++i){
                double cb = (j + 0.5)/nb;
                double cd = (i + 0.5)/nd;
                double w  = acc[j][i];
                if(w < 0) w = 0;
                if(inject) w *= injW(cb) * injEps(cd);
                s += w;
                cum[(size_t)j*nd + i] = s;
            }
        }
        total = s;
    }

    void sample(TRandom3& rnd, double& cb, double& cd) const {
        double u = rnd.Uniform(0.0, total);
        size_t k = std::upper_bound(cum.begin(), cum.end(), u) - cum.begin();
        if(k >= cum.size()) k = cum.size() - 1;
        int j = (int)(k / nd);
        int i = (int)(k % nd);
        cb = (j + rnd.Uniform())/nb;   // posicion dentro de la celda fina
        cd = (i + rnd.Uniform())/nd;
    }
};

// muestreador "plano en celda" (forma angular distinta -> accidentales/fondo)
struct FlatSampler {
    int nb, nd;
    FlatSampler(int nb_, int nd_) : nb(nb_), nd(nd_) {}
    void sample(TRandom3& rnd, double& cb, double& cd) const {
        cb = rnd.Uniform();
        cd = rnd.Uniform();
    }
};

// ================================================================
//  Llenado directo de los Vec3D de cuentas (TEST A)
// ================================================================
static void generateCounts(const CellSampler& s,
                           int nbins_e, long N_per_ebin,
                           Vec3D& counts, Vec3D& u_counts,
                           TRandom3& rnd)
{
    counts.assign(nbins_e,
        Vec2D(nbins_beam, std::vector<double>(nbins_det, 0.0)));
    u_counts = counts;

    for(int e = 0; e < nbins_e; ++e){
        for(long k = 0; k < N_per_ebin; ++k){
            double cb, cd;
            s.sample(rnd, cb, cd);
            int j = std::min(nbins_beam - 1, (int)(cb * nbins_beam));
            int i = std::min(nbins_det  - 1, (int)(cd * nbins_det ));
            counts[e][j][i] += 1.0;
        }
        for(int j = 0; j < nbins_beam; ++j)
            for(int i = 0; i < nbins_det; ++i)
                u_counts[e][j][i] = std::sqrt(counts[e][j][i]);
    }
}

// ================================================================
//  TEST A — cadena angular pura
// ================================================================
void toy_closure()
{
    TRandom3 rnd(toy::seed);

    // --- aceptancia (la misma que usa el analisis real) ---
    Vec2D acc_fine;
    if(!loadAcceptanceCSV(toy::acc_csv, acc_fine)){
        std::cerr << "No se pudo cargar " << toy::acc_csv << "\n";
        return;
    }
    Vec2D acceptance = rebin(acc_fine);
    CellSampler sampler(acc_fine, true);

    std::cout << "CSV fino: " << acc_fine.size() << " x " << acc_fine[0].size()
              << "   rebin: " << acceptance.size() << " x " << acceptance[0].size()
              << "   (nbins_beam=" << nbins_beam
              << ", nbins_det=" << nbins_det << ")\n";

    // --- binning identico al analisis real ---
    const std::vector<double> energy_bins_eff = {1, 10, 100, 500, 1000, 2000};
    const int nbins_eff = (int)energy_bins_eff.size() - 1;

    const int nbins_aniso = 50;
    std::vector<double> energy_bins_aniso = buildLogBins(nbins_aniso, 1.6, 2000.0);

    AnalysisConfig cfg = makeUraniumConfig(energy_bins_aniso, "aniso");

    // --- eficiencia ---
    Vec3D c_eff, u_c_eff;
    generateCounts(sampler, nbins_eff, toy::N_per_ebin_eff, c_eff, u_c_eff, rnd);

    std::vector<EfficiencyResult> eff(nbins_eff);
    for(int e = 0; e < nbins_eff; ++e)
        eff[e] = computeEfficiency(nbins_det - 1, nbins_beam, nbins_det,
                                   c_eff, u_c_eff, acceptance, e);

    std::cout << "\n=== Eficiencia: recuperada vs inyectada (ebin 0) ===\n";
    double eps_ref = injEps((nbins_det - 0.5)/nbins_det);
    for(int i = 0; i < nbins_det; ++i){
        double cd  = (i + 0.5)/nbins_det;
        double tru = injEps(cd)/eps_ref;   // normalizada al bin de referencia
        std::cout << "  cd=" << cd
                  << "  eps_rec=" << eff[0].eps[i] << " +/- " << eff[0].u_eps[i]
                  << "   eps_true=" << tru
                  << "   pull=" << (eff[0].u_eps[i] > 0
                                    ? (eff[0].eps[i]-tru)/eff[0].u_eps[i] : 0.0)
                  << "\n";
    }

    // --- anisotropia ---
    Vec3D c_a, u_c_a;
    generateCounts(sampler, nbins_aniso, toy::N_per_ebin_aniso, c_a, u_c_a, rnd);

    std::vector<AnisotropyResult> aniso(nbins_aniso);
    for(int e = 0; e < nbins_aniso; ++e){
        double Ec = std::sqrt(energy_bins_aniso[e]*energy_bins_aniso[e+1]);
        int e_eff = findBin(energy_bins_eff, Ec);
        if(e_eff < 0)          e_eff = 0;
        if(e_eff >= nbins_eff) e_eff = nbins_eff - 1;

        aniso[e] = anisotropy(nbins_beam, nbins_det,
                              c_a, u_c_a, acceptance, e,
                              eff[e_eff].eps, eff[e_eff].u_eps, cfg);
    }

    // --- comparacion con la verdad ---
    double cb0   = 0.5*dcos_beam;
    double Wnorm = injW(cb0);

    double chi2 = 0.0;
    int    ndf  = 0;
    std::vector<double> mean_w(nbins_beam, 0.0), mean_u(nbins_beam, 0.0);
    std::vector<int>    n_ok  (nbins_beam, 0);

    for(int e = 0; e < nbins_aniso; ++e){
        for(int j = 0; j < nbins_beam; ++j){
            double w  = aniso[e].w[j];
            double uw = aniso[e].u_w[j];
            if(!(uw > 0) || !std::isfinite(w)) continue;
            double cb  = (j + 0.5)*dcos_beam;
            double tru = injW(cb)/Wnorm;
            double p   = (w - tru)/uw;
            chi2 += p*p;
            ++ndf;
            mean_w[j] += w/(uw*uw);
            mean_u[j] += 1.0/(uw*uw);
            ++n_ok[j];
        }
    }

    std::cout << "\n=== W(theta)/W(90): promedio ponderado sobre las "
              << nbins_aniso << " bins de energia ===\n";
    for(int j = 0; j < nbins_beam; ++j){
        if(!n_ok[j]) continue;
        double w  = mean_w[j]/mean_u[j];
        double uw = 1.0/std::sqrt(mean_u[j]);
        double cb = (j + 0.5)*dcos_beam;
        std::cout << "  cos_beam=" << cb
                  << "   W=" << w << " +/- " << uw
                  << "   true=" << injW(cb)/Wnorm
                  << "   desv=" << (w - injW(cb)/Wnorm)/uw << " sigma\n";
    }
    std::cout << "\nchi2/ndf global frente a la verdad = "
              << chi2/std::max(1,ndf) << "  (" << chi2 << "/" << ndf << ")\n";
    std::cout << (std::fabs(chi2/std::max(1,ndf) - 1.0) < 0.3
                  ? "  -> el procedimiento cierra\n"
                  : "  -> hay sesgo o los errores no estan bien propagados\n");

    // --- salida grafica ---
    gSystem->mkdir(toy::outdir.c_str(), kTRUE);

    std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0), ey(nbins_beam);
    for(int j = 0; j < nbins_beam; ++j){
        x[j]  = (j + 0.5)*dcos_beam;
        y[j]  = n_ok[j] ? mean_w[j]/mean_u[j] : 0.0;
        ey[j] = n_ok[j] ? 1.0/std::sqrt(mean_u[j]) : 0.0;
    }
    TCanvas* c = new TCanvas("c_toy", "closure", 900, 600);
    TGraphErrors* g = new TGraphErrors(nbins_beam, x.data(), y.data(),
                                       ex.data(), ey.data());
    g->SetName("w_closure");
    g->SetTitle("Closure test;cos(#theta_{beam});W(#theta)/W(90)");
    g->SetMarkerStyle(20);
    g->SetMinimum(0.8); g->SetMaximum(1.2);
    g->Draw("AP");
    TLine* l = new TLine(0.0, 1.0, 1.0, 1.0);
    l->SetLineStyle(2); l->SetLineColor(kRed+1); l->Draw();
    c->SaveAs((toy::outdir + "closure_isotropic.pdf").c_str());

    plotAnisotropy(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                   toy::outdir + "anisotropy_toy.pdf");
    plotAnisotropyRatio(aniso, nbins_aniso, nbins_beam, energy_bins_aniso,
                        toy::outdir + "anisotropy_ratio_toy");
}

// ================================================================
//  TEST B — arbol sintetico para la cadena completa
//
//  >>>>>>>>>>  ADAPTAR ESTE BLOQUE A TUS RAMAS  <<<<<<<<<<
//  Los nombres/tipos tienen que coincidir EXACTAMENTE con los que
//  lee fillHistograms() y con los cortes de getCuts().
// ================================================================
struct ToyEvent {
    double en       = 0;   // energia del neutron (MeV)
    double tof      = 0;   // tiempo de coincidencia (ns)
    double cos_beam = 0;   // |cos| respecto al haz
    double cos_det  = 0;   // |cos| respecto a la normal del detector
    double x = 0, y = 0;   // posicion (para la ventana mx0..my1)
};

static void bindBranches(TTree* t, ToyEvent& ev){
    t->Branch("en",       &ev.en,       "en/D");
    t->Branch("tof",      &ev.tof,      "tof/D");
    t->Branch("cos_beam", &ev.cos_beam, "cos_beam/D");
    t->Branch("cos_det",  &ev.cos_det,  "cos_det/D");
    t->Branch("x",        &ev.x,        "x/D");
    t->Branch("y",        &ev.y,        "y/D");
}
// >>>>>>>>>>  FIN DEL BLOQUE A ADAPTAR  <<<<<<<<<<

void toy_tree()
{
    TRandom3 rnd(toy::seed + 1);

    Vec2D acc_fine;
    if(!loadAcceptanceCSV(toy::acc_csv, acc_fine)){
        std::cerr << "No se pudo cargar " << toy::acc_csv << "\n";
        return;
    }
    CellSampler sig(acc_fine, true);
    FlatSampler flat((int)acc_fine.size(), (int)acc_fine[0].size());

    gSystem->mkdir(toy::outdir.c_str(), kTRUE);

    // ---- arbol de datos ----
    TFile* fout = TFile::Open((toy::outdir + "toy_data.root").c_str(), "RECREATE");
    TTree* t    = new TTree("toy", "toy isotropic");
    ToyEvent ev;
    bindBranches(t, ev);

    const double lgmin = std::log(toy::emin_gen), lgmax = std::log(toy::emax_gen);

    for(long k = 0; k < toy::N_tree; ++k){
        ev.en = std::exp(rnd.Uniform(lgmin, lgmax));   // plano en log(E)
        ev.x  = rnd.Uniform(-3.0, 3.0);
        ev.y  = rnd.Uniform(-3.0, 3.0);

        double u = rnd.Uniform();
        if(u < toy::f_acc + toy::f_bkg){
            // fondo / accidentales: plano en todo el rango temporal,
            // y con una forma angular DISTINTA (para que la resta importe)
            flat.sample(rnd, ev.cos_beam, ev.cos_det);
            ev.tof = rnd.Uniform(-toy::tof_full, toy::tof_full);
        } else {
            // senal: PLANA en tiempo entre -5 y 5, isotropa en cos_theta
            sig.sample(rnd, ev.cos_beam, ev.cos_det);
            ev.tof = rnd.Uniform(toy::tof_lo, toy::tof_hi);
        }
        t->Fill();
    }
    t->Write();
    fout->Close();

    // ---- arbol de accidentales (plantilla de forma) ----
    TFile* facc = TFile::Open((toy::outdir + "toy_acc.root").c_str(), "RECREATE");
    TTree* ta   = new TTree("toy_acc", "toy accidentals");
    ToyEvent eva;
    bindBranches(ta, eva);
    for(long k = 0; k < toy::N_tree/2; ++k){
        eva.en  = std::exp(rnd.Uniform(lgmin, lgmax));
        eva.x   = rnd.Uniform(-3.0, 3.0);
        eva.y   = rnd.Uniform(-3.0, 3.0);
        eva.tof = rnd.Uniform(-toy::tof_full, toy::tof_full);
        flat.sample(rnd, eva.cos_beam, eva.cos_det);
        ta->Fill();
    }
    ta->Write();
    facc->Close();

    std::cout << "Escritos " << toy::outdir << "toy_data.root y toy_acc.root\n"
              << "Ahora crea un makeToyConfig() que apunte a estos ficheros y\n"
              << "corre la misma secuencia que uranium_analysis().\n";
}

// ================================================================
void test_isotropy(const char* which = "A")
{
    if(std::string(which) == "A") toy_closure();
    else                          toy_tree();
}