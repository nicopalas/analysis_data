#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TROOT.h"

void copiar_flux() {

    // 1. Abrir archivo
    TFile *f = TFile::Open("/Users/nico/Desktop/Tese/Macros/Macros/n_tof_cerium/Energy_EAR1_FLUKA_500bdp.root");
    if (!f || f->IsZombie()) {
        printf("No puedo abrir el archivo.\n");
        return;
    }

    // 2. Leer h_Flux
    TH1 *h_Flux = (TH1*)f->Get("h_Flux");
    if (!h_Flux) {
        printf("No se encontró el histograma h_Flux.\n");
        return;
    }

    // 3. Obtener binning logarítmico del histograma original
    Int_t nbins = h_Flux->GetNbinsX();

    // Crear array dinámico de bordes (size nbins+1)
    double *edges = new double[nbins+1];

    for (int i = 0; i <= nbins; i++) {
        edges[i] = h_Flux->GetXaxis()->GetBinLowEdge(i+1);
    }

    // 4. Crear histograma nuevo con binning variable
    TH1F *flux = new TH1F("flux", "Neutron Flux", nbins, edges);

    // 5. Copiar contenido y eliminar errores
    for (int i = 1; i <= nbins; i++) {
        flux->SetBinContent(i, h_Flux->GetBinContent(i));
        flux->SetBinError(i, 0);
    }

    // 6. Títulos de ejes
    flux->GetXaxis()->SetTitle("Neutron energy (eV)");
    flux->GetYaxis()->SetTitle("Ed#phi/dE/pulse");

    // 7. Estilo
    gStyle->SetCanvasPreferGL(1);

    flux->SetFillColor(kRed-3);
    flux->SetFillStyle(1001);
    flux->SetFillColorAlpha(kRed-3, 0.3);
    flux->SetLineColor(kRed);
    flux->SetLineWidth(2);

    // 8. Dibujar usando HIST (sin barras de error)
    TCanvas *c = new TCanvas("c", "Flux", 900, 700);
    flux->Draw("HIST");

    // Limpieza
    delete[] edges;
}
