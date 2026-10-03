#pragma once
#include "types.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TLegend.h"
#include "TFile.h"
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "../include/fit_anisotropy.h"

// ========================================================================
//  ESTILO DE PUBLICACION
//
//  Todas las figuras comparten:
//    - una sola funcion de estilo (setPubStyle) y unos pocos helpers de
//      marco, leyenda y cabecera, para que tamanos y margenes coincidan;
//    - una paleta sobria: tinta casi negra para ejes, azul profundo para
//      los datos, terracota para los ajustes, verde azulado para los
//      residuos, granate para W(0)/W(90), grises claros para referencias;
//    - una paleta secuencial (azul -> verde azulado -> ambar -> terracota
//      -> granate) cuando las series son bines de energia ORDENADOS, y una
//      categorica (Paul Tol "muted", apta para daltonicos) para conjuntos
//      sin orden, como los datos de EXFOR.
//  Las leyendas se colocan en un hueco reservado en la parte superior del
//  marco (el rango en y se amplia lo necesario), asi nunca tapan datos.
//
//  Novedades de esta version:
//    - etiqueta de figura opcional (setFigureTag) arriba a la derecha de
//      las figuras de un solo panel, p. ej. "^{197}Au(n,f), MC";
//    - ejes log con etiquetas en {1,2,5}x10^n cuando el rango cubre
//      menos de ~2 decadas (10-900 MeV ya no queda con solo "10, 100");
//    - series superpuestas desplazadas ligeramente en x (dodge) y con
//      marcadores distintos, legibles tambien en blanco y negro;
//    - nombres de canvas unicos: sin avisos "Deleting canvas with same name".
// ========================================================================

static int hexColor(const char* hex) { return TColor::GetColor(hex); }

// se conserva por compatibilidad con codigo que la use
static int okabeIto(double r, double g, double b)
{
    return TColor::GetColor((Float_t)(r/255.), (Float_t)(g/255.), (Float_t)(b/255.));
}

// --- colores de la casa --------------------------------------------------
static const int kInk      = hexColor("#1B1B1B");   // ejes y texto
static const int kRefGray  = hexColor("#8C8C8C");   // lineas de referencia
static const int kBandGray = hexColor("#E6E6E6");   // bandas de referencia

static const int kAnisoColor    = hexColor("#1F4E79");  // azul profundo   -> datos W(theta)
static const int kThisWorkColor = hexColor("#1F4E79");  // azul profundo   -> "este trabajo"
static const int kAnisoFill     = hexColor("#CFDDEA");  // azul muy claro  -> bandas de dispersion
static const int kFitColor      = hexColor("#D1603D");  // terracota       -> ajustes
static const int kResidColor    = hexColor("#2A9D8F");  // verde azulado   -> residuos
static const int kRatioColor    = hexColor("#7A1F4F");  // granate         -> W(0)/W(90) vs E
static const int kBkgColor      = hexColor("#3A3A3A");  // gris carbon     -> espectros crudos
static const int kBkgFitColor   = hexColor("#D1603D");  // terracota       -> ajuste de fondo
static const int kSubColor      = hexColor("#2E7D4F");  // verde bosque    -> espectros restados
static const int kBkgFill       = hexColor("#DCDCDC");  // relleno del espectro crudo
static const int kSubFill       = hexColor("#CFE5D8");  // relleno del espectro restado

// --- paletas -------------------------------------------------------------
// categorica (Paul Tol "muted"): series sin orden
static std::vector<int> categoricalPalette()
{
    return { hexColor("#332288"), hexColor("#CC6677"), hexColor("#117733"),
             hexColor("#882255"), hexColor("#44AA99"), hexColor("#999933"),
             hexColor("#AA4499") };
}
// compatibilidad con el nombre antiguo
static std::vector<int> okabeItoPalette() { return categoricalPalette(); }

// EXFOR: la categorica sin el indigo, que se confundiria con "este trabajo"
static std::vector<int> exforPalette()
{
    return { hexColor("#CC6677"), hexColor("#117733"), hexColor("#882255"),
             hexColor("#44AA99"), hexColor("#999933"), hexColor("#AA4499"),
             hexColor("#555555") };
}

// secuencial para bines de energia: el color codifica el orden
static std::vector<int> energyPalette(int n)
{
    static const double stops[5][3] = {
        { 31,  78, 121},   // azul profundo
        { 42, 157, 143},   // verde azulado
        {214, 160,  60},   // ambar
        {209,  96,  61},   // terracota
        {122,  31,  79}    // granate
    };
    std::vector<int> out;
    for(int i = 0; i < n; ++i){
        const double t = (n > 1) ? double(i) / (n - 1) : 0.0;
        const double s = t * 4.0;
        const int    k = std::min(3, (int)s);
        const double u = s - k;
        double rgb[3];
        for(int j = 0; j < 3; ++j) rgb[j] = stops[k][j] + (stops[k+1][j] - stops[k][j]) * u;
        out.push_back(TColor::GetColor((Float_t)(rgb[0]/255.), (Float_t)(rgb[1]/255.),
                                       (Float_t)(rgb[2]/255.)));
    }
    return out;
}

// marcadores abiertos para conjuntos externos (sin triangulos)
static const int kOpenMarkers[6] = {24, 25, 27, 28, 42, 46};

// marcadores llenos para series propias ordenadas (circulo, cuadrado,
// rombo, cruz, estrella, aspa): distinguibles sin color
static const int kSeriesMarkers[6] = {20, 21, 33, 34, 29, 47};

// rombos y estrellas se ven mas pequenos a igual tamano nominal
static double markerScale(int m)
{
    if(m == 33 || m == 29) return 1.35;
    if(m == 34 || m == 47) return 1.10;
    return 1.0;
}

// desplazamiento horizontal de la serie i de n, repartidas en 'width'
static double dodge(int i, int n, double width)
{
    return (n > 1) ? (i - 0.5 * (n - 1)) * width / n : 0.0;
}

// nombres unicos para canvas y pads
static std::string uniqueName(const std::string& base)
{
    static int counter = 0;
    return base + "_" + std::to_string(counter++);
}

// etiqueta de figura (vacia = no se dibuja)
static std::string& figureTag()
{
    static std::string tag;
    return tag;
}
static void setFigureTag(const std::string& s) { figureTag() = s; }

// ========================================================================
//  helpers de estilo
// ========================================================================
static void setPubStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);

    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameLineColor(kInk);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.030, "X");
    gStyle->SetTickLength(0.020, "Y");
    gStyle->SetNdivisions(510, "X");
    gStyle->SetNdivisions(505, "Y");

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetTitleSize(0.050, "xyz");
    gStyle->SetLabelOffset(0.010, "xyz");
    gStyle->SetAxisColor(kInk, "xyz");
    gStyle->SetLabelColor(kInk, "xyz");
    gStyle->SetTitleColor(kInk, "xyz");

    gStyle->SetEndErrorSize(0);      // barras de error sin remates
    gStyle->SetErrorX(0.);           // sin barra horizontal en histogramas
    gStyle->SetMarkerSize(1.0);

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.038);

    TGaxis::SetMaxDigits(4);
}

// cuadricula de pads aproximadamente cuadrada
static void gridLayout(int n, int& ncol, int& nrow)
{
    ncol = (int)std::ceil(std::sqrt((double)n));
    nrow = (int)std::ceil((double)n / ncol);
}

static void stylePad(TVirtualPad* pad)
{
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.04);
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.14);
    pad->SetTickx(1);
    pad->SetTicky(1);
    pad->SetFillColor(kWhite);
}

// Ejes de un marco o histograma. s escala tamanos cuando el pad es mas bajo
// que la celda (paneles de residuos); los offsets de titulo no se escalan
// porque ya son proporcionales al tamano de letra.
static void styleFrame(TH1* fr, const char* xtitle, const char* ytitle,
                       double s = 1.0, double yoff = 1.45, double xoff = 1.10)
{
    TAxis* ax = fr->GetXaxis();
    TAxis* ay = fr->GetYaxis();
    ax->SetTitle(xtitle);
    ay->SetTitle(ytitle);
    for(TAxis* a : {ax, ay}){
        a->SetLabelFont(42);
        a->SetTitleFont(42);
        a->SetLabelSize(0.045 * s);
        a->SetTitleSize(0.050 * s);
        a->SetAxisColor(kInk);
        a->SetLabelColor(kInk);
        a->SetTitleColor(kInk);
    }
    ax->SetLabelOffset(0.010 * s);
    ay->SetLabelOffset(0.010);
    ax->SetTitleOffset(xoff);
    ay->SetTitleOffset(yoff);
    ax->SetTickLength(0.030 * s);
    ay->SetTickLength(0.020);
}

// eje log estandar: solo decadas (1, 10, 100, 1000), sin notacion 10^n
// (se conserva por compatibilidad; las figuras nuevas usan logXLabels)
static void logLabels(TH1* fr)
{
    fr->GetXaxis()->SetMoreLogLabels(kFALSE);
    fr->GetXaxis()->SetNoExponent(kTRUE);
}

static std::string fmtTick(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof buf, "%g", v);
    return buf;
}

// Etiquetas de un eje X logaritmico. Si el rango cubre mas de ~2.3 decadas
// basta con las decadas de ROOT; si no, se dibujan a mano en {1,2,5}x10^n,
// que es lo habitual en revistas (p. ej. 10, 20, 50, 100, 200, 500).
// Llamar despues de styleFrame y con los margenes del pad ya fijados.
static void logXLabels(TVirtualPad* pad, TH1* fr, double xmin, double xmax,
                       double size = 0.045, double offset = 0.012)
{
    TAxis* ax = fr->GetXaxis();
    ax->SetMoreLogLabels(kFALSE);
    ax->SetNoExponent(kTRUE);
    if(xmin <= 0.0 || xmax <= xmin) return;

    const double span = std::log10(xmax / xmin);
    if(span > 2.3) return;

    ax->SetLabelSize(0.0);
    pad->cd();
    const double L = pad->GetLeftMargin();
    const double R = pad->GetRightMargin();
    const double B = pad->GetBottomMargin();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    t.SetTextColor(kInk);
    t.SetTextAlign(23);

    const double mult[3] = {1.0, 2.0, 5.0};
    const int d0 = (int)std::floor(std::log10(xmin));
    const int d1 = (int)std::ceil (std::log10(xmax));
    for(int d = d0; d <= d1; ++d)
        for(double m : mult){
            const double v = m * std::pow(10.0, d);
            if(v < xmin * (1.0 - 1e-9) || v > xmax * (1.0 + 1e-9)) continue;
            const double u = L + (1.0 - L - R) * std::log10(v / xmin) / span;
            t.DrawLatex(u, B - offset, fmtTick(v).c_str());
        }
}

static void styleGraph(TGraph* g, int color, int marker, double size, int lw = 1)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(lw);
    g->SetFillStyle(0);
}

static TLegend* makeLegend(double x1, double y1, double x2, double y2, double tsize)
{
    TLegend* lg = new TLegend(x1, y1, x2, y2);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->SetTextFont(42);
    lg->SetTextSize(tsize);
    lg->SetTextColor(kInk);
    lg->SetMargin(0.22);
    return lg;
}

// Leyenda en la franja superior del marco. xFrom/xTo en fraccion del ancho
// util, para poder compartir la franja con otro texto.
static TLegend* topLegend(TVirtualPad* pad, int rows, double rowH, double tsize,
                          double xFrom = 0.0, double xTo = 1.0)
{
    const double L  = pad->GetLeftMargin() + 0.03;
    const double R  = 1.0 - pad->GetRightMargin() - 0.02;
    const double yT = 1.0 - pad->GetTopMargin() - 0.015;
    return makeLegend(L + xFrom * (R - L), yT - rows * rowH, L + xTo * (R - L), yT, tsize);
}

// fraccion del marco que ocupa un bloque de 'rows' lineas de altura rowH (NDC)
static double legendFrac(TVirtualPad* pad, int rows, double rowH)
{
    const double frameH = 1.0 - pad->GetTopMargin() - pad->GetBottomMargin();
    return std::min(0.65, (rows * rowH + 0.035) / frameH);
}

// ymax tal que [ymin, ymax_datos] ocupa la parte inferior (1-frac) del marco
static double withHeadroom(double ymin, double ymax, double frac)
{
    frac = std::max(0.0, std::min(frac, 0.8));
    return ymin + (ymax - ymin) / (1.0 - frac);
}
static double withHeadroomLog(double ymin, double ymax, double frac)
{
    frac = std::max(0.0, std::min(frac, 0.8));
    return ymin * std::pow(ymax / ymin, 1.0 / (1.0 - frac));
}

static TLine* drawHLine(double x1, double x2, double y,
                        int style = 2, int color = kRefGray, int width = 1)
{
    TLine* l = new TLine(x1, y, x2, y);
    l->SetLineStyle(style);
    l->SetLineColor(color);
    l->SetLineWidth(width);
    l->Draw();
    return l;
}

// banda rellena solida (valida en ejes logaritmicos y en cualquier formato)
static TGraph* drawBand(double x1, double x2, double ylo, double yhi, int color = kBandGray)
{
    TGraph* b = new TGraph(4);
    b->SetPoint(0, x1, ylo);
    b->SetPoint(1, x2, ylo);
    b->SetPoint(2, x2, yhi);
    b->SetPoint(3, x1, yhi);
    b->SetFillColor(color);
    b->SetFillStyle(1001);
    b->SetLineColor(color);
    b->SetLineWidth(0);
    b->Draw("F SAME");
    return b;
}

// texto en el margen superior del pad: izquierda y derecha
static void padHeader(const std::string& left, const std::string& right = "",
                      double size = 0.050, int leftColor = kInk)
{
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(size);
    const double y = 1.0 - 0.5 * gPad->GetTopMargin();
    if(!left.empty()){
        t.SetTextColor(leftColor);
        t.SetTextAlign(12);
        t.DrawLatex(gPad->GetLeftMargin(), y, left.c_str());
    }
    if(!right.empty()){
        t.SetTextColor(kInk);
        t.SetTextAlign(32);
        t.DrawLatex(1.0 - gPad->GetRightMargin(), y, right.c_str());
    }
}

// hueco en el margen superior para la etiqueta de figura (si la hay)
static void reserveTag(TVirtualPad* pad, double margin = 0.075)
{
    if(!figureTag().empty())
        pad->SetTopMargin(std::max<double>(pad->GetTopMargin(), margin));
}
static void drawTag(double size = 0.040)
{
    if(!figureTag().empty()) padHeader("", figureTag(), size);
}

// --- etiquetas de energia -----------------------------------------------
static std::string fmtE(double v)
{
    char buf[32];
    if(v >= 100.) std::snprintf(buf, sizeof buf, "%.0f", v);
    else          std::snprintf(buf, sizeof buf, "%.3g", v);
    return buf;
}
// cabeceras de panel:  1.46 < E_n < 2.14 MeV
static std::string energyHeader(double lo, double hi)
{
    return fmtE(lo) + " < E_{n} < " + fmtE(hi) + " MeV";
}
// entradas de leyenda:  [1.46, 2.14] MeV
static std::string energyRange(double lo, double hi)
{
    return "[" + fmtE(lo) + ", " + fmtE(hi) + "] MeV";
}

// ------------------------------------------------------------------------
// Resta de fondo: un panel por bin de energia, en escala logaritmica.
//   espectro crudo    -> area gris claro con contorno gris
//   espectro restado  -> linea de color, sin relleno
//   ajuste del fondo  -> linea discontinua, solo si drawFits = true
// Para el uranio (el ajuste no describe el fondo) llamar con drawFits =
// false; en ese caso el vector fits puede ir vacio.
// En log, los bines <= 0 del espectro restado no aparecen.
// ------------------------------------------------------------------------
static void plotBackgroundFits(
    std::vector<TH1D*>& hists_tof,
    std::vector<TH1D*>& hists_sub,
    std::vector<TF1*>&  fits,
    int   nbins,
    const std::string& outname,
    const std::vector<double>& energy_bins = {},
    bool  drawFits = true)
{
    setPubStyle();

    int ncol, nrow;
    gridLayout(nbins, ncol, nrow);

    TCanvas* c = new TCanvas(uniqueName("c_bkg").c_str(), "Background subtraction",
                             520*ncol, 440*nrow);
    c->Divide(ncol, nrow, 0.001, 0.001);

    const bool haveLabels = (int)energy_bins.size() == nbins + 1;

    for(int i = 0; i < nbins; ++i){
        const std::string label = haveLabels
            ? energyHeader(energy_bins[i], energy_bins[i+1])
            : std::string(Form("bin %d", i));
        const bool fitHere = drawFits && i < (int)fits.size() && fits[i];

        TVirtualPad* pad = c->cd(i + 1);
        stylePad(pad);
        pad->SetLeftMargin(0.20);
        pad->SetBottomMargin(0.16);
        pad->SetTopMargin(0.10);
        pad->SetLogy();

        // copias: no se tocan los histogramas del analisis
        TH1D* hr = (TH1D*)hists_tof[i]->Clone(Form("%s_pub", hists_tof[i]->GetName()));
        TH1D* hs = (TH1D*)hists_sub[i]->Clone(Form("%s_pub", hists_sub[i]->GetName()));
        hr->SetDirectory(nullptr);
        hs->SetDirectory(nullptr);

        // rango en y a partir de los contenidos positivos de ambos
        double yPosMin = 1e30, yMax = 0.0;
        for(TH1D* h : {hr, hs})
            for(int b = 1; b <= h->GetNbinsX(); ++b){
                const double v = h->GetBinContent(b);
                if(v > 0.0){ yPosMin = std::min(yPosMin, v); yMax = std::max(yMax, v); }
            }
        if(yMax <= 0.0){ padHeader(label, "", 0.052); continue; }

        const bool   withLeg = (i == 0);
        const int    nLeg    = fitHere ? 3 : 2;
        const double rowH    = 0.062;
        const double ylo     = 0.5 * yPosMin;
        const double yTop    = withHeadroomLog(ylo, 1.5 * yMax,
                                               withLeg ? legendFrac(pad, nLeg, rowH) : 0.0);

        TH1F* fr = pad->DrawFrame(hr->GetXaxis()->GetXmin(), ylo,
                                  hr->GetXaxis()->GetXmax(), yTop);
        styleFrame(fr, "#Deltat (ns)", "Counts", 1.0, 1.35);

        hr->SetLineColor(kRefGray);
        hr->SetLineWidth(1);
        hr->SetFillColor(kBkgFill);
        hr->SetFillStyle(1001);
        hr->SetMarkerSize(0);
        hr->Draw("HIST SAME");

        hs->SetLineColor(kSubColor);
        hs->SetLineWidth(2);
        hs->SetFillStyle(0);
        hs->SetMarkerSize(0);
        hs->Draw("HIST SAME");

        if(fitHere){
            fits[i]->SetLineColor(kBkgFitColor);
            fits[i]->SetLineStyle(2);
            fits[i]->SetLineWidth(2);
            fits[i]->SetNpx(500);
            fits[i]->Draw("SAME");
        }

        if(withLeg){
            TLegend* lg = topLegend(pad, nLeg, rowH, 0.046, 0.0, 0.85);
            lg->AddEntry(hr, "Raw spectrum",          "f");
            lg->AddEntry(hs, "Background subtracted", "l");
            if(fitHere) lg->AddEntry(fits[i], "Background fit", "l");
            lg->Draw();
        }

        padHeader(label, "", 0.052);
        pad->RedrawAxis();
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Eficiencia frente a cos(theta'), una serie por bin de energia. Colores en
// orden de energia (paleta secuencial), marcadores distintos por serie y
// un pequeno desplazamiento en x para que las barras no se solapen.
// ------------------------------------------------------------------------
static void plotEfficiency(
    const std::vector<EfficiencyResult>& eff,
    int nbins,
    int nbins_det,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    setPubStyle();

    const double w = 1.0 / nbins_det;
    const std::vector<int> col = energyPalette(nbins);

    double ydata = 1.0;
    for(int e = 0; e < nbins; ++e)
        for(int i = 0; i < nbins_det; ++i){
            const double v = eff[e].eps[i] + eff[e].u_eps[i];
            if(std::isfinite(v)) ydata = std::max(ydata, v);
        }

    TCanvas* c = new TCanvas(uniqueName("c_eff").c_str(), "Efficiency", 820, 640);
    stylePad(c);
    reserveTag(c);

    const int    ncols = nbins > 12 ? 3 : (nbins > 4 ? 2 : 1);
    const int    rows  = (nbins + ncols - 1) / ncols;
    const double rowH  = 0.046;
    const double yTop  = withHeadroom(0.0, 1.03 * ydata, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(0.0, 0.0, 1.0, yTop);
    styleFrame(fr, "cos#theta'", "#varepsilon(cos#theta')");
    drawHLine(0.0, 1.0, 1.0);

    TLegend* lg = topLegend(c, rows, rowH, 0.034, 0.0, ncols == 1 ? 0.55 : 1.0);
    lg->SetNColumns(ncols);

    for(int e = 0; e < nbins; ++e){
        std::vector<double> x(nbins_det);
        for(int i = 0; i < nbins_det; ++i)
            x[i] = (i + 0.5) * w + dodge(e, nbins, 0.45 * w);
        TGraphErrors* gr = new TGraphErrors(
            nbins_det, x.data(), eff[e].eps.data(), nullptr, eff[e].u_eps.data());
        const int mk = kSeriesMarkers[e % 6];
        styleGraph(gr, col[e], mk, 1.0 * markerScale(mk), 2);
        gr->Draw("P");
        lg->AddEntry(gr, energyRange(energy_bins[e], energy_bins[e+1]).c_str(), "pe");
    }
    lg->Draw();
    drawTag();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Resolucion relativa de la eficiencia, sigma(eps)/eps. Mismo formato,
// colores y marcadores que plotEfficiency, para leerlas como pareja.
// ------------------------------------------------------------------------
static void plotEfficiencyResolution(
    const std::vector<EfficiencyResult>& eff,
    int nbins,
    int nbins_det,
    const std::vector<double>& energy_bins,
    const std::string& outname,
    bool logy = true)
{
    setPubStyle();

    std::vector<double> centers(nbins_det);
    for(int i = 0; i < nbins_det; ++i) centers[i] = (i + 0.5) * (1.0 / nbins_det);

    const std::vector<int> col = energyPalette(nbins);

    std::vector<TGraph*> graphs;
    double ymax = 0.0, yminPos = 1e30;
    for(int e = 0; e < nbins; ++e){
        std::vector<double> x, y;
        for(int i = 0; i < nbins_det; ++i){
            const double eps = eff[e].eps[i];
            if(eps <= 0.0) continue;
            const double rel = eff[e].u_eps[i] / eps;
            if(!std::isfinite(rel)) continue;
            x.push_back(centers[i]);
            y.push_back(rel);
            ymax = std::max(ymax, rel);
            if(rel > 0.0) yminPos = std::min(yminPos, rel);
        }
        TGraph* gr = new TGraph((int)x.size(), x.data(), y.data());
        const int mk = kSeriesMarkers[e % 6];
        styleGraph(gr, col[e], mk, 1.0 * markerScale(mk), 2);
        graphs.push_back(gr);
    }
    if(ymax <= 0.0){
        std::cerr << "[WARN] plotEfficiencyResolution: nothing to draw\n";
        return;
    }

    TCanvas* c = new TCanvas(uniqueName("c_eff_res").c_str(), "Efficiency resolution", 820, 640);
    stylePad(c);
    reserveTag(c);
    if(logy) c->SetLogy();

    const int    ncols = nbins > 12 ? 3 : (nbins > 4 ? 2 : 1);
    const int    rows  = (nbins + ncols - 1) / ncols;
    const double rowH  = 0.046;
    const double frac  = legendFrac(c, rows, rowH);

    double ylo, yTop;
    if(logy){
        ylo  = 0.5 * yminPos;
        yTop = withHeadroomLog(ylo, 1.5 * ymax, frac);
    } else {
        ylo  = 0.0;
        yTop = withHeadroom(0.0, 1.05 * ymax, frac);
    }

    TH1F* fr = c->DrawFrame(0.0, ylo, 1.0, yTop);
    styleFrame(fr, "cos#theta'", "#sigma_{#varepsilon} / #varepsilon");
    if(logy) fr->GetYaxis()->SetMoreLogLabels(kTRUE);

    TLegend* lg = topLegend(c, rows, rowH, 0.034, 0.0, ncols == 1 ? 0.55 : 1.0);
    lg->SetNColumns(ncols);
    for(int e = 0; e < nbins; ++e){
        if(graphs[e]->GetN() == 0) continue;
        graphs[e]->Draw("LP");
        lg->AddEntry(graphs[e], energyRange(energy_bins[e], energy_bins[e+1]).c_str(), "lp");
    }
    lg->Draw();
    drawTag();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// W(theta) por bin de energia. Todos los paneles comparten el rango en y,
// para que se comparen de un vistazo.
// ------------------------------------------------------------------------
static void plotAnisotropy(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    int nbins_beam,
    const std::vector<double>& energy_bins,
    const std::string& outname)
{
    setPubStyle();

    double gmin = 1.0, gmax = 1.0;
    for(int e = 0; e < nbins; ++e)
        for(int i = 0; i < nbins_beam; ++i){
            const double v = aniso[e].w[i], u = aniso[e].u_w[i];
            if(!std::isfinite(v) || !std::isfinite(u)) continue;
            gmin = std::min(gmin, v - u);
            gmax = std::max(gmax, v + u);
        }
    const double m   = 0.08 * (gmax - gmin);
    const double ylo = std::max(0.0, gmin - m);
    const double yhi = gmax + m;

    int ncol, nrow;
    gridLayout(nbins, ncol, nrow);

    TCanvas* c = new TCanvas(uniqueName("c_aniso").c_str(), "Anisotropy", 500*ncol, 440*nrow);
    c->Divide(ncol, nrow, 0.001, 0.001);

    for(int e = 0; e < nbins; ++e){
        std::vector<double> x(nbins_beam), y(nbins_beam), ex(nbins_beam, 0.0);
        for(int i = 0; i < nbins_beam; ++i){
            x[i] = (i + 0.5) * dcos_beam;
            y[i] = aniso[e].w[i];
        }

        TVirtualPad* pad = c->cd(e + 1);
        stylePad(pad);
        pad->SetLeftMargin(0.20);
        pad->SetBottomMargin(0.16);
        pad->SetTopMargin(0.10);

        TH1F* fr = pad->DrawFrame(0.0, ylo, 1.0, yhi);
        styleFrame(fr, "cos#theta_{beam}", "W(#theta)/W(90^{#circ})", 1.0, 1.35);
        drawHLine(0.0, 1.0, 1.0);

        TGraphErrors* g = new TGraphErrors(
            nbins_beam, x.data(), y.data(), ex.data(), aniso[e].u_w.data());
        styleGraph(g, kAnisoColor, 20, 1.0, 2);
        g->Draw("P");

        padHeader(energyHeader(energy_bins[e], energy_bins[e+1]), "", 0.056);
        pad->RedrawAxis();
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// W(0)/W(90) frente a la energia. En la figura, las barras horizontales
// cubren el bin real (asimetricas en escala log); en el .root se guarda el
// mismo TGraphErrors que antes para no romper nada aguas abajo.
// ------------------------------------------------------------------------
static void plotAnisotropyRatio(
    const std::vector<AnisotropyResult>& aniso,
    int nbins,
    int nbins_beam,
    const std::vector<double>& energy_bins,
    const std::string& outname,
    const std::string& reaction_label = "^{238}U(n,f)")
{
    setPubStyle();

    const int bin_0 = nbins_beam - 1;

    std::vector<double> E_centers(nbins), ratio(nbins), u_ratio(nbins), ex(nbins);
    for(int e = 0; e < nbins; ++e){
        E_centers[e] = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        ratio[e]     = aniso[e].w[bin_0];
        u_ratio[e]   = aniso[e].u_w[bin_0];
        ex[e]        = 0.5 * (energy_bins[e+1] - energy_bins[e]);
    }
    TGraphErrors* g = new TGraphErrors(
        nbins, E_centers.data(), ratio.data(), ex.data(), u_ratio.data());

    // grafico para dibujar: barra horizontal = anchura real del bin
    TGraphAsymmErrors* gd = new TGraphAsymmErrors(nbins);
    double ymin = 0.95, ymax = 1.05;
    for(int e = 0; e < nbins; ++e){
        gd->SetPoint(e, E_centers[e], ratio[e]);
        gd->SetPointError(e, E_centers[e] - energy_bins[e], energy_bins[e+1] - E_centers[e],
                          u_ratio[e], u_ratio[e]);
        ymin = std::min(ymin, ratio[e] - u_ratio[e]);
        ymax = std::max(ymax, ratio[e] + u_ratio[e]);
    }
    styleGraph(gd, kRatioColor, 20, 1.1, 2);
    styleGraph(g,  kRatioColor, 20, 1.1);

    const double xmin = energy_bins.front(), xmax = energy_bins.back();

    TCanvas* c = new TCanvas(uniqueName("c_ratio").c_str(), "", 820, 620);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const double m    = 0.08 * (ymax - ymin);
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, 3, 0.052));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", "W(0^{#circ}) / W(90^{#circ})");
    logXLabels(c, fr, xmin, xmax);

    TGraph* band = drawBand(xmin, xmax, 0.95, 1.05);
    TLine*  line = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);
    gd->Draw("P");

    TLegend* lg = topLegend(c, 3, 0.052, 0.038, 0.0, 0.62);
    lg->AddEntry(gd,   (reaction_label + "  W(0^{#circ})/W(90^{#circ})").c_str(), "pe");
    lg->AddEntry(line, "Isotropic",  "l");
    lg->AddEntry(band, "#pm5% band", "f");
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    c->SaveAs((outname + ".pdf").c_str());

    TFile* fout = TFile::Open((outname + ".root").c_str(), "RECREATE");
    if(!fout || fout->IsZombie()){
        std::cerr << "[ERROR] Cannot create output ROOT file: " << outname << ".root\n";
        return;
    }
    g->Write("anisotropy_ratio");
    line->Write("isotropic_line");
    band->Write("five_percent_band");
    fout->Close();

    std::cout << "[INFO] Saved: " << outname << ".pdf\n";
    std::cout << "[INFO] Saved: " << outname << ".root\n";
}

// ========================================================================
//  Comparacion con EXFOR
//
//  Los CSV de EXFOR para magnitudes DA no tienen columnas fijas (el error
//  puede llamarse "DATA-ERR (NO-DIM)" o "ERR-T (NO-DIM)", la resolucion
//  puede ser semi-anchura o anchura completa): el lector busca por PREFIJO.
// ========================================================================

// CSV con campos entre comillas que contienen comas (columna Reacode)
static std::vector<std::string> parseCsvLine(const std::string& line)
{
    std::vector<std::string> fields;
    std::string cur;
    bool inQuotes = false;
    for(char ch : line){
        if(ch == '"'){ inQuotes = !inQuotes; continue; }
        if(ch == ',' && !inQuotes){ fields.push_back(cur); cur.clear(); continue; }
        cur += ch;
    }
    fields.push_back(cur);
    return fields;
}

struct ExforSource {
    std::string path;
    std::string label = "";   // vacio: se construye como "Autor (anio)"
};

struct ExforData {
    TGraphErrors* graph = nullptr;
    std::string   label;
};

static ExforData loadExforAniso(const ExforSource& src)
{
    ExforData out;

    std::ifstream f(src.path);
    if(!f.is_open()){
        std::cerr << "[ERROR] Cannot open " << src.path << "\n";
        return out;
    }

    std::string header_line;
    std::getline(f, header_line);
    std::vector<std::string> header = parseCsvLine(header_line);

    auto findCol = [&](std::initializer_list<const char*> candidates) -> int {
        for(auto cand : candidates)
            for(size_t i = 0; i < header.size(); ++i)
                if(header[i].rfind(cand, 0) == 0)
                    return (int)i;
        return -1;
    };

    int iID    = findCol({"DatasetID"});
    int iYear  = findCol({"year1"});
    int iAuth  = findCol({"author1"});
    int iEN    = findCol({"EN (EV)"});
    int iDATA  = findCol({"DATA (NO-DIM)"});
    int iERR_abs = findCol({"DATA-ERR (NO-DIM)", "ERR-T (NO-DIM)", "ERR (NO-DIM)"});
    int iERR_pct = findCol({"DATA-ERR (PER-CENT)", "ERR-T (PER-CENT)", "ERR (PER-CENT)"});
    int iERR     = (iERR_abs >= 0) ? iERR_abs : iERR_pct;
    bool errIsPercent = (iERR_abs < 0) && (iERR_pct >= 0);
    int iRSL   = findCol({"EN-RSL-HW (EV)", "EN-RSL (EV)"});
    bool rslIsHalfWidth = (iRSL >= 0) && header[iRSL].rfind("EN-RSL-HW", 0) == 0;

    if(iEN < 0 || iDATA < 0 || iERR < 0){
        std::cerr << "[ERROR] " << src.path << ": required columns not found\n";
        std::cerr << "        looked for EN (EV): " << (iEN>=0?"found":"MISSING")
                  << " | DATA (NO-DIM): " << (iDATA>=0?"found":"MISSING")
                  << " | error column: "  << (iERR>=0?"found":"MISSING") << "\n";
        std::cerr << "        header columns actually present in this file:\n";
        for(size_t i = 0; i < header.size(); ++i)
            std::cerr << "          [" << i << "] '" << header[i] << "'\n";
        return out;
    }

    std::vector<double> ex, ey, ex_err, ey_err;
    std::string line, ds_id, author, year;

    while(std::getline(f, line)){
        if(line.empty()) continue;
        auto fld = parseCsvLine(line);
        if((int)fld.size() <= std::max({iEN, iDATA, iERR})) continue;

        if(ds_id.empty()  && iID   >= 0 && (int)fld.size() > iID)   ds_id  = fld[iID];
        if(author.empty() && iAuth >= 0 && (int)fld.size() > iAuth) author = fld[iAuth];
        if(year.empty()   && iYear >= 0 && (int)fld.size() > iYear) year   = fld[iYear];

        try{
            double E   = std::stod(fld[iEN]) * 1e-6;    // eV -> MeV
            double y   = std::stod(fld[iDATA]);
            double eyv = std::stod(fld[iERR]);
            if(errIsPercent) eyv = y * eyv / 100.0;
            double exv = 0.0;
            if(iRSL >= 0 && (int)fld.size() > iRSL && !fld[iRSL].empty()){
                double rsl = std::stod(fld[iRSL]) * 1e-6;
                exv = rslIsHalfWidth ? rsl : rsl / 2.0;
            }
            ex.push_back(E); ey.push_back(y);
            ex_err.push_back(exv); ey_err.push_back(eyv);
        } catch(...){ continue; }
    }

    if(ex.empty()){
        std::cerr << "[WARN] " << src.path << ": no valid data rows parsed\n";
        return out;
    }

    out.graph = new TGraphErrors((int)ex.size(),
        ex.data(), ey.data(), ex_err.data(), ey_err.data());

    if(!src.label.empty()){
        out.label = src.label;
    } else {
        out.label = !author.empty() ? author
                  : (!ds_id.empty() ? "EXFOR " + ds_id : std::string("EXFOR"));
        if(year.size() >= 4) out.label += " (" + year.substr(0, 4) + ")";
    }
    return out;
}

// Metrica rapida de compatibilidad (solo el error de EXFOR): ver comentario
// de computePulls para la version con ambos errores.
static bool chi2PerPointVsThis(TGraphErrors* ext, TGraphErrors* mine,
                               double& chi2_per_n, int& n_used)
{
    int nMine = mine->GetN();
    if(nMine < 2) return false;
    double* mx = mine->GetX();
    double* my = mine->GetY();

    double sum = 0.0; int n = 0;
    for(int i = 0; i < ext->GetN(); ++i){
        double x = ext->GetX()[i];
        if(x < mx[0] || x > mx[nMine-1]) continue;

        int k = (int)(std::lower_bound(mx, mx+nMine, x) - mx);
        if(k == 0) k = 1;
        double x0 = mx[k-1], x1 = mx[k], y0 = my[k-1], y1 = my[k];
        double t = (std::log(x) - std::log(x0)) / (std::log(x1) - std::log(x0));
        double yInterp = y0 + t * (y1 - y0);

        double ey = ext->GetEY()[i];
        if(ey <= 0.0) continue;
        double pull = (ext->GetY()[i] - yInterp) / ey;
        sum += pull*pull;
        ++n;
    }
    if(n == 0) return false;
    chi2_per_n = sum / n;
    n_used = n;
    return true;
}

// ------------------------------------------------------------------------
// Este trabajo frente a varios conjuntos de EXFOR, en una figura.
// ------------------------------------------------------------------------
static void plotAnisoVsExfor(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.0, double xmax = 2000.0)
{
    setPubStyle();

    if(!g_this){
        std::cerr << "[ERROR] plotAnisoVsExfor: this-work graph is null\n";
        return;
    }
    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }

    double ymin = 1.0, ymax = 1.0;
    auto updateRange = [&](TGraphErrors* g){
        for(int i = 0; i < g->GetN(); ++i){
            const double x = g->GetX()[i];
            if(x < xmin || x > xmax) continue;
            const double v = g->GetY()[i], e = g->GetEY()[i];
            ymin = std::min(ymin, v - e);
            ymax = std::max(ymax, v + e);
        }
    };
    updateRange(g_this);
    for(auto& d : exfor) updateRange(d.graph);
    const double m = 0.08 * (ymax - ymin);

    const int nEntries = (int)exfor.size() + 2;
    const int ncols    = nEntries > 4 ? 2 : 1;
    const int rows     = (nEntries + ncols - 1) / ncols;
    const double rowH  = 0.050;

    TCanvas* c = new TCanvas(uniqueName("c_aniso_exfor").c_str(), "Anisotropy vs EXFOR", 880, 640);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", "W(0^{#circ}) / W(90^{#circ})");
    logXLabels(c, fr, xmin, xmax);

    TLine* line = drawHLine(xmin, xmax, 1.0);

    TLegend* lg = topLegend(c, rows, rowH, 0.036);
    lg->SetNColumns(ncols);

    // este trabajo: estilo fijado ya, se dibuja al final (primer plano)
    styleGraph(g_this, kThisWorkColor, 21, 1.3, 2);
    lg->AddEntry(g_this, "This work", "pe");

    const std::vector<int> pal = exforPalette();
    for(size_t k = 0; k < exfor.size(); ++k){
        styleGraph(exfor[k].graph, pal[k % pal.size()], kOpenMarkers[k % 6], 1.0, 1);
        exfor[k].graph->Draw("P");
        lg->AddEntry(exfor[k].graph, exfor[k].label.c_str(), "pe");
    }
    g_this->Draw("P");

    lg->AddEntry(line, "Isotropic", "l");
    lg->Draw();
    drawTag();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Pull por punto:  (y_ext - y_this_interp) / sqrt(sigma_ext^2 + sigma_this^2)
// Valor e incertidumbre de este trabajo se interpolan linealmente en log E.
// Un conjunto compatible da pulls centrados en 0 con RMS cercano a 1.
// ------------------------------------------------------------------------
struct PullSeries {
    TGraph* graph = nullptr;   // x = E_n (MeV), y = pull
    double  mean  = 0.0;
    double  rms   = 0.0;
    int     n     = 0;
};

static PullSeries computePulls(TGraphErrors* ext, TGraphErrors* mine)
{
    PullSeries out;
    int nMine = mine->GetN();
    if(nMine < 2) return out;

    double* mx  = mine->GetX();
    double* my  = mine->GetY();
    double* mey = mine->GetEY();

    std::vector<double> px, py;
    double sum = 0.0, sum2 = 0.0;

    for(int i = 0; i < ext->GetN(); ++i){
        double x = ext->GetX()[i];
        if(x < mx[0] || x > mx[nMine-1]) continue;

        int k = (int)(std::lower_bound(mx, mx+nMine, x) - mx);
        if(k == 0) k = 1;
        double x0 = mx[k-1], x1 = mx[k];
        double t  = (std::log(x) - std::log(x0)) / (std::log(x1) - std::log(x0));

        double yInterp  = my[k-1]  + t * (my[k]  - my[k-1]);
        double eyInterp = mey[k-1] + t * (mey[k] - mey[k-1]);

        double ey_ext = ext->GetEY()[i];
        double sigma  = std::sqrt(ey_ext*ey_ext + eyInterp*eyInterp);
        if(sigma <= 0.0) continue;

        double pull = (ext->GetY()[i] - yInterp) / sigma;
        px.push_back(x);
        py.push_back(pull);
        sum  += pull;
        sum2 += pull*pull;
    }

    if(px.empty()) return out;

    out.n     = (int)px.size();
    out.mean  = sum / out.n;
    out.rms   = std::sqrt(sum2 / out.n);
    out.graph = new TGraph(out.n, px.data(), py.data());
    return out;
}

// dibuja el fondo comun de un panel de pulls: banda +-1, lineas 0 y +-2
static void drawPullGuides(double xmin, double xmax)
{
    drawBand(xmin, xmax, -1.0, 1.0);
    drawHLine(xmin, xmax,  0.0, 1, kRefGray);
    drawHLine(xmin, xmax,  2.0, 3, kRefGray);
    drawHLine(xmin, xmax, -2.0, 3, kRefGray);
}

// ------------------------------------------------------------------------
// Pulls, un panel por conjunto de EXFOR.
// ------------------------------------------------------------------------
static void plotPullsVsExfor(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 2000.0)
{
    setPubStyle();

    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;

    int ncol, nrow;
    gridLayout((int)exfor.size(), ncol, nrow);

    TCanvas* c = new TCanvas(uniqueName("c_pulls").c_str(), "Pulls vs EXFOR", 520*ncol, 440*nrow);
    c->Divide(ncol, nrow, 0.001, 0.001);

    const std::vector<int> pal = exforPalette();

    for(size_t k = 0; k < exfor.size(); ++k){
        const int color = pal[k % pal.size()];
        PullSeries p = computePulls(exfor[k].graph, g_this);

        TVirtualPad* pad = c->cd((int)k + 1);
        stylePad(pad);
        pad->SetLeftMargin(0.20);
        pad->SetBottomMargin(0.16);
        pad->SetTopMargin(0.10);
        pad->SetLogx();

        double ymax = 3.0;
        if(p.graph)
            for(int i = 0; i < p.graph->GetN(); ++i)
                ymax = std::max(ymax, 1.2 * std::abs(p.graph->GetY()[i]));
        // hueco arriba para la linea de estadisticos
        const double yTop = withHeadroom(-ymax, ymax, legendFrac(pad, 1, 0.06));

        TH1F* fr = pad->DrawFrame(xmin, -ymax, xmax, yTop);
        styleFrame(fr, "E_{n} (MeV)", "(data #minus this work) / #sigma", 1.0, 1.35);
        logXLabels(pad, fr, xmin, xmax);
        drawPullGuides(xmin, xmax);

        if(p.graph){
            styleGraph(p.graph, color, 20, 0.9);
            p.graph->Draw("P");

            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.046);
            t.SetTextColor(kInk); t.SetTextAlign(13);
            t.DrawLatex(pad->GetLeftMargin() + 0.03, 1.0 - pad->GetTopMargin() - 0.02,
                Form("#LTpull#GT = %.2f    RMS = %.2f    N = %d", p.mean, p.rms, p.n));
        }
        padHeader(exfor[k].label, "", 0.055, color);
        pad->RedrawAxis();
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Pulls de todos los conjuntos superpuestos.
// ------------------------------------------------------------------------
static void plotPullsOverlay(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 2000.0)
{
    setPubStyle();

    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;

    std::vector<PullSeries> pulls;
    double ymax = 3.0;
    for(auto& d : exfor){
        PullSeries p = computePulls(d.graph, g_this);
        if(p.graph)
            for(int i = 0; i < p.graph->GetN(); ++i)
                ymax = std::max(ymax, 1.2 * std::abs(p.graph->GetY()[i]));
        pulls.push_back(p);
    }

    TCanvas* c = new TCanvas(uniqueName("c_pulls_overlay").c_str(), "Pulls overlay", 880, 640);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const int    rows = (int)exfor.size();
    const double rowH = 0.048;
    const double yTop = withHeadroom(-ymax, ymax, legendFrac(c, rows, rowH));

    TH1F* fr = c->DrawFrame(xmin, -ymax, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", "(data #minus this work) / #sigma");
    logXLabels(c, fr, xmin, xmax);
    drawPullGuides(xmin, xmax);

    const std::vector<int> pal = exforPalette();
    TLegend* lg = topLegend(c, rows, rowH, 0.034);

    for(size_t k = 0; k < exfor.size(); ++k){
        if(!pulls[k].graph) continue;
        styleGraph(pulls[k].graph, pal[k % pal.size()], kOpenMarkers[k % 6], 1.0, 1);
        pulls[k].graph->Draw("P");
        lg->AddEntry(pulls[k].graph,
            Form("%s   #LTpull#GT = %.2f,  RMS = %.2f",
                 exfor[k].label.c_str(), pulls[k].mean, pulls[k].rms), "p");
    }
    lg->Draw();
    drawTag();
    c->RedrawAxis();
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Este trabajo frente a cada conjunto de EXFOR, un panel por conjunto.
// ------------------------------------------------------------------------
static void plotAnisoVsExforIndividual(
    TGraphErrors* g_this,
    const std::vector<ExforSource>& sources,
    const std::string& outname,
    double xmin = 1.1, double xmax = 1000.0)
{
    setPubStyle();

    std::vector<ExforData> exfor;
    for(auto& s : sources){
        ExforData d = loadExforAniso(s);
        if(d.graph) exfor.push_back(d);
    }
    if(exfor.empty() || !g_this) return;

    int ncol, nrow;
    gridLayout((int)exfor.size(), ncol, nrow);

    TCanvas* c = new TCanvas(uniqueName("c_aniso_exfor_grid").c_str(),
                             "This work vs each EXFOR set", 520*ncol, 440*nrow);
    c->Divide(ncol, nrow, 0.001, 0.001);

    const std::vector<int> pal = exforPalette();

    for(size_t k = 0; k < exfor.size(); ++k){
        const int color = pal[k % pal.size()];

        double ymin = 1.0, ymax = 1.0;
        auto updateRange = [&](TGraphErrors* g){
            for(int i = 0; i < g->GetN(); ++i){
                const double x = g->GetX()[i];
                if(x < xmin || x > xmax) continue;
                const double v = g->GetY()[i], e = g->GetEY()[i];
                ymin = std::min(ymin, v - e);
                ymax = std::max(ymax, v + e);
            }
        };
        updateRange(g_this);
        updateRange(exfor[k].graph);
        const double m = 0.08 * (ymax - ymin);

        TVirtualPad* pad = c->cd((int)k + 1);
        stylePad(pad);
        pad->SetLeftMargin(0.20);
        pad->SetBottomMargin(0.16);
        pad->SetTopMargin(0.10);
        pad->SetLogx();

        const double ylo  = ymin - m;
        const double yTop = withHeadroom(ylo, ymax + m, legendFrac(pad, 1, 0.065));

        TH1F* fr = pad->DrawFrame(xmin, ylo, xmax, yTop);
        styleFrame(fr, "E_{n} (MeV)", "W(0^{#circ})/W(90^{#circ})", 1.0, 1.35);
        logXLabels(pad, fr, xmin, xmax);
        drawHLine(xmin, xmax, 1.0);

        TGraphErrors* ge = (TGraphErrors*)exfor[k].graph->Clone();
        styleGraph(ge, color, 24, 0.9, 1);
        ge->Draw("P");

        TGraphErrors* gt = (TGraphErrors*)g_this->Clone();
        styleGraph(gt, kThisWorkColor, 21, 1.0, 2);
        gt->Draw("P");

        TLegend* lg = topLegend(pad, 1, 0.065, 0.046);
        lg->SetNColumns(2);
        lg->AddEntry(gt, "This work", "pe");
        lg->AddEntry(ge, exfor[k].label.c_str(), "pe");
        lg->Draw();

        double chi2n; int nUsed;
        const std::string right = chi2PerPointVsThis(exfor[k].graph, g_this, chi2n, nUsed)
            ? std::string(Form("#chi^{2}/N #approx %.2f  (N = %d)", chi2n, nUsed))
            : std::string("");
        padHeader(exfor[k].label, right, 0.050, color);
        pad->RedrawAxis();
    }
    c->SaveAs(outname.c_str());
}

// ------------------------------------------------------------------------
// Curva normalizada W(cos theta)/W(90) de la serie de Legendre
// ------------------------------------------------------------------------
static double legendreRatioCurve(double* x, double* p)
{
    double c = x[0];
    double N = 1. + p[0]*legP2(c)  + p[1]*legP4(c);
    double D = 1. + p[0]*legP2(0.) + p[1]*legP4(0.);
    return N / D;
}

// ------------------------------------------------------------------------
// W(0)/W(90) del ajuste de Legendre frente a la energia.
// Guarda <outname>.pdf y <outname>.root (grafico "anisotropy_ratio") y
// devuelve el grafico para plotAnisoVsExfor / plotPulls*.
// ------------------------------------------------------------------------
static TGraphErrors* plotAnisotropyRatioFit(
    const std::vector<LegendreResult>& leg,
    const std::vector<double>& energy_bins,
    const std::string& outname,
    const std::string& reaction_label = "W(0^{#circ})/W(90^{#circ})")
{
    setPubStyle();

    std::vector<double> x, y, ex, ey, xl, xh;
    for(int e = 0; e < (int)leg.size(); ++e){
        if(!leg[e].valid) continue;
        const double xc = std::sqrt(energy_bins[e] * energy_bins[e+1]);
        x.push_back(xc);
        ex.push_back(0.0);
        xl.push_back(xc - energy_bins[e]);
        xh.push_back(energy_bins[e+1] - xc);
        y.push_back(leg[e].anisotropy);
        ey.push_back(leg[e].u_anisotropy);
    }
    if(x.empty()){
        std::cerr << "[WARN] plotAnisotropyRatioFit: no valid fits\n";
        return nullptr;
    }

    // resultado (se devuelve y se guarda como siempre)
    TGraphErrors* g = new TGraphErrors(
        (int)x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetName("anisotropy_ratio");
    styleGraph(g, kRatioColor, 20, 1.1);

    // para dibujar: barra horizontal = bin de energia
    TGraphAsymmErrors* gd = new TGraphAsymmErrors(
        (int)x.size(), x.data(), y.data(), xl.data(), xh.data(), ey.data(), ey.data());
    styleGraph(gd, kRatioColor, 20, 1.15, 2);

    double ymin = 0.95, ymax = 1.05;
    for(size_t i = 0; i < y.size(); ++i){
        ymin = std::min(ymin, y[i] - ey[i]);
        ymax = std::max(ymax, y[i] + ey[i]);
    }

    const double xmin = energy_bins.front(), xmax = energy_bins.back();

    TCanvas* c = new TCanvas(uniqueName("c_ratio_fit").c_str(), "", 820, 620);
    stylePad(c);
    reserveTag(c);
    c->SetLogx();

    const double m    = 0.08 * (ymax - ymin);
    const double ylo  = ymin - m;
    const double yTop = withHeadroom(ylo, ymax + m, legendFrac(c, 3, 0.052));

    TH1F* fr = c->DrawFrame(xmin, ylo, xmax, yTop);
    styleFrame(fr, "E_{n} (MeV)", "W(0^{#circ}) / W(90^{#circ})");
    logXLabels(c, fr, xmin, xmax);

    TGraph* band = drawBand(xmin, xmax, 0.95, 1.05);
    TLine*  line = drawHLine(xmin, xmax, 1.0, 7, kRefGray, 2);
    gd->Draw("P");

    TLegend* lg = topLegend(c, 3, 0.052, 0.038, 0.0, 0.72);
    lg->AddEntry(gd,   (reaction_label + "  (Legendre fit)").c_str(), "pe");
    lg->AddEntry(line, "Isotropic",  "l");
    lg->AddEntry(band, "#pm5% band", "f");
    lg->Draw();

    drawTag();
    c->RedrawAxis();
    c->SaveAs((outname + ".pdf").c_str());

    TFile* fout = TFile::Open((outname + ".root").c_str(), "RECREATE");
    if(fout && !fout->IsZombie()){
        g->Write("anisotropy_ratio");
        fout->Close();
    } else {
        std::cerr << "[ERROR] Cannot create " << outname << ".root\n";
    }
    return g;
}

// ------------------------------------------------------------------------
// W(cos theta)/W(90) por bin de energia: puntos + ajuste de Legendre y,
// debajo, los residuos normalizados del ajuste
//
//      pull_k = (y_k - f(cos theta_k)) / sigma_k ,
//
// con y_k = w_k/W90 y sigma_k = u_w_k/W90. Dividir por W90 no cambia el
// cociente, asi que la suma de pull^2 es el chi^2 del ajuste.
//
// Colores: datos en azul, ajuste en terracota, residuos en verde azulado.
// showResiduals = false da la figura sin panel de residuos.
// La leyenda (solo panel 1) ocupa una fila en dos columnas y los
// parametros del ajuste van con letra algo menor: el bloque de texto
// quita menos altura a los datos que antes.
// ------------------------------------------------------------------------
static void plotAnisotropyFit(
    const std::vector<LegendreResult>& leg,
    const std::vector<double>& energy_bins,
    const std::string& outname,
    int  step          = 1,
    bool showResiduals = true)
{
    setPubStyle();

    std::vector<int> sel;
    for(int e = 0; e < (int)leg.size(); e += std::max(step, 1))
        if(leg[e].valid && leg[e].W90 > 0.) sel.push_back(e);
    if(sel.empty()){
        std::cerr << "[WARN] plotAnisotropyFit: no valid fits to draw\n";
        return;
    }

    int ncol, nrow;
    gridLayout((int)sel.size(), ncol, nrow);

    const int cellW = 520;
    const int cellH = showResiduals ? 590 : 450;
    TCanvas* c = new TCanvas(uniqueName("c_aniso_fit").c_str(), "Anisotropy fits",
                             cellW*ncol, cellH*nrow);
    c->Divide(ncol, nrow, 0.001, 0.001);

    // reparto de la celda y factores de escala de la letra
    const double fUp = showResiduals ? 0.70 : 1.0;
    const double fLo = 1.0 - fUp;
    const double sUp = 1.0 / fUp;
    const double sLo = showResiduals ? 1.0 / fLo : 1.0;

    // margenes en fraccion de la altura de CELDA
    const double kTopCell = 0.085;
    const double kBotCell = 0.135;
    const double kLeft    = 0.21;
    const double kRight   = 0.04;

    for(size_t p = 0; p < sel.size(); ++p){
        const int e = sel[p];
        const LegendreResult& L = leg[e];

        // --- puntos normalizados -------------------------------------------
        const int n = (int)L.w.size();
        std::vector<double> y(n), ey(n), ex(n, 0.0);
        double ymin = 1.0, ymax = std::max(1.0, L.anisotropy);
        for(int k = 0; k < n; ++k){
            y[k]  = L.w[k]   / L.W90;
            ey[k] = L.u_w[k] / L.W90;
            ymin = std::min(ymin, y[k] - ey[k]);
            ymax = std::max(ymax, y[k] + ey[k]);
        }
        const double m   = 0.08 * (ymax - ymin);
        const double ylo = std::max(0.0, ymin - m);
        const double yhi = ymax + m;

        // --- residuos normalizados -----------------------------------------
        std::vector<double> px, py;
        double pmax = 0.0;
        double par[2] = {L.a2, L.a4};
        for(int k = 0; k < n; ++k){
            if(ey[k] <= 0.0) continue;
            double xk = L.cos_theta[k];
            const double fk = legendreRatioCurve(&xk, par);
            const double pk = (y[k] - fk) / ey[k];
            px.push_back(xk);
            py.push_back(pk);
            pmax = std::max(pmax, std::abs(pk));
        }

        // --- sub-pads -------------------------------------------------------
        TVirtualPad* cell = c->cd((int)p + 1);
        TVirtualPad* pUp  = cell;
        TPad*        pLo  = nullptr;
        if(showResiduals){
            cell->cd();
            TPad* up = new TPad(uniqueName(Form("pad_up_%d", e)).c_str(), "", 0.0, fLo, 1.0, 1.0);
            TPad* lo = new TPad(uniqueName(Form("pad_lo_%d", e)).c_str(), "", 0.0, 0.0, 1.0, fLo);
            up->Draw();
            lo->Draw();
            pUp = up;
            pLo = lo;
        }

        // ===================== panel principal ==============================
        pUp->cd();
        stylePad(pUp);
        pUp->SetLeftMargin(kLeft);
        pUp->SetRightMargin(kRight);
        pUp->SetTopMargin(kTopCell * sUp);
        pUp->SetBottomMargin(showResiduals ? 0.015 : kBotCell);

        // bloque de texto: [leyenda en una fila, solo panel 1] + a2, (a4),
        // W(0)/W(90), chi2
        const bool   showA4  = L.u_a4 > 0.0;
        const bool   withLeg = (p == 0);
        const int    nText   = (showA4 ? 4 : 3) + (withLeg ? 1 : 0);
        const double tsz     = 0.038 * sUp;
        const double lineH   = 0.052 * sUp;
        const double yTop    = withHeadroom(ylo, yhi, legendFrac(pUp, nText, lineH));

        TH1F* frame = pUp->DrawFrame(0.0, ylo, 1.0, yTop);
        styleFrame(frame, "cos#theta_{beam}", "W(#theta)/W(90^{#circ})", 0.92 * sUp, 1.30);
        if(showResiduals){
            frame->GetXaxis()->SetLabelSize(0.0);
            frame->GetXaxis()->SetTitleSize(0.0);
        }
        drawHLine(0.0, 1.0, 1.0);

        TF1* f = new TF1(uniqueName(Form("f_legfit_%d", e)).c_str(), legendreRatioCurve, 0., 1., 2);
        f->SetParameters(L.a2, L.a4);
        f->SetLineColor(kFitColor);
        f->SetLineWidth(2);
        f->SetNpx(300);
        f->Draw("SAME");

        TGraphErrors* g = new TGraphErrors(
            n, L.cos_theta.data(), y.data(), ex.data(), ey.data());
        styleGraph(g, kAnisoColor, 20, 1.1, 2);
        g->Draw("P");

        const double xT = kLeft + 0.035;
        double yT = 1.0 - pUp->GetTopMargin() - 0.02;

        // leyenda solo en el primer panel, una fila, dos columnas
        if(withLeg){
            TLegend* lg = makeLegend(xT - 0.01, yT - lineH, 0.78, yT, tsz);
            lg->SetNColumns(2);
            lg->SetMargin(0.30);
            lg->AddEntry(g, "Data",         "pe");
            lg->AddEntry(f, "Legendre fit", "l");
            lg->Draw();
            yT -= lineH + 0.006;
        }

        // parametros del ajuste, debajo
        TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
        lat.SetTextColor(kInk); lat.SetTextAlign(13);
        lat.SetTextSize(tsz);
        lat.DrawLatex(xT, yT, Form("a_{2} = %.3f #pm %.3f", L.a2, L.u_a2));             yT -= lineH;
        if(showA4){
            lat.DrawLatex(xT, yT, Form("a_{4} = %.3f #pm %.3f", L.a4, L.u_a4));         yT -= lineH;
        }
        lat.DrawLatex(xT, yT, Form("W(0^{#circ})/W(90^{#circ}) = %.3f #pm %.3f",
                                   L.anisotropy, L.u_anisotropy));                      yT -= lineH;
        lat.DrawLatex(xT, yT, Form("#chi^{2}/ndf = %.2f", L.chi2ndf));

        padHeader(energyHeader(energy_bins[e], energy_bins[e+1]), "", 0.048 * sUp);
        pUp->RedrawAxis();

        if(!showResiduals) continue;

        // ===================== panel de residuos ============================
        pLo->cd();
        pLo->SetLeftMargin(kLeft);
        pLo->SetRightMargin(kRight);
        pLo->SetTopMargin(0.035);
        pLo->SetBottomMargin(kBotCell * sLo);
        pLo->SetTickx(1);
        pLo->SetTicky(1);
        pLo->SetFillColor(kWhite);

        const double pr = std::max(2.8, 1.25 * pmax);
        TH1F* fr2 = pLo->DrawFrame(0.0, -pr, 1.0, pr);
        styleFrame(fr2, "cos#theta_{beam}", "Pull", 0.92 * sLo, 1.30);
        fr2->GetYaxis()->SetNdivisions(503);
        fr2->GetYaxis()->CenterTitle();
        drawPullGuides(0.0, 1.0);

        if(!px.empty()){
            TGraph* gp = new TGraph((int)px.size(), px.data(), py.data());
            styleGraph(gp, kResidColor, 20, 0.75);
            gp->Draw("P");
        }
        pLo->RedrawAxis();
    }
    c->SaveAs(outname.c_str());
}