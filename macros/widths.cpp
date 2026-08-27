#include "/Users/nico/Desktop/Tese/Analysis/cross_section/include/utils.h"

static int okabeIto(double r, double g, double b)
{
    return TColor::GetColor((Float_t)(r/255.), (Float_t)(g/255.), (Float_t)(b/255.));
}
static std::vector<int> okabeItoPalette()
{
    return {
        okabeIto(  0, 114, 178),  // navy blue      (primary accent 1)
        okabeIto(230, 159,   0),  // burnt orange   (primary accent 2)
        okabeIto(  0, 158, 115),  // bluish green
        okabeIto(204, 121, 167),  // reddish purple
        okabeIto(  0,   0,   0),  // black
        okabeIto(240, 228,  66)   // yellow (use sparingly — low contrast on white)
    };
}
static const int kAnisoColor  = okabeIto(  0, 114, 178);
const int nbins = 20;
const double emin = 1.;
const double emax = 2000.;
void widths(){
    TFile *fin = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root");
    TTree *tin = (TTree*) fin->Get("events_uranium");

    double tof1, tof0, neutron_energy;
    float amp0, amp1;
    double cos_theta, cos_theta_det;

    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tin->SetBranchAddress("tof1", &tof1);
    tin->SetBranchAddress("tof0", &tof0);
    tin->SetBranchAddress("amp1", &amp1);
    tin->SetBranchAddress("amp0", &amp0);
    tin->SetBranchAddress("neutron_energy", &neutron_energy);

    Long64_t nentries = tin->GetEntries();

    std::vector<double> neutron_energy_bins = buildLogBins(nbins, emin, emax);
    std::vector<double> ecenters(nbins);
    for (int i=0;i<nbins;i++) ecenters[i] = std::sqrt(neutron_energy_bins[i]*neutron_energy_bins[i+1]);
    std::vector<TH1D*> hists_tof(nbins);
    int nbins_hist = 100;
    double min_tof = -20;
    double max_tof = +20;
    for (int i = 0; i<nbins; i++){
        hists_tof[i] = new TH1D(Form("hists_ebin%d", i), "", nbins_hist, min_tof, max_tof);
        hists_tof[i]->SetDirectory(nullptr);
    }

    std::vector<double> fwtm_vec(nbins);
    for (int i = 0 ; i<nentries ; i++){
        tin->GetEntry(i);
        double rat_amp = (amp1-amp0)/(amp0+amp1);
        double sum_amps = amp1+amp0;
        if (sum_amps<9000 || rat_amp>0.3) continue;
        int bin = findBin(neutron_energy_bins, neutron_energy);
        if (bin < 0 || bin >= nbins) continue;
        hists_tof[bin]->Fill(tof1-tof0);
    }
    fin->Close();
    for (int ebin = 0; ebin<nbins; ebin++){
        double plateau = 0;
        double left_edge = 0;
        double right_edge = 0;
        int peak_bin = -1;
        int plateau_left = 0;
        int plateau_right = 0;
        double max_peak = -100;
        for (int i = 1; i<=nbins_hist; i++){
            if (hists_tof[ebin]->GetBinContent(i)>max_peak){
                peak_bin = i;
                max_peak = hists_tof[ebin]->GetBinContent(i);
            }
        }
        int bin_right = hists_tof[ebin]->FindBin(max_tof);
        int bin_left = hists_tof[ebin]->FindBin(min_tof);
        for (int i = 1 ; i<bin_right; i++)
        for (int i = 1; i<peak_bin; i++){
            if(hists_tof[ebin]->GetBinContent(i)>max_peak/2){
                left_edge = hists_tof[ebin]->GetBinCenter(i);
                break;
            }
        }
        for (int i = nbins_hist; i>peak_bin; i--){
            if(hists_tof[ebin]->GetBinContent(i)>max_peak/2){
                right_edge = hists_tof[ebin]->GetBinCenter(i);
                break;
            }
        }
        double fwtm = right_edge-left_edge;
        fwtm_vec[ebin] = fwtm;
    }
    TCanvas *c = new TCanvas("c","",800,600);
    c->SetLogx();
    TGraphErrors *gr = new TGraphErrors(nbins, ecenters.data(), fwtm_vec.data(), nullptr, nullptr);
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kAnisoColor);
    gr->SetMarkerSize(1.);
    gr->SetTitle("; E_{n} [MeV]; FWTM (ns)");
    gr->GetXaxis()->SetTitleOffset(1.1);
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->Draw("AP");

    TFile *fout = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/output/widths.root", "RECREATE");
    gr->Write();
    fout->Close();


}