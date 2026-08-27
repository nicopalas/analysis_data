void intrinsic_efficiency(){

    TFile *mc = TFile::Open("/Users/nico/Desktop/Tese/Analysis/mc_acceptance.root");
    TH1D *h_mc = (TH1D*) mc->Get("cos_theta_det");
    int nbins = h_mc->GetNbinsX();

    TFile *eff = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/output/U-238/output_efficiency_uranium.root");
    TH1D *h_eff_ebin_0 = (TH1D*) eff->Get("heff_ebin0");
    TH1D *h_eff_ebin_1 = (TH1D*) eff->Get("heff_ebin1");
    TH1D *h_eff_ebin_2 = (TH1D*) eff->Get("heff_ebin2");
    TH1D *h_eff_ebin_3 = (TH1D*) eff->Get("heff_ebin3");

    if (!h_mc || !h_eff_ebin_0 || !h_eff_ebin_1) {
        std::cerr << "[ERROR] falta algún histograma de entrada (revisa nombres)\n";
        return;
    }
    if (h_eff_ebin_0->GetNbinsX() != nbins || h_eff_ebin_1->GetNbinsX() != nbins) {
        std::cerr << "[ERROR] binning incompatible entre eff y mc\n";
        return;
    }

    TH1D *h_intrinsic_ebin0 = (TH1D*) h_eff_ebin_0->Clone("h_intrinsic_uranium_ebin0");
    h_intrinsic_ebin0->Divide(h_mc);

    TH1D *h_intrinsic_ebin1 = (TH1D*) h_eff_ebin_1->Clone("h_intrinsic_uranium_ebin1");
    h_intrinsic_ebin1->Divide(h_mc);

    TH1D *h_intrinsic_ebin2 = (TH1D*) h_eff_ebin_2->Clone("h_intrinsic_uranium_ebin2");
    h_intrinsic_ebin2->Divide(h_mc);

    TH1D *h_intrinsic_ebin3 = (TH1D*) h_eff_ebin_3->Clone("h_intrinsic_uranium_ebin3");
    h_intrinsic_ebin3->Divide(h_mc);

    TFile *fout = new TFile("intrinsic_uranium.root", "RECREATE");
    h_intrinsic_ebin0->Write();
    h_intrinsic_ebin1->Write();
    h_intrinsic_ebin2->Write();
    h_intrinsic_ebin3->Write();
    fout->Close();
}