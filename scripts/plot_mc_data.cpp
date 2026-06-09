void plot_mc_data(){
    TFile* f0 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root");
    TFile* f1 = TFile::Open("/Users/nico/Desktop/Tese/Analysis/montecarlo/output/output_mc.root");

    TTree* t_data = (TTree*)f0->Get("events_gold");
    TTree* t_mc   = (TTree*)f1->Get("CoincTree");

    t_data->Draw("tof1-tof0>>h_data(100,-20,20)",
                 "neutron_energy>100 && neutron_energy<800  ", "goff");
    t_mc->Draw("time_j-time_i>>h_mc(100,-20,20)",
               "ppac_i==7 && ppac_j==8 && neutronE>100 && neutronE<800", "goff");

    TH1D* h_data = (TH1D*)gDirectory->Get("h_data");
    TH1D* h_mc   = (TH1D*)gDirectory->Get("h_mc");

    h_data->Scale(1./h_data->Integral());
    h_mc->Scale(1./h_mc->Integral());

    h_data->SetLineColor(kBlack);  h_data->SetLineWidth(2);
    h_mc->SetLineColor(kRed);      h_mc->SetLineWidth(2);

    h_data->GetXaxis()->SetTitle("#Deltat (ns)");
    h_data->GetYaxis()->SetTitle("Normalized counts");
    h_data->SetTitle("Data vs MC: #Deltat between adjacent PPACs");

    TCanvas* c = new TCanvas("c","",800,600);
    h_data->Draw("HIST");
    h_mc->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.7, 0.88, 0.88);
    leg->AddEntry(h_data, "Data (U-238)", "l");
    leg->AddEntry(h_mc,   "MC",        "l");
    leg->Draw();

    c->Update();
    c->WaitPrimitive();
}