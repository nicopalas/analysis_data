void efficiency_toy() {

    std::vector<double> E    = {46.3, 66.6, 73.9, 94.1, 132.9, 144.6, 173.3};
    std::vector<double> sig  = {0.103, 0.81, 1.20, 2.81, 6.1,   8.1,   10.3};
    std::vector<double> usig = {0.019, 0.12, 0.17, 0.39, 0.9,   1.2,   1.6};
    std::vector<double> ex(E.size(), 0.0);

    // convert mb to barn
    for (auto &s : sig)  s  *= 1e-3;
    for (auto &u : usig) u  *= 1e-3;

    TGraphErrors *gr = new TGraphErrors(
        (int)E.size(),
        E.data(), sig.data(),
        ex.data(), usig.data());
    gr->SetName("g_gold_ref");
    gr->SetTitle("Au-197 reference;E_{n} (MeV);#sigma (barn)");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kRed+1);
    gr->SetLineColor(kRed+1);

    TFile *fout = TFile::Open(
        "/Users/nico/Desktop/Tese/Analysis/cross_section/output/cs_toy_gold.root",
        "UPDATE");
    gr->Write();
    fout->Close();
    std::cout << "[DONE] gold_reference.root\n";
}