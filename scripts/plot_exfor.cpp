void plot_exfor(const char* csvfile = "/Users/nico/Downloads/417560032.csv"){
    std::ifstream f(csvfile);
    std::string line;
    std::getline(f, line);

    std::vector<double> x, y, ey;
    while (std::getline(f,line)){
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')){
            cols.push_back(token);
        }
        if (cols.size()<11) continue;
        double cs = std::stod(cols[3])*1e3;
        double en = std::stod(cols[11])/1e6;
        double scs = cs*std::stod(cols[4])/100;
        x.push_back(en);
        y.push_back(cs);
        ey.push_back(scs);
    }
    TGraph *g = new TGraphErrors(x.size(), x.data(), y.data(), nullptr ,ey.data());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlue+1);
    g->SetLineColor(kBlue+1);
    g->SetTitle("EXFOR 41756;E_{n} (MeV);#sigma (mb)");
    TCanvas* c = new TCanvas("c","EXFOR",800,600);
    c->SetLogx();
    c->SetLogy();
    g->Draw("AP");
}
    