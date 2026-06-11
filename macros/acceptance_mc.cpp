// code to compute the acceptance
// toy simulation in which we emit isotropically from a target
//  then check if the fragments arrive at the ppacs active area.
// according to D. Tarrío, the acceptance should be the same for tilted and perpendicular setups
// to do the tilted, the easiest thing we can do is work on the detector frame, where target and ppacs are not rotated
// then, after checking whether it hit or not, convert to beam axis.


void acceptance_mc(){

    const int nevents = 1e7;
    const double target_radius = 4; // cm
    const double target_thickness = 0.4e-4; //cm
    const double ppac_half = 10; 
    int missed = 0;
    int coincidence = 0;
    int counts_forward = 0;
    int counts_backward = 0;
    double z_front = 2.5+0.32;
    double z_back = -2.5-0.32;
    const int nbins_beam = 40;
    const int nbins_det = 40;
    const double dcos_beam = double(1./nbins_beam);
    const double dcos_det = double(1./nbins_det);

    std::vector<std::vector<double>>cell_counts(nbins_beam, std::vector<double>(nbins_det,0.0));


    TRandom3 rng(0);

    TH2D* hist_thetas =  new TH2D ("theta_det_beam", "", nbins_beam,0,1,nbins_det,0,1);
    TH2D* eff_beam =  new TH2D ("eff_beam", "",100,0,1, 100,-TMath::Pi(),TMath::Pi());
    TH2D* eff_det =  new TH2D ("eff_det", "",100,0,1, 100,-180,180);
    TH2D* hist_emitted =  new TH2D ("emitted", "", nbins_beam,0,1,nbins_det,0,1);
    TH1D* hist_cos_theta_det = new TH1D("cos_theta_det", "", nbins_det, 0,1);
    TH1D* hist_cos_theta = new TH1D("cos_theta", "", nbins_det, 0,1);
    TH1D* hist_x0 = new TH1D("x0", "", 100, -15,15);
    TH1D* hist_cos_theta_emitted = new TH1D("cos_theta_emitted", "", nbins_det, 0,1);

    std::vector<TH1D*> solid_angle_curves(nbins_det, nullptr);

    for (int i = 0 ; i < nbins_det; i++){
        TH1D* solid_angle = new TH1D (Form("solid_angle_cell_%d", i), "", nbins_beam, 0, 1);
        solid_angle_curves[i] = solid_angle;
    }

    for (int i = 0 ; i < nevents; i++){
        double origin_r =std::sqrt(rng.Uniform(0,1))*target_radius;
        double phi_origin = 2*TMath::Pi() * rng.Uniform(0,1);
        double x0 = origin_r * cos(phi_origin);
        double y0 = origin_r * sin(phi_origin);

        double origin_z = (0.5-rng.Uniform(0,1))*target_thickness;

        std::vector<double> origin(3, 0.0);

        origin = {x0, y0, origin_z};

        // once we have the point in the target->emit isotropically
        double cos_theta_dir = rng.Uniform(0,1);
        double sin_theta_dir = sqrt(1-cos_theta_dir*cos_theta_dir);
        double phi_dir = TMath::Pi()*2*rng.Uniform(0,1);
        std::vector<double> dir_forward(3,0.0);
        std::vector<double> dir_backward(3,0.0);
        dir_forward = {std::cos(phi_dir)*sin_theta_dir, std::sin(phi_dir)*sin_theta_dir, cos_theta_dir}; //forward emitted
        dir_backward =  {-std::cos(phi_dir)*sin_theta_dir, -std::sin(phi_dir)*sin_theta_dir, -cos_theta_dir}; // assuming back-to-back

        double t_front = (z_front-origin[2])/dir_forward[2]; 
        double t_back = (z_back-origin[2])/dir_backward[2];

        if (t_front<0 || t_back<0) continue;

        double x_intersect_forward = origin[0]+t_front * dir_forward[0];
        double y_intersect_forward = origin[1]+t_front * dir_forward[1];
        double x_intersect_backward = origin[0]+t_back * dir_backward[0];
        double y_intersect_backward = origin[1]+t_back * dir_backward[1];

        bool hit_front = x_intersect_forward>-ppac_half-2.5 && x_intersect_forward<ppac_half-2.5 && std::fabs(y_intersect_forward)<=ppac_half;
        bool hit_back = x_intersect_backward>-ppac_half+2.5 && x_intersect_backward<ppac_half+2.5 && std::fabs(y_intersect_backward)<=ppac_half;

        double dz = z_front - z_back;
        double dx = x_intersect_forward-x_intersect_backward;
        double dy = y_intersect_forward-y_intersect_backward;
        double norm = sqrt(dx*dx+dy*dy+dz*dz);

        double cos_theta_det = dz/norm;
        double sin_theta_det = std::sqrt(1-cos_theta_det*cos_theta_det);
        double phi_det = TMath::ATan2(dy,dx);

        double cos_theta = 1/sqrt(2) * (-sin_theta_det*std::cos(phi_det)+cos_theta_det);
        double phi = TMath::ATan2(sqrt(2)*sin_theta_det*std::sin(phi_det), sin_theta_det*std::cos(phi_det)+cos_theta_det);
        hist_emitted->Fill(std::fabs(cos_theta), cos_theta_det);
        hist_cos_theta_emitted->Fill(std::fabs(cos_theta));

        if (hit_front){
            counts_forward++;
        }
        if (hit_back){
            counts_backward++;
        }
        if (!hit_back || !hit_front){
            missed++;
        }
        if (hit_front && hit_back){
            coincidence++;
            int bin_beam = int(std::fabs(cos_theta)/dcos_beam);
            int bin_det = int(std::fabs(cos_theta_det)/dcos_det);
            hist_x0->Fill(x_intersect_backward);

            hist_cos_theta_det->Fill(cos_theta_det);
            hist_cos_theta->Fill(abs(cos_theta));

            cell_counts[bin_beam][bin_det]++;
            solid_angle_curves[bin_det]->Fill(cos_theta);

            hist_thetas->Fill(std::fabs(cos_theta), cos_theta_det);
            eff_beam->Fill(cos_theta, phi);
            eff_det->Fill(cos_theta_det, phi_det*180/TMath::Pi());
        }
    }

    std::cout << "----------RESULTS-------------" << std::endl;
    std::cout << " % hits front = " << 1.0*counts_forward/nevents*100. << std::endl;
    std::cout << " % hits back = " << 1.0*counts_backward/nevents*100. << std::endl;
    std::cout << " % hits in both = " << 1.0*coincidence/nevents*100. << std::endl;
    std::cout << " % missed = " << 1.0*missed/nevents*100. << std::endl;

    TH1D *efficiency = new TH1D ("eff_solid_angle", "", nbins_det, 0, 1);
    TH1D *acceptance = new TH1D ("acceptance", "", nbins_det, 0, 1);

    TFile *fout  = new TFile ("mc_acceptance.root", "RECREATE");
    eff_det->Write();
    eff_beam->Write();
    hist_thetas->Divide(hist_emitted);
    hist_emitted->Write();
    hist_thetas->Write();

    for (int i = 0 ; i < nbins_det; i++){
        double acc = 0;
        for (int j = 0; j < nbins_beam; j++){
            acc += hist_emitted->GetBinContent(j+1, i+1);  // cuentas emitidas antes del divide
        }
        double hits = 0;
        for (int j = 0; j < nbins_beam; j++){
            hits += cell_counts[j][i];  // cuentas que dieron coincidencia
        }
    double eff = (acc > 0) ? hits / acc : 0;
    acceptance->SetBinContent(i+1, eff);
    acceptance->SetBinError(i+1, sqrt(eff * (1 - eff) / acc));
}
    acceptance->Scale(1/acceptance->Integral());
    acceptance->Write();
    hist_cos_theta->Divide(hist_cos_theta_emitted);
    hist_cos_theta_det->Write();
    hist_cos_theta->Write();
    hist_x0->Write();
    for (int i = 0 ; i < nbins_det; i++){
        solid_angle_curves[i]->Write();
        double eff = solid_angle_curves[i]->Integral() / coincidence;
        efficiency->SetBinContent(i+1, eff);
        efficiency->SetBinError(i+1, sqrt(eff * (1 - eff) / coincidence));
    }

    efficiency->Divide(acceptance); 
    std::cout << "Efficiency bin content (0) after divide: " << efficiency->GetBinContent(30) << std::endl;
    efficiency->Write();    
    fout->Close();


}