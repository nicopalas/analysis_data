// what we want is to run this macro with the raw data and to see if there is any 

static const double NEUTRON_MEV = 939.56542;
static const double C_NS        = 0.299792458;   // m/ns
static const double AMP_CUT = 1000;

double tof_to_energy(double tof_ns, double t_offset_ns, double flight_m)
{
    double t_light_ns = (flight_m / C_NS);
    double t_s = (tof_ns - t_offset_ns + t_light_ns) * 1e-9;
    if (t_s <= 0) return -1.0;
    double c = 299792458.0;
    double b2 = (flight_m / (c * t_s));
    b2 *= b2;
    if (b2 >= 1.0) return -1.0;
    return NEUTRON_MEV * (1.0 / std::sqrt(1.0 - b2) - 1.0);
}
double flight_path(int detn){
    return 183.5+sqrt(2)*5.0*detn;
}
void run_stability(int run_number){

    const int ndets_anode  = 10; 
    const int ndets_cathode = 40; 
    TString *fname = "/nucl_lustre/n_tof_INTC_P_665/Analysis/time_calibrations/time_calibrated_run"+std::to_str(run_number)+".root";
    TFile *fin = TFile::Open(fname.to_str(), "RECREATE");

    TTree *ppan = (TTree*) fin->Get("calibrated_ppan");

    Long64_t nentries_ppan = ppan->GetEntries();

    std::cout << "Number of entries in PPAN tree : " << nentries_ppan << ". " << std::endl; 
    
    std::vector<TH1D> hists_ppan;
    for (int i = 0 ; i<ndets_anode; i++){
        TH1D *h_ppan = new TH1D(Form("hist_det_%d", i), "", 400, 0, 2000);
        hists_ppan.push_back(hists_ppan);
    }

    double tof, amp, gamma_flash;
    float PulseIntensity;
    int BunchNumber, RunNumber, detn;

    ppan->SetBranchAddress("tof", &tof);
    ppan->SetBranchAddress("amp", &amp);
    ppan->SetBranchAddress("gamma_flash", &gamma_flash);
    ppan->SetBranchAddress("PulseIntensity", &PulseIntensity);
    ppan->SetBranchAddress("BunchNumber", &BunchNumber);
    ppan->SetBranchAddress("RunNumber", &RunNumber);
    ppan->SetBranchAddress("detn", &detn);

    ppan->GetEntry(nentries_ppan-1);
    int bunches = BunchNumber+1;

    std::vector<double> counts_bunch;
    
    for (int i = 0; i < nentries_ppan; i++){
        ppan->GetEntry(i);
        if (amp<AMP_CUT) continue;
        double flight_path_m = flight_path(detn);
        double neutron_energy = tof_to_energy(tof, gamma_flash, flight_path_m);
        hists_ppan[detn]->Fill(neutron_energy);
    }
}