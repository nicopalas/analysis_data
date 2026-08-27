#include "../include/config.h"
#include "../include/cuts.h"
#include "../include/histograms.h"


void run_number_analysis(){
    TFile *fin = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/events_selection.root");
    TFile *fold = TFile::Open("/Users/nico/Desktop/Tese/Analysis/cross_section/data/coincidences.root")
    TTree *tin = (TTree*) fin->Get("events_uranium");
    TTree *told = (TTree*) fold->Get("events_uranium");

    double cos_theta, cos_theta_det;
    double tof1, tof0;
    double amp0, amp1;
    int RunNumber;

    tin->SetBranchAddress("amp0", &amp0);
    tin->SetBranchAddress("amp1", &amp1);
    tin->SetBranchAddress("tof0", &tof0);
    tin->SetBranchAddress("tof1", &tof1);
    tin->SetBranchAddress("cos_theta", &cos_theta);
    tin->SetBranchAddress("cos_theta_det", &cos_theta_det);
    tin->SetBranchAddress("RunNumber", &RunNumber);

    double cos_theta_old, cos_theta_det_old;
    double tof1_old, tof0_old;
    double amp0_old, amp1_old;
    int RunNumber_old;

    //what we want to see is that the calibration was performed correctly:
    // in order to do this, we can see the amplitude of the FFs along the different runs and see if their distributions are constant
    // we can compare with the non aligned case to have a reference






    
}