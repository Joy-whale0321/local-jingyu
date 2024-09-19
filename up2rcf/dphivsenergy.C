#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <TChain.h>
#include <fstream>
#include <string>
#include <TGraphErrors.h>

void f_guassfit(TH1D *h_fit, double &peak, double &mean, double &sigma, double fvalue);

void dphivsenergy() {
    TFile *newf = new TFile("/home/jingyu/sphenix/fromxvdong/PhotonEMC-main/output/photon_dphidz_energy_10K.root","recreate");

    TFile *file = TFile::Open("/home/jingyu/sphenix/fromxvdong/PhotonEMC-main/output/TrackCalo_0_ana.root");
    TTree *tree = (TTree*)file->Get("tree");
    // TChain *chain = new TChain("tree");
    
    // std::ifstream infile("list.txt");
    // std::string filename;

    // while (std::getline(infile, filename)) {
    //     if (!filename.empty()) {
    //         chain->Add(filename.c_str());
    //         std::cout << "Added file: " << filename << std::endl;
    //     }
    // }

    TH2D *h_deltaphi_energy = new TH2D("h_deltaphi_energy", "h_deltaphi_energy", 10000, -0.5, 0.5, 12, -1, 11);
    TH2D *h_deltaZ_energy  = new TH2D("h_deltaZ_energy",  "h_deltaZ_energy", 10000, -10, 10, 12, -1, 11);
    TH1D *h_dfactor_energy = new TH1D("h_dfactor_energy", "h_dfactor_energy", 10001, -0.005, 0.005);

    TH1D *h_deltaZ_1D = new TH1D("h_deltaZ_1D", "h_deltaZ_1D", 10000, 50, 50);

    // 创建指向vector的指针，来接收分支的数据
    std::vector<double> *_emcal_phi           = nullptr;
    std::vector<double> *_emcal_z             = nullptr;
    std::vector<double> *_truth_phi           = nullptr;
    std::vector<double> *_truth_e             = nullptr;
    std::vector<double> *_truth_px            = nullptr;
    std::vector<double> *_truth_py            = nullptr;
    std::vector<double> *_truth_pz            = nullptr;
    std::vector<double> *_truth_pt            = nullptr;

    // 设置分支地址
    tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
    tree->SetBranchAddress("_emcal_z", &_emcal_z);
    tree->SetBranchAddress("_truth_phi", &_truth_phi);
    tree->SetBranchAddress("_truth_e", &_truth_e);
    tree->SetBranchAddress("_truth_px", &_truth_px);
    tree->SetBranchAddress("_truth_py", &_truth_py);
    tree->SetBranchAddress("_truth_pz", &_truth_pz);
    tree->SetBranchAddress("_truth_pt", &_truth_pt);

    // 获取条目数量
    Long64_t nentries = tree->GetEntries();

    // 按顺序读取每个条目中的vector
    for (Long64_t i = 0; i < nentries; i++) 
    {
        tree->GetEntry(i);

        // map 1 - truth; 2 - g4hit; 3 - cluster reco; 4 - g4particle;
        if (_emcal_phi->size() == 1) 
        {
            double energy_truth = _truth_e->at(0);
            double phi_truth = _truth_phi->at(0);

            double px_truth = _truth_px->at(0);
            double py_truth = _truth_py->at(0);
            double pz_truth = _truth_pz->at(0);
            double pt_truth = _truth_pt->at(0);

            double phi_emcal_cluster = _emcal_phi->at(0);
            double Z_emcal_cluster   = _emcal_z->at(0);

            // get delta phi
            double delta_phi_13 = phi_truth - phi_emcal_cluster;            
            h_deltaphi_energy->Fill(delta_phi_13, energy_truth);

            // get delta Z
            double factor1 = pz_truth * sqrt( 1 / (energy_truth*energy_truth - pz_truth*pz_truth));
            double factor2 = pz_truth / pt_truth;           
            double delta_factor = factor1 - factor2;

            h_dfactor_energy->Fill(delta_factor);

            double R = 93.5;
            double delta_Z = factor1*R - Z_emcal_cluster;
            h_deltaZ_energy->Fill(delta_Z,energy_truth);

            h_deltaZ_1D->Fill(delta_Z);
        }
    }

    std::vector<double> *x_energy  = new std::vector<double>;
    std::vector<double> *y_dphi    = new std::vector<double>;
    std::vector<double> *y_dz    = new std::vector<double>;
    std::vector<double> *ex_energy = new std::vector<double>;
    std::vector<double> *ey_dphi   = new std::vector<double>;
    std::vector<double> *ey_dz   = new std::vector<double>;

    int nYBins = h_deltaphi_energy->GetNbinsY();
    for (int i = 2; i <= (nYBins-1); ++i) {

        double peak =0;
        double mean =0;
        double sigma=0;

        TH1D *h_delta_phi_13 = h_deltaphi_energy->ProjectionX(Form("h_delta_phi_bin_%d", i), i, i);
        double fitrange = 0.01;
        f_guassfit(h_delta_phi_13, peak, mean, sigma, fitrange);

        TString histName = Form("h_delta_phi_bin_%d", i);
        newf->WriteTObject(h_delta_phi_13, histName);
    
        x_energy ->push_back(i-1.5);
        y_dphi   ->push_back(mean);
        ex_energy->push_back(0.5);
        ey_dphi  ->push_back(sigma);
    }
    
    for (int i = 2; i <= (nYBins-1); ++i) {

        double peak =0;
        double mean =0;
        double sigma=0;

        TH1D *h_delta_Z = h_deltaZ_energy->ProjectionX(Form("h_delta_Z_bin_%d", i), i, i);
        double fitrange = 1;
        f_guassfit(h_delta_Z, peak, mean, sigma, fitrange);

        TString histName = Form("h_delta_Z_bin_%d", i);
        newf->WriteTObject(h_delta_Z, histName);
    
        y_dz ->push_back(mean);
        ey_dz->push_back(sigma);
    }

    TGraphErrors *g_deltaphi_energy= new TGraphErrors(x_energy->size(), x_energy->data(), y_dphi->data(), ex_energy->data(), ey_dphi->data());
    g_deltaphi_energy->SetTitle("Delta Phi vs Energy;Energy (GeV);#Phi_{truth - cluster}");
    g_deltaphi_energy->SetMarkerStyle(20);
    g_deltaphi_energy->SetMarkerSize(1.0);
    g_deltaphi_energy->SetMarkerColor(kBlack); 

    TGraphErrors *g_deltaZ_energy= new TGraphErrors(x_energy->size(), x_energy->data(), y_dz->data(), ex_energy->data(), ey_dz->data());
    g_deltaZ_energy->SetTitle("Delta Z vs Energy; Energy (GeV); Z_{truth - cluster} (cm)");
    g_deltaZ_energy->SetMarkerStyle(20);
    g_deltaZ_energy->SetMarkerSize(1.0);
    g_deltaZ_energy->SetMarkerColor(kBlack); 

    newf->WriteTObject(g_deltaphi_energy, "g_deltaphi_energy");
    newf->WriteTObject(h_deltaphi_energy, "h_deltaphi_energy");
    newf->WriteTObject(g_deltaZ_energy, "g_deltaZ_energy");
    newf->WriteTObject(h_deltaZ_energy, "h_deltaZ_energy");
    newf->WriteTObject(h_dfactor_energy,"h_dfactor_energy");

    newf->WriteTObject(h_deltaZ_1D,"h_deltaZ_1D");

    file->Close();
}

void f_guassfit(TH1D *h_fit, double &peak, double &mean, double &sigma, double fvalue)
{
    gStyle->SetOptStat(1111);

    double xmin = -1*fvalue;
    double xmax = fvalue;
    TF1 *gausFit = new TF1("gausFit", "gaus", xmin, xmax);

    h_fit->Fit(gausFit, "R"); // "R" 表示只在指定的范围内进行拟合

    peak = gausFit->GetParameter(0);
    mean = gausFit->GetParameter(1);
    sigma= gausFit->GetParameter(2);

    cout<<"sigma is: "<<sigma<<endl;
}