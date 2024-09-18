#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <TChain.h>
#include <fstream>
#include <string>
#include <TGraphErrors.h>

void f_guassfit(TH1D *h_fit, double &peak, double &mean, double &sigma);

void dphivsenergy() {
    TFile *newf = new TFile("/sphenix/user/jzhang1/testcode4all/PhotonEMC/macro/photon_deltaphi_energy_1M.root","recreate");

    TFile *file = TFile::Open("/sphenix/user/jzhang1/testcode4all/PhotonEMC/macro/Reconstructed/17/truthcluster_1M.root");
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

    // 创建指向vector的指针，来接收分支的数据
    std::vector<double> *_emcal_phi           = nullptr;
    std::vector<double> *_truth_phi           = nullptr;
    std::vector<double> *_truth_eta           = nullptr;
    std::vector<double> *_truth_e             = nullptr;

    // 设置分支地址
    tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
    tree->SetBranchAddress("_truth_phi", &_truth_phi);
    tree->SetBranchAddress("_truth_eta", &_truth_eta);
    tree->SetBranchAddress("_truth_e", &_truth_e);

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

            double phi_emcal_cluster = _emcal_phi->at(0);

            double delta_phi_13 = phi_truth - phi_emcal_cluster;
            
            h_deltaphi_energy->Fill(delta_phi_13, energy_truth);
        }
    }

    std::vector<double> *x_energy  = nullptr;
    std::vector<double> *y_dphi    = nullptr;
    std::vector<double> *ex_energy = nullptr;
    std::vector<double> *ey_dphi   = nullptr;

    int nYBins = h_deltaphi_energy->GetNbinsY();
    for (int i = 2; i <= (nYBins-1); ++i) {

        double peak =0;
        double mean =0;
        double sigma=0;

        TH1D *h_delta_phi_13 = h_deltaphi_energy->ProjectionX(Form("h_delta_phi_bin_%d", i), i, i);
        f_guassfit(h_delta_phi_13, peak, mean, sigma);

        TString histName = Form("h_delta_phi_bin_%d", i);
        newf->WriteTObject(h_delta_phi_13, histName);
    
        x_energy ->push_back((i-1.5));
        y_dphi   ->push_back(mean);
        ex_energy->push_back(0.5);
        ey_dphi  ->push_back(sigma);
    }

    TGraphErrors *g_deltaphi_energy= new TGraphErrors(x_energy->size(), x_energy->data(), y_dphi->data(), ex_energy->data(), ey_dphi->data());

    newf->WriteTObject(h_deltaphi_energy, "h_deltaphi_energy");
    newf->WriteTObject(g_deltaphi_energy, "g_deltaphi_energy");

    // file->Close();
}

void f_guassfit(TH1D *h_fit, double &peak, double &mean, double &sigma)
{
    gStyle->SetOptStat(1111);

    TF1 *gausFit = new TF1("gausFit", "gaus", -0.01, 0.01);

    h_fit->Fit(gausFit, "R"); // "R" 表示只在指定的范围内进行拟合

    peak = gausFit->GetParameter(0);
    mean = gausFit->GetParameter(1);
    sigma= gausFit->GetParameter(2);

    cout<<"sigma is: "<<sigma<<endl;
}