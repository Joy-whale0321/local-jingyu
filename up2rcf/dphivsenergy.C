#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <TChain.h>
#include <fstream>
#include <string>

void f_guassfit(TH1D *h_fit);

void dphivsenergy() {
    TFile *newf = new TFile("/home/jingyu/sphenix/fromxvdong/PhotonEMC-main/output/photon_deltaphi_energy.root","recreate");

    // TFile *file = TFile::Open("/home/jingyu/sphenix/fromxvdong/PhotonEMC-main/output/TrackCalo_0_ana_truthcluster.root");
    // TTree *tree = (TTree*)file->Get("tree");
    TChain *chain("tree");
    
    std::ifstream infile("list.txt");
    std::string filename;

    while (std::getline(infile, filename)) {
        if (!filename.empty()) {
            chain->Add(filename.c_str());
            std::cout << "Added file: " << filename << std::endl;
        }
    }
    // 检查总共添加了多少事件
    Long64_t nentries = chain->GetEntries();
    std::cout << "Total number of entries in the chain: " << nentries << std::endl;

    TH2D *h_deltaphi_energy = new TH2D("h_deltaphi_energy", "h_deltaphi_energy", 10000, -0.5, 0.5, 12, -1, 11);

    // 创建指向vector的指针，来接收分支的数据
    std::vector<double> *_emcal_phi           = nullptr;
    std::vector<double> *_truth_phi           = nullptr;
    std::vector<double> *_truth_eta           = nullptr;
    std::vector<double> *_truth_e             = nullptr;

    // 设置分支地址
    chain->SetBranchAddress("_emcal_phi", &_emcal_phi);
    chain->SetBranchAddress("_truth_phi", &_truth_phi);
    chain->SetBranchAddress("_truth_eta", &_truth_eta);
    chain->SetBranchAddress("_truth_e", &_truth_e);

    // 获取条目数量
    Long64_t nentries = chain->GetEntries();

    // 按顺序读取每个条目中的vector
    for (Long64_t i = 0; i < nentries; i++) 
    {
        chain->GetEntry(i);

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

    int nYBins = h_deltaphi_energy->GetNbinsY();
    for (int i = 1; i <= nYBins; ++i) {
        TH1D *h_delta_phi_13 = h_deltaphi_energy->ProjectionX(Form("h_delta_phi_bin_%d", i), i, i);
        f_guassfit(h_delta_phi_13);

        TString histName = Form("h_delta_phi_bin_%d", i);
        newf->WriteTObject(h_delta_phi_13, histName);
    }

    newf->WriteTObject(h_deltaphi_energy, "h_deltaphi_energy");

    file->Close();
}

void f_guassfit(TH1D *h_fit)
{
    gStyle->SetOptStat(1111);

    TF1 *gausFit = new TF1("gausFit", "gaus", -0.05, 0.05);

    h_fit->Fit(gausFit, "R"); // "R" 表示只在指定的范围内进行拟合

    double peak = gausFit->GetParameter(0);
    double mean = gausFit->GetParameter(1);
    double sigma = gausFit->GetParameter(2);

    cout<<"sigma is: "<<sigma<<endl;
}