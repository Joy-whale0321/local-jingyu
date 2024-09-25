#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include "stdlib.h"
#include <iostream>
#include "Riostream.h"
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <complex>
#include <vector>
#include <algorithm>

#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TGraph.h"
#endif

//#include "test_ampt.h"
#include "AMPT.h"

using namespace std;

#define PI 3.1415926

//ifstream inputFile;

vector<TVector3> Dm;
vector<TVector3> Hm;
vector<TVector3> uds_quark;
vector<TVector3> c_quark;

double eta_cut=1.0;//MLCC
double delta_eta_cut=0.5;//MLCC                              
double pt_low_cut=0.3;//MLCC
double pt_high_cut=5.0;//MLCC
//double npair;

int nmult1,nmult2;
double sum_r2sin1phi, sum_r2cos1phi, sum_r2sin2phi,sum_r2cos2phi,sum_r2sin3phi,sum_r2cos3phi,sum_r2sin4phi,sum_r2cos4phi,sum_r2,sum_r3,sum_r4;
double sum_r2sin1phi_charm, sum_r2cos1phi_charm, sum_r2sin2phi_charm,sum_r2cos2phi_charm,sum_r2sin3phi_charm,sum_r2cos3phi_charm,sum_r2sin4phi_charm,sum_r2cos4phi_charm,sum_r2_charm,sum_r3_charm,sum_r4_charm;

double centrality,psi2,psi3,psi4,epsilon2,epsilon3,epsilon4;

double sum_r2sin1phi_ep, sum_r2cos1phi_ep, sum_r2sin2phi_ep,sum_r2cos2phi_ep,sum_r2sin3phi_ep,sum_r2cos3phi_ep,sum_r2sin4phi_ep,sum_r2cos4phi_ep,sum_r2_ep,sum_r3_ep,sum_r4_ep;
double psiEP2, psiEP3, psiEP4,v2EP,v3EP,v4EP,v2EP_s,v3EP_s,v4EP_s,v2EP_h,v3EP_h,v4EP_h;

double res_x[10]={0, 4.42, 6.25, 7.65, 8.83, 9.88, 10.82, 11.68, 12.50, 15.0};
double res_y[10]={0.61, 0.75, 0.82, 0.8, 0.77, 0.70, 0.59, 0.46, 0.31, 0.2};

double qx2,qy2,qx3,qy3,qx4,qy4,q2,q3,q4,q2pt,q3pt,q4pt,v22c,v32c,v42c,vn20,vn21,vn22;
int multiq, hmu, dmu;

int gen, ebin1, ebin2, ebin3, ebin4, ebin5, ebin6, ebin7, ebin8, ebin9, ebin10, ebin11, ebin12, ebin13, ebin14, ebin15, ebin16, ebin17, ebin18, ebin19, ebin20, ebin21, ebin22, ebin23, ebin24, ebin25, ebin26, ebin27, ebin28, ebin29, ebin30, ebin31, ebin32, ebin33, ebin34, ebin35, ebin36, ebin37, ebin38, ebin39, ebin40, ebin41, ebin42, ebin43, ebin44, ebin45, ebin46, ebin47, ebin48;

int numv2h[10]={0};
int numv2d[10]={0};
int numhs=0;
int numds=0;
double arrv2h[10]={0};
double arrv2d[10]={0};

int numv2_uds[10]={0};
int numv2_c[10]={0};
double arrv2_uds[10]={0};
double arrv2_c[10]={0};

int noverd = 0;
int noverh = 0;
double oneoverd = 0;
double oneoverh = 0;
int nassd0 = 0;
int nassh0 = 0;

int charm_mu;
int Dmeson_mu;

bool passEvent(AMPT* event1,AMPT* event2);//逐事件提取，这里event1对应同一个事件的initial数据，event2对应同一个事件的final数据
bool parton_loop(AMPT* track,int i,AMPT* track2);//initial数据中逐个径迹提取
bool hadron_loop(AMPT* track,int i,TGraph *tgr);//final数据中逐个径迹提取
void writeHistograms(char *outFile);//填好的histogram写入文件
void deleteHistograms();

bool event_loop(AMPT* track,int i);

void makeDmcorrelation(vector<TVector3> Dm,vector<TVector3> Hm,TH2D* dphi_Dm_pt);//计算D介子和强子的方位角关联
void makeHmcorrelation(vector<TVector3> Hm,TH2D* dphi_Hm_pt);//计算强子和强子的方位角关联

bool udsc_loop(AMPT* track,int i);

void make_c_correlation(vector<TVector3> c_quark,vector<TVector3> uds_quark,TH2D* dphi_c_pt);
void make_uds_correlation(vector<TVector3> uds_quark,TH2D* dphi_uds_pt);

bool bad_event(AMPT* track,int i);

TH1D *spect_D0_pt;

TH2D* dphi_Dm_pt;
TH2D* dphi_Hm_pt;

TH2D* dphi_c_pt;
TH2D* dphi_uds_pt;

TH2D* epsilon2_centrality;
TH2D* epsilon3_centrality;
TH2D* epsilon4_centrality;

TProfile* v2part_Dm_pt;//PP v-pt
TProfile* v3part_Dm_pt;
TProfile* v4part_Dm_pt;

TProfile* v2part_Hm_pt;//PP v-pt
TProfile* v3part_Hm_pt;
TProfile* v4part_Hm_pt;

TProfile* v2EP_Dm_pt;//EP v-pt
TProfile* v3EP_Dm_pt;
TProfile* v4EP_Dm_pt;

TGraph* event_plane_resolution; 

TH1D *hbf; // v2 without event plane resolution
TH1D *haf; // v2 after event plane resolution

TH1D *hv2pp; // v2 cal by pp distribution
TH1D *hv2pd; // D0 v2 cal by pp distribution

TH1D *hv3pp;
TH1D *hv3pd;

TProfile* TPr2; // v2 cal by pp average value

TH2D *he2q2;
TH2D *he3q3;
TH2D *he4q4;

TH2D *he2qpt;
TH2D *he3qpt;
TH2D *he4qpt;

TH2D *De2v2_pp;
TH2D *He2v2_pp;
TH2D *De3v3_pp;
TH2D *He3v3_pp;

TH2D *De2v2_2c;
TH2D *He2v2_2c;
TH2D *De3v3_2c;
TH2D *He3v3_2c;

TH1D *hHad0;
TH1D *hHadmp;
TH1D *hDad0;
TH1D *hDadmp;

TH1D *hHadm3;
TH1D *hDadm3;

TH1D *hpsi2ad;
TH1D *hpsi3ad;

TProfile* H_b_v2e2;
TProfile* H_b_v3e3;
TProfile* D_b_v2e2;
TProfile* D_b_v3e3;

TH2D* H_Hb_v2e2;
TH2D* H_Hb_v3e3;
TH2D* D_Hb_v2e2;
TH2D* D_Hb_v3e3;

TH1D *hadbf; //a~angel
TH1D *hadaf;


//nass flow
TH2D* D_flow_nass;
TH1D* D_f_n_test;

//pp
TProfile* H_b_v2_pp;
TProfile* H_b_v3_pp;
TProfile* D_b_v2_pp;
TProfile* D_b_v3_pp;

TProfile* H_b_e2_pp;
TProfile* H_b_e3_pp;
TProfile* D_b_e2_pp;
TProfile* D_b_e3_pp;

//2c
TProfile* H_b_v2_2c;
TProfile* H_b_v3_2c;
TProfile* D_b_v2_2c;
TProfile* D_b_v3_2c;

TProfile* H_b_e2_2c;
TProfile* H_b_e3_2c;
TProfile* D_b_e2_2c;
TProfile* D_b_e3_2c;

TProfile* H_b_v2_p_2c;
TProfile* H_b_v3_p_2c;
TProfile* D_b_v2_p_2c;
TProfile* D_b_v3_p_2c;

//event by event 2c
TProfile* Dv2_e_s;
TProfile* Dv3_e_s;
TProfile* Hv2_e_s;
TProfile* Hv3_e_s;

TProfile* Dv2_2c_ebe;
TProfile* Dv3_2c_ebe;
TProfile* Hv2_2c_ebe;
TProfile* Hv3_2c_ebe;

// TH1D *ndna;
// TH1D *nhna;

double VD2,VH2,ptd,pth,VD2_s,VH2_s,VD3,VH3,VD3_s,VH3_s,VD4,VH4,VD4_s,VH4_s,VD2_2c_s,VH2_2c_s,VD3_2c_s,VH3_2c_s;   
double VH2_2c_b,VH3_2c_b,VD2_2c_b,VD3_2c_b;
int NassDH, NassHH;

double az_VD2,az_VH2,az_ptd,az_pth,az_VD3,az_VH3;
int az_NassDH, az_NassHH;

TProfile* VDH2_pt;
TProfile* VHH2_pt;

TProfile* VDH3_pt;
TProfile* VHH3_pt;

TProfile* az_VDH2_pt;
TProfile* az_VHH2_pt;
TProfile* az_VDH3_pt;
TProfile* az_VHH3_pt;

TProfile* V2H_ebes;
TProfile* V3H_ebes;

TProfile* V2H_ps;
TProfile* V3H_ps;


// test something
TH1D *h_numd;
TH1D *h_numh;
TH2D *h_rati;
TH2D *h_arrh;
TH2D *h_arrd;

// bad event
TH2D *h_bad_mult;
TH2D *h_bad_e2e3;
TH2D *h_bad_dhmu;

TH2D *h_bad_aft;
TH2D *h_bad_zpc;

double r_phi_bad;
double r_bads;
double rsin_bads;
double rcos_bads;

int main(int argc, char **argv)
{
    TChain *chain1 = new TChain("AMPT");
    TChain *chain2 = new TChain("AMPT");	
    TChain *chain3 = new TChain("AMPT");

    if(argc!=5 && argc!=1) return 0;
    
    char *FileInput1=0;
    char *FileInput2=0;
    char *FileInput3=0;
    char *FileOutput=0;

    if(argc==1){
      	//FileInput  = "example.list";
		FileOutput = "example.root";
    }

    if(argc==5){
		FileInput1 = argv[1];
		FileInput2 = argv[2];
		FileInput3 = argv[3];
		FileOutput = argv[4];
    }  

    int fileNumber = 0;
    char FileList[512];

    ifstream* inputStream1 = new ifstream;
    ifstream* inputStream2 = new ifstream;
    ifstream* inputStream3 = new ifstream;

    inputStream1->open(FileInput1);
    if (!(inputStream1))
    {
		printf("can not open zpct0 list file\n");
		return 0;
    }

    inputStream2->open(FileInput2);
    if (!(inputStream2))
    {
		printf("can not open ampt list file\n");
		return 0;
    }

	inputStream3->open(FileInput3);
    if (!(inputStream3))
    {
		printf("can not open zpc list file\n");
		return 0;
    }

    for (;inputStream1->good();)
    {
		inputStream1->getline(FileList,512);
		if  ( inputStream1->good() )
		{
	    	TFile *ftmp1 = new TFile(FileList);
	    	if(!ftmp1||!(ftmp1->IsOpen())||!(ftmp1->GetNkeys()))
	    	{
				printf("input1 zpct0 file %s error in opening!!!\n",FileList);
	    	}
	    	else
	    	{
				printf("input1 zpct0 read in file %s\n",FileList);
				chain1->Add(FileList);
				fileNumber++;
	    	}
	    	delete ftmp1;
		}
    }

   	for (;inputStream2->good();)
    {
		inputStream2->getline(FileList,512);
		if  ( inputStream2->good() )
		{
		    TFile *ftmp2 = new TFile(FileList);
		    if(!ftmp2||!(ftmp2->IsOpen())||!(ftmp2->GetNkeys()))
		    {
				printf("input2 afterART file %s error in opening!!!\n",FileList);
		    }
		    else
		    {
				printf("input2 afterART read in file %s\n",FileList);
				chain2->Add(FileList);
				fileNumber++;
		    }
		    delete ftmp2;
		}
    }

	for (;inputStream3->good();)
    {
		inputStream3->getline(FileList,512);
		if  ( inputStream3->good() )
		{
	    	TFile *ftmp3 = new TFile(FileList);
	    	if(!ftmp3||!(ftmp3->IsOpen())||!(ftmp3->GetNkeys()))
	    	{
				printf("input3 after zpc file %s error in opening!!!\n",FileList);
	    	}
	    	else
	    	{
				printf("input1 after zpc read in file %s\n",FileList);
				chain3->Add(FileList);
				fileNumber++;
	    	}
	    	delete ftmp3;
		}
    }

    printf(" files read in %d\n",fileNumber/3);

    AMPT* ampt1 = new AMPT(chain1);//zpct0 chain
    AMPT* ampt2 = new AMPT(chain2);//afterART chain
    AMPT* ampt3 = new AMPT(chain3);//afterzpc chain

	int mNEvents1 = (int)chain1->GetEntries();//number of events
	int mNEvents2 = (int)chain2->GetEntries();//number of events
	int mNEvents3 = (int)chain3->GetEntries();//number of events

 	TH1::SetDefaultSumw2(1);

 	cout<<"No. zpct0="<<mNEvents1<<" No. afterART="<<mNEvents2<<" No. after zpc="<<mNEvents3<<endl;

 	if(mNEvents1!=mNEvents2)return 0;
 	if(mNEvents1!=mNEvents3)return 0;

	hbf = new TH1D("hbf","hbf",40,-2.0,2.0);
	haf = new TH1D("haf","haf",40,-2.0,2.0);

	hv2pp = new TH1D("hv2pp","hv2pp",200,-1.0,1.0);
	hv2pd = new TH1D("hv2pd","hv2pd",200,-1.0,1.0);

	hv3pp = new TH1D("hv3pp","hv3pp",200,-1.0,1.0);
	hv3pd = new TH1D("hv3pd","hv3pd",200,-1.0,1.0);

	TPr2 = new TProfile("TPr2","TPr2",40,0,4.0,-2.0,2.0);

	double adp=2*PI;
	hadbf = new TH1D("hadbf","hadbf",40,-adp,adp);
	hadaf = new TH1D("hadaf","hadaf",40,-adp,adp);

	// ndna = new TH1D("ndna","ndna",10,0,5.0);
	// nhna = new TH1D("nhna","nhna",10,0,5.0);

	VDH2_pt = new TProfile("VDH2_pt","VDH2_pt",10,0,5.0,-1.0,1.0);
	VHH2_pt = new TProfile("VHH2_pt","VHH2_pt",10,0,5.0,-1.0,1.0);

	VDH3_pt = new TProfile("VDH3_pt","VDH3_pt",10,0,5.0,-1.0,1.0);
	VHH3_pt = new TProfile("VHH3_pt","VHH3_pt",10,0,5.0,-1.0,1.0);

	az_VDH2_pt = new TProfile("az_VDH2_pt","az_VDH2_pt",10,0,5.0,-1.0,1.0);
	az_VHH2_pt = new TProfile("az_VHH2_pt","az_VHH2_pt",10,0,5.0,-1.0,1.0);
	az_VDH3_pt = new TProfile("az_VDH3_pt","az_VDH3_pt",10,0,5.0,-1.0,1.0);
	az_VHH3_pt = new TProfile("az_VHH3_pt","az_VHH3_pt",10,0,5.0,-1.0,1.0);
	
	V2H_ebes = new TProfile("V2H_ebes","V2H_ebes",10,0,5.0,-1.0,1.0);
	V3H_ebes = new TProfile("V3H_ebes","V3H_ebes",10,0,5.0,-1.0,1.0);

	V2H_ps = new TProfile("V2H_ps","V2H_ps",10,0,5.0,-1.0,1.0);
	V3H_ps = new TProfile("V3H_ps","V3H_ps",10,0,5.0,-1.0,1.0);

	v2part_Dm_pt=new TProfile("v2part_Dm_pt","v2part_Dm_pt",10,0,5.0,-1.0,1.0);
	v3part_Dm_pt=new TProfile("v3part_Dm_pt","v3part_Dm_pt",10,0,5.0,-1.0,1.0);
	v4part_Dm_pt=new TProfile("v4part_Dm_pt","v4part_Dm_pt",10,0,5.0,-1.0,1.0);

	v2part_Hm_pt=new TProfile("v2part_Hm_pt","v2part_Hm_pt",10,0,5.0,-1.0,1.0);
	v3part_Hm_pt=new TProfile("v3part_Hm_pt","v3part_Hm_pt",10,0,5.0,-1.0,1.0);
	v4part_Hm_pt=new TProfile("v4part_Hm_pt","v4part_Hm_pt",10,0,5.0,-1.0,1.0);

	v2EP_Dm_pt=new TProfile("v2EP_Dm_pt","v2EP_Dm_pt",10,0,5.0,-5.0,5.0);
	v3EP_Dm_pt=new TProfile("v3EP_Dm_pt","v3EP_Dm_pt",10,0,5.0,-1.0,1.0);
	v4EP_Dm_pt=new TProfile("v4EP_Dm_pt","v4EP_Dm_pt",10,0,5.0,-1.0,1.0);

	H_b_v2e2 = new TProfile("H_b_v2e2","H_b_v2e2",20,0,20.0,-100,100);
	H_b_v3e3 = new TProfile("H_b_v3e3","H_b_v3e3",20,0,20.0,-100,100);
	D_b_v2e2 = new TProfile("D_b_v2e2","D_b_v2e2",20,0,20.0,-100,100);
	D_b_v3e3 = new TProfile("D_b_v3e3","D_b_v3e3",20,0,20.0,-100,100);

	H_Hb_v2e2 = new TH2D("H_Hb_v2e2","H_Hb_v2e2",20,0,20.0,200,-100,100);
	H_Hb_v3e3 = new TH2D("H_Hb_v3e3","H_Hb_v3e3",20,0,20.0,200,-100,100);
	D_Hb_v2e2 = new TH2D("D_Hb_v2e2","D_Hb_v2e2",20,0,20.0,200,-100,100);
	D_Hb_v3e3 = new TH2D("D_Hb_v3e3","D_Hb_v3e3",20,0,20.0,200,-100,100);

	// pp and 2c v(e) vs b
	H_b_v2_pp = new TProfile("H_b_v2_pp","H_b_v2_pp",20,0,20.0,-1,1);
	H_b_v3_pp = new TProfile("H_b_v3_pp","H_b_v3_pp",20,0,20.0,-1,1);
	D_b_v2_pp = new TProfile("D_b_v2_pp","D_b_v2_pp",20,0,20.0,-1,1);
	D_b_v3_pp = new TProfile("D_b_v3_pp","D_b_v3_pp",20,0,20.0,-1,1);

	H_b_e2_pp = new TProfile("H_b_e2_pp","H_b_e2_pp",20,0,20.0,-1,1);
	H_b_e3_pp = new TProfile("H_b_e3_pp","H_b_e3_pp",20,0,20.0,-1,1);
	D_b_e2_pp = new TProfile("D_b_e2_pp","D_b_e2_pp",20,0,20.0,-1,1);
	D_b_e3_pp = new TProfile("D_b_e3_pp","D_b_e3_pp",20,0,20.0,-1,1);

	H_b_v2_2c = new TProfile("H_b_v2_2c","H_b_v2_2c",20,0,20.0,-1,1);
	H_b_v3_2c = new TProfile("H_b_v3_2c","H_b_v3_2c",20,0,20.0,-1,1);
	D_b_v2_2c = new TProfile("D_b_v2_2c","D_b_v2_2c",20,0,20.0,-1,1);
	D_b_v3_2c = new TProfile("D_b_v3_2c","D_b_v3_2c",20,0,20.0,-1,1);

	H_b_e2_2c = new TProfile("H_b_e2_2c","H_b_e2_2c",20,0,20.0,-1,1);
	H_b_e3_2c = new TProfile("H_b_e3_2c","H_b_e3_2c",20,0,20.0,-1,1);
	D_b_e2_2c = new TProfile("D_b_e2_2c","D_b_e2_2c",20,0,20.0,-1,1);
	D_b_e3_2c = new TProfile("D_b_e3_2c","D_b_e3_2c",20,0,20.0,-1,1);

	H_b_v2_p_2c = new TProfile("H_b_v2_p_2c","H_b_v2_p_2c",20,0,20.0,-1,1);
	H_b_v3_p_2c = new TProfile("H_b_v3_p_2c","H_b_v3_p_2c",20,0,20.0,-1,1);
	D_b_v2_p_2c = new TProfile("D_b_v2_p_2c","D_b_v2_p_2c",20,0,20.0,-1,1);
	D_b_v3_p_2c = new TProfile("D_b_v3_p_2c","D_b_v3_p_2c",20,0,20.0,-1,1);

	//event by event average
	Dv2_2c_ebe = new TProfile("Dv2_2c_ebe","Dv2_2c_ebe",10,0,5,-1,1);
	Dv3_2c_ebe = new TProfile("Dv3_2c_ebe","Dv3_2c_ebe",10,0,5,-1,1);
	Hv2_2c_ebe = new TProfile("Hv2_2c_ebe","Hv2_2c_ebe",10,0,5,-1,1);
	Hv3_2c_ebe = new TProfile("Hv3_2c_ebe","Hv3_2c_ebe",10,0,5,-1,1);

	// H_b_v2e2->SetMinimum(-INFINITY);
	// H_b_v2e2->SetMaximum(INFINITY);
	// H_b_v3e3->SetMinimum(-INFINITY);
	// H_b_v3e3->SetMaximum(INFINITY);
	// D_b_v2e2->SetMinimum(-INFINITY);
	// D_b_v2e2->SetMaximum(INFINITY);
	// D_b_v3e3->SetMinimum(-INFINITY);
	// D_b_v3e3->SetMaximum(INFINITY);

	epsilon2_centrality=new TH2D("epsilon2_centrality","epsilon2_centrality",30,0,15.0,100,0,1.0);
	epsilon3_centrality=new TH2D("epsilon3_centrality","epsilon3_centrality",30,0,15.0,100,0,1.0);
	epsilon4_centrality=new TH2D("epsilon4_centrality","epsilon4_centrality",30,0,15.0,100,0,1.0);

	dphi_Dm_pt=new TH2D("dphi_Dm_pt","dphi_Dm_pt",50,0,5.0,40,0,PI);
	dphi_Hm_pt=new TH2D("dphi_Hm_pt","dphi_Hm_pt",50,0,5.0,40,0,PI);

	dphi_c_pt=new TH2D("dphi_c_pt","dphi_c_pt",50,0,5.0,40,0,PI);
	dphi_uds_pt=new TH2D("dphi_uds_pt","dphi_uds_pt",50,0,5.0,40,0,PI);

	spect_D0_pt=new TH1D("spect_D0_pt","spect_D0_pt",50,0,10.0);

	he2q2=new TH2D("he2q2","he2q2",100,0,1.0,100,0,1.0);
	he3q3=new TH2D("he3q3","he3q3",100,0,1.0,100,0,1.0);
	he4q4=new TH2D("he4q4","he4q4",100,0,1.0,100,0,1.0);

	he2qpt=new TH2D("he2qpt","he2qpt",100,0,1.0,100,0,1.0);
	he3qpt=new TH2D("he3qpt","he3qpt",100,0,1.0,100,0,1.0);
	he4qpt=new TH2D("he4qpt","he4qpt",100,0,1.0,100,0,1.0);
	
	//envn pp
	De2v2_pp=new TH2D("De2v2_pp","De2v2_pp",200,0,1.0,300,-1.5,1.5);
	He2v2_pp=new TH2D("He2v2_pp","He2v2_pp",200,0,1.0,300,-1.5,1.5);
	De3v3_pp=new TH2D("De3v3_pp","De3v3_pp",200,0,1.0,300,-1.5,1.5);
	He3v3_pp=new TH2D("He3v3_pp","He3v3_pp",200,0,1.0,300,-1.5,1.5);

	//envn 2c
	De2v2_2c=new TH2D("De2v2_2c","De2v2_2c",200,0,1.0,2000,-10,10);
	He2v2_2c=new TH2D("He2v2_2c","He2v2_2c",200,0,1.0,300,-1.5,1.5);
	De3v3_2c=new TH2D("De3v3_2c","De3v3_2c",200,0,1.0,300,-1.5,1.5);
	He3v3_2c=new TH2D("He3v3_2c","He3v3_2c",200,0,1.0,300,-1.5,1.5);

	hHad0  = new TH1D("hHad0","hHad0",100,-adp,adp);
	hHadmp = new TH1D("hHadmp","hHadmp",100,-adp,adp); 
	hDad0  = new TH1D("hDad0","hDad0",100,-adp,adp);
	hDadmp = new TH1D("hDadmp","hDadmp",100,-adp,adp);

	hHadm3 = new TH1D("hHadm3","hHadm3",100,-adp,adp); 
	hDadm3 = new TH1D("hDadm3","hDadm3",100,-adp,adp);

	hpsi2ad = new TH1D("hpsi2ad","hpsi2ad",100,-adp,adp);
	hpsi3ad = new TH1D("hpsi3ad","hpsi3ad",100,-adp,adp);

	// test something
	h_numd = new TH1D("h_numd","h_numd",10000,0,1000000);
	h_numh = new TH1D("h_numh","h_numh",10000,0,1000000);
	h_rati = new TH2D("h_rati","h_rati",10,0,5.0,1000,0,10000);
	h_arrh = new TH2D("h_arrh","h_arrh",10,0,5.0,10000,-10000,10000);
	h_arrd = new TH2D("h_arrd","h_arrd",10,0,5.0,10000,-10000,10000);

	D_flow_nass = new TH2D("D_flow_nass","D_flow_nass",2000,-1000,1000,2000,-1000,1000);
	D_f_n_test = new TH1D("D_f_n_test","D_f_n_test",2000,-10,10);

	// bad event
	h_bad_mult = new TH2D("h_bad_mult","h_bad_mult",500,-2000,2000,500,-2000,2000);
	h_bad_e2e3 = new TH2D("h_bad_e2e3","h_bad_e2e3",500,-2000,2000,500,-2000,2000);
	h_bad_dhmu = new TH2D("h_bad_dhmu","h_bad_dhmu",500,-5000,5000,500,-5000,5000);

	h_bad_aft = new TH2D("h_bad_aft","h_bad_aft",200,-10,10,200,-10,10);
	h_bad_zpc = new TH2D("h_bad_zpc","h_bad_zpc",400,-20,20,400,-20,20);

	gen = 0;
	ebin1  = 0; 
	ebin2  = 0; 
	ebin3  = 0; 
	ebin4  = 0; 
	ebin5  = 0; 
	ebin6  = 0; 
	ebin7  = 0; 
	ebin8  = 0; 
	ebin9  = 0; 
	ebin10 = 0; 
	ebin11 = 0; 
	ebin12 = 0; 
	ebin13 = 0; 
	ebin14 = 0; 
	ebin15 = 0; 
	ebin16 = 0; 
	ebin17 = 0; 
	ebin18 = 0; 
	ebin19 = 0; 
	ebin20 = 0; 
	ebin21 = 0; 
	ebin22 = 0; 
	ebin23 = 0; 
	ebin24 = 0; 
	ebin25 = 0; 
	ebin26 = 0; 
	ebin27 = 0; 
	ebin28 = 0; 
	ebin29 = 0; 
	ebin30 = 0; 
	ebin31 = 0; 
	ebin32 = 0; 
	ebin33 = 0; 
	ebin34 = 0; 
	ebin35 = 0; 
	ebin36 = 0; 
	ebin37 = 0; 
	ebin38 = 0; 
	ebin39 = 0; 
	ebin40 = 0; 
	ebin41 = 0; 
	ebin42 = 0; 
	ebin43 = 0; 
	ebin44 = 0; 
	ebin45 = 0; 
	ebin46 = 0; 
	ebin47 = 0; 
	ebin48 = 0; 

	//	event plane resolution
	event_plane_resolution=new TGraph(10,res_x,res_y);

	for(int mNEventCount = 0; mNEventCount < mNEvents1; mNEventCount++)
	{
		if(mNEventCount%1000==0) cout << "Working on event #" << mNEventCount << endl;	

		chain1->GetEntry(mNEventCount);
		chain2->GetEntry(mNEventCount);
		chain3->GetEntry(mNEventCount);
		
	  	sum_r2sin2phi=0.;
	  	sum_r2cos2phi=0.;			
	  	sum_r2sin3phi=0.;
	  	sum_r2cos3phi=0.;
	  	sum_r2sin4phi=0.;
	  	sum_r2cos4phi=0.;

	  	sum_r2=0;

		sum_r2sin2phi_charm=0.;
	  	sum_r2cos2phi_charm=0.;			
	  	sum_r2sin3phi_charm=0.;
	  	sum_r2cos3phi_charm=0.;
	  	sum_r2sin4phi_charm=0.;
	  	sum_r2cos4phi_charm=0.;

	  	sum_r2_charm=0;

		sum_r2sin2phi_ep=0.;
	  	sum_r2cos2phi_ep=0.;			
	  	sum_r2sin3phi_ep=0.;
	  	sum_r2cos3phi_ep=0.;
	  	sum_r2sin4phi_ep=0.;
	  	sum_r2cos4phi_ep=0.;

	  	sum_r2_ep=0;

	  	//sum_v2=0;
	  	//sum_v3=0;
	  	//sum_v4=0;
	  	//ntrack=0;
	  	//sum_num=0;

		qx2=0;//EP
	  	qx3=0;
	  	qx4=0;
	  	qy2=0;
	  	qy3=0;
	  	qy4=0;
		multiq=0;

	  	centrality=0;
      	nmult1 = 0;
	  	nmult2 = 0;

		charm_mu = 0;
		Dmeson_mu = 0;

	  	if(passEvent(ampt1,ampt2))
	  	{
			//event by event 2c
			// Hv2_e_s = new TProfile("Hv2_e_s","Hv2_e_s",10,0,5,-1,1);
			// Hv3_e_s = new TProfile("Hv3_e_s","Hv3_e_s",10,0,5,-1,1);
			// Dv2_e_s = new TProfile("Dv2_e_s","Dv2_e_s",10,0,5,-1,1);
			// Dv3_e_s = new TProfile("Dv3_e_s","Dv3_e_s",10,0,5,-1,1);

	  	  	for(int i=0; i<nmult1; i++){
	  	    	parton_loop(ampt1, i, ampt2);//zpct0					  
	  	  	}
			//是不是忘记取平均了？nt 比值已经消除取平均了
			// sum_r2sin2phi = sum_r2sin2phi / nmult1;
			// sum_r2cos2phi = sum_r2cos2phi / nmult1;
			// sum_r2 = sum_r2 / nmult1;
			// sum_r2sin3phi = sum_r2sin3phi / nmult1;
			// sum_r2cos3phi = sum_r2cos3phi / nmult1;
			// sum_r2sin4phi = sum_r2sin4phi / nmult1;
			// sum_r2cos4phi = sum_r2cos4phi / nmult1;

			// if (charm_mu<4)
			// {
			// 	continue;	
			// }

			if (Dmeson_mu<4)
			{
				continue;	
			}

			//par-plane
	    	psi2=(atan2(sum_r2sin2phi,sum_r2cos2phi)+PI)/2.;
	    	psi3=(atan2(sum_r2sin3phi,sum_r2cos3phi)+PI)/3.;
	    	psi4=(atan2(sum_r2sin4phi,sum_r2cos4phi)+PI)/4.;

			if(psi2 > PI/2)
			{
				psi2 = psi2 - PI;
			}
			hpsi2ad->Fill(psi2);
			if(psi3 > PI/3)
			{
				psi3 = psi3 - 2*PI/3;
			}
			hpsi3ad->Fill(psi3);

			// epsilon2=sqrt((sum_r2sin2phi_charm*sum_r2sin2phi_charm + sum_r2cos2phi_charm*sum_r2cos2phi_charm)/(sum_r2_charm*sum_r2_charm));
	    	// epsilon3=sqrt((sum_r2sin3phi_charm*sum_r2sin3phi_charm + sum_r2cos3phi_charm*sum_r2cos3phi_charm)/(sum_r2_charm*sum_r2_charm));
	    	// epsilon4=sqrt((sum_r2sin4phi_charm*sum_r2sin4phi_charm + sum_r2cos4phi_charm*sum_r2cos4phi_charm)/(sum_r2_charm*sum_r2_charm));

	    	epsilon2=sqrt((sum_r2sin2phi*sum_r2sin2phi+sum_r2cos2phi*sum_r2cos2phi)/(sum_r2*sum_r2));
	    	epsilon3=sqrt((sum_r2sin3phi*sum_r2sin3phi+sum_r2cos3phi*sum_r2cos3phi)/(sum_r2*sum_r2));
	    	epsilon4=sqrt((sum_r2sin4phi*sum_r2sin4phi+sum_r2cos4phi*sum_r2cos4phi)/(sum_r2*sum_r2));

			if (epsilon2>0.99999)
			{
				h_bad_mult->Fill(sum_r2sin2phi,sum_r2cos2phi);
				h_bad_e2e3->Fill(sum_r2sin3phi,sum_r2cos3phi);
				h_bad_dhmu->Fill(sum_r2,sum_r2);
	  	  	
				r_bads = 0.;
				rsin_bads = 0.;
				rcos_bads = 0.;

				for(int i=0; i<nmult1; i++){
					bad_event(ampt1,i);
				}

				// cout<<"r_bads sum is: "<<r_bads<<", rsin_bads sum is:"<<rsin_bads<<", rcos_bads sum is: "<<rcos_bads<<endl;

				// cout<<"sum_r2sin2phi is: "<<sum_r2sin2phi<<", sum_r2cos2phi is: "<<sum_r2cos2phi<<", sum_r2 is: "<<sum_r2<<", r_phi_bad is:"<<r_phi_bad<<endl;
			}

			double cce2 = epsilon2 * epsilon2;
			double cce3 = epsilon3 * epsilon3;

			// extract event below guassian distribution from all data
			// if(gen>=99976) {cout<<"gen is "<<gen<<endl; break;}

			// extract event below guassian distribution from all data
			// if     ( (epsilon2>=0.210) && (epsilon2<0.215) )    { ebin1  = ebin1  + 4;   if(ebin1  > 20    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.215) && (epsilon2<0.220) )	{ ebin2  = ebin2  + 4;   if(ebin2  > 112   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.220) && (epsilon2<0.225) )	{ ebin3  = ebin3  + 4;   if(ebin3  > 486   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.225) && (epsilon2<0.230) )	{ ebin4  = ebin4  + 4;   if(ebin4  > 1654  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.230) && (epsilon2<0.235) )	{ ebin5  = ebin5  + 4;   if(ebin5  > 4406  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.235) && (epsilon2<0.240) )	{ ebin6  = ebin6  + 4;   if(ebin6  > 9185  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.240) && (epsilon2<0.245) )	{ ebin7  = ebin7  + 4;   if(ebin7  > 14988 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.245) && (epsilon2<0.250) )	{ ebin8  = ebin8  + 4;   if(ebin8  > 19146 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.250) && (epsilon2<0.255) )	{ ebin9  = ebin9  + 4;   if(ebin9  > 19146 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.255) && (epsilon2<0.260) )	{ ebin10 = ebin10 + 4;   if(ebin10 > 14988 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.260) && (epsilon2<0.265) )	{ ebin11 = ebin11 + 4;   if(ebin11 > 9185  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.265) && (epsilon2<0.270) )	{ ebin12 = ebin12 + 4;   if(ebin12 > 4406  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.270) && (epsilon2<0.275) )	{ ebin13 = ebin13 + 4;   if(ebin13 > 1654  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.275) && (epsilon2<0.280) )	{ ebin14 = ebin14 + 4;   if(ebin14 > 486   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.280) && (epsilon2<0.285) )	{ ebin15 = ebin15 + 4;   if(ebin15 > 112   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.285) && (epsilon2<0.290) )	{ ebin16 = ebin16 + 4;   if(ebin16 > 20    ) { continue; } gen = gen + 4;}

			// if     ( (epsilon2>=0.01) && (epsilon2<0.02) )	{ ebin1  = ebin1  + 4;   if(ebin1  > 3    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.02) && (epsilon2<0.03) )	{ ebin2  = ebin2  + 4;   if(ebin2  > 6    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.03) && (epsilon2<0.04) )	{ ebin3  = ebin3  + 4;   if(ebin3  > 11   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.04) && (epsilon2<0.05) )	{ ebin4  = ebin4  + 4;   if(ebin4  > 20   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.05) && (epsilon2<0.06) )	{ ebin5  = ebin5  + 4;   if(ebin5  > 34   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.06) && (epsilon2<0.07) )	{ ebin6  = ebin6  + 4;   if(ebin6  > 58   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.07) && (epsilon2<0.08) )	{ ebin7  = ebin7  + 4;   if(ebin7  > 95   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.08) && (epsilon2<0.09) )	{ ebin8  = ebin8  + 4;   if(ebin8  > 153  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.09) && (epsilon2<0.10) )	{ ebin9  = ebin9  + 4;   if(ebin9  > 6    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.10) && (epsilon2<0.11) )	{ ebin10 = ebin10 + 4;   if(ebin10 > 14   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.11) && (epsilon2<0.12) )	{ ebin11 = ebin11 + 4;   if(ebin11 > 34   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.12) && (epsilon2<0.13) )	{ ebin12 = ebin12 + 4;   if(ebin12 > 77   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.13) && (epsilon2<0.14) )	{ ebin13 = ebin13 + 4;   if(ebin13 > 163  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.14) && (epsilon2<0.15) )	{ ebin14 = ebin14 + 4;   if(ebin14 > 323  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.15) && (epsilon2<0.16) )	{ ebin15 = ebin15 + 4;   if(ebin15 > 601  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.16) && (epsilon2<0.17) )	{ ebin16 = ebin16 + 4;   if(ebin16 > 1053 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.17) && (epsilon2<0.18) )	{ ebin17 = ebin17 + 4;   if(ebin17 > 20   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.18) && (epsilon2<0.19) )	{ ebin18 = ebin18 + 4;   if(ebin18 > 112  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.19) && (epsilon2<0.20) )	{ ebin19 = ebin19 + 4;   if(ebin19 > 486  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.20) && (epsilon2<0.21) )	{ ebin20 = ebin20 + 4;   if(ebin20 > 1654 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.21) && (epsilon2<0.22) )	{ ebin21 = ebin21 + 4;   if(ebin21 > 4406 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.22) && (epsilon2<0.23) )	{ ebin22 = ebin22 + 4;   if(ebin22 > 9185 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.23) && (epsilon2<0.24) )	{ ebin23 = ebin23 + 4;   if(ebin23 > 14988) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.24) && (epsilon2<0.25) )	{ ebin24 = ebin24 + 4;   if(ebin24 > 19146) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.25) && (epsilon2<0.26) )	{ ebin25 = ebin25 + 4;   if(ebin25 > 19146) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.26) && (epsilon2<0.27) )	{ ebin26 = ebin26 + 4;   if(ebin26 > 14988) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.27) && (epsilon2<0.28) )	{ ebin27 = ebin27 + 4;   if(ebin27 > 9185 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.28) && (epsilon2<0.29) )	{ ebin28 = ebin28 + 4;   if(ebin28 > 4406 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.29) && (epsilon2<0.30) )	{ ebin29 = ebin29 + 4;   if(ebin29 > 1654 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.30) && (epsilon2<0.31) )	{ ebin30 = ebin30 + 4;   if(ebin30 > 486  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.31) && (epsilon2<0.32) )	{ ebin31 = ebin31 + 4;   if(ebin31 > 112  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.32) && (epsilon2<0.33) )	{ ebin32 = ebin32 + 4;   if(ebin32 > 20   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.33) && (epsilon2<0.34) )	{ ebin33 = ebin33 + 4;   if(ebin33 > 1053 ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.34) && (epsilon2<0.35) )	{ ebin34 = ebin34 + 4;   if(ebin34 > 601  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.35) && (epsilon2<0.36) )	{ ebin35 = ebin35 + 4;   if(ebin35 > 323  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.36) && (epsilon2<0.37) )	{ ebin36 = ebin36 + 4;   if(ebin36 > 163  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.37) && (epsilon2<0.38) )	{ ebin37 = ebin37 + 4;   if(ebin37 > 77   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.38) && (epsilon2<0.39) )	{ ebin38 = ebin38 + 4;   if(ebin38 > 34   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.39) && (epsilon2<0.40) )	{ ebin39 = ebin39 + 4;   if(ebin39 > 14   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.40) && (epsilon2<0.41) )	{ ebin40 = ebin40 + 4;   if(ebin40 > 6    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.41) && (epsilon2<0.42) )	{ ebin41 = ebin41 + 4;   if(ebin41 > 153  ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.42) && (epsilon2<0.43) )	{ ebin42 = ebin42 + 4;   if(ebin42 > 95   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.43) && (epsilon2<0.44) )	{ ebin43 = ebin43 + 4;   if(ebin43 > 58   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.44) && (epsilon2<0.45) )	{ ebin44 = ebin44 + 4;   if(ebin44 > 34   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.45) && (epsilon2<0.46) )	{ ebin45 = ebin45 + 4;   if(ebin45 > 20   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.46) && (epsilon2<0.47) )	{ ebin46 = ebin46 + 4;   if(ebin46 > 11   ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.47) && (epsilon2<0.48) )	{ ebin47 = ebin47 + 4;   if(ebin47 > 6    ) { continue; } gen = gen + 4;}
			// else if( (epsilon2>=0.48) && (epsilon2<0.49) )	{ ebin48 = ebin48 + 4;   if(ebin48 > 3    ) { continue; } gen = gen + 4;}
			// else   {continue;}  

			// az - after zpc
			// for(int az=0; az<nmult2; az++){
	    	// 	udsc_loop(ampt3, az);
        	// }	
			// make_c_correlation(c_quark,uds_quark,dphi_c_pt);
 			// make_uds_correlation(uds_quark,dphi_uds_pt);

			if(epsilon2>0.99)
			{
				continue;
			}

			for(int i=0; i<nmult2; i++){
	  	    	event_loop(ampt2, i);					  
	  	  	}
			
			q2=sqrt(qx2*qx2+qy2*qy2)/multiq;
			q3=sqrt(qx3*qx3+qy3*qy3)/multiq;
			q4=sqrt(qx4*qx4+qy4*qy4)/multiq;

			q2pt=sqrt(sum_r2sin2phi_ep*sum_r2sin2phi_ep+sum_r2cos2phi_ep*sum_r2cos2phi_ep)/sum_r2_ep;
			q3pt=sqrt(sum_r2sin3phi_ep*sum_r2sin3phi_ep+sum_r2cos3phi_ep*sum_r2cos3phi_ep)/sum_r2_ep;
			q4pt=sqrt(sum_r2sin4phi_ep*sum_r2sin4phi_ep+sum_r2cos4phi_ep*sum_r2cos4phi_ep)/sum_r2_ep;

			he2q2->Fill(epsilon2,q2);
			he3q3->Fill(epsilon3,q3);
			he4q4->Fill(epsilon4,q4);

			psiEP2=atan2(sum_r2sin2phi_ep,sum_r2cos2phi_ep)/2.0;
			psiEP3=atan2(sum_r2sin3phi_ep,sum_r2cos3phi_ep)/3.0;
			psiEP4=atan2(sum_r2sin4phi_ep,sum_r2cos4phi_ep)/4.0;

			v2EP_s = 0;
			v3EP_s = 0;
			v4EP_s = 0;

			VH2_s = 0;
			VH3_s = 0;
			VH4_s = 0;

			VD2_s = 0;
			VD3_s = 0;
			VD4_s = 0;

			hmu = 0;
			dmu = 0;

			VD2_2c_s = 0;
			VD3_2c_s = 0;
			VH2_2c_s = 0;
			VH3_2c_s = 0;

			numhs = 0;
			numds = 0;

			fill(arrv2d, arrv2d + 10, 0);
			fill(arrv2h, arrv2h + 10, 0);

			fill(numv2d, numv2d + 10, 0);
			fill(numv2h, numv2h + 10, 0);

        	for(int j=0; j<nmult2; j++){
	    		hadron_loop(ampt2, j, event_plane_resolution);
        	}

			makeDmcorrelation(Dm,Hm,dphi_Dm_pt);
 			makeHmcorrelation(Hm,dphi_Hm_pt);	

			v2EP_s = v2EP_s/hmu;
			v3EP_s = v3EP_s/hmu;
			v4EP_s = v4EP_s/hmu;

			if (dmu>0.5)
			{
				VD2_s = VD2_s / dmu;
				hv2pd->Fill(VD2_s);
			}

			VH2_s = VH2_s / hmu;
			hv2pp->Fill(VH2_s); // TH1D

			if (dmu>0.5)
			{
				VD3_s = VD3_s / dmu;
				hv3pd->Fill(VD3_s);
			}

			VH3_s = VH3_s / hmu;
			hv3pp->Fill(VH3_s); // TH1D

			// v2part_Dm_pt->Fill(pt,VD2_s);//填入histogram，得到vn的pt关系 不好分pt bin
			// v3part_Dm_pt->Fill(pt,VD3_s);
			// v4part_Dm_pt->Fill(pt,v4);

			TPr2->Fill(1.95,VH2_s); // TProfile

			if(nmult2<5) continue;

			if(epsilon2>0&&epsilon2<1.0){epsilon2_centrality->Fill(centrality,epsilon2);}  //计算出的各阶偏心率填到histogram
			if(epsilon3>0&&epsilon3<1.0){epsilon3_centrality->Fill(centrality,epsilon3);}
			if(epsilon4>0&&epsilon4<1.0){epsilon4_centrality->Fill(centrality,epsilon4);}
						
			D_flow_nass->Fill(NassHH,VH2_2c_s);
			VH2_2c_s = VH2_2c_s / NassHH;
			VH3_2c_s = VH3_2c_s / NassHH;

			VH2_2c_b = sqrt(VH2_2c_s);
			VH3_2c_b = sqrt(VH3_2c_s);

			if (dmu>0.5)
			{
				if(NassDH!=0){
					noverd = noverd + 1;
					oneoverd = oneoverd + 1/NassDH;

					h_numd->Fill(NassDH);
				}
				
				double D_v2e2 = VD2_s / epsilon2;			
				double D_v3e3 = VD3_s / epsilon3;
				D_b_v2e2->Fill(centrality,D_v2e2);
				D_b_v3e3->Fill(centrality,D_v3e3);
				
				D_Hb_v2e2->Fill(centrality,D_v2e2);
				D_Hb_v3e3->Fill(centrality,D_v3e3);

				//pp d
				D_b_v2_pp->Fill(centrality,VD2_s);
				D_b_v3_pp->Fill(centrality,VD3_s);
				D_b_e2_pp->Fill(centrality,epsilon2);
				D_b_e3_pp->Fill(centrality,epsilon3);

				// D_flow_nass->Fill(NassDH,VD2_2c_s);
				//2c d
				VD2_2c_s = VD2_2c_s / NassDH ;
				VD3_2c_s = VD3_2c_s / NassDH ;

				D_f_n_test->Fill(VD2_2c_s);

				VD2_2c_b = VD2_2c_s / sqrt(VH2_2c_s);
				VD3_2c_b = VD3_2c_s / sqrt(VH3_2c_s);
				
				if (NassDH!=0){
					D_b_v2_2c->Fill(centrality,VD2_2c_s);
					D_b_v3_2c->Fill(centrality,VD3_2c_s);
					D_b_e2_2c->Fill(centrality,cce2);
					D_b_e3_2c->Fill(centrality,cce3);
				}

				De2v2_pp->Fill(epsilon2,VD2_s); // long
				De3v3_pp->Fill(epsilon3,VD3_s); 
				De2v2_2c->Fill(epsilon2,VD2_2c_b); // long
				De3v3_2c->Fill(epsilon3,VD3_2c_b); 
			}
			if(NassHH!=0){
				noverh = noverh + 1;
				oneoverh = oneoverh + 10000/NassHH;

				h_numh->Fill(NassHH);
			}


			if(dmu==0){
				nassd0 = nassd0 + 1;
			}
			if(hmu==0){
				nassh0 = nassh0 + 1;
			}

			He2v2_pp->Fill(epsilon2,VH2_s);	// long
			He3v3_pp->Fill(epsilon3,VH3_s);
			He2v2_2c->Fill(epsilon2,VH2_2c_b);	// long
			He3v3_2c->Fill(epsilon3,VH3_2c_b);

			he2qpt->Fill(epsilon2,q2pt);
			he3qpt->Fill(epsilon3,q3pt);
			he4qpt->Fill(epsilon4,q4pt);

			double H_v2e2 = VH2_s / epsilon2;
			double H_v3e3 = VH3_s / epsilon3;
			
			H_b_v2e2->Fill(centrality,H_v2e2);
			H_b_v3e3->Fill(centrality,H_v3e3);

			H_Hb_v2e2->Fill(centrality,H_v2e2);
			H_Hb_v3e3->Fill(centrality,H_v3e3);

			//pp h
			H_b_v2_pp->Fill(centrality,VH2_s);
			H_b_v3_pp->Fill(centrality,VH3_s);
			H_b_e2_pp->Fill(centrality,epsilon2);
			H_b_e3_pp->Fill(centrality,epsilon3);

			//2c h
			if (NassHH!=0){
				H_b_v2_2c->Fill(centrality,VH2_2c_s);
				H_b_v3_2c->Fill(centrality,VH3_2c_s);
				H_b_e2_2c->Fill(centrality,cce2);
				H_b_e3_2c->Fill(centrality,cce3);
			}
			
			// //event by event average
			// for(int ebe = 1; ebe<=10; ebe++)
			// {
			// 	double x_dv2 = 0.5*ebe - 0.25;
			// 	double x_dv3 = 0.5*ebe - 0.25;
			// 	double x_hv2 = 0.5*ebe - 0.25;
			// 	double x_hv3 = 0.5*ebe - 0.25;				

			// 	int aebe = ebe - 1;

			// 	if (dmu>0.5)
			// 	{
			// 		double y_dv2 = arrv2d[aebe]/numv2d[aebe]; // Dv2_e_s->GetBinContent(ebe);
			// 		double y_dv3 = Dv3_e_s->GetBinContent(ebe);	
										
			// 		// Dv2_2c_ebe->Fill(x_dv2,y_dv2);						
			// 		Dv3_2c_ebe->Fill(x_dv3,y_dv3);
			// 	}
			// 	double y_hv2 = arrv2h[aebe]/numv2h[aebe]; // Hv2_e_s->GetBinContent(ebe);
			// 	double y_hv3 = Hv3_e_s->GetBinContent(ebe);

			// 	// Hv2_2c_ebe->Fill(x_hv2,y_hv2);
			// 	Hv3_2c_ebe->Fill(x_hv3,y_hv3);

			// 	// test something				
			// 	h_arrh->Fill(x_dv2,arrv2h[aebe]);
			// 	h_arrd->Fill(x_dv2,arrv2d[aebe]);
			// }
			// h_rati->Fill(1.1,Dm.size());

			// V2H_ebes->Fill(1.1,VH2_2c_s);
			// V3H_ebes->Fill(1.1,VH3_2c_s);

			// test something
			// h_numd->Fill(NassDH);
			// h_numh->Fill(NassHH);

			// cout<<"numhs is "<<numhs<<endl;
			
			// delete Hv2_e_s;
			// delete Hv3_e_s;
			// delete Dv2_e_s;
			// delete Dv3_e_s;
		}

		c_quark.clear();
		uds_quark.clear();

		Dm.clear();
		Hm.clear();
	}//end event loop

	writeHistograms(FileOutput);//histogram write to file
	deleteHistograms();//释放内存
	delete chain1;
	delete chain2;
	delete chain3;

	cout<<"oneoverh is "<<oneoverh<<" ; oneoverd is "<<oneoverd<<endl;
	cout<<"noverh is "<<noverh<<" ; noverd is "<<noverd<<endl;

	cout<<"nassh0 is "<<nassh0<<" ; nassd0 is "<<nassd0<<endl;

	return 0;
}


bool passEvent(AMPT* event1,AMPT* event2) //判断两文件是否是同一事件初末态
{
	//in >> nevent>>nrun>>nls>>impactpar>>NpartP>>NpartT>>NELP>>NINP>>NELT>>NINT;//input event information
	double impactpar1 = (double)event1->Event_impactpar;
	double impactpar2 = (double)event2->Event_impactpar;

	centrality=impactpar1;

	nmult1=(int)event1->Event_multi;
	nmult2=(int)event2->Event_multi;

 	if(centrality>15||fabs(impactpar1-impactpar2)>0.01){return 0;}//make sure initial and final data from the same event

 	if(nmult2<5){return 0;}//reject bad event

	// if((centrality<=7)||(centrality>=8)) {return 0;} // cut centrality
	//if(!(impactpar1>0.0&&impactpar1<13.5)){return 0;}//MLCC

    return 1;
}


bool hadron_loop(AMPT* track,int i,TGraph *tgr)
{
	//in >> id >> px >> py >> pz >> am >> others1 >>others2 >>others3 >>others4;//input track(particle) information
	
	int id = (int)track->ID[i];
	//double index = (double)track->Indx[i];
	double px = (double)track->Px[i];
	double py = (double)track->Py[i];
	double pz = (double)track->Pz[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];
	double am = (double)track->Mass[i];
	TVector3 mP(px,py,pz);

	if(mP.Pt()<0.3 || mP.Pt()>99.0) return 0;      

	double theta=acos(pz/sqrt(px*px+py*py+pz*pz));
	double pseorap=-log(tan(theta/2.));
	double pt=sqrt(px*px+py*py);
	double phi=atan2(py,px); //momentum coordinate MLCC!!
		
	// double phiaf = phi-psi2;
	// hsdbf->Fill(phi);
	// hsdaf->Fill(phiaf);

	if(fabs(pseorap)<eta_cut && (abs(id)==2212||abs(id)==321||abs(id)==211) && pt<5.0)//选出普通强子存到数组
	{
		TVector3 hHv(pseorap,phi,pt);//eta/phi/pt
		Hm.push_back(hHv);

		double phimp2 = phi - psi2;		
		if(phimp2 > PI)
		{
			phimp2 = phimp2 - 2*PI;	
		}
		else if(phimp2 < -PI)
		{
			phimp2 = phimp2 + 2*PI;	
		}		
		hHad0->Fill(phi);
		hHadmp->Fill(phimp2);

		double phimp3 = phi - psi3;	
		if(phimp3 > PI)
		{
			phimp3 = phimp3 - 2*PI;	
		}
		else if(phimp3 < -PI)
		{
			phimp3 = phimp3 + 2*PI;	
		}	
		hHadm3->Fill(phimp3);

		double v2_h=cos(2.0*phimp2);//PP
		double v3_h=cos(3.0*(phi-psi3));
		double v4_h=cos(4.0*(phi-psi4));

		v2part_Hm_pt->Fill(pt,v2_h);//填入histogram，得到vn的pt关系
		v3part_Hm_pt->Fill(pt,v3_h);
		v4part_Hm_pt->Fill(pt,v4_h);

		VH2_s = VH2_s + v2_h;
		VH3_s = VH3_s + v3_h;
		VH4_s = VH4_s + v4_h;

		v2EP_h=cos(2.0*(phi-psiEP2));//EP
		v3EP_h=cos(3.0*(phi-psiEP3));
		v4EP_h=cos(4.0*(phi-psiEP4));

		v2EP_s = v2EP_s + v2EP_h;
		v3EP_s = v3EP_s + v3EP_h;
		v4EP_s = v4EP_s + v4EP_h;
	
		hmu = hmu + 1;

		hbf->Fill(v3_h);
	}

	if(fabs(pseorap)<eta_cut && (abs(id)==421||abs(id)==411||abs(id)==413) && pt<5.0)//选出D介子存到数组
	{
		TVector3 hDv(pseorap,phi,pt);//eta/phi/pt
		Dm.push_back(hDv);

		{spect_D0_pt->Fill(pt,1.0/(2.0*PI*pt*0.2*2.0));}//后一项是填谱的权重，2PIptdptdy
	
		double phimp2 = phi - psi2;
		if(phimp2 > PI)
		{
			phimp2 = phimp2 - 2*PI;	
		}
		else if(phimp2 < -PI)
		{
			phimp2 = phimp2 + 2*PI;	
		}
		hDad0->Fill(phi);
		hDadmp->Fill(phimp2);

		double phimp3 = phi - psi3;	
		if(phimp3 > PI)
		{
			phimp3 = phimp3 - 2*PI;	
		}
		else if(phimp3 < -PI)
		{
			phimp3 = phimp3 + 2*PI;	
		}		
		hDadm3->Fill(phimp3);

		double v2=cos(2.0*phimp2);//计算D的各阶flow,用的是参与平面法
		double v3=cos(3.0*(phi-psi3));
		double v4=cos(4.0*(phi-psi4));

		VD2_s = VD2_s + v2;
		VD3_s = VD3_s + v3;

		v2part_Dm_pt->Fill(pt,v2);//填入histogram，得到vn的pt关系
		v3part_Dm_pt->Fill(pt,v3);
		v4part_Dm_pt->Fill(pt,v4);

		v2EP=cos(2.0*(phi-psiEP2));//EP
		v3EP=cos(3.0*(phi-psiEP3));
		v4EP=cos(4.0*(phi-psiEP4));

		double e_p_r=event_plane_resolution->Eval(centrality);
		double v2EP_aft=v2EP/e_p_r;

		v2EP_Dm_pt->Fill(pt,v2EP_aft);
		v3EP_Dm_pt->Fill(pt,v3EP);
		v4EP_Dm_pt->Fill(pt,v4EP);

		haf->Fill(v3);
		// hbf->Fill(v2EP);
		// haf->Fill(v2EP_aft);

		dmu = dmu + 1;			
	}

    double teval=event_plane_resolution->Eval(centrality);
	if((teval)<0||(teval>1)){
		cout<<"event_plane_resolution is "<<teval<<" ;  with centrality is "<<centrality<<endl;
	}

	return 1;
}


bool parton_loop(AMPT* track,int i,AMPT* track2)
{
	int id = (int)track->ID[i];
	//double index = (double)track->Indx[i];
	double px = (double)track->Px[i];
	double py = (double)track->Py[i];
	double pz = (double)track->Pz[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];
	double am = (double)track->Mass[i];
	TVector3 mP(px,py,pz);

	if((x*x+y*y)<1.0e-7) return 0;
	//if(mP.Pt()<10e-7 || mP.Mag()<10e-7 || mP.Mag()>90 || mP.Pt()>90 || fabs(x)<=0.0001 || fabs(y)<=0.0001) return 0;
	if(sqrt(px*px+py*py)<1.0e-7||sqrt(px*px+py*py)>99.0) return 0;

	//energy=pow(px*px+py*py+pz*pz+am*am,0.5);
	//inrap=((energy-pz)>1e-6)?((energy+pz)/(energy-pz)):((energy+pz)/1e-6);                      
	//theta=acos(pz/sqrt(px*px+py*py+pz*pz));
	//pseorap=-log(tan(theta/2.));
	//pt=sqrt(px*px+py*py);
	//phi=atan2(py,px);
	//if(fabs(pseorap)>5.0){return 0;}

	if(abs(id)==4)//选出charm quark存到数组
	{
		charm_mu = charm_mu + 1;		

		double r2_charm = x*x+y*y;   
		double r_phi_charm = atan2(y,x);//spatial coordinate MLCC!! (MLCC?)

		sum_r2sin2phi_charm = sum_r2sin2phi_charm + r2_charm*sin(2.*r_phi_charm);//计算各阶偏心率epsilon的分子和分母项
		sum_r2cos2phi_charm = sum_r2cos2phi_charm + r2_charm*cos(2.*r_phi_charm);		
		sum_r2sin3phi_charm = sum_r2sin3phi_charm + r2_charm*sin(3.*r_phi_charm);
		sum_r2cos3phi_charm = sum_r2cos3phi_charm + r2_charm*cos(3.*r_phi_charm);
		sum_r2sin4phi_charm = sum_r2sin4phi_charm + r2_charm*sin(4.*r_phi_charm);
		sum_r2cos4phi_charm = sum_r2cos4phi_charm + r2_charm*cos(4.*r_phi_charm);

		sum_r2_charm = sum_r2_charm+r2_charm;	
	}

	double r2 = x*x+y*y;   
	double r_phi = atan2(y,x);//spatial coordinate MLCC!! (MLCC?)

	r_phi_bad = r_phi;

	sum_r2sin2phi = sum_r2sin2phi + r2*sin(2.*r_phi);//计算各阶偏心率epsilon的分子和分母项
	sum_r2cos2phi = sum_r2cos2phi + r2*cos(2.*r_phi);		
	sum_r2sin3phi = sum_r2sin3phi + r2*sin(3.*r_phi);
	sum_r2cos3phi = sum_r2cos3phi + r2*cos(3.*r_phi);
	sum_r2sin4phi = sum_r2sin4phi + r2*sin(4.*r_phi);
	sum_r2cos4phi = sum_r2cos4phi + r2*cos(4.*r_phi);

	sum_r2 = sum_r2+r2;

	if(i>nmult2) return 1;

	int id2 = (int)track2->ID[i];
	double px2 = (double)track2->Px[i];
	double py2 = (double)track2->Py[i];
	double pz2 = (double)track2->Pz[i];
	double theta2=acos(pz2/sqrt(px2*px2+py2*py2+pz2*pz2));
	double pseorap2=-log(tan(theta2/2.));
	double pt2=sqrt(px2*px2+py2*py2);

	if(fabs(pseorap2)<eta_cut && (abs(id2)==421||abs(id2)==411||abs(id2)==413) && pt2<5.0 && i<nmult2)//选出D介子存到数组
	{
		Dmeson_mu = Dmeson_mu +1;
	}

	return 1;
}


void makeDmcorrelation(vector<TVector3> Dm,vector<TVector3> Hm,TH2D* dphi_Dm_pt)
{
  	if((int)(Dm.size())==0) return;

	NassDH = 0;

	for(int i=0;i<(int)(Dm.size());i++)
  	{
		for(int j=0;j<(int)(Hm.size());j++)
  		{
    		if(fabs(Dm[i].X()-Hm[j].X())>delta_eta_cut && Dm[i].Z()>pt_low_cut && Dm[i].Z()<pt_high_cut && Hm[j].Z()>pt_low_cut && Hm[j].Z()<pt_high_cut)
			{
				NassDH = NassDH + 1;
			}
		}
	}

	for(int i=0;i<(int)(Dm.size());i++)
  	{
		for(int j=0;j<(int)(Hm.size());j++)
  		{
    		if(fabs(Dm[i].X()-Hm[j].X())>delta_eta_cut && Dm[i].Z()>pt_low_cut && Dm[i].Z()<pt_high_cut && Hm[j].Z()>pt_low_cut && Hm[j].Z()<pt_high_cut)
			{
				if(fabs(Dm[i].Y()-Hm[j].Y())>PI)
				{
					dphi_Dm_pt->Fill(Dm[i].Z(),2.0*PI-fabs(Dm[i].Y()-Hm[j].Y()),1.0/NassDH); //why 1.0/Dm.size() ???					
				 
					VD2 = cos(2*(2.0*PI-fabs(Dm[i].Y()-Hm[j].Y())));
					VD3 = cos(3*(2.0*PI-fabs(Dm[i].Y()-Hm[j].Y())));

					D_b_v2_p_2c->Fill(centrality,VD2);
					D_b_v3_p_2c->Fill(centrality,VD3);

					VD2_2c_s = VD2_2c_s + VD2;
					VD3_2c_s = VD3_2c_s + VD3;

					ptd = Dm[i].Z();
					VDH2_pt->Fill(ptd,VD2);
					VDH3_pt->Fill(ptd,VD3);

					// Dv2_e_s->Fill(ptd,VD2);
					// Dv3_e_s->Fill(ptd,VD3);

					numds = numds + 1;

					double jyd2 = 0;

					if     ( (ptd>=0.0) && (ptd<0.5) ) {arrv2d[0] = arrv2d[0] + VD2; numv2d[0] = numv2d[0] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=0.5) && (ptd<1.0) ) {arrv2d[1] = arrv2d[1] + VD2; numv2d[1] = numv2d[1] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=1.0) && (ptd<1.5) ) {arrv2d[2] = arrv2d[2] + VD2; numv2d[2] = numv2d[2] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=1.5) && (ptd<2.0) ) {arrv2d[3] = arrv2d[3] + VD2; numv2d[3] = numv2d[3] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=2.0) && (ptd<2.5) ) {arrv2d[4] = arrv2d[4] + VD2; numv2d[4] = numv2d[4] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=2.5) && (ptd<3.0) ) {arrv2d[5] = arrv2d[5] + VD2; numv2d[5] = numv2d[5] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=3.0) && (ptd<3.5) ) {arrv2d[6] = arrv2d[6] + VD2; numv2d[6] = numv2d[6] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=3.5) && (ptd<4.0) ) {arrv2d[7] = arrv2d[7] + VD2; numv2d[7] = numv2d[7] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=4.0) && (ptd<4.5) ) {arrv2d[8] = arrv2d[8] + VD2; numv2d[8] = numv2d[8] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=4.5) && (ptd<5.0) ) {arrv2d[9] = arrv2d[9] + VD2; numv2d[9] = numv2d[9] + 1; jyd2 = VD2 / NassDH;}

					double dhad = 2.0*PI-fabs(Dm[i].Y()-Hm[j].Y());
					hadbf->Fill(dhad);

					// ndna->Fill(ptd);

					if (jyd2 != 0){
						Dv2_2c_ebe->Fill(ptd,jyd2);
					}						
				}
				else
				{
					dphi_Dm_pt->Fill(Dm[i].Z(),fabs(Dm[i].Y()-Hm[j].Y()),1.0/NassDH);

					VD2 = cos(2*(fabs(Dm[i].Y()-Hm[j].Y())));
					VD3 = cos(3*(fabs(Dm[i].Y()-Hm[j].Y())));
					
					VD2_2c_s = VD2_2c_s + VD2;
					VD3_2c_s = VD3_2c_s + VD3;

					D_b_v2_p_2c->Fill(centrality,VD2);
					D_b_v3_p_2c->Fill(centrality,VD3);

					ptd = Dm[i].Z();
					VDH2_pt->Fill(ptd,VD2);
					VDH3_pt->Fill(ptd,VD3);

					// Dv2_e_s->Fill(ptd,VD2);
					// Dv3_e_s->Fill(ptd,VD3);

					numds = numds + 1;

					double jyd2 = 0;

					if     ( (ptd>=0.0) && (ptd<0.5) ) { arrv2d[0] = arrv2d[0] + VD2; numv2d[0] = numv2d[0] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=0.5) && (ptd<1.0) ) { arrv2d[1] = arrv2d[1] + VD2; numv2d[1] = numv2d[1] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=1.0) && (ptd<1.5) ) { arrv2d[2] = arrv2d[2] + VD2; numv2d[2] = numv2d[2] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=1.5) && (ptd<2.0) ) { arrv2d[3] = arrv2d[3] + VD2; numv2d[3] = numv2d[3] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=2.0) && (ptd<2.5) ) { arrv2d[4] = arrv2d[4] + VD2; numv2d[4] = numv2d[4] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=2.5) && (ptd<3.0) ) { arrv2d[5] = arrv2d[5] + VD2; numv2d[5] = numv2d[5] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=3.0) && (ptd<3.5) ) { arrv2d[6] = arrv2d[6] + VD2; numv2d[6] = numv2d[6] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=3.5) && (ptd<4.0) ) { arrv2d[7] = arrv2d[7] + VD2; numv2d[7] = numv2d[7] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=4.0) && (ptd<4.5) ) { arrv2d[8] = arrv2d[8] + VD2; numv2d[8] = numv2d[8] + 1; jyd2 = VD2 / NassDH;}
					else if( (ptd>=4.5) && (ptd<5.0) ) { arrv2d[9] = arrv2d[9] + VD2; numv2d[9] = numv2d[9] + 1; jyd2 = VD2 / NassDH;}

					double dhad = fabs(Dm[i].Y()-Hm[j].Y());
					hadbf->Fill(dhad);

					// ndna->Fill(ptd);

					if (jyd2 != 0){
						Dv2_2c_ebe->Fill(ptd,jyd2);
					}		
				}
			}
  		}
  	}

  	return;
}


void makeHmcorrelation(vector<TVector3> Hm,TH2D* dphi_Hm_pt)
{
  	if((int)(Hm.size())==0) return;

	NassHH = 0;

	for(int i=0;i<(int)(Hm.size());i++)
  	{
		for(int j=0;j<(int)(Hm.size());j++)
		{
		    if(i!=j && fabs(Hm[i].X()-Hm[j].X())>delta_eta_cut && Hm[i].Z()>pt_low_cut && Hm[i].Z()<pt_high_cut && Hm[j].Z()>pt_low_cut && Hm[j].Z()<pt_high_cut)
			{
				NassHH = NassHH +1;
			}
		}
	}

	for(int i=0;i<(int)(Hm.size());i++)
  	{
		for(int j=0;j<(int)(Hm.size());j++)
		{
		    if(i!=j && fabs(Hm[i].X()-Hm[j].X())>delta_eta_cut && Hm[i].Z()>pt_low_cut && Hm[i].Z()<pt_high_cut && Hm[j].Z()>pt_low_cut && Hm[j].Z()<pt_high_cut)
			{
				if(fabs(Hm[i].Y()-Hm[j].Y())>PI)
				{
					dphi_Hm_pt->Fill(Hm[i].Z(),2.0*PI-fabs(Hm[i].Y()-Hm[j].Y()),1.0/NassHH);
				
					VH2 = cos(2*(2.0*PI-fabs(Hm[i].Y()-Hm[j].Y())));
					VH3 = cos(3*(2.0*PI-fabs(Hm[i].Y()-Hm[j].Y())));

					H_b_v2_p_2c->Fill(centrality,VH2);
					H_b_v3_p_2c->Fill(centrality,VH3);

					VH2_2c_s = VH2_2c_s + VH2;
					VH3_2c_s = VH3_2c_s + VH3;

					pth = Hm[i].Z();
					VHH2_pt->Fill(pth,VH2);
					VHH3_pt->Fill(pth,VH3);

					V2H_ps->Fill(1.1,VH2); 
					V3H_ps->Fill(1.1,VH3); 
					
					// Hv2_e_s->Fill(pth,VH2);
					// Hv3_e_s->Fill(pth,VH3);

					numhs = numhs + 1;

					double jyh2 = 0;

					if     ( (pth>=0.0) && (pth<0.5) ) { arrv2h[0] = arrv2h[0] + VH2; numv2h[0] = numv2h[0] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=0.5) && (pth<1.0) ) { arrv2h[1] = arrv2h[1] + VH2; numv2h[1] = numv2h[1] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=1.0) && (pth<1.5) ) { arrv2h[2] = arrv2h[2] + VH2; numv2h[2] = numv2h[2] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=1.5) && (pth<2.0) ) { arrv2h[3] = arrv2h[3] + VH2; numv2h[3] = numv2h[3] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=2.0) && (pth<2.5) ) { arrv2h[4] = arrv2h[4] + VH2; numv2h[4] = numv2h[4] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=2.5) && (pth<3.0) ) { arrv2h[5] = arrv2h[5] + VH2; numv2h[5] = numv2h[5] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=3.0) && (pth<3.5) ) { arrv2h[6] = arrv2h[6] + VH2; numv2h[6] = numv2h[6] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=3.5) && (pth<4.0) ) { arrv2h[7] = arrv2h[7] + VH2; numv2h[7] = numv2h[7] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=4.0) && (pth<4.5) ) { arrv2h[8] = arrv2h[8] + VH2; numv2h[8] = numv2h[8] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=4.5) && (pth<5.0) ) { arrv2h[9] = arrv2h[9] + VH2; numv2h[9] = numv2h[9] + 1; jyh2 = VH2 / NassHH;}

					double hhad = 2.0*PI-fabs(Hm[i].Y()-Hm[j].Y());
					hadaf->Fill(hhad);

					// nhna->Fill(pth);

					if (jyh2 != 0){
						Hv2_2c_ebe->Fill(ptd,jyh2);
					}	
				}
				else
				{
					dphi_Hm_pt->Fill(Hm[i].Z(),fabs(Hm[i].Y()-Hm[j].Y()),1.0/NassHH);
				
					VH2 = cos(2*(fabs(Hm[i].Y()-Hm[j].Y())));
					VH3 = cos(3*(fabs(Hm[i].Y()-Hm[j].Y())));

					H_b_v2_p_2c->Fill(centrality,VH2);
					H_b_v3_p_2c->Fill(centrality,VH3);

					VH2_2c_s = VH2_2c_s + VH2;
					VH3_2c_s = VH3_2c_s + VH3;

					pth = Hm[i].Z();
					VHH2_pt->Fill(pth,VH2);
					VHH3_pt->Fill(pth,VH3);

					V2H_ps->Fill(1.1,VH2); 
					V3H_ps->Fill(1.1,VH3); 

					// Hv2_e_s->Fill(pth,VH2);
					// Hv3_e_s->Fill(pth,VH3);

					numhs = numhs + 1;

					double jyh2 = 0;

					if     ( (pth>=0.0) && (pth<0.5) ) { arrv2h[0] = arrv2h[0] + VH2; numv2h[0] = numv2h[0] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=0.5) && (pth<1.0) ) { arrv2h[1] = arrv2h[1] + VH2; numv2h[1] = numv2h[1] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=1.0) && (pth<1.5) ) { arrv2h[2] = arrv2h[2] + VH2; numv2h[2] = numv2h[2] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=1.5) && (pth<2.0) ) { arrv2h[3] = arrv2h[3] + VH2; numv2h[3] = numv2h[3] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=2.0) && (pth<2.5) ) { arrv2h[4] = arrv2h[4] + VH2; numv2h[4] = numv2h[4] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=2.5) && (pth<3.0) ) { arrv2h[5] = arrv2h[5] + VH2; numv2h[5] = numv2h[5] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=3.0) && (pth<3.5) ) { arrv2h[6] = arrv2h[6] + VH2; numv2h[6] = numv2h[6] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=3.5) && (pth<4.0) ) { arrv2h[7] = arrv2h[7] + VH2; numv2h[7] = numv2h[7] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=4.0) && (pth<4.5) ) { arrv2h[8] = arrv2h[8] + VH2; numv2h[8] = numv2h[8] + 1; jyh2 = VH2 / NassHH;}
					else if( (pth>=4.5) && (pth<5.0) ) { arrv2h[9] = arrv2h[9] + VH2; numv2h[9] = numv2h[9] + 1; jyh2 = VH2 / NassHH;}

					double hhad = fabs(Hm[i].Y()-Hm[j].Y());
					hadaf->Fill(hhad);

					// nhna->Fill(pth);

					if (jyh2 != 0){
						Hv2_2c_ebe->Fill(ptd,jyh2);
					}				
				}
			}
		}
  	}

  	return;
}


void writeHistograms(char *outFile)//histogram write to file
{
	char file[256];
	sprintf(file,"%s",outFile);
	TFile *f = new TFile(file,"RECREATE");

	//b_ntrack->Write();

	hbf->Write();
	haf->Write();

	hv2pp->Write();
	hv2pd->Write();

	hv3pp->Write();
	hv3pd->Write();

	TPr2->Write();

	he2q2->Write();
	he3q3->Write();
	he4q4->Write();

	he2qpt->Write();
	he3qpt->Write();
	he4qpt->Write();

	De2v2_pp->Write();
	He2v2_pp->Write();
	De3v3_pp->Write();
	He3v3_pp->Write();

	De2v2_2c->Write();
	He2v2_2c->Write();
	De3v3_2c->Write();
	He3v3_2c->Write();

	hHad0->Write();  
	hHadmp->Write(); 
	hDad0->Write();  
	hDadmp->Write(); 

	hHadm3->Write(); 
	hDadm3->Write(); 

	hpsi2ad->Write();
	hpsi3ad->Write();

	hadbf->Write();
	hadaf->Write();

	// ndna->Write();
	// nhna->Write();

	VDH2_pt->Write();
	VHH2_pt->Write();
	VDH3_pt->Write();
	VHH3_pt->Write();

 	v2part_Dm_pt->Write();
 	v3part_Dm_pt->Write();
 	v4part_Dm_pt->Write();

	v2part_Hm_pt->Write();
 	v3part_Hm_pt->Write();
 	v4part_Hm_pt->Write();

	epsilon2_centrality->Write();
	epsilon3_centrality->Write();
	epsilon4_centrality->Write();

	v2EP_Dm_pt->Write();
 	v3EP_Dm_pt->Write();
 	v4EP_Dm_pt->Write();

	dphi_Dm_pt->Write();
	dphi_Hm_pt->Write();

	dphi_c_pt->Write();
	dphi_uds_pt->Write();

	spect_D0_pt->Write();

	H_b_v2e2->Write();
	H_b_v3e3->Write();
	D_b_v2e2->Write();
	D_b_v3e3->Write();

	H_Hb_v2e2->Write();
	H_Hb_v3e3->Write();
	D_Hb_v2e2->Write();
	D_Hb_v3e3->Write();

	// pp and 2c __ v and e TProfile 
	H_b_v2_pp->Write();
	H_b_v3_pp->Write();
	D_b_v2_pp->Write();
	D_b_v3_pp->Write();

	H_b_e2_pp->Write();
	H_b_e3_pp->Write();
	D_b_e2_pp->Write();
	D_b_e3_pp->Write();

	H_b_v2_2c->Write();
	H_b_v3_2c->Write();
	D_b_v2_2c->Write();
	D_b_v3_2c->Write();

	H_b_e2_2c->Write();
	H_b_e3_2c->Write();
	D_b_e2_2c->Write();
	D_b_e3_2c->Write();

	H_b_v2_p_2c->Write();
	H_b_v3_p_2c->Write();
	D_b_v2_p_2c->Write();
	D_b_v3_p_2c->Write();

	Dv2_2c_ebe->Write();
	Dv3_2c_ebe->Write();
	Hv2_2c_ebe->Write();
	Hv3_2c_ebe->Write();

	V2H_ebes->Write();
	V3H_ebes->Write();

	V2H_ps->Write();
	V3H_ps->Write();

	h_numd->Write();
	h_numh->Write();
	h_rati->Write();
	h_arrh->Write();
	h_arrd->Write();

	D_flow_nass->Write();
	D_f_n_test->Write();

	az_VDH2_pt->Write();
	az_VHH2_pt->Write();
	az_VDH3_pt->Write();
	az_VHH3_pt->Write();

	// bad event
	h_bad_mult->Write();
	h_bad_e2e3->Write();
	h_bad_dhmu->Write();

	h_bad_aft->Write();
	h_bad_zpc->Write();

	f->Write(); ////!!!!!!!!
	f->Close();
	delete f;
	return;

}


void deleteHistograms()
{	
	delete epsilon2_centrality;
	delete epsilon3_centrality;
	delete epsilon4_centrality;

 	delete v2part_Dm_pt;
 	delete v3part_Dm_pt;
 	delete v4part_Dm_pt;

	delete v2part_Hm_pt;
 	delete v3part_Hm_pt;
 	delete v4part_Hm_pt;

	delete v2EP_Dm_pt;
 	delete v3EP_Dm_pt;
 	delete v4EP_Dm_pt;

 	delete dphi_Dm_pt;
	delete dphi_Hm_pt;
	delete dphi_c_pt;
	delete dphi_uds_pt;
 	delete spect_D0_pt;

 	return;
}

bool event_loop(AMPT* track,int i)
{
	//in >> id >> px >> py >> pz >> am >> others1 >>others2 >>others3 >>others4;//input track(particle) information
	
	int id = (int)track->ID[i];
	//double index = (double)track->Indx[i];
	double px = (double)track->Px[i];
	double py = (double)track->Py[i];
	double pz = (double)track->Pz[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];
	double am = (double)track->Mass[i];
	TVector3 mP(px,py,pz);

	if(mP.Pt()<0.3 || mP.Pt()>99.0) return 0;      

	double theta=acos(pz/sqrt(px*px+py*py+pz*pz));
	double pseorap=-log(tan(theta/2.));
	double pt=sqrt(px*px+py*py);
	double phi=atan2(py,px); //momentum coordinate MLCC!!

	double r_phi = atan2(py,px);

	if(fabs(pseorap)<1.0 && (abs(id)==2212||abs(id)==321||abs(id)==211) && pt<5.0)//选出普通强子存到数组
	{
		sum_r2sin2phi_ep = sum_r2sin2phi_ep + pt*sin(2.*r_phi);
		sum_r2cos2phi_ep = sum_r2cos2phi_ep + pt*cos(2.*r_phi);
		sum_r2sin3phi_ep = sum_r2sin3phi_ep + pt*sin(3.*r_phi);
		sum_r2cos3phi_ep = sum_r2cos3phi_ep + pt*cos(3.*r_phi);
		sum_r2sin4phi_ep = sum_r2sin4phi_ep + pt*sin(4.*r_phi);
		sum_r2cos4phi_ep = sum_r2cos4phi_ep + pt*cos(4.*r_phi);

		sum_r2_ep = sum_r2_ep + pt;

		qx2 = qx2 + cos(2.*r_phi);
		qy2 = qy2 + sin(2.*r_phi);
		qx3 = qx3 + cos(3.*r_phi);
		qy3 = qy3 + sin(3.*r_phi);
		qx4 = qx4 + cos(4.*r_phi);
		qy4 = qy4 + sin(4.*r_phi);

		multiq = multiq + 1;	
	}

	// if(fabs(pseorap)<eta_cut && (abs(id)==421||abs(id)==411||abs(id)==413) && pt<5.0)//选出D介子存到数组
	// {
	// 	Dmeson_mu = Dmeson_mu +1;
	// }

  	return 1;
}

bool udsc_loop(AMPT* track,int i)
{
	//in >> id >> px >> py >> pz >> am >> others1 >>others2 >>others3 >>others4;//input track(particle) information
	
	int id = (int)track->ID[i];
	//double index = (double)track->Indx[i];
	double px = (double)track->Px[i];
	double py = (double)track->Py[i];
	double pz = (double)track->Pz[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];
	double am = (double)track->Mass[i];
	TVector3 mP(px,py,pz);

	if(mP.Pt()<0.3 || mP.Pt()>99.0) return 0;      

	double theta=acos(pz/sqrt(px*px+py*py+pz*pz));
	double pseorap=-log(tan(theta/2.));
	double pt=sqrt(px*px+py*py);
	double phi=atan2(py,px); //momentum coordinate MLCC!!

	if(fabs(pseorap)<eta_cut && (abs(id)==1||abs(id)==2||abs(id)==3) && pt<5.0)//选出普通强子存到数组
	{
		TVector3 uds_v(pseorap,phi,pt);//eta/phi/pt
		uds_quark.push_back(uds_v);
	}

	if(fabs(pseorap)<eta_cut && (abs(id)==4) && pt<5.0)//选出普通强子存到数组
	{
		TVector3 c_v(pseorap,phi,pt);//eta/phi/pt
		c_quark.push_back(c_v);
	}
}

void make_c_correlation(vector<TVector3> c_quark,vector<TVector3> uds_quark,TH2D* dphi_c_pt)
{
  	if((int)(c_quark.size())==0) return;

	for(int i=0;i<(int)(c_quark.size());i++)
  	{
		for(int j=0;j<(int)(uds_quark.size());j++)
		{
		    if(i!=j && fabs(c_quark[i].X()-uds_quark[j].X())>delta_eta_cut && c_quark[i].Z()>pt_low_cut && c_quark[i].Z()<pt_high_cut && uds_quark[j].Z()>pt_low_cut && uds_quark[j].Z()<pt_high_cut)
			{
				if(fabs(c_quark[i].Y()-uds_quark[j].Y())>PI)
				{
					dphi_c_pt->Fill(c_quark[i].Z(),2.0*PI-fabs(c_quark[i].Y()-uds_quark[j].Y()));
				
					az_VD2 = cos(2*(2.0*PI-fabs(c_quark[i].Y()-uds_quark[j].Y())));
					az_VD3 = cos(3*(2.0*PI-fabs(c_quark[i].Y()-uds_quark[j].Y())));

					az_ptd = c_quark[i].Z();

					// arrv2_c[]=arrv2_c[]+az_VD2;
					az_VDH2_pt->Fill(az_ptd,az_VD2);
					az_VDH3_pt->Fill(az_ptd,az_VD3);
				}
				else
				{
					dphi_c_pt->Fill(c_quark[i].Z(),fabs(c_quark[i].Y()-uds_quark[j].Y()));
				
					az_VD2 = cos(2*(fabs(c_quark[i].Y()-uds_quark[j].Y())));
					az_VD3 = cos(3*(fabs(c_quark[i].Y()-uds_quark[j].Y())));

					az_ptd = c_quark[i].Z();

					// arrv2_c[]=arrv2_c[]+az_VD2;
					az_VDH2_pt->Fill(az_ptd,az_VD2);
					az_VDH3_pt->Fill(az_ptd,az_VD3);
				}
			}
		}
  	}

  	return;
}


void make_uds_correlation(vector<TVector3> uds_quark,TH2D* dphi_uds_pt)
{
   	if((int)(uds_quark.size())==0) return;

	for(int i=0;i<(int)(uds_quark.size());i++)
  	{
		for(int j=0;j<(int)(uds_quark.size());j++)
		{
		    if(i!=j && fabs(uds_quark[i].X()-uds_quark[j].X())>delta_eta_cut && uds_quark[i].Z()>pt_low_cut && uds_quark[i].Z()<pt_high_cut && uds_quark[j].Z()>pt_low_cut && uds_quark[j].Z()<pt_high_cut)
			{
				if(fabs(uds_quark[i].Y()-uds_quark[j].Y())>PI)
				{
					dphi_uds_pt->Fill(uds_quark[i].Z(),2.0*PI-fabs(uds_quark[i].Y()-uds_quark[j].Y()));
				
					// az_VH2 = cos(2*(2.0*PI-fabs(uds_quark[i].Y()-uds_quark[j].Y())));
					// az_VH3 = cos(3*(2.0*PI-fabs(uds_quark[i].Y()-uds_quark[j].Y())));

					// az_pth = uds_quark[i].Z();
					// az_VHH2_pt->Fill(az_pth,az_VH2);
					// az_VHH3_pt->Fill(az_pth,az_VH3);
				}
				else
				{
					dphi_uds_pt->Fill(uds_quark[i].Z(),fabs(uds_quark[i].Y()-uds_quark[j].Y()));
				
					// az_VH2 = cos(2*(fabs(uds_quark[i].Y()-uds_quark[j].Y())));
					// az_VH3 = cos(3*(fabs(uds_quark[i].Y()-uds_quark[j].Y())));

					// az_pth = uds_quark[i].Z();
					// az_VHH2_pt->Fill(az_pth,az_VH2);
					// az_VHH3_pt->Fill(az_pth,az_VH3);
				}
			}
		}
  	}

  	return;
}


bool bad_event(AMPT* track, int i)
{
	int id = (int)track->ID[i];
	double x=(double)track->X[i];
	double y=(double)track->Y[i];
	double z=(double)track->Z[i];

	double r2 = x*x+y*y;   
	double r_phi = atan2(y,x);

	double rsin = r2*sin(2.*r_phi);
	double rcos = r2*cos(2.*r_phi);

	double delta_phi = abs(r_phi-0.785398);

	if(delta_phi<0.00001){
		cout<<"r2 is: "<<r2<<", r_phi is: "<<r_phi<<", rsin is: "<<rsin<<", rcos is: "<<rcos<<", id is: "<<id<<endl;
		cout<<"x is: "<<x<<", y is: "<<y<<", z is: "<<z<<endl;
	}

	r_bads = r_bads + r2;
	rsin_bads = rsin_bads + rsin;
	rcos_bads = rcos_bads + rcos;

	// cout<<"r_bads is: "<<r_bads<<", rsin_bads is:"<<rsin_bads<<", rcos_bads is: "<<rcos_bads<<endl;

	// h_bad_aft->Fill(px,py);dwdw
	// h_bad_zpc->Fill(x,y);

  	return 1;
}