#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TLegend.h"

#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include <TPad.h>
#include <TLine.h>
#include <TRandom.h>

#include "TUnfold.h"
#include "TUnfoldDensity.h"


using namespace std;

bool muMinus = false;
bool mlApplied = false;

bool var_p = true;

int version = 108;

double chisquare_p = 0;
double chisquare_phi = 0;

//theta
int nthetabins_gen = 40;
int nthetabins_reco = 80;
double theta_low = 0.;
double theta_high = 90.;

//phi
int nphibins_gen = 40;
int nphibins_reco = 80;
double phi_low = -180.;
double phi_high = 180.;

//Momentum
int npbins_gen_init=8;
int npbins_reco = 16;
double p_low = 0.8;
double p_high = 3.0;


double binwidth_gen = (p_high-p_low)/npbins_gen_init;
double binwidth_reco = (p_high-p_low)/npbins_reco;

int extra_bins = (p_low - 0.5)/binwidth_gen;
int npbins_gen = npbins_gen_init;//+2*extra_bins;

void DivideHistogramByBinWidth(TH1D *histogram) {
    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        double binContent = histogram->GetBinContent(i);
        double binWidth = histogram->GetBinWidth(i);
        histogram->SetBinContent(i, binContent / binWidth);
        histogram->SetBinError(i, histogram->GetBinError(i) / binWidth);
    }
}

//int main()
void unfolding(){

	string filename = "Mc_Data_Magnetic.root";

	char name[100];

	double p_low_temp = p_low;

	if(muMinus) {p_low = -1*p_high; p_high = -1*p_low_temp;}

	double pival = acos(-1.);
    const int numThetaRanges = 5;
    double thetaRanges[numThetaRanges][2] = {{0, 17}, {17, 26}, {26, 34}, {34, 44}, {44,90}};

	//double reco_bin_sch[23] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1};
	//double gen_bin_sch[13] = {0,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.1};

	//double reco_bin_sch[35] = {0, 1, 1.06149, 1.12298, 1.18446, 1.24595, 1.30744, 1.36893, 1.43041, 1.4919, 1.55339, 1.61488, 1.67636, 1.73785, 1.79934, 1.86083, 1.92231, 1.9838, 2.04529, 2.10678, 2.16826, 2.22975, 2.29124, 2.35273, 2.41421, 2.4757, 2.53719, 2.59868, 2.66016, 2.72165, 2.78314, 2.84463, 2.90611, 2.9676, 3.1};
	//double gen_bin_sch[17] = {0, 1, 1.14142, 1.28284, 1.42426, 1.56569, 1.70711, 1.84853, 1.98995, 2.13137, 2.27279, 2.41421, 2.55563, 2.69706, 2.83848, 2.9799, 3.1};

	//double reco_bin_sch[20] = {0, 1, 1.11314, 1.22627, 1.33941, 1.45255, 1.56569, 1.67882, 1.79196, 1.9051, 2.01823, 2.13137, 2.24451, 2.35765, 2.47078, 2.58392, 2.69706, 2.81019, 2.92333, 3.1};
	//double gen_bin_sch[10] = {0, 1, 1.28284, 1.56569, 1.84853, 2.13137, 2.41421, 2.69706, 2.9799, 3.1};

	//Variable bin width based on resolution MuPlus
	//double reco_bin_sch_muplus[18] = {0, 1, 1.11314, 1.22819, 1.3452, 1.4642, 1.58521, 1.70828, 1.83344, 1.96072, 2.09016, 2.2218, 2.35567, 2.49181, 2.63027, 2.77107, 2.91427, 3.1};
	//double gen_bin_sch_muplus[9] = {0, 1, 1.28284, 1.57769, 1.88504, 2.20543, 2.53941, 2.88757, 3.1};
	// double reco_bin_sch_muplus[47] = {0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.95, 0.98, 1.01, 1.04, 1.07, 1.1, 1.13, 1.16, 1.2, 1.23, 1.27, 1.3, 1.34, 1.38, 1.42, 1.46, 1.51, 1.55, 1.59, 1.64, 1.69, 1.74, 1.79, 1.84, 1.89, 1.95, 2.01, 2.06, 2.13, 2.19, 2.25, 2.32, 2.38, 2.45, 2.52, 2.6, 2.67, 2.75, 2.83, 2.92, 3};
	// double gen_bin_sch_muplus[23] = {0.8, 0.85, 0.9, 0.96, 1.02, 1.08, 1.15, 1.22, 1.29, 1.37, 1.46, 1.55, 1.65, 1.75, 1.86, 1.97, 2.09, 2.22, 2.36, 2.51, 2.66, 2.83, 3};
  	double reco_bin_sch_muplus[17] = {0.8, 0.87, 0.94, 1.02, 1.11, 1.21, 1.31, 1.43, 1.55, 1.68, 1.83, 1.98, 2.16, 2.34, 2.54, 2.76, 3};
	double gen_bin_sch_muplus[9] = {0.8, 0.94, 1.11, 1.31, 1.55, 1.83, 2.16, 2.54, 3};

	//Start from 0.6 - 3.5
	//0.7-0.85, 0.85-1.0

	//suppose 0.8-3.0 its 20 bins total bins =22, number of points =23
	//suppose 0.8-3.0 its 8 bins total bins =10 , number of points =11


	//Variable bin width based on resolution MuMinus
	//double reco_bin_sch_muminus[47] = {-3, -2.92, -2.83, -2.75, -2.67, -2.6, -2.52, -2.45, -2.38, -2.32, -2.25, -2.19, -2.13, -2.06, -2.01, -1.95, -1.89, -1.84, -1.79, -1.74, -1.69, -1.64, -1.59, -1.55, -1.51, -1.46, -1.42, -1.38, -1.34, -1.3, -1.27, -1.23, -1.2, -1.16, -1.13, -1.1, -1.07, -1.04, -1.01, -0.98, -0.95, -0.92, -0.9, -0.87, -0.85, -0.82, -0.8};
	//double gen_bin_sch_muminus[23] = {-3, -2.83, -2.66, -2.51, -2.36, -2.22, -2.09, -1.97, -1.86, -1.75, -1.65, -1.55, -1.46, -1.37, -1.29, -1.22, -1.15, -1.08, -1.02, -0.96, -0.9, -0.85, -0.8,};
	double reco_bin_sch_muminus[17] = {-3, -2.76, -2.54, -2.34, -2.16, -1.98, -1.83, -1.68, -1.55, -1.43, -1.31, -1.21, -1.11, -1.02, -0.94, -0.87, -0.8};
	double gen_bin_sch_muminus[9] = {-3, -2.54, -2.16, -1.83, -1.55, -1.31, -1.11, -0.94, -0.8};

	double reco_bin_sch[17];
	double gen_bin_sch[9];

	if(!muMinus) {
		std::copy(std::begin(reco_bin_sch_muplus), std::end(reco_bin_sch_muplus), std::begin(reco_bin_sch));
		std::copy(std::begin(gen_bin_sch_muplus), std::end(gen_bin_sch_muplus), std::begin(gen_bin_sch));
	} else {
		std::copy(std::begin(reco_bin_sch_muminus), std::end(reco_bin_sch_muminus), std::begin(reco_bin_sch));
		std::copy(std::begin(gen_bin_sch_muminus), std::end(gen_bin_sch_muminus), std::begin(gen_bin_sch));
	}

	//Theta
	TH1D *hist_th_reco = new TH1D("h_theta_mc_reco","",nthetabins_reco,theta_low,theta_high); hist_th_reco->Sumw2();
	TH1D *hist_th_data = new TH1D("h_theta_data_reco","",nthetabins_reco,theta_low,theta_high); hist_th_data->Sumw2();
	TH1D *hist_th_gen = new TH1D("h_theta_gen","",nthetabins_gen,theta_low,theta_high);  hist_th_gen->Sumw2();
	TH2D *mat_th_rm = new TH2D("h_theta_gen_reco","",nthetabins_reco,theta_low,theta_high, nthetabins_gen,theta_low,theta_high); mat_th_rm->Sumw2();

	//Phi
	TH1D *hist_phi_reco[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_mc_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_reco[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_reco[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_mc_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_reco[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_reco[j]->Sumw2();}
	}
	TH1D *hist_phi_data[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_data_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_data[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_data[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_data_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_data[j] = new TH1D(name,name,nphibins_reco, phi_low, phi_high); hist_phi_data[j]->Sumw2();}
	}
	TH1D *hist_phi_gen[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_gen_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_gen[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_gen[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_gen_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_gen[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_gen[j]->Sumw2();}
	}
	if(!muMinus) {sprintf(name,"h_phi_gen_reco_muplus_v%d",version);}
	if(muMinus)  {sprintf(name,"h_phi_gen_reco_muminus_v%d",version);}
	TH2D *mat_phi_rm[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_gen_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		mat_phi_rm[j] = new TH2D(name,name,nphibins_reco, phi_low, phi_high, nphibins_gen, phi_low, phi_high); mat_phi_rm[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_gen_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		mat_phi_rm[j] = new TH2D(name,name,nphibins_reco, phi_low, phi_high, nphibins_gen, phi_low, phi_high); mat_phi_rm[j]->Sumw2();}
	}

	//Momentum
	//TH1D *hist_p_reco = new TH1D(name,"",npbins_reco,reco_bin_sch); hist_p_reco->Sumw2();
	TH1D *hist_p_reco[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_p_mc_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_reco[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_reco[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_p_mc_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_reco[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_reco[j]->Sumw2();}
	}


	//TH1D *hist_p_data = new TH1D(name,"",npbins_reco,reco_bin_sch); hist_p_data->Sumw2();
	TH1D *hist_p_data[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_p_data_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_data[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_data[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_p_data_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_data[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_data[j]->Sumw2();}
	}


	//TH1D *hist_p_gen = new TH1D(name,"",npbins_gen,gen_bin_sch); hist_p_gen->Sumw2();
	TH1D *hist_p_gen[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_p_gen_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_gen[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_p_gen[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_p_gen_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_p_gen[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_p_gen[j]->Sumw2();}
	}


	if(!muMinus) {sprintf(name,"h_p_gen_reco_muplus_v%d",version);}
	if(muMinus)  {sprintf(name,"h_p_gen_reco_muminus_v%d",version);}
	//TH2D *mat_p_rm = new TH2D(name,"",npbins_reco,reco_bin_sch, npbins_gen,gen_bin_sch); mat_p_rm->Sumw2();
	TH2D *mat_p_rm[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_p_gen_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		mat_p_rm[j] = new TH2D(name,name,npbins_reco,reco_bin_sch, npbins_gen,gen_bin_sch); mat_p_rm[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_p_gen_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		mat_p_rm[j] = new TH2D(name,name,npbins_reco,reco_bin_sch, npbins_gen,gen_bin_sch); mat_p_rm[j]->Sumw2();}
	}


	TH1D *hist_data[numThetaRanges], *hist_reco[numThetaRanges], *hist_gen[numThetaRanges];
	TH2D *mat_rm[numThetaRanges];

	TH1D *fake_rate = new TH1D("fake_rate","",npbins_reco,reco_bin_sch);
	TH1D *efficiency = new TH1D("efficiency","", npbins_gen,gen_bin_sch);

  // hist_p_reco->Rebin(2);
  // hist_p_reco->Rebin(2);
  // hist_p_gen->Rebin(2);
  // mat_p_rm->Rebin(2,2);
  // fake_rate->Rebin(2);
  // efficiency->Rebin(2);

	int nbins_gen, nbins_reco;
	double lowedge, upedge;

	double thgen, phgen, pgen;
	double threco, phreco, preco;
	double chi2reco;
	int ndfreco;

	//In tree
	Float_t momin[1], thein[1], phiin[1];
	Float_t momrf[1], therf[1], phirf[1];
	float chisquare[1];
	UInt_t ndof[1];
  UInt_t fitfailed;
	Float_t learnedVal, learnedErr;

	Float_t momrfd[1], therfd[1], phirfd[1];
	float chisquared[1];
	UInt_t ndofd[1];
	UInt_t fitfailedd;
	Float_t learnedVald, learnedErrd;

	sprintf(name,"%s",filename.c_str());
	TFile *file1 = new TFile(name,"read");
	file1->cd();
	TTree *T1;
	T1 = (TTree*)file1->Get("T3");
	T1->SetBranchAddress("thein",thein);
	T1->SetBranchAddress("phiin",phiin);
	T1->SetBranchAddress("momin",momin);
	T1->SetBranchAddress("therf",therf);
	T1->SetBranchAddress("phirf",phirf);
	T1->SetBranchAddress("momrf",momrf);
	T1->SetBranchAddress("chisquare",chisquare);
	T1->SetBranchAddress("ndof",ndof);
	T1->SetBranchAddress("fitfailed",&fitfailed);
	T1->SetBranchAddress("learnedVal",&learnedVal);
	T1->SetBranchAddress("learnedErr",&learnedErr);

	TTree *T2;
	T2 = (TTree*)file1->Get("T5");
	T2->SetBranchAddress("therf",therfd);
	T2->SetBranchAddress("phirf",phirfd);
	T2->SetBranchAddress("momrf",momrfd);
	T2->SetBranchAddress("chisquare",chisquared);
	T2->SetBranchAddress("ndof",ndofd);
	T2->SetBranchAddress("fitfailed",&fitfailedd);
	T2->SetBranchAddress("learnedVal",&learnedVald);
	T2->SetBranchAddress("learnedErr",&learnedErrd);

	int totalevents= 0;


	


	for(int ientry=0; ientry<T1->GetEntries(); ientry++) {

		T1->GetEntry(ientry);

		//if(abs(momin[0])>3) {continue;}   To check to see it has any effect on assymettry

		if(!mlApplied && !muMinus) {learnedVal=1.;}
		if(!mlApplied && muMinus)  {learnedVal=0.;}

		thgen = thein[0]*180./pival; phgen = phiin[0]*180./pival; pgen = momin[0];
		threco = therf[0]*180./pival; phreco = phirf[0]*180./pival;
    	if(momrf[0]> 0){preco = -0.014765 + 1.60720*momrf[0];}//pol1 fit in range 0.6 to 2.0
    	//if(momrf[0]< 0){preco = -0.0398542 + 1.46176*momrf[0];}//pol1 fit in range -2.0 to -0.6
    	//if(momrf[0]< 0){preco = -0.0162093 + 1.59788*momrf[0];}//pol1 fit in range -2.0 to -0.6
    	if(momrf[0]< 0){preco = -0.185643 + 1.21464*momrf[0];}//pol1 fit in range -2.0 to -0.6
    	//if(momrf[0]< 0){preco = momrf[0]-0.2;}
		chi2reco = chisquare[0]; ndfreco = ndof[0];
		bool selected =  (fitfailed==1 && ndfreco>=5 && chi2reco/ndfreco<2)? true: false;

		if(fitfailed==1 || fitfailed !=1){
			//if(totalevents>1000000) break;

//////////////////MuPlus//////////////////////////////////////////////////////////////////////////////
			//if(pgen>p_low && pgen<p_high) {   //pgen should be outside
			if(!muMinus) {
				hist_th_gen->Fill(thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_gen[j]->Fill(pgen);
						hist_phi_gen[j]->Fill(phgen);
					}
				}
			}

			if(selected && learnedVal>0.9 && !muMinus) {
				hist_th_reco->Fill(threco);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_reco[j]->Fill(preco);
						hist_phi_reco[j]->Fill(phreco);
					}
				}
					totalevents++;
		  	} else if(!muMinus) {
				hist_th_reco->Fill(-100.);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_reco[j]->Fill(-100);
						hist_phi_reco[j]->Fill(-200.);
					}
				}
			}

			// if(selected && preco>p_low && preco<p_high && ientry<T1->GetEntries()/2.) {
			// 		hist_p_data->Fill(preco);
			// } else if(ientry<T1->GetEntries()/2.){hist_p_data->Fill(-100);}

			if(selected && learnedVal>0.9 && !muMinus) {
				mat_th_rm->Fill(threco,thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						mat_p_rm[j]->Fill(preco,pgen);                                    //Selected events
						mat_phi_rm[j]->Fill(phreco,phgen);
					}
				}
				
			} else if(!muMinus) {
				mat_th_rm->Fill(-100., thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						mat_p_rm[j]->Fill(-100.,pgen);                                    //Selected events
						mat_phi_rm[j]->Fill(-200.,phgen);
					}
				}
			}

//////////////MuMinus----------------------------------------------------------
			//Gen
			if(muMinus) {
				hist_th_gen->Fill(thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_gen[j]->Fill(pgen);
						hist_phi_gen[j]->Fill(phgen);
					}
				}
			}
			//Reco
			if(selected && learnedVal<0.1 && muMinus) {
				hist_th_reco->Fill(threco);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_reco[j]->Fill(preco);
						hist_phi_reco[j]->Fill(phreco);
					}
				}
					totalevents++;
		  	} else if(muMinus) {
				hist_th_reco->Fill(-100.);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_reco[j]->Fill(100);
						hist_phi_reco[j]->Fill(-200.);
					}
				}
			}
			//Response Matrix
			if(selected && learnedVal<0.1 && muMinus) {
				mat_th_rm->Fill(threco, thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						mat_p_rm[j]->Fill(preco,pgen);                                    //Selected events
						mat_phi_rm[j]->Fill(phreco,phgen);
					}
				}
				
			} else if(muMinus) {
				mat_th_rm->Fill(-100.,thgen);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						mat_p_rm[j]->Fill(100.,pgen);                                    //Selected events
						mat_phi_rm[j]->Fill(-200.,phgen);
					}
				}
			}

/*
			if(pgen>p_low && pgen<p_high && preco>p_low && preco<p_high && selected) {
				mat_p_rm->Fill(preco,pgen);
			}
			if((pgen<p_low || pgen>p_high) && preco>p_low && preco<p_high && selected) {   //Fakes
				mat_p_rm->Fill(preco,-100);
			}
			// if(pgen>p_high  && preco>p_low && preco<p_high) {
			// 	mat_p_rm->Fill(preco,900);
			// }
			if(pgen>p_low && pgen<p_high && (preco<p_low || preco>p_high || (!selected))) {   //Misses
				mat_p_rm->Fill(-100,pgen);
			}
*/
		}
	}

	//Fake Rate Calcualtion
	int nbinx = mat_p_rm[0]->GetNbinsX();
  	int nbiny = mat_p_rm[0]->GetNbinsY();

	static const int nbinmx = 120; //max(hResp->GetNbinsY(),hResp->GetNbinsX()) + 5;
	double totalgen[nbinmx]={0.};
  double totalreco[nbinmx]={0.};

	for(int ib=0; ib< nbinmx; ib++){
		totalgen[ib] = 0;
		totalreco[ib] = 0;
	}

	double fakerate[nbinmx];
	double effi[nbinmx];

	//TH2D* mat_p_rm_initial = (TH2D*)mat_p_rm->Clone();

	for (int ix=0; ix<nbinx+1; ix++) {
		for (int iy=0; iy<nbiny+1; iy++) {
			if(ix==0&&iy==0) continue;
			if(mat_p_rm[0]->GetBinContent(ix,iy) < -0.1)	mat_p_rm[0]->SetBinContent(ix,iy,0);
			totalreco[ix] += mat_p_rm[0]->GetBinContent(ix, iy);          // Total number of reco events in bin ix
			if (iy==0) fakerate[ix] = mat_p_rm[0]->GetBinContent(ix,iy);  //Number of muMinus events accepted in bin ix
			totalgen[iy] +=mat_p_rm[0]->GetBinContent(ix, iy);            //Total number of gen events in bin iy
			if (ix==0) effi[iy] = mat_p_rm[0]->GetBinContent(ix, iy) ;     //Number events not reconstructed but generated in bin iy   // Actually it should be ineffi
			if (ix==nbinx) effi[iy] += mat_p_rm[0]->GetBinContent(ix, iy) ;
			if (ix==0 || iy==0) {
				//mat_p_rm->SetBinContent(ix, iy, 0.0);
				//mat_p_rm->SetBinError(ix, iy, 0.0);
			}
		}//iy
	}//ix

	for (int iy=1; iy<nbiny+1; iy++) {
    	effi[iy] = (totalgen[iy] - effi[iy])/max(1.e-10, totalgen[iy]);      //Now inefficiency changed to efficiency
      if (abs(effi[iy]) > 1) effi[iy] = 1;
      else if ( effi[iy] < 0) effi[iy] = 0.0000001;
			efficiency->SetBinContent(iy, effi[iy]);
  	} //iy

  	for (int ix=1; ix<nbinx+1; ix++) {
		cout<<"jim, ix="<<ix<<"\tfake="<<fakerate[ix]<<"\ttotal="<<totalreco[ix]<<endl;//Actually fakes are too much for the same width bins
      fakerate[ix] = fakerate[ix] / max(1.e-10, totalreco[ix]);
      if(abs(fakerate[ix]) > 1) fakerate[ix] = 0.99999999;
      else if ( fakerate[ix] < 0) fakerate[ix] = 0.0000001;
			fake_rate->SetBinContent(ix, fakerate[ix]);
  	}//ixfakefake

	for(int ix=0; ix <((hist_p_reco[0]->GetNbinsX())); ix++){
		cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_p_reco[0]->GetBinContent(ix+1)<<endl;
		//hist_p_reco->SetBinContent(ix+1,(1-fakerate[ix+1])*(hist_p_reco->GetBinContent(ix+1)));
		//hist_p_reco->SetBinError(ix+1,(1-fakerate[ix+1])*(hist_p_reco->GetBinError(ix+1)));
	}
	//hist_p_reco->SetBinContent(0,0);
	//hist_p_reco->SetBinError(0,0);




/*
	for(int ix=0; ix<(hist_p_gen->GetNbinsX()); ix++){
		cout<<"gen "<<ix+1<<" = "<<effi[ix+1]<<" * "<<hist_p_gen->GetBinContent(ix+1)<<endl;
		hist_p_gen->SetBinContent(ix+1,(hist_p_gen->GetBinContent(ix+1))*effi[ix+1]);
		hist_p_gen->SetBinError(ix+1,(hist_p_gen->GetBinError(ix+1))*effi[ix+1]);
	}

	hist_p_gen->SetBinContent(0,0);
	hist_p_gen->SetBinError(0,0);
*/
	cout<<"MC TotalEvents Reco good"<<totalevents<<endl;

	totalevents=0;

	for(int ientry=0; ientry<T2->GetEntries(); ientry++){

		T2->GetEntry(ientry);

		if(!mlApplied && !muMinus) {learnedVald=1.;}
		if(!mlApplied && muMinus)  {learnedVald=0.;}

		threco = therfd[0]*180./pival; phreco = phirfd[0]*180./pival;
		if(momrfd[0]>0) {preco = -0.014765 + 1.60720*momrfd[0];}//pol1 fit in range 0.6 to 2.0
		//if(momrfd[0]< 0){preco = -0.0398542 + 1.46176*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    	//if(momrfd[0]< 0){preco = -0.0162093 + 1.59788*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    	if(momrfd[0]< 0){preco = -0.185643 + 1.21464*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    	//if(momrfd[0]< 0){preco = momrfd[0]-0.2;}
		chi2reco = chisquared[0]; ndfreco = ndofd[0];

		bool selected = (fitfailedd==1 && ndfreco>=5 && chi2reco/ndfreco<2)? true: false;

		if(fitfailedd==1 || fitfailedd!=1){
			//if(totalevents>1000000) break;
			if(selected && learnedVald>0.9 && !muMinus) {
				hist_th_data->Fill(threco);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_data[j]->Fill(preco);
						hist_phi_data[j]->Fill(phreco);
					}
				}
				totalevents++;
			} else if(!muMinus) {
				hist_th_data->Fill(-100.);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_data[j]->Fill(-100.);
						hist_phi_data[j]->Fill(-200.);
					}
				}
			}

			if(selected && learnedVald<0.1 && muMinus) {
				hist_th_data->Fill(threco);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_data[j]->Fill(preco);
						hist_phi_data[j]->Fill(phreco);
					}
				}
				totalevents++;
			} else if(muMinus) {
				hist_th_data->Fill(-100.);
				for (int j = 0; j < numThetaRanges; ++j) {
					if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
						hist_p_data[j]->Fill(100.);
						hist_phi_data[j]->Fill(-200.);
					}
				}
			}


		}
	}

	cout<<"Data TotalEvents Reco good"<<totalevents<<endl;

	for(int ix=0; ix <((hist_p_data[0]->GetNbinsX())); ix++){
		//cout<<"reco "<<ix+1<<" = "<<1-fakerate[ix+1]<<" * "<<hist_p_data->GetBinContent(ix+1)<<endl;
		//hist_p_data->SetBinContent(ix+1,(1-fakerate[ix+1])*(hist_p_data->GetBinContent(ix+1)));
		//hist_p_data->SetBinError(ix+1,(1-fakerate[ix+1])*(hist_p_data->GetBinError(ix+1)));
	}

	//hist_p_data->SetBinContent(0,0);
	//hist_p_data->SetBinError(0,0);

	hist_th_data->Scale(0.8*T1->GetEntries()/T2->GetEntries());
	for (int j = 0; j < numThetaRanges; ++j) {
		hist_p_data[j]->Scale(0.8*T1->GetEntries()/T2->GetEntries());
		hist_phi_data[j]->Scale(0.8*T1->GetEntries()/T2->GetEntries());
	}


	// hist_p_data = (TH1D*)hist_p_reco->Clone();
	// hist_p_data->SetName("h_p_data");

	if(var_p){
		for (int j = 0; j < numThetaRanges; ++j) {
			hist_data[j] = (TH1D*)hist_p_data[j]->Clone();
			hist_reco[j] = (TH1D*)hist_p_reco[j]->Clone();
			hist_gen[j] = (TH1D*)hist_p_gen[j]->Clone();
			mat_rm[j] = (TH2D*)mat_p_rm[j]->Clone();
		}
		nbins_gen = npbins_gen;
		nbins_reco = npbins_reco;
		lowedge = p_low;
		upedge =  p_high;
	}
	else{
		hist_data[0] = (TH1D*)hist_th_data->Clone();
		hist_reco[0] = (TH1D*)hist_th_reco->Clone();
		hist_gen[0] = (TH1D*)hist_th_gen->Clone();
		mat_rm[0] = (TH2D*)mat_th_rm->Clone();
		nbins_gen = nthetabins_gen;
		nbins_reco = nthetabins_reco;
		lowedge = theta_low;
		upedge =  theta_high;
	}

	TFile *fileout ;
	if(!muMinus)sprintf(name,"TUnfold_momentum_muplus_v%d.root",version);
	if(muMinus)sprintf(name,"TUnfold_momentum_muminus_v%d.root",version);
	fileout = new TFile(name,"recreate") ;

	//Theta
	TH1D *hist_th_unf_reco;
	if(!muMinus) {sprintf(name,"h_phi_unfolded_mc_muplus_v%d",version);}
	if(muMinus) {sprintf(name,"h_phi_unfolded_mc_muminus_v%d",version);}
	hist_th_unf_reco = new TH1D(name, name, nthetabins_gen, theta_low, theta_high);

	TH1D *hist_th_unf_data;
	if(!muMinus) {sprintf(name,"h_phi_unfolded_data_muplus_v%d",version);}
	if(muMinus) {sprintf(name,"h_phi_unfolded_data_muminus_v%d",version);}
	hist_th_unf_data = new TH1D(name, name, nthetabins_gen, theta_low, theta_high);

	//Phi
	TH1D *hist_phi_unf_reco[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_unfolded_mc_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_unf_reco[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_unf_reco[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_unfolded_mc_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_unf_reco[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_unf_reco[j]->Sumw2();}
	}

	TH1D *hist_phi_unf_data[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_phi_unfolded_data_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_unf_data[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_unf_data[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_phi_unfolded_data_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_phi_unf_data[j] = new TH1D(name,name,nphibins_gen, phi_low, phi_high); hist_phi_unf_data[j]->Sumw2();}
	}
	//Momentum
	TH1D *hist_unf_reco[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_unfolded_mc_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_unf_reco[j] = new TH1D(name,name,nbins_gen,gen_bin_sch); hist_unf_reco[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_unfolded_mc_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_unf_reco[j] = new TH1D(name,name,nbins_gen,gen_bin_sch); hist_unf_reco[j]->Sumw2();}
	}

	TH1D *hist_unf_data[numThetaRanges];
	for (int j = 0; j < numThetaRanges; ++j) {
		if(!muMinus) {sprintf(name,"h_unfolded_data_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_unf_data[j] = new TH1D(name,name,nbins_gen,gen_bin_sch); hist_unf_data[j]->Sumw2();}
		if(muMinus)  {sprintf(name,"h_unfolded_data_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
		hist_unf_data[j] = new TH1D(name,name,nbins_gen,gen_bin_sch); hist_unf_data[j]->Sumw2();}
	}

	// regularize curvature
	TUnfold::ERegMode regMode=TUnfold::kRegModeCurvature;
	// preserve the area
   	TUnfold::EConstraint constraintMode=TUnfold::kEConstraintArea;
	//TUnfold::EConstraint constraintMode=TUnfold::kEConstraintNone;
  	// bin content is divided by the bin width
  	TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;


	TUnfoldDensity unfoldBbBreco_th(mat_th_rm,TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	TUnfoldDensity unfoldBbBdata_th(mat_th_rm,TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
	unfoldBbBreco_th.SetInput(hist_th_reco, 1.0);
	unfoldBbBdata_th.SetInput(hist_th_reco, 1.0);	
	unfoldBbBreco_th.RegularizeBins(20,1,20, regMode);
	unfoldBbBdata_th.RegularizeBins(20,1,20, regMode);
	Double_t tauMin_th=pow(10,-1);
	Double_t tauMax_th=pow(10,0);
	Int_t nScan_th=50;
	Int_t iBest_th;
	TSpline *logTauX_th,*logTauY_th;
	TGraph *lCurvem_th;
	TGraph *lCurved_th;

	iBest_th = unfoldBbBreco_th.ScanLcurve(nScan_th,tauMin_th,tauMax_th,&lCurvem_th,&logTauX_th,&logTauY_th);
	Double_t t_th[1],x_th[1],y_th[1];
	logTauX_th->GetKnot(iBest_th,t_th[0],x_th[0]);
	logTauY_th->GetKnot(iBest_th,t_th[0],y_th[0]);
	hist_th_unf_reco = (TH1D*)unfoldBbBreco_th.GetOutput("Unfolded");

	iBest_th = unfoldBbBdata_th.ScanLcurve(nScan_th,tauMin_th,tauMax_th,&lCurvem_th,&logTauX_th,&logTauY_th);
	logTauX_th->GetKnot(iBest_th,t_th[0],x_th[0]);
	logTauY_th->GetKnot(iBest_th,t_th[0],y_th[0]);
	hist_th_unf_data = (TH1D*)unfoldBbBdata_th.GetOutput("Unfolded");


  	hist_th_data->SetLineColor(kBlack);
	hist_th_reco->SetLineColor(kBlue);
	hist_th_gen->SetLineColor(kRed);
	hist_th_unf_reco->SetLineColor(kMagenta);
	hist_th_unf_data->SetLineColor(kGreen);

	hist_th_reco->SetStats(0);
	hist_th_reco->SetStats(0);
	hist_th_gen->SetStats(1);
	hist_th_unf_reco->SetStats(0);
	hist_th_unf_data->SetStats(0);
	hist_th_gen->SetLineStyle(kDashed);

	for (int j = 0; j < numThetaRanges; ++j) {
		TUnfoldDensity unfoldBbBreco(mat_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
		TUnfoldDensity unfoldBbBdata(mat_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
		unfoldBbBreco.SetInput(hist_reco[j],1.0);
		unfoldBbBdata.SetInput(hist_data[j],1.0);	
		//Regulatisation
		Double_t estimatedPeakPosition=1.2;
		Int_t nPeek=3;
		Int_t iPeek=(Int_t)(npbins_gen*(estimatedPeakPosition-p_low)/(p_high-p_low)
                      // offset 1.5
                      // accounts for start bin 1
                      // and rounding errors +0.5
                      +1.5);

		// regularize output bins 1..iPeek-nPeek
		//unfoldBbB.RegularizeBins(iPeek+10,1,npbins_gen,regMode);
		unfoldBbBreco.RegularizeBins(3,1,9, regMode);
		unfoldBbBdata.RegularizeBins(3,1,9, regMode);

		// unfoldBbBreco.DoUnfold(pow(10,-1));
	 	// unfoldBbBdata.DoUnfold(pow(10,-1));

		// hist_unf_reco[j] = (TH1D*)unfoldBbBreco.GetOutput("Unfolded");
		// hist_unf_data[j] = (TH1D*)unfoldBbBdata.GetOutput("Unfolded");

		Double_t tauMin=pow(10,-1);
		Double_t tauMax=pow(10,3);
		Int_t nScan=50;
		Int_t iBest;
		TSpline *logTauX,*logTauY;
		TGraph *lCurvem;
		TGraph *lCurved;

		iBest = unfoldBbBreco.ScanLcurve(nScan,tauMin,tauMax,&lCurvem,&logTauX,&logTauY);
		std::cout<<"chi**2_reco="<<unfoldBbBreco.GetChi2A()<<"+"<<unfoldBbBreco.GetChi2L()<<" / "<<unfoldBbBreco.GetNdf()<<"\n";
		Double_t t[1],x[1],y[1];
		logTauX->GetKnot(iBest,t[0],x[0]);
		logTauY->GetKnot(iBest,t[0],y[0]);
		cout<<"tau = "<<t[0] <<"\tx = "<<x[0]<<"\ty = "<<y[0]<<endl;
		hist_unf_reco[j] = (TH1D*)unfoldBbBreco.GetOutput("Unfolded");

		iBest = unfoldBbBdata.ScanLcurve(nScan,tauMin,tauMax,&lCurvem,&logTauX,&logTauY);
		std::cout<<"chi**2_reco="<<unfoldBbBdata.GetChi2A()<<"+"<<unfoldBbBdata.GetChi2L()<<" / "<<unfoldBbBdata.GetNdf()<<"\n";
		logTauX->GetKnot(iBest,t[0],x[0]);
		logTauY->GetKnot(iBest,t[0],y[0]);
		cout<<"tau = "<<t[0] <<"\tx = "<<x[0]<<"\ty = "<<y[0]<<endl;
		hist_unf_data[j] = (TH1D*)unfoldBbBdata.GetOutput("Unfolded");

	
  		hist_data[j]->SetLineColor(kBlack);
		hist_reco[j]->SetLineColor(kBlue);
		hist_gen[j]->SetLineColor(kRed);
		hist_unf_reco[j]->SetLineColor(kMagenta);
		hist_unf_data[j]->SetLineColor(kGreen);

		hist_reco[j]->SetStats(0);
		hist_reco[j]->SetStats(0);
		hist_gen[j]->SetStats(1);
		hist_unf_reco[j]->SetStats(0);
		hist_unf_data[j]->SetStats(0);
		hist_gen[j]->SetLineStyle(kDashed);
  
		DivideHistogramByBinWidth(hist_data[j]);
		DivideHistogramByBinWidth(hist_reco[j]);
		DivideHistogramByBinWidth(hist_gen[j]);
		DivideHistogramByBinWidth(hist_unf_reco[j]);
		DivideHistogramByBinWidth(hist_unf_data[j]);

		unfoldBbBreco.Delete();
		unfoldBbBdata.Delete();

		//Phi
		TUnfoldDensity unfoldBbBreco_phi(mat_phi_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
		TUnfoldDensity unfoldBbBdata_phi(mat_phi_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
		unfoldBbBreco_phi.SetInput(hist_phi_reco[j],1.0);
		unfoldBbBdata_phi.SetInput(hist_phi_data[j],1.0);	
		//Regulatisation
		unfoldBbBreco_phi.RegularizeBins(20,1,39, regMode);
		unfoldBbBdata_phi.RegularizeBins(20,1,39, regMode);

		Double_t tauMin_phi=pow(10,-1);
		Double_t tauMax_phi=pow(10,0);
		Int_t nScan_phi=50;
		Int_t iBest_phi;
		TSpline *logTauX_phi,*logTauY_phi;
		TGraph *lCurvem_phi;
		TGraph *lCurved_phi;

		iBest_phi = unfoldBbBreco_phi.ScanLcurve(nScan_phi,tauMin_phi,tauMax_phi,&lCurvem_phi,&logTauX_phi,&logTauY_phi);
		std::cout<<"chi**2_reco="<<unfoldBbBreco_phi.GetChi2A()<<"+"<<unfoldBbBreco_phi.GetChi2L()<<" / "<<unfoldBbBreco_phi.GetNdf()<<"\n";
		Double_t t_phi[1],x_phi[1],y_phi[1];
		logTauX_phi->GetKnot(iBest_phi,t_phi[0],x_phi[0]);
		logTauY_phi->GetKnot(iBest_phi,t_phi[0],y_phi[0]);
		cout<<"tau = "<<t_phi[0] <<"\tx = "<<x_phi[0]<<"\ty = "<<y_phi[0]<<endl;
		hist_phi_unf_reco[j] = (TH1D*)unfoldBbBreco_phi.GetOutput("Unfolded");

		iBest_phi = unfoldBbBdata_phi.ScanLcurve(nScan_phi,tauMin_phi,tauMax_phi,&lCurvem_phi,&logTauX_phi,&logTauY_phi);
		std::cout<<"chi**2_reco="<<unfoldBbBdata_phi.GetChi2A()<<"+"<<unfoldBbBdata_phi.GetChi2L()<<" / "<<unfoldBbBdata_phi.GetNdf()<<"\n";
		logTauX_phi->GetKnot(iBest_phi,t_phi[0],x_phi[0]);
		logTauY_phi->GetKnot(iBest_phi,t_phi[0],y_phi[0]);
		cout<<"tau = "<<t_phi[0] <<"\tx = "<<x_phi[0]<<"\ty = "<<y_phi[0]<<endl;
		hist_phi_unf_data[j] = (TH1D*)unfoldBbBdata_phi.GetOutput("Unfolded");

	

	
  		hist_phi_data[j]->SetLineColor(kBlack);
		hist_phi_reco[j]->SetLineColor(kBlue);
		hist_phi_gen[j]->SetLineColor(kRed);
		hist_phi_unf_reco[j]->SetLineColor(kMagenta);
		hist_phi_unf_data[j]->SetLineColor(kGreen);

		hist_phi_reco[j]->SetStats(0);
		hist_phi_reco[j]->SetStats(0);
		hist_phi_gen[j]->SetStats(1);
		hist_phi_unf_reco[j]->SetStats(0);
		hist_phi_unf_data[j]->SetStats(0);
		hist_phi_gen[j]->SetLineStyle(kDashed);

		unfoldBbBreco_phi.Delete();
		unfoldBbBdata_phi.Delete();


	}

	if(muMinus){
		gStyle->SetStatX(0.3);
		gStyle->SetStatY(0.3);
		gStyle->SetStatW(0.25);
		gStyle->SetStatH(0.25);
	}
	//Momentum
	TCanvas *c1 = new TCanvas("Tunfold_output_mc_reco_momentum","Tunfold_output_mc_reco_momentum",800,1200);
	c1->Divide(2,3);

	for (int j = 0; j < numThetaRanges; ++j) {
		TLegend *leg1 = new TLegend(0.4,0.2,0.6,0.4);
		leg1->SetBorderSize(0);
		c1->cd(j+1);

		//gPad->SetLogx(1);
		gPad->SetLogy(1);

		hist_gen[j]->GetYaxis()->SetRangeUser(1000,100000);
		//hist_gen[j]->SetTitle("MC Reco Unfolding Without ML");

		hist_gen[j]->Draw("hist");  leg1->AddEntry(hist_gen[j],"GEN","l");
		hist_data[j]->Draw("hist:same");  leg1->AddEntry(hist_data[j],"Data Reco","l");
		hist_reco[j]->Draw("hist:sames"); leg1->AddEntry(hist_reco[j],"MC Reco","l");
		hist_unf_reco[j]->Draw("sames"); leg1->AddEntry(hist_unf_reco[j],"Unfolded","l");

		leg1->Draw();
		
	}

	sprintf(name,"%s.png",c1->GetName());
	c1->SaveAs(name);

	TCanvas *c2 = new TCanvas("Tunfold_output_data_reco_momentum","Tunfold_output_data_reco_momentum",800,1200);
	c2->Divide(2,3);
	for (int j = 0; j < numThetaRanges; ++j) {
		TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
		leg2->SetBorderSize(0);
		c2->cd(j+1);
		gPad->SetLogy(1);

		hist_gen[j]->GetYaxis()->SetRangeUser(1000,100000);
		//hist_gen[j]->SetTitle("Data Reco Unfolding Without ML");

		hist_gen[j]->Draw("hist");  leg2->AddEntry(hist_gen[j],"GEN","l");
		hist_data[j]->Draw("hist:same");  leg2->AddEntry(hist_data[j],"Data Reco","l");
		hist_reco[j]->Draw("hist:sames"); leg2->AddEntry(hist_reco[j],"MC Reco","l");
		hist_unf_data[j]->Draw("sames"); leg2->AddEntry(hist_unf_data[j],"Unfolded","l");

		leg2->Draw();
	}

	sprintf(name,"%s.png",c2->GetName());
	c2->SaveAs(name);


	//Phi
	TCanvas *c3 = new TCanvas("Tunfold_output_mc_reco_phi","Tunfold_output_mc_reco_phi",800,1200);
	c3->Divide(2,3);

	for (int j = 0; j < numThetaRanges; ++j) {
		TLegend *leg1 = new TLegend(0.4,0.2,0.6,0.4);
		leg1->SetBorderSize(0);
		c3->cd(j+1);

		//gPad->SetLogx(1);
		gPad->SetLogy(0);

		hist_phi_gen[j]->GetYaxis()->SetRangeUser(1000,20000);
		//hist_phi_gen[j]->SetTitle("MC Reco Unfolding Without ML");

		hist_phi_gen[j]->Draw("hist");  leg1->AddEntry(hist_phi_gen[j],"GEN","l");
		hist_phi_data[j]->Draw("hist:same");  leg1->AddEntry(hist_phi_data[j],"Data Reco","l");
		hist_phi_reco[j]->Draw("hist:sames"); leg1->AddEntry(hist_phi_reco[j],"MC Reco","l");
		hist_phi_unf_reco[j]->Draw("sames"); leg1->AddEntry(hist_phi_unf_reco[j],"Unfolded","l");

		leg1->Draw();
		
	}

	sprintf(name,"%s.png",c3->GetName());
	c3->SaveAs(name);

	TCanvas *c4 = new TCanvas("Tunfold_output_data_reco_phi","Tunfold_output_data_reco_phi",800,1200);
	c4->Divide(2,3);
	for (int j = 0; j < numThetaRanges; ++j) {
		TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
		leg2->SetBorderSize(0);
		c4->cd(j+1);
		gPad->SetLogy(0);

		hist_phi_gen[j]->GetYaxis()->SetRangeUser(1000,20000);
		//hist_phi_gen[j]->SetTitle("Data Reco Unfolding Without ML");

		hist_phi_gen[j]->Draw("hist");  leg2->AddEntry(hist_phi_gen[j],"GEN","l");
		hist_phi_data[j]->Draw("hist:same");  leg2->AddEntry(hist_phi_data[j],"Data Reco","l");
		hist_phi_reco[j]->Draw("hist:sames"); leg2->AddEntry(hist_phi_reco[j],"MC Reco","l");
		hist_phi_unf_data[j]->Draw("sames"); leg2->AddEntry(hist_phi_unf_data[j],"Unfolded","l");

		leg2->Draw();
	}

	sprintf(name,"%s.png",c4->GetName());
	c4->SaveAs(name);

	//Theta

	//Phi
	TCanvas *c5 = new TCanvas("Tunfold_output_mc_reco_theta","Tunfold_output_mc_reco_theta",800,1200);
	for (int j = 0; j < 1; ++j) {
		TLegend *leg1 = new TLegend(0.4,0.2,0.6,0.4);
		leg1->SetBorderSize(0);
		c5->cd(j+1);

		//gPad->SetLogx(1);
		gPad->SetLogy(0);

		//hist_th_gen[j]->GetYaxis()->SetRangeUser(1000,100000);
		//hist_th_gen[j]->SetTitle("MC Reco Unfolding Without ML");

		hist_th_gen->Draw("hist");  leg1->AddEntry(hist_th_gen,"GEN","l");
		hist_th_data->Draw("hist:same");  leg1->AddEntry(hist_th_data,"Data Reco","l");
		hist_th_reco->Draw("hist:sames"); leg1->AddEntry(hist_th_reco,"MC Reco","l");
		hist_th_unf_reco->Draw("sames"); leg1->AddEntry(hist_th_unf_reco,"Unfolded","l");

		leg1->Draw();
		
	}

	sprintf(name,"%s.png",c5->GetName());
	c5->SaveAs(name);

	TCanvas *c6 = new TCanvas("Tunfold_output_data_reco_theta","Tunfold_output_data_reco_theta",800,1200);
	for (int j = 0; j < 1; ++j) {
		TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
		leg2->SetBorderSize(0);
		c6->cd(j+1);
		gPad->SetLogy(0);

		//hist_th_gen[j]->GetYaxis()->SetRangeUser(1000,100000);
		//hist_th_gen[j]->SetTitle("Data Reco Unfolding Without ML");

		hist_th_gen->Draw("hist");  leg2->AddEntry(hist_th_gen,"GEN","l");
		hist_th_data->Draw("hist:same");  leg2->AddEntry(hist_th_data,"Data Reco","l");
		hist_th_reco->Draw("hist:sames"); leg2->AddEntry(hist_th_reco,"MC Reco","l");
		hist_th_unf_data->Draw("sames"); leg2->AddEntry(hist_th_unf_data,"Unfolded","l");

		leg2->Draw();
	}

	sprintf(name,"%s.png",c6->GetName());
	c6->SaveAs(name);


	fileout->cd();
	
	hist_th_gen->Write();
	hist_th_reco->Write();
	hist_th_data->Write();
	hist_th_unf_reco->Write();
	hist_th_unf_data->Write();
	mat_th_rm->Write();

	for (int j = 0; j < numThetaRanges; ++j) {
		hist_gen[j]->Write();
		hist_reco[j]->Write();
		hist_data[j]->Write();
		hist_unf_reco[j]->Write();
		hist_unf_data[j]->Write();
		mat_rm[j]->Write();

		hist_phi_gen[j]->Write();
		hist_phi_reco[j]->Write();
		hist_phi_data[j]->Write();
		hist_phi_unf_reco[j]->Write();
		hist_phi_unf_data[j]->Write();
		mat_phi_rm[j]->Write();
	}
	fileout->Write();
	fileout->Close();


}
