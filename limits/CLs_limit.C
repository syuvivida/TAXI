#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <iostream>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TLeaf.h>
#include "compute.C"

using namespace std;

const double kB_constant=1.38e-23; //1.380649e-23; // in J*K-1
const double noiseT =3 ; // in K
const double binSize = 1000; // in Hz
const double nspec = 52443000;
const double mean_noise_exp = kB_constant*noiseT*binSize;
const double sigmaN_exp= mean_noise_exp/sqrt(nspec);

// for 95% CL limits, TF1::GetX(0.05), for Gaussian distribution
// TF1* f1 = new TF1("f1", CLs,-10,10,1);
// f1->SetParameter(0,xxx) where xxx is (observed power - noise_mean)/sigma
// f1->GetX(0.05) gives the 95% CL upper limit on signal power in terms of sigma
Double_t CLs_Gaus(Double_t* xx, Double_t* par){

  Double_t mu_s=xx[0]; // mu_s, now set mu_b=0; sigma_b=1;
  Double_t obs =par[0]; // in terms of sigma_b, can range from -inf to inf

  Double_t CLspb = (1+TMath::Erf((obs-mu_s)/sqrt(2)))*0.5;
  Double_t CLb = (1+TMath::Erf((obs)/sqrt(2)))*0.5;
  Double_t CLs = CLb>1e-16? CLspb/CLb:-1;
  return CLs;
}

void CLs_limit(std::string inputFile, std::string treeName="data")
{



  TChain* thisTree = new TChain(treeName.data());

  thisTree->Add(inputFile.data());


  Long64_t nentries = (Long64_t)thisTree->GetEntries();
  TH1F* hlimit = new TH1F("hlimit","95% CLs limit on signal power",200,-10,10);
  TH1F* hobs_sigma = new TH1F("hobs_sigma","Observed power difference from noise in terms of #sigma_{N}",100,-5,5);
  TF1* f1 = new TF1("f1",CLs_Gaus,-10,10,1);
  f1->SetNpx(5000);
  const unsigned int nBins = nentries;
  double x[nBins], y[nBins];
  double y_gaAA[nBins];
  double y_ratio[nBins];
  cout << "Expected mean_noise = " <<  mean_noise_exp << endl;
  cout << "Expected #sigma_N = " << sigmaN_exp << endl;

  double mean_obs = 0;
  
  // loop over event entries
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // this line is very important!!
    thisTree->GetEntry(jentry);
    double frequency=-9999, power_spectrum=-9999;
    frequency      = thisTree->GetBranch("frequency")->GetLeaf("frequency")->GetValue();
    power_spectrum = thisTree->GetBranch("power_spectrum")->GetLeaf("power_spectrum")->GetValue();

    mean_obs += power_spectrum;
    
    double obs_diffFromNoise_sigma = (power_spectrum - mean_noise_exp)/sigmaN_exp;
    f1->SetParameter(0, obs_diffFromNoise_sigma);
    double limit =f1->GetX(0.05);
    x[jentry] = frequency/1e9;
    y[jentry] = limit; // in units of sigmaN
    hlimit->Fill(limit);
    hobs_sigma->Fill(obs_diffFromNoise_sigma);

    double limit_on_gaAA =  computeLimit(limit,frequency);
    y_gaAA[jentry] = limit_on_gaAA;

    double limit_KSVZ = computegAgg(frequency);
    double ratio = limit_on_gaAA/limit_KSVZ;
    y_ratio[jentry] = ratio;

  } // end of loop over entries

  mean_obs /= nentries;
  cout << "Measured mean_noise = " <<  mean_obs << endl;

  hlimit->Draw();
  hobs_sigma->Draw();
  
  TGraph* g_limit = new TGraph(nBins, x,y);
  g_limit->GetXaxis()->SetTitle("Frequency [GHz]");
  g_limit->GetYaxis()->SetTitle("95% CLS limit on signal power in terms of #sigma_{N}");
  g_limit->Draw("AC*");


  TGraph* gaAA_limit = new TGraph(nBins, x,y_gaAA);
  gaAA_limit->GetXaxis()->SetTitle("Frequency [GHz]");
  gaAA_limit->GetYaxis()->SetTitle("95% CLS limit on g_{a#gamma#gamma} [GeV^{-1}]");
  gaAA_limit->Draw("AC*");

  TGraph* g_ratio = new TGraph(nBins, x,y_ratio);
  g_ratio->GetXaxis()->SetTitle("Frequency [GHz]");
  g_ratio->GetYaxis()->SetTitle("95% CLS limit on g_{a#gamma#gamma}/ g_{a#gamma#gamma}^{KSVZ}");
  g_ratio->Draw("AC*");

  
}
