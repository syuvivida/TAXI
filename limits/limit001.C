#include <TRandom3.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <iostream>
#include <string>
#include <fstream>

// a program to compute upper limits on signal power considering a very narrow
// axion width (less than the binwidth 125 Hz)

using namespace std;

// Pout = 5.7e-27 W [ (g_gamma/0.97)^2 (rho_a/0.45) ] [B^2 V (f/1 GHz)(Q_L/10000)(Cmnp/0.5)(1-2S11)/(1-S11)]
// corr Pout =
// 5.7e-27 W [ (g_gamma/0.97)^2 (rho_a/0.45) * (axion_width/digitizer resolution)*rho] [B^2 V (f/5 GHz)(Q_L/2000)(Cmnp/0.5)(1-2S11)/(1-S11)*1/(1+4Q^2(f/f0-1)^2)]

const Double_t QL=2000;  // 1/QL = 1/Q_antenna + 1/Q_cavity
const Double_t B=8.0;
const Double_t Cmnp=0.5; // form factor
const Double_t S11=0;
const Double_t g_gamma=0.97;
const Double_t eta=1;
const Double_t scale=1; // due to the difference of axion width and digitizer resolution
const Double_t baseline=1.14; // (*1e-26)
const Double_t Tsys=4; // in Kelvin
const Double_t bandwidth=5000; // 5 kHz
const Double_t kB=1.38e3; // ( *1e-26)

const Double_t N95=1.64485; // corresponds to 95% confidence level


#ifndef FATAL
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)
#endif


Double_t P_exp(Double_t* x, Double_t* par)
{
  Double_t f_vary = x[0];
  Double_t f0 = par[0];
  Double_t Power =  baseline * (f0/5.) * (QL/2000.) * B * B * (Cmnp/0.5) * (1-2*S11)/(1-S11) *
    (1/(1+4*QL*QL*pow(f_vary/f0-1,2))) * pow(g_gamma/0.97,2) * eta * scale;

  return Power;

}


Double_t Int_G(Double_t* x, Double_t* par)
{
  Double_t xx = x[0]; // excluded power in terms of noise sigma
  Double_t meas_sig = par[0]; // measured power in terms of noise sigma
  Double_t integral = 0.5 + 0.5*TMath::Erf( (xx-meas_sig)/sqrt(2));
  return integral;

}


// frequency discussed here is in GHz
void limit001(const unsigned long int nTrials = 1000, const Double_t freq_lo = 4.0, const Double_t freq_hi = 6.0, const Double_t f0 = 5.0)
{

  if(S11 > 1 || fabs(S11-1)<1e-6 || S11 < 0) FATAL("S11 should be between 0 and 1");

  TF1* fsig = new TF1("fsig", P_exp, freq_lo, freq_hi, 1);
  fsig->SetNpx(2500);
  fsig->SetParameter(0,f0);

  
  TF1* fint = new TF1("fint", Int_G,-10,10,1);
  fint->SetNpx(2500);
  fint->SetParameter(0,0);

  
  const unsigned int nBins = (freq_hi-freq_lo)*1e9/bandwidth;

  TH1D* htemp = new TH1D("htemp","", 100, -5,5);
  htemp->SetYTitle("Number of trials");

  TH1D* hnoise = (TH1D*)htemp->Clone("hnoise");
  hnoise->SetTitle("Measured power in terms of #sigma_{noise}");

  TProfile* ptemp = new TProfile("ptemp", "",nBins, freq_lo, freq_hi);
  ptemp->SetXTitle("Frequency [Hz]");
  ptemp->SetYTitle(Form("Average Power per %d Hz bin",(int)bandwidth));

  // assuming no signal is observed
  TProfile* pMeas = (TProfile*)ptemp->Clone("pMeas");
  pMeas->SetYTitle("Average measured power / #sigma_{noise}");

  TProfile* pLimit= (TProfile*)ptemp->Clone("pLimit");
  pLimit->SetYTitle("95% Upper limit on signal power / #sigma_{noise}");

  TProfile* pCoupling= (TProfile*)ptemp->Clone("pCoupling");
  pCoupling->SetTitle("95% Upper limit on g_{a#gamma#gamma} vs. frequency");
  pCoupling->SetYTitle(" 10^{-12} GeV^{-1} ");

  
  
  

  TRandom3* gRandom= new TRandom3();
  Double_t sigma_noise = kB*Tsys*bandwidth;

  for(unsigned long int it=0; it < nTrials; it++)
    {

      Double_t f_any = freq_lo + (freq_hi-freq_lo)*gRandom->Rndm();
      

      Double_t mea_sig = gRandom->Gaus(0,1);
      Double_t p_measure = mea_sig*sigma_noise;
      // fint -> SetParameter(0, mea_sig);
      // Double_t upper_sig = fint->GetX(0.95);

      // //      cout << "measured significance = " << mea_sig << "\t upper_sig = " << upper_sig << endl;
      // cout << "Diff = " << upper_sig - mea_sig << endl;
      
      // Double_t p_limit = upper_sig*sigma_noise;

      Double_t upper_sig = (mea_sig + N95);
      
      Double_t p_limit = upper_sig*sigma_noise;


      hnoise->Fill(mea_sig);

      pMeas->Fill(f_any, mea_sig);
      pLimit->Fill(f_any, upper_sig);

      if(p_limit<0)continue;
      Double_t p_coupling = 5*sqrt(p_limit/fsig->Eval(f_any));
      pCoupling->Fill(f_any, p_coupling);


    } // end of loop over nTrials
  
  

  // writing example output file
  TFile* outFile = new TFile("test_taxi_limit001.root","recreate");
  fsig->Write();
  fint->Write();
  hnoise->Write();
  pMeas->Write();
  pLimit->Write();
  pCoupling->Write();
  outFile->Close();

  
} // end of program









