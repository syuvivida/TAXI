#include <TCanvas.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

// a program to search for axion using numbers of Ed Daw 1807.09369
// assuming a maxwell-boltzmann distribution for signal by default
// If the "narrow" option is turned on, the signal power is deposited in one bin

using namespace std;

// Pout_original = 3.0e-22 W [ V/220L * B/7.6T * B/7.6T* Cnlm *
//                 g_gamma/0.97*g_gamma/0.97 * rhoa /0.45 * f0/750MHz * Q/70000]
// corr Pout =
// Pout_original * (1-2S11)/(1-S11)/(1+4Q^2(f/f0-1)^2)] *eta

const double c = 3e5; // in km/s
const double eta = 1;
const double baseline = 3.0; // (*1e-22)
const double QL = 70000;
const double S11 = 0;
const double Tsys = 5.6; // in Kelvin
const double bandwidth = 125e-6; //  in MHz, 125 Hz
const double kB = 1.38e-1; // ( *1e-22)
//const double integration_time = 80; // in seconds
const int N_integral = 10000;
const double sigma_noise = kB*Tsys*bandwidth*1e6/sqrt(N_integral); // for each bin
const double rangeSpec = 50e-3; // 50 kHz
const int nBins = rangeSpec/bandwidth;
const double step_size = 2e-3; // in MHz, 2 kHz

#ifndef FATAL
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)
#endif


// frequency in MHz
double P_exp(double* x, double* par)
{
  double f_vary = x[0];
  double f0 = par[0];
  double Power =  baseline * (f0/750.) * (1-2*S11)/(1-S11) *
    (1/(1+4*QL*QL*pow(f_vary/f0-1,2))) * eta;
  return Power;
}


double maxwell_boltzmann(double x)
{
  double v = 226;
  double func = 4*TMath::Pi()*x*x/pow(v*v*TMath::Pi(),1.5)*TMath::Exp(-x*x/v/v);
  return func;
}


double P_exp_boltzmann(double* x, double* par)
{
  double f_vary = x[0];
  double f0 = par[0]; // resonance frequency
  double faxion = par[1];
  double deltaf = fabs(f_vary-faxion);
  double v = c*sqrt(2*deltaf/faxion);
  
  // the factor c*c/v/faxion comes from the transformation of variables
  
  double Power =  baseline * (f0/750.) * (1-2*S11)/(1-S11) *
    (1/(1+4*QL*QL*pow(f_vary/f0-1,2))) * eta *
    maxwell_boltzmann(v)*(c*c/v/faxion);
  return Power;
}


void generate_background(TH1F* hbkg, const double sigma)
{
  const unsigned int nbinsx = hbkg->GetNbinsX();
  for(unsigned int ib=1; ib <= nbinsx; ib++){
    double noise = gRandom->Gaus(0, sigma);
    hbkg->SetBinContent(ib, noise);
  }
  return;
}

void generate_signal(bool narrow, bool debug, TH1F* hsig,
		     TF1* fsig, TF1* fsig_wide,
		     const double f_axion,
		     const double resFreq, const double end_freq)
{
  // if axion signal frequency is higher than the highest frequency of
  // this spectrum, skip to the next resonance frequency
  
  if(f_axion > end_freq)return;

  fsig->SetParameter(0,resFreq);

  fsig_wide->SetParameter(0,resFreq);
  fsig_wide->SetParameter(1,f_axion);
  
  // first generate signal
  // find the first bin of axion signal in this spectrum
  // binStart could be equal to zero
  const unsigned int binStart = hsig->FindBin(f_axion);
  if(debug) cout << "binStart = " << binStart << endl;

  // if assuming narrow distribution
  if(narrow){
    // signal power for a given axion frequency and cavity resonance frequency
    double power_sig = fsig->Eval(f_axion); 
    hsig->SetBinContent(binStart,power_sig);
    // if(debug) cout << "power signal " <<  power_sig << endl;      
  }
  // if assuming a maxwell distribution
  else{
    
    double sum_signal_power = 0;
    
    for(unsigned int ib=binStart; ib <= nBins; ib++){

      double binlo = (ib==binStart)? f_axion: hsig->GetBinLowEdge(ib);
      double binhi = hsig->GetBinLowEdge(ib+1);
      // if(debug)
      // 	cout << "binlo = " << setprecision(10) << binlo <<
      // 	  "\t binhi = " << setprecision(10) << binhi << endl;
    	
      // set signal power for every bin, assuming maxwell distribution
      double power_for_this_bin = fsig_wide->Integral(binlo,binhi);
      hsig->SetBinContent(ib, power_for_this_bin);
      sum_signal_power += power_for_this_bin;	

      // if(debug) cout << "power for this bin = " << power_for_this_bin  << endl;
    	
    } // end of loop over nBins

    if(debug) cout << "sum of signal power is " << sum_signal_power << endl;
  } // if use a wider distribution
  return;
}
		     

// frequency discussed here is in GHz, and the power is expressed
// with a unit of 1e-22 W
//
// 0. Generate noise power spectrum with sigma= kT*b/sqrt(N)
// 1. Generate fake data that contains an axion signal with some
// random frequency between 749 and 751 MHz
// 2. Overlay the signal with background from thermal noise
// 3. Sum of signal and background is the measured spectrum
// 4. Weight each spectrum
// 5. Add to get signal/noise ratio and see if we could find a signal
//
// The full scanning range is lo-25 kHz to hi + 25 kHz

void search001(bool narrow=false,
	       bool debug=false,
	       bool changeSeed=false,
	       const unsigned int nTrials=1, 
	       const double lo = 749, const double hi = 751)
{

  if(S11 > 1 || fabs(S11-1)<1e-6 || S11 < 0) FATAL("S11 should be between 0 and 1");

  double freq_lo = lo - 0.5*rangeSpec;
  double freq_hi = hi + 0.5*rangeSpec;
  const unsigned int nSteps = (hi-lo)/step_size;

  cout << "Preparing a study with " << nTrials << " trials and ";
  cout << nSteps << " steps of frequency changes" << endl;
  cout << "sigma of noise is " << sigma_noise << endl;
  cout << "Grand spectrum frequency range is " <<
    freq_lo << " -- " << freq_hi << " MHz" << endl;
  
  TF1* fsig = new TF1("fsig", P_exp, freq_lo, freq_hi, 1);
  fsig->SetNpx(2500);

  TF1* fsig_wide = new TF1("fsig_wide", P_exp_boltzmann,
			   freq_lo, freq_hi, 2);
  fsig_wide->SetNpx(2500);


  double sigmaN[nSteps];
  for(unsigned int is=0; is<nSteps; is++)sigmaN[is]=sigma_noise;

  
  const unsigned int nTotalBins = (freq_hi-freq_lo)/bandwidth;
  TH1F* hGrand = new TH1F("hGrand","Grand spectrum",
			  nTotalBins, freq_lo, freq_hi);
  hGrand->SetXTitle("Frequency [MHz]");
  hGrand->SetYTitle("Power [10^{-22} Watts]");

  TH1F* htotal_delta[nTrials];
  TH1F* htotal_sigma[nTrials];
  TH1F* htotal_SoN[nTrials];
  
  TH1F* hsig[nTrials][nSteps];
  TH1F* hbkg[nTrials][nSteps];
  TH1F* hmea[nTrials][nSteps];
  TH1F* hSoN[nTrials][nSteps];

  TCanvas* c1 = new TCanvas("c1");
  // random number generators
  TString number = gSystem->GetFromPipe("./random.sh");
  UInt_t seed=number.Atoll();
  if(debug) cout << number << "\t" << seed << endl;
  TRandom3* gRandom= changeSeed? new TRandom3(seed): new TRandom3();

  // writing example output file
  string endstring = narrow? "_narrow":"_wide";
  TFile* outFile = new TFile(Form("sim_search001%s.root",endstring.data()),
			     "recreate");

    
  for(unsigned int itrial=0; itrial < nTrials; itrial++){

    htotal_delta[itrial]=(TH1F*)hGrand->Clone(Form("htotal_delta%03d",itrial));
    htotal_delta[itrial]->SetTitle("Weighted Spectrum");
    
    htotal_sigma[itrial]=(TH1F*)hGrand->Clone(Form("htotal_sigma%03d",itrial));
    htotal_sigma[itrial]->SetTitle("Standard Deviation of Weighted Sum");
    
    htotal_SoN[itrial]=(TH1F*)hGrand->Clone(Form("htotal_SoN%03d",itrial));
    htotal_SoN[itrial]->SetTitle("Signal/Noise from the Weighted Spectrum");
    htotal_SoN[itrial]->SetYTitle("");

    double delta_sum[nTotalBins];
    double weight_sum[nTotalBins];
    double sigma2_sum[nTotalBins];
    for(unsigned int ib=0;ib<nTotalBins;ib++)
      {delta_sum[ib]=0; weight_sum[ib]=0; sigma2_sum[ib]=0;}

    // first random peak a signal frequency
    double f_axion = lo + (hi-lo)*gRandom->Rndm();
    cout << "search for a signal with frequency = " << f_axion
	 << " MHz" << endl;
    // now start stepping frequency in size of step_size

    unsigned int nSpectra=0;
    
    for(unsigned int istep=0; (istep<nSteps) || (debug && nSpectra<2) ;
	istep++){

      double start_freq = freq_lo + istep*step_size;
      double end_freq   = start_freq + rangeSpec;

      if(debug)
	cout << "start:end frequencies = " << start_freq << "\t"
	     << end_freq << endl;

      TH1F* htemp = new TH1F(Form("htemp%03d%04d",itrial,istep),
			     "template of frequency",
			     nBins, start_freq, end_freq);
      htemp->SetXTitle("Frequency [MHz]");
      htemp->SetYTitle("Power [10^{-22} W]");

      hsig[itrial][istep] = (TH1F*)htemp->Clone(Form("hsig%03d%04d",
						     itrial,istep));
      hsig[itrial][istep] -> SetTitle(Form("Signal for trial %03d and"
					   " frequency step %04d",
					   itrial,istep));

      hbkg[itrial][istep] = (TH1F*)htemp->Clone(Form("hbkg%03d%04d",
						     itrial,istep));
      hbkg[itrial][istep] -> SetTitle(Form("Background for trial %03d and "
					   "frequency step %04d",itrial,istep));

      hmea[itrial][istep] = (TH1F*)htemp->Clone(Form("hmea%03d%04d",
						     itrial,istep));
      hmea[itrial][istep] -> SetTitle(Form("Measured for trial %03d and "
					   "frequency step %04d",itrial,istep));

      hSoN[itrial][istep] = (TH1F*)htemp->Clone(Form("hSoN%03d%04d",
						     itrial,istep));
      hSoN[itrial][istep] -> SetTitle(Form("Measured S/N for trial %03d and "
					   "frequency step %04d",itrial,istep));
      hSoN[itrial][istep] -> SetYTitle("");


      // now generate background spectrum
      generate_background(hbkg[itrial][istep], sigmaN[istep]);

      
      double resFreq = start_freq + 0.5*rangeSpec;
      if(debug) cout << "resonance frequency = " <<  resFreq << endl;
      
      generate_signal(narrow, debug, hsig[itrial][istep], fsig, fsig_wide,
		      f_axion, resFreq, end_freq);
      
      // add signal and background generation to get measured spectrum
      hmea[itrial][istep]->Add(hsig[itrial][istep],
			       hbkg[itrial][istep]);

      // compute the RMS over all frequency bins in each spectrum
      double RMS_for_this_spectra=0;
      for(int ib=1; ib <= nBins; ib++)
	RMS_for_this_spectra += pow(hmea[itrial][istep]->GetBinContent(ib),2);
      RMS_for_this_spectra /= nBins;
      RMS_for_this_spectra = sqrt(RMS_for_this_spectra);

      if(debug)
	cout << "RMS for this spectra = " << RMS_for_this_spectra << "\t"
	     << "input sigma = " << sigmaN[istep] << endl;

      for(int ib=1; ib<= nBins; ib++){
	if(RMS_for_this_spectra<1e-10)continue;
	hSoN[itrial][istep]->SetBinContent(ib,
	  hmea[itrial][istep]->GetBinContent(ib)/RMS_for_this_spectra);
      }
      
      
      if(f_axion < end_freq){
	nSpectra++;
	if(debug){
	  hsig[itrial][istep]->Write();
	  hsig[itrial][istep]->Draw();
	  c1->Print(Form("hsig%03d%04d.png",itrial,istep));

	  hbkg[itrial][istep]->Write();
	  hbkg[itrial][istep]->Draw();
	  c1->Print(Form("hbkg%03d%04d.png",itrial,istep));

	  hmea[itrial][istep]->Write();
	  hmea[itrial][istep]->Draw();
	  c1->Print(Form("hmea%03d%04d.png",itrial,istep));

	  hSoN[itrial][istep]->Write();
	  hSoN[itrial][istep]->Draw();
	  c1->Print(Form("hSoN%03d%04d.png",itrial,istep));

	}
      }


      // compute grand spectrum by adding spectra that overlapped frequencies
      
      fsig->SetParameter(0,resFreq);

      
      // then assign this RMS to compute weight      
      for(int ib=1; ib <= nBins; ib++){
	double binCenter = hmea[itrial][istep]->GetBinCenter(ib);
	int binGrand = hGrand->FindBin(binCenter);

	// double w = fsig->Eval(binCenter)/pow(sigmaN[istep],2);
     	
	double w = RMS_for_this_spectra < 1e-10? 0:
		   fsig->Eval(binCenter)/pow(RMS_for_this_spectra,2);

	double sigma2 = RMS_for_this_spectra < 1e-10? 0:
	  pow(fsig->Eval(binCenter)/(RMS_for_this_spectra),2);

	weight_sum[binGrand-1] += w;
	delta_sum[binGrand-1] += w*hmea[itrial][istep]->GetBinContent(ib);
	sigma2_sum[binGrand-1] += sigma2;
      }
      
	
    } // end of loop over steps of frequency

    // fill the measured grand spectrum with weighted average
    for(unsigned int ib=1; ib<= nTotalBins; ib++){
      // if(debug) cout << "delta_sum " << ib << " = " <<  delta_sum[ib-1] << endl;
      // if(debug) cout << "weight_sum " << ib << " = " <<  weight_sum[ib-1] << endl;
      if(weight_sum[ib]<1e-10)continue;
      double normalized = delta_sum[ib-1]/weight_sum[ib-1];
      htotal_delta[itrial]->SetBinContent(ib, normalized);
      // if(debug) cout << "normalized " << ib << " = " << normalized << endl;

      normalized =sqrt(sigma2_sum[ib-1])/weight_sum[ib-1];
      htotal_sigma[itrial]->SetBinContent(ib, normalized);

      double SN = delta_sum[ib-1]/sqrt(sigma2_sum[ib-1]);
      htotal_SoN[itrial]->SetBinContent(ib, SN);

    }
    htotal_delta[itrial]->Write();
    htotal_delta[itrial]->Draw();
    c1->Print(Form("htotal_delta%03d.png",itrial));

    htotal_sigma[itrial]->Write();
    htotal_sigma[itrial]->Draw();
    c1->Print(Form("htotal_sigma%03d.png",itrial));

    htotal_SoN[itrial]->Write();
    htotal_SoN[itrial]->Draw();
    c1->Print(Form("htotal_SoN%03d.png",itrial));

    
  } // end of loop over trials


  
   outFile->Close();


  
} // end of program









