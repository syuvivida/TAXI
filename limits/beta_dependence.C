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
const double fa=5e9;

Double_t betaFunc(Double_t* xx, Double_t* par){

  Double_t beta=xx[0]; // starts from 0
  Double_t Q0 =par[0]; // intrinsic Q
  Double_t QL = Q0/(1+beta);
  Double_t step = fa/QL*0.5;
  UInt_t    nSteps = par[1]; // from -nSteps to nSteps;
  Double_t factor = pow(beta/(1+beta)/(1+beta),2);
  
  Double_t sum = 0;
  for(int i=0; i< nSteps+1; i++)
    {
      Double_t thisFrequency = fa+ (i)*step;
      sum += factor*pow( 1/(1+pow(2*QL*(thisFrequency/fa-1),2)),2);

      if(i>0)
	{
	  thisFrequency = fa- (i)*step;
	  sum += factor*pow( 1/(1+pow(2*QL*(thisFrequency/fa-1),2)),2);

	}
    }
  
  return sum;
}

