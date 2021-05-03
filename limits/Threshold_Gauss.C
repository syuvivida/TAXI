#include <TMath.h>
#include <TF1.h>

Double_t Threshold(Double_t* xx, Double_t* par){

  Double_t cut=xx[0]; // in terms of sigma
  Double_t nBins =par[0]; 

  Double_t nCan = nBins*0.5*(1-TMath::Erf((cut)/sqrt(2)));
  return nCan;
}
