#include <TMath.h>
#include <TF1.h>

Double_t CLEdge_Gauss(Double_t* xx, Double_t* par){

  Double_t chi=xx[0]; // in terms of sigma
  Double_t cut =par[0]; // in terms of sigma
  Double_t power=par[1];
  Double_t x= chi*TMath::Pi()/2.0; 
  Double_t weightPower = power/pow(x,2)*
    pow(TMath::Sin(x),2);
    
  Double_t cl=0.5*(1+TMath::Erf((weightPower-cut)/sqrt(2)));
  return cl;
}
