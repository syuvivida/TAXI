#include <TMath.h>
#include <TF1.h>

// for 95% CL limits, TF1::GetX(0.05), for Poisson distribution
Double_t CLs(Double_t* xx, Double_t* par){

  Double_t mu_s=xx[0]; // mu_s
  Double_t mu_b=par[0];
  Double_t obs=par[1];

  Double_t prob_deno = TMath::Prob(2*mu_b,2*(obs+1));
  Double_t prob_numr = TMath::Prob(2*(mu_b+mu_s),2*(obs+1));

  return prob_numr/prob_deno;

}
