#include <TMath.h>
#include <TF1.h>

// for 95% CL limits, TF1::GetX(0.05), for Gaussian distribution
// TF1* f1 = new TF1("f1", CLs,-10,10,1);
// f1->SetParameter(0,xxx) where xxx is (observed power - noise_mean)/sigma
// f1->GetX(0.05) gives the 95% CL upper limit on signal power in terms of sigma
Double_t CLs(Double_t* xx, Double_t* par){

  Double_t mu_s=xx[0]; // mu_s, now set mu_b=0; sigma_b=1;
  Double_t obs =par[0]; // in terms of sigma_b, can range from -inf to inf

  Double_t CLspb = (1+TMath::Erf((obs-mu_s)/sqrt(2)))*0.5;
  Double_t CLb = (1+TMath::Erf((obs)/sqrt(2)))*0.5;
  Double_t CLs = CLb>1e-16? CLspb/CLb:-1;
  return CLs;
}
