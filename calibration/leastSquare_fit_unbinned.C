#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <TMinuit.h>
#include <vector>
#include <TMath.h>
#include "TVirtualFitter.h"
#include <TPaveText.h>
#include "TFile.h"

using namespace std;

struct myData
{
  Double_t x;
  Double_t y;
  Double_t z;
};



vector<myData> dataColl;
vector<Double_t> info;
vector<Double_t> info_err;

#define NPAR 5

Double_t pol1(Double_t *v, Double_t *par)
{
  Double_t x = v[0];
  Double_t y = v[1];
  Double_t z = v[2];
  Double_t value = pow(par[0]+par[1]*x+par[2]*y-z,2);
  return value;
}

Double_t pol2(Double_t *v, Double_t *par)
{
  Double_t x = v[0];
  Double_t y = v[1];
  Double_t z = v[2];
  Double_t value = pow(par[0]+par[1]*x+par[2]*y+par[3]*pow(x,2)+par[4]*pow(y,2)-z,2);
  return value;
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  Double_t Nevt=dataColl.size();

  
  for ( int i=0; i<dataColl.size(); i++ ) {
    //PDF for signal and background

    Double_t data[]
      ={
	dataColl[i].x,
	dataColl[i].y,
	dataColl[i].z
    };
    
    //Get sum of least square
    Lsum += pol2(data,par);
  }
  f=Lsum;
}



//___________________________________________________________________________
Double_t* Ifit(std::string dataText="toy.txt")
{

  dataColl.clear();
  Double_t* fitted = new Double_t[NPAR*2];
  for(int i=0; i<NPAR*2;i++)fitted[i]=0.0;
  
  FILE *infile =  fopen(dataText.data(),"r");
  myData tmpData;
  Double_t x, y, z;
  int flag = 1;
  while (flag!=-1){
    flag =fscanf(infile,"%lf %lf %lf",&x, &y, &z);

    tmpData.x = x;
    tmpData.y = y;
    tmpData.z = z;
    
    dataColl.push_back(tmpData);
  }

  long int ndata = dataColl.size();
  cout << "There are " << ndata << " data points " << endl;
  
  if(ndata==0) {
    printf(" ---  no evetns in the fit \n");
    return fitted;
  }
    
  int fit_status;


  //--------------------------------------------------
  //init parameters for fit
  Double_t vstart[NPAR];
  for(unsigned int i=0; i<NPAR; i++)vstart[i]=1.1;
 
  TMinuit *gMinuit = new TMinuit(NPAR);  
  gMinuit->Command("SET STR 1");
  gMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  arglist[0] = 1;
  gMinuit->mnexcm("SET PRINT", arglist ,1,ierflg);

  Double_t step[NPAR];
  for(unsigned int i=0; i<NPAR; i++)step[i]=0.1;

  for(int i=0; i<NPAR; i++){
    gMinuit->mnparm(i,  Form("p%d",i), vstart[i],  step[i], -9999,9999 , ierflg);
  }
  

  printf(" --------------------------------------------------------- \n");
  printf(" Setting stragety = 2 \n ----------------------\n");

  arglist[0] = 2;
  gMinuit->mnexcm("SET STRAT", arglist ,1,ierflg);

  printf(" --------------------------------------------------------- \n");
  printf(" Now ready for minimization step \n ----------------------\n");

  arglist[0] = 2000; // number of iteration
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  printf (" -------------------------------------------- \n");
  printf("Finished.  ierr = %d \n", ierflg);

  info.clear();
  info_err.clear();

  double para[NPAR],errpara[NPAR];
  if ( ierflg == 0 ) 
    {
      for(unsigned int j=0; j<NPAR;j++) {
        gMinuit->GetParameter(j, para[j],errpara[j]);
        info.push_back(para[j]);
        info_err.push_back(errpara[j]);
        printf("Parameter  %d = %f +- %f\n",j,para[j],errpara[j]);
	
      }
    }
  else {
    printf(" *********** Fit failed! ************\n");
  }

  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(1,amin);  
  gMinuit->mnmatu(1);
  printf(" ========= happy ending !? =========================== \n");  
  printf("FCN =  %3.3f \n", amin);


  for(unsigned int i=0; i< NPAR; i++)
    {
      fitted[i*2] = para[i];
      fitted[i*2+1] = errpara[i];
    }

  return fitted;
}


