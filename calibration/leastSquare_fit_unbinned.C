#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <TMinuit.h>
#include <vector>
#include <TMath.h>
#include "TVirtualFitter.h"
#include <TPaveText.h>
#include "TFile.h"
#include "TGraph.h"

using namespace std;


// global variables

vector<vector<Double_t>> dataColl;
vector<Double_t> info;
vector<Double_t> info_err;

// number of variables, including the target variable
unsigned int NVAR=0;
unsigned int NORDER=0;

Double_t polN(Double_t *v, Double_t *par)
{
  Double_t funvalue = par[0];

  for(unsigned int i=0; i < NVAR; i++)
    {
      for(unsigned int j=1; j <= NORDER; j++)
	{
	  funvalue += par[j+NORDER*i]*pow(v[i+1],j);
	}
    }
      
  return funvalue;
}

Double_t DistPolN(Double_t *v, Double_t *par)
{
  
  Double_t value = pow(polN(v,par) - v[0],2);
  return value;
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  
  for (unsigned int i=0; i<dataColl.size(); i++ ) {
    //PDF for signal and background

    Double_t* data = new Double_t[NVAR];
    for(unsigned int j=0; j< NVAR; j++)
      data[j] = dataColl[i][j];
    //Get sum of least square
    Lsum += DistPolN(data,par);
  }
  f=Lsum;
}



//___________________________________________________________________________
Double_t* Ifit(std::string dataText="toy.txt", const unsigned int order=1, bool DEBUG=false)
{

  NORDER=order;
  dataColl.clear();

  unsigned int rows=0, cols = 0;
  string line, item;

  // count the number of variables by reading the first line

  ifstream file(dataText.data());
  while ( getline( file, line ) )
   {
     rows++;
     if ( rows == 1 )                 // First row only: determine the number of columns
       {
         stringstream ss( line );      // Set up up a stream from this line
         while ( ss >> item ) cols++;  // Each item delineated by spaces adds one to cols
      }
   }
   file.close();

   cout << "\n File has " << rows << " rows and " << cols << " columns" << endl;

   NVAR= cols;

   const unsigned int NPAR = (NVAR-1)*NORDER+1;
   const unsigned int NFIT = NPAR*2;
   // return fit result with mean and error
   Double_t* fitted = new Double_t[NFIT];
   for(unsigned int i=0; i<NFIT;i++)fitted[i]=0.0;

   if(NORDER<1)return fitted;
   vector<Double_t> myData;
   ifstream fin(dataText.data());
   const unsigned int nDataPoints=rows;
   Double_t xData[nDataPoints];
   Double_t yData[nDataPoints];
   
   for(unsigned int line=0; line < nDataPoints; line++)
     {
       myData.clear();
       xData[line]=0.0;
       yData[line]=0.0;
       
       for(unsigned int i=0; i < NVAR; i++)
	 {
	   Double_t tempVar;
	   fin >> tempVar;
	   myData.push_back(tempVar);
	   if(i==0)
	     {
	       xData[line]=line;       
	       yData[line]=tempVar;
	     }
	 }
       dataColl.push_back(myData);
     }
   

  long int ndata = dataColl.size();
  cout << "There are " << ndata << " data points " << endl;

  if(DEBUG){
    for(unsigned int i=0; i < ndata; i++)
      {
	for(unsigned int j=0; j< NVAR; j++)
	  cout << dataColl[i][j] << " ";
	cout << endl;
      }
  }

  
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

  for(unsigned int i=0; i<NPAR; i++){
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

  // plotting

  TGraph* gData = new TGraph(nDataPoints,xData,yData);
  gData->SetMarkerColor(2);
  gData->SetLineColor(2);
  // fill the function value
  Double_t funData[nDataPoints];
  for (unsigned int i=0; i<dataColl.size(); i++ ) {
    Double_t* data = new Double_t[NVAR];
    for(unsigned int j=0; j< NVAR; j++)
      data[j] = dataColl[i][j];
    //Get sum of least square
    funData[i]=polN(data,para);
  }

  TGraph* gFunc = new TGraph(nDataPoints,xData,funData);  
  gFunc->SetMarkerColor(4);
  gFunc->SetLineColor(4);

  gFunc->Draw("AL");
  gData->Draw("L");
  
  return fitted;
}


