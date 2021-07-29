#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "setTDRStyle.C"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>  


using namespace std;


// global variables

struct myData
{
  // the dimension of the vector depends on the number of variables
  vector<Double_t> value; 
  vector<Double_t> error;
};

vector<myData> dataColl;
vector<vector<Double_t>> originalData;

// number of variables, including the target variable
unsigned int NVAR=0;
unsigned int NORDER=0;
unsigned int NPOINTS=0;

// return multi-variable polynomial functional values
Double_t polN(Double_t *v, Double_t *par)
{
  Double_t funvalue = par[0];

  for(unsigned int i=0; i < NVAR; i++)
    {
      for(unsigned int j=1; j <= NORDER; j++)
	{
	  funvalue += (par[j+NORDER*i]*pow(v[i+1],j));
	}
    }
      
  return funvalue;
}

// return multi-variable polynomial functional errors
Double_t polNError2(Double_t *v, Double_t *verr, Double_t *par)
{
  Double_t funvalue = 0;

  for(unsigned int i=0; i < NVAR; i++)
    {
      Double_t errTemp=0;
      for(unsigned int j=1; j <= NORDER; j++)
	{
	  errTemp += (par[j+NORDER*i]*pow(v[i+1],j-1)*j);
	}

      funvalue += (pow(errTemp,2)*pow(verr[i+1],2));
    }
      
  return funvalue;
}


// chi-square for each bin
Double_t Chi2PolN(Double_t *v, Double_t* verr, Double_t *par)
{
  Double_t dataError = verr[0];
  Double_t funcError2 = polNError2(v, verr, par);
  //  Double_t totalError2 = pow(dataError,2)+ funcError2;
  Double_t totalError2 = pow(dataError,2);
  Double_t value = (totalError2<1e-12 || NPOINTS<2) ?
				pow(polN(v,par) - v[0],2) : pow(polN(v,par) - v[0],2)/totalError2;
  return value;
}



// summing up the chi-square of all bins
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  
  for (unsigned int i=0; i<dataColl.size(); i++ ) {
    //PDF for signal and background

    Double_t* dataValue = new Double_t[NVAR];
    Double_t* dataError = new Double_t[NVAR];
    
    for(unsigned int j=0; j< NVAR; j++)
      {
	dataValue[j] = dataColl[i].value[j];
	dataError[j] = dataColl[i].error[j];
      }
    //Get sum of least square
    Lsum += Chi2PolN(dataValue,dataError,par);
  }
  f=Lsum;
}



// Major function to call the fit and display results,
// when nPoints = 1, the fit is not a chi-square fit, the FCN becomes the distance.
//___________________________________________________________________________
Double_t* Ifit(std::string dataText="toy.txt", const unsigned int nPoints=10, const unsigned int order=1, bool DEBUG=false)
{
  setTDRStyle();
  NORDER=order;
  NPOINTS=nPoints;
  originalData.clear();
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
   vector<Double_t> myDataValue;
   vector<Double_t> myDataError;
   ifstream fin(dataText.data());
   const unsigned int nDataPoints=rows;
   Double_t xData[nDataPoints];
   Double_t yData[nDataPoints];

   Double_t sum[NVAR], sumSquare[NVAR];

   
   unsigned int ibin=0;
   const unsigned int NBINS= nDataPoints/nPoints + ceil( double(nDataPoints%nPoints)/double(nPoints));
   Double_t xBinData[NBINS];
   Double_t yBinData[NBINS];
   Double_t yBinErrorData[NBINS];


   //-----------------------------------------------------------------------
   // Organizing data
   //-----------------------------------------------------------------------

   unsigned int nCount=0;
   for(unsigned int line=0; line < nDataPoints; line++)
     {
       myDataValue.clear();
       myDataError.clear();
       xData[line]=0.0;
       yData[line]=0.0;

       if(line%nPoints==0)nCount=0;
       nCount++;

       // dimension = NVAR
       vector<Double_t> tempOriginalData; 
       vector<Double_t> tempValueVec;
       vector<Double_t> tempErrorVec;
       
       
       for(unsigned int i=0; i < NVAR; i++)
	 {
	   // Reset every nPoints
	   if(line%nPoints==0){
	     sum[i]=0;
	     sumSquare[i]=0;
	   }

	   Double_t tempVar;
	   fin >> tempVar;
	   tempOriginalData.push_back(tempVar);

	   if(DEBUG)
	     cout << "before  line" << line << "\t i = "  << i << " : " <<
	       setprecision(12) << tempVar << "\t" << 
	       setprecision(12) << sumSquare[i] << "\t" << setprecision(12)
		  << sum[i] << endl;
	   
	   sum[i] += tempVar;
	   sumSquare[i] += pow(tempVar,2);	   
	   if(i==0)
	     {
	       xData[line]=line;       
	       yData[line]=tempVar;
	     }

	   if(DEBUG)
	     cout << "debug  line" << line << "\t i = "  << i << " : " <<
	       setprecision(12) << tempVar << "\t" << 
	       setprecision(12) << sumSquare[i] << "\t" << setprecision(12)
		  << sum[i] << endl;

	   // compute the error on the mean when we reach nPoints for average
	   if( line%nPoints==(nPoints-1) ||
	      line== (nDataPoints-1) )
	     {
	       Double_t mean = sum[i]/(Double_t)nCount;
	       // taking into account that the number of data points we have is small
	       Double_t error = nCount<2? 1:
		 sqrt((sumSquare[i] - 2*sum[i]*mean + mean*mean*nCount)/(Double_t)(nCount-1));
	       error = nCount<1? 1: error/sqrt(nCount);
	       // // dump numbers;
	       if(i==0)
		 {
		   xBinData[ibin] = ibin;
		   yBinData[ibin] = mean;		   
		   yBinErrorData[ibin] = error;
		   if(DEBUG)
		     cout << xBinData[ibin] << "\t" << yBinData[ibin] << " +- " << yBinErrorData[ibin] << endl;
		 }
		
	       tempValueVec.push_back(mean);
	       tempErrorVec.push_back(error);
	     
	     } // when we reach nPoints cycle

	 } // end of loop over variables
       
       if(line%nPoints==(nPoints-1) || line== (nDataPoints-1))
	 {
	   if(DEBUG)
	     for(unsigned int iv=0; iv < tempValueVec.size(); iv++)
		 cout << line << ": " << iv << "\t " << tempValueVec[iv] << " +- " << tempErrorVec[iv] << endl;
	   myData tempDataStruct;
	   tempDataStruct.value = tempValueVec;
	   tempDataStruct.error = tempErrorVec;
	   
	   dataColl.push_back(tempDataStruct);
	   ibin++;
	   
	 }

       originalData.push_back(tempOriginalData);
     } // looping over input file lines


   long int ndata = originalData.size();
   cout << "There are " << ndata << " data points " << endl;

  long int nbins = dataColl.size();
  cout << "There are " << nbins << " bins " << endl;

  if(DEBUG){
    for(unsigned int i=0; i < ndata; i++)
      {
	for(unsigned int j=0; j< NVAR; j++)
	  cout << originalData[i][j] << " ";
	cout << endl;
      }
  }
  

  if(DEBUG){
    for(unsigned int i=0; i < nbins; i++)
      {
	for(unsigned int j=0; j< NVAR; j++)
	  cout << dataColl[i].value[j] << "+-" << dataColl[i].error[j] << " ";
	cout << endl;
      }
  }

  
  if(ndata==0) {
    printf(" ---  no evetns in the fit \n");
    return fitted;
  }

  //--------------------------------------------------
  // Fitting starts

  int fit_status;


  //--------------------------------------------------
  //init parameters for fit
  Double_t vstart[NPAR];
  for(unsigned int i=0; i<NPAR; i++)vstart[i]=0.0;
 
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


  double para[NPAR],errpara[NPAR];
  if ( ierflg == 0 ) 
    {
      for(unsigned int j=0; j<NPAR;j++) {
        gMinuit->GetParameter(j, para[j],errpara[j]);
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

  //-----------------------------------------------------------------------

  // plotting the unbinned data first 

  TGraph* gData = new TGraph(nDataPoints,xData,yData);
  gData->SetMarkerColor(4);
  gData->SetLineColor(4);
  // fill the function value
  Double_t funData[nDataPoints];
  Double_t ratioData[nDataPoints];
  for (unsigned int i=0; i< originalData.size(); i++ ) {
    Double_t* data = new Double_t[NVAR];
    for(unsigned int j=0; j< NVAR; j++)
      data[j] = originalData[i][j];
    //Get sum of least square
    funData[i]=polN(data,para);
    ratioData[i]=funData[i]/yData[i];
  }

  TGraph* gFunc = new TGraph(nDataPoints,xData,funData);  
  gFunc->SetMarkerColor(2);
  gFunc->SetLineColor(2);
  gData->SetTitle("Step 18, 19");
  gData->GetXaxis()->SetTitle("Time");
  gData->GetYaxis()->SetTitle("Temperature (K)");
  gData->GetXaxis()->SetDecimals();
  gData->GetYaxis()->SetDecimals();
  gData->GetYaxis()->SetRangeUser(2.41,2.51);
  gData->GetYaxis()->SetNdivisions(5);

  const double LABELSIZE = 20.0;

  TCanvas* c1 = new TCanvas("c1","",1400,1000);
  c1->Divide(2,2,0.01,0);
  c1->cd(1);
  double temp1_pad = gPad->GetWh()*gPad->GetAbsHNDC();
  double label1_size = LABELSIZE/temp1_pad;
  gData->GetXaxis()->SetLabelSize(label1_size);
  gData->GetYaxis()->SetLabelSize(label1_size);
  
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  gData->Draw("AL*");
  gFunc->Draw("L");

  TLegend* leg = new TLegend(0.193,0.763,0.470,0.936);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.040);
  leg->SetBorderSize(0);
  leg->SetHeader("Fit with order 1");
  leg->AddEntry(gData,"Unbinned Data");
  leg->AddEntry(gFunc,"Fitted Function");
  leg->Draw("same");


  // plotting the binned data second
  Double_t xBinErrorData[NBINS];

  for(unsigned int ibin=0; ibin<NBINS; ibin++)
    xBinErrorData[ibin] = 0.25;

  TGraphErrors* gBinData = new TGraphErrors(NBINS,xBinData,yBinData, xBinErrorData,yBinErrorData);
  gBinData->SetMarkerColor(4);
  gBinData->SetLineColor(4);
  // fill the function value
  Double_t funBinData[NBINS];
  Double_t ratioBinData[NBINS];
  for (unsigned int i=0; i< dataColl.size(); i++ ) {
    Double_t* data = new Double_t[NVAR];
    for(unsigned int j=0; j< NVAR; j++)
      data[j] = dataColl[i].value[j];
    //Get sum of least square
    funBinData[i]=polN(data,para);
    ratioBinData[i]=funBinData[i]/yBinData[i];
  }

  TGraph* gBinFunc = new TGraph(NBINS,xBinData,funBinData);  
  gBinFunc->SetMarkerColor(2);
  gBinFunc->SetLineColor(2);
  gBinData->SetTitle("Step 18, 19");
  gBinData->GetXaxis()->SetTitle("Time");
  gBinData->GetYaxis()->SetTitle("Temperature (K)");
  gBinData->GetXaxis()->SetLabelSize(label1_size);
  gBinData->GetYaxis()->SetLabelSize(label1_size);
  gBinData->GetXaxis()->SetDecimals();
  gBinData->GetYaxis()->SetDecimals();
  gBinData->GetYaxis()->SetRangeUser(2.41,2.51);
  gBinData->GetYaxis()->SetNdivisions(5);

  c1->cd(2);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.08);

  
  gBinData->Draw("AL*");
  gBinFunc->Draw("L");

 
  TLegend* leg2 = new TLegend(0.129,0.711,0.408,0.885);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.045);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.040);
  leg2->SetBorderSize(0);
  leg2->SetHeader("Fit with order 1");
  leg2->AddEntry(gBinData,"binned Data");
  leg2->AddEntry(gBinFunc,"Fitted Function");
  leg2->Draw("same");


  // plotting ratios
  
  c1->cd(3);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetTickx();
  double temp2_pad = gPad->GetWh()*gPad->GetAbsHNDC();
  double label2_size = LABELSIZE/temp2_pad;
  
  TGraph* gRatio = new TGraph(nDataPoints,xData,ratioData);
  gRatio->SetMarkerColor(1);
  gRatio->SetLineColor(1);
  gRatio->SetTitle("");
  gRatio->GetXaxis()->SetTitle("Time");
  gRatio->GetYaxis()->SetTitle("Fit/Data");
  gRatio->GetXaxis()->SetLabelSize(label2_size);
  gRatio->GetYaxis()->SetLabelSize(label2_size);
  gRatio->GetXaxis()->SetDecimals();
  gRatio->GetYaxis()->SetDecimals();
  gRatio->Draw("AL");

  
  c1->cd(4);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.08);
  gPad->SetTickx();
  
  TGraph* gRatioBin = new TGraph(NBINS,xBinData,ratioBinData);
  gRatioBin->SetMarkerColor(1);
  gRatioBin->SetLineColor(1);
  gRatioBin->SetTitle("");
  gRatioBin->GetXaxis()->SetTitle("Time");
  gRatioBin->GetYaxis()->SetTitle("Fit/Data");
  gRatioBin->GetXaxis()->SetLabelSize(label2_size);
  gRatioBin->GetYaxis()->SetLabelSize(label2_size);
  gRatioBin->GetXaxis()->SetDecimals();
  gRatioBin->GetYaxis()->SetDecimals();
  gRatioBin->Draw("AL");

  
  return fitted;
}


