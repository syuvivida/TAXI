#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>

using namespace std;

void toy_correlation(unsigned int nExp=10000, bool correlated=false){

  TRandom3* gRandom = new TRandom3();
  vector<double> xdata;
  vector<double> ydata;
  double mean_x=0;
  double mean_y=0;
  double mean_x2=0;
  double mean_y2=0;
  
  for(unsigned int i=0; i< nExp; i++){

    double x = gRandom->Gaus(1,2);
    xdata.push_back(x);
    mean_x  += x;
    mean_x2 += x*x;
    
    double y= correlated? -5*x: gRandom->Gaus(2,0.5);
    ydata.push_back(y);

    mean_y  += y;
    mean_y2 += y*y;

  }

  mean_x /= (double)nExp;
  mean_y /= (double)nExp;
  
  mean_x2 /= (double)nExp;
  mean_y2 /= (double)nExp;
  

  double sigma_x = sqrt(mean_x2-mean_x*mean_x);
  double sigma_y = sqrt(mean_y2-mean_y*mean_y);

  cout << "Mean x = " << mean_x << "\t sigma_x = " << sigma_x << endl;
  cout << "Mean y = " << mean_y << "\t sigma_y = " << sigma_y << endl;

  double mean_cov=0;
  for(unsigned int i=0; i<nExp; i++)
    {
      double cov = (xdata[i]-mean_x)*(ydata[i]-mean_y);
      mean_cov += cov;
    }

  mean_cov /= (double)nExp;

  cout << "covariance = " << mean_cov << endl;

  double correlation = mean_cov/sigma_x/sigma_y;

  
  cout << "correlation= " << correlation << endl;

}
