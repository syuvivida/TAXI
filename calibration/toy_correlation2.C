#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <vector>
#include <TProfile.h>
#include <fstream>

using namespace std;

void toy_correlation2(unsigned int nExp=10000){

  TRandom3* gRandom = new TRandom3();
  ofstream fout;
  fout.open("toy.txt");
  
  for(unsigned int i=0; i< nExp; i++){

    double x = gRandom->Gaus(1,2);    
    double y = gRandom->Gaus(2,0.5);
    double z = 2*x -3*y;
    fout << z << " " << x << " " << y << endl;
  }

  fout.close();

}
