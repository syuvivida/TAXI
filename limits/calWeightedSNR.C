{

  float sum =0;
  float sumerr=0;
  float x0=0; // sigma
  float x1=1; // power
  float x2=0; // weight
  cin >> x0 >> x1 >> x2;
  while(x2 > 0)
    {
      sum += x1*x2; 
      sumerr += x0*x0*x2*x2;
      cin >> x0 >> x1 >> x2;
    }


  sumerr = sqrt(sumerr);
  float SNR = sum/sumerr;
  printf("SNR = %1.6lf \n", SNR);

}
