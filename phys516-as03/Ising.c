/* Monte Carlo simulation of Ising spings */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define L 20
int s[L][L];
double exp_dV[2][5];
double JdivT; // J/kBT
double HdivT; // H/kBT
int hist[2*L*L+1];

void table_set() {
  for (int k=0; k<2; k++) {
    int s_new = 2*k-1;
    for (int l=0; l<5; l++){
      int S = 2*l-4;
      exp_dV[k][l] = exp(2*s_new*(JdivT*S+HdivT));
    }
  }
}
int main() {
  double runM;
  double sumM = 0.0, sumM2 = 0.0;
  int Sta_step;
  double avgM, sigM;

  printf("Input JdivT HdivT Sta_step\n");
  scanf("%le %le %d",&JdivT, &HdivT, &Sta_step);

  table_set();

  //initialize the spins, s[i][j] (0 ≤ i,j ≤ L-1)
  for (int i=0; i<L; i++){
    for (int j=0; j<L; j++){
      s[i][j] = 1;
    }
  }
  runM = 1.0*L*L;
  for (int i=0; i<2*L*L+1; i++ ){
    hist[i] = 0;
  }

  srand((unsigned)time((long *)0));
//for (try=0; try<ntry; try++) {
  for (int step=0; step<Sta_step; step++) {
//  x = rand()/(double)RAND_MAX;
//  sum += 4.0/(1.0 + x*x);
    int i = rand()%L;
    int j = rand()%L;
    int s_new = -s[i][j];
    int k = (1+s_new)/2;
    int im = (i + L - 1) % L;
    int ip = (i + 1) % L;
    int jm = (j + L - 1) % L;
    int jp = (j + 1) % L;
    int S = s[im][j]+s[ip][j]+s[i][jm]+s[i][jp];
    int l = (4+S)/2;
    double exp_val = exp_dV[k][l];
    if (exp_val > 1.0) {
      s[i][j] = s_new;
      runM += 2*s_new;
      
    }
    else if (rand()/(double)RAND_MAX <= exp_val) {
      s[i][j] = s_new;
      runM += 2*s_new;
    }
    sumM += runM;
    sumM2 += runM*runM;
    hist[(int)runM+L*L] += 1;
  }
  //pi = sum/ntry;
  avgM = sumM/Sta_step;
  sigM = sqrt(sumM2/Sta_step-avgM*avgM);

  printf("M and StdV = %le\t%le\n", avgM, sigM);
  for (int i=0; i<2*L*L+1; i++){
    printf("%d\t%d\n", i-L*L, hist[i]);
  }

  return 0;
}
