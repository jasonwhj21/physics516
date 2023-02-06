/* Monte Carlo integration of PI by sample mean */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define Nseed 100
int main() {
  double x, pi, sum = 0.0, pi_av = 0.0, pi2_av = 0.0, stdv;;
  int try, ntry, outer;
  printf("Input the number of MC trials\n");
  scanf("%d",&ntry);
  srand((unsigned)time((long *)0));
  for (outer=0; outer<Nseed; outer++) {
    sum = 0.0;
    for (try=0; try<ntry; try++) {
      x = rand()/(double)RAND_MAX;
      sum += 4.0/(1.0 + x*x);
    }
    pi = sum/ntry;
    pi_av += pi;
    pi2_av += pi*pi;
  }
  pi_av /= Nseed;
  pi2_av /= Nseed;
  stdv = sqrt(pi2_av-pi_av*pi_av);
  printf("Ntry & stdv: %d = %f\n", ntry,stdv);
  return 0;
}
