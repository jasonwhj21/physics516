#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAX_DAY 365
#define DT (1.0/MAX_DAY)
#define MU 0.14
#define SIGMA 0.2
#define N_WALKER 1000
#define S_INIT 20
#define N_HIST 50

double rand_normal() {
  double r1 = rand()/(double)RAND_MAX;
  double r2 = rand()/(double)RAND_MAX;
  return sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2); 
}

int main() {
  //int Max_step; // Maximum number of random-walk steps
  //int N_walker; // Number of walkers
  int day,walker,k;
  double s; // Drunkard's position
  int hist[N_HIST];

  double mu_dt = MU*DT;
  double sigma_rtdt = SIGMA*sqrt(DT);

  FILE *pst;
  pst = fopen("st.out", "w");



  for (k=0; k<=N_HIST; k++)
    hist[k] = 0;

  srand((unsigned)time((long *)0)); // Initialize the random-number sequence

  for (walker=1; walker<=N_WALKER; walker++) {
    s = S_INIT;
    for (day=1; day<=MAX_DAY; day++) {
      s += s*(mu_dt+sigma_rtdt*rand_normal());
      if (s<0.0) break;
      if (walker == 1) fprintf(pst, "%d\t%le\n", day, s);

    } // Endfor step
    s = s>0.0 ? s:0.0;
    ++hist[(int)s];
  } // Endfor walker

  for (k=0; k<N_HIST; k++)
    printf("%d\t%d\n",k,hist[k]);

  fclose(pst);
  
  return 0;
}
