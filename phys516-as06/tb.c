#include "tb.h"
#include "math.h"
#include "stdio.h"


/*----------------------------------------------------------------------------*/
void InitConf() {
/*------------------------------------------------------------------------------
   r[][] is initialized to diamond lattice positions.
------------------------------------------------------------------------------*/
   double gap[3];  // Unit cell size
   double c[3];
   /* Atom positions in a unit diamond crystalline unit cell */
   double origAtom[NAUC][3] = {{0.0, 0.0, 0.0 }, {0.0, 0.5, 0.5 },
                                {0.5, 0.0, 0.5 }, {0.5, 0.5, 0.0 },
                                {0.25,0.25,0.25}, {0.25,0.75,0.75},
                                {0.75,0.25,0.75}, {0.75,0.75,0.25}};

   /* Read the # of unit cells in the x, y & z directions */
   scanf("%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);

   /* Sets up a diamond lattice */
   for (int k=0; k<3; k++) gap[k] = LCNS;
   nAtom = 0;
   for (int nZ=0; nZ<InitUcell[2]; nZ++) {
      c[2] = nZ*gap[2];
      for (int nY=0; nY<InitUcell[1]; nY++) {
         c[1] = nY*gap[1];
         for (int nX=0; nX<InitUcell[0]; nX++) {
            c[0] = nX*gap[0];
            for (int j=0; j<NAUC; j++) {
               for (int k=0; k<3; k++)
                  r[nAtom][k] = c[k] + gap[k]*origAtom[j][k];
               ++nAtom;
            }
         }
      }
   }
}

void htb() {

   double RegionH[3];
   double dr[3];

   for (int a=0; a<3; a++)
      RegionH[a] = 0.5*InitUcell[a]*LCNS;

   for (int i=1; i<=n4; i++)
      for (int j=1; j<=n4; j++)
         h[i][j] = 0.0;

   for (int i=0; i<nAtom; i++){
      int i40 = 4*i+1;
      int i41 = 4*i+2;
      int i42 = 4*i+3;
      int i43 = 4*i+4;

      h[i40][i40] = ES;
      h[i41][i41] = EP;
      h[i42][i42] = EP;
      h[i43][i43] = EP;

      for (int j=i+1; j<nAtom; j++) {
         int j40 = 4*j+1;
         int j41 = 4*j+2;
         int j42 = 4*j+3;
         int j43 = 4*j+4;

         double r2 = 0.0;
         for (int a=0; a<3; a++){
            dr[a] = r[i][a] - r[j][a];
            dr[a] = dr[a]-SignR(RegionH[a],dr[a] - RegionH[a]) - SignR(RegionH[a], dr[a] + RegionH[a]);
            r2 += dr[a]*dr[a];
         }
         double r1 = sqrt(r2);
         for (int a=0; a<3; a++) dr[a] /= r1;

         double sss = HSSS*pow(R0/r1,N)* exp( N* (-pow(r1/RSSS,NSSS)+pow(R0/RSSS, NSSS)) );
         double sps = HSPS*pow(R0/r1,N)* exp( N* (-pow(r1/RSPS,NSPS)+pow(R0/RSPS, NSPS)) );
         double pps = HPPS*pow(R0/r1,N)* exp( N* (-pow(r1/RPPS,NPPS)+pow(R0/RPPS, NPPS)) );
         double ppp = HPPP*pow(R0/r1,N)* exp( N* (-pow(r1/RPPP,NPPP)+pow(R0/RPPP, NPPP)) );


         h[i40][j40] = sss;
         h[i40][j41] = dr[0]*sps;
         h[i40][j42] = dr[1]*sps;
         h[i40][j43] = dr[2]*sps;

         h[i41][j40] = -h[i40][j41];
         h[i42][j40] = -h[i40][j42];
         h[i43][j40] = -h[i40][j43];

         h[i41][j41] = dr[0]*dr[0]*pps+(1.0-dr[0]*dr[0]) * ppp;
         h[i41][j42] = dr[0]*dr[1]*(pps-ppp);
         h[i41][j43] = dr[0]*dr[2]*(pps-ppp);
         h[i42][j42] = dr[1]*dr[1]*pps+(1.0-dr[1]*dr[1]) * ppp;
         h[i42][j43] = dr[1]*dr[2]*(pps-ppp);
         h[i43][j43] = dr[2]*dr[2]*pps+(1.0-dr[2]*dr[2]) * ppp;

         h[i42][j41] = h[i41][j42];
         h[i43][j41] = h[i41][j43];
         h[i43][j42] = h[i42][j43];

         for (int idummy=4*i+1; idummy<=4*i+4; idummy++)
            for (int jdummy=4*j+1; jdummy<=4*j+4; jdummy++)
               h[jdummy][idummy] = h[idummy][jdummy]; 


      }
   }
}





int main(){
   double *e;

   InitConf();

   n4 = 4*nAtom; // Hamiltonian matrix size with s-p basis h = dmatrix(1,n4,1,n4); // Use h[1:n4][1:n4]
   h = dmatrix(1, n4, 1 ,n4);
   d = dvector(1,n4); // d[1:n4]
   e = dvector(1,n4); // e[1:n4]
   /* Set up the Hamiltonian matrix elements h here */
   /* Diagonalize the Hamiltonian matrix */
   htb();

   tred2(h,n4,d,e);
   tqli(d,e,n4,h);

   for (int i=1; i<n4; i++)
      for (int j=i+1; j<=n4; j++)
         if (d[i] < d[j]){
            double dummy = d[i];
            d[i] = d[j];
            d[j] = dummy;
            for (int k=1; k<=n4; k++) e[k] = h[k][i];
            for (int k=1; k<=n4; k++) h[k][i] = h[k][j];
            for (int k=1; k<=n4; k++) h[k][j] = e[k];
         }
   for (int i=1; i<= n4; i++)
      printf("%d\t%le\n", i, d[i]);

   double de = (EMAX - EMIN)/NBIN;
   double fac = 1.0/(sqrt(M_PI)*SIGMA);
   for (int i=0; i<NBIN; i++) {
      eng[i] = EMIN + de*i;
      dos[i] = 0.0;
      for (int j=1; j<=n4; j++)
         dos[i] += fac*exp(-pow((eng[i]-d[j])/SIGMA,2));
      printf("%le\t%le\n", eng[i], dos[i]);
   }

   return 0;

}
   