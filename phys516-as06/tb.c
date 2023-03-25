

int main(){
   double *e

   InitConf()

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
            for (int k=1; k<=n4; k++) h[k][j] = e[k]
         }
   }
