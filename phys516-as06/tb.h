#define NMAX 100
#define NAUC 8
#define LCNS (1.0*5.43)
#define NBIN 500
#define EMIN -15.0
#define EMAX 10.0
#define SIGMA 0.1

#define R0 2.360352
#define N 2
#define ES (-5.25)
#define EP 1.2
#define HSSS (-2.038)
#define NSSS 9.5
#define RSSS 3.4
#define HSPS 1.745
#define NSPS 8.5
#define RSPS 3.55
#define HPPS 2.75
#define NPPS 7.5
#define RPPS 3.7
#define HPPP (-1.075)
#define NPPP 7.5
#define RPPP 3.7


double SignR(double v, double x) {if (x>0) return v; else return -v;}
double **dmatrix(int, int, int, int);
double *dvector(int, int);
void tred2(double **, int, double *, double *);
void tqli(double *, double *, int, double **);

int InitUcell[3];

int nAtom;
int n4;
double r[NMAX][3];
double **h;
double *d;
double eng[NBIN], dos[NBIN];



