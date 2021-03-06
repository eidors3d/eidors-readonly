// Jacobian.cpp : Defines the entry point for the console application.
//
#if defined(_OPENMP)
#include <omp.h>
#endif

#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <sstream>
using std::istringstream;
using std::ostringstream;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void original( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{

/*
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x,*y;
  mwSize mrows,ncols;

  if(nrhs!=3) {
    mexErrMsgTxt("Three inputs required.");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mxGetM(prhs[0])==1)|| !(mxGetN(prhs[0])==4) ) 
  {  
    mexErrMsgTxt("1x4 expected.");
  }
*/

  double * sz = S_sz;
  
  int dims[] = {sz[0],sz[1],sz[2]};
//  double DE [DE_SZ]; // test eq to S_DE
  // printf("DE [%d x %d x %d]\n", dims[0],dims[1],dims[2]);
  
//   d1 = sz(4)-1;
//   for k= 1: sz(3);
//       idx= d1*(k-1) + (1:d1);
//       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
//       DE(:,:,k)= dq;
//   end

  int d1 = (int) sz[3] - 1;
  double * M1 = S_zi2E_FCt;
  double * M2 = S_FC_sv;
  int R = GetM_S_zi2E_FCt;
  int L = GetN_S_zi2E_FCt;
//printf("M1[%d x %d] M2[%d x %d] d1=%d\n",
//       GetM_S_zi2E_FCt, GetN_S_zi2E_FCt,
//       GetM_S_FC_sv,    GetN_S_FC_sv,     d1);

  int E = (int) sz[0]; // assume sz[1] == sz[0]
  int E2= E*E;
  for (int l=0; l< sz[2]; l++) {
     int ofs = d1*l;
     for (int i=0; i< R; i++) {
     for (int j=0; j< R; j++) {
        double sum= 0.0;
        for (int k=0; k< d1; k++) {
// DE(i,j,l) = sum_k M1(i, l*d1 + k) * M2(l*d1+k, j)
          int m1o=  i + R*(ofs + k);
          int m2o=  ofs + k + L*j;
//if (i<=1 && l==0 && j<=1)
//printf("i=%d j=%d k=%d m1[%d]=%f m2[%d]=%f\n",i,j,k, m1o, M1[m1o], m2o, M2[m2o]);
          sum+= M1[m1o] * M2[m2o];
        }
        DE[i + E*j + E2*l] = sum;
     }
     }
  }

  //double diff = 0.0;
  //for(int i=0; i<DE_SZ; i++) {
    //diff += fabs( DE[i] - S_DE[i]); 
  //}

  //printf("fini. Diff (should be 0) is %g\n",diff);
}

void alpha( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{

    const double * sz = S_sz;
    const int dims[] = {sz[0],sz[1],sz[2]};

    // cout << "DE [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << endl;

    const int d1 = (int) sz[3] - 1;
    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    const int R = GetM_S_zi2E_FCt;
    const int L = GetN_S_zi2E_FCt;

    // cout << "M1[" << GetM_S_zi2E_FCt << " x " << GetN_S_zi2E_FCt << "]"
    //     << " M2[" << GetM_S_FC_sv << " x " << GetN_S_FC_sv << "] d1=" << d1 << endl;

    const int E = (int) sz[0]; // assume sz[1] == sz[0]
    const int E2= E*E;
    const int Rd1 = R * d1;

    // long long mult = 0;
        for (int l=0, ofs=0, E2l=0, Rofs = 0; l< sz[2]; l++, ofs += d1, E2l += E2, Rofs += Rd1 )
        {
            for (int i=0, m1o = Rofs; i< R; i++, m1o++)
            {
                for (int j=0, Lj = 0, Ej = 0, m2o = ofs; j< R; j++, m2o += L, Ej += E )
                {
                    double sum= 0.0;
                    for (int k=0, rk=0; k< d1; k++, rk += R )
                    {
                        //                    if (i<=1 && l==0 && j<=1)
                        //                        printf("i=%d j=%d k=%d m1[%d]=%f m2[%d]=%f\n",i,j,k, m1o+rk, M1[m1o+rk], m2o+k, M2[m2o+k]);
                        sum += M1[m1o+rk] * M2[m2o+k];
                        // ++mult;
                    }
                    DE[i + Ej + E2l] = sum;
                }
            }
        }

    // cout << "Mult: " << mult << endl;
}

void beta( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{

    const double * sz = S_sz;
    const int dims[] = {sz[0],sz[1],sz[2]};

#if defined(_OPENMP)
    // cout << "Found " << omp_get_num_procs() << " processors" << endl;
#endif
    // cout << "DE [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << endl;

    const int d1 = (int) sz[3] - 1;
    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    const int R = GetM_S_zi2E_FCt;
    const int L = GetN_S_zi2E_FCt;

    // cout << "M1[" << GetM_S_zi2E_FCt << " x " << GetN_S_zi2E_FCt << "]"
    //     << " M2[" << GetM_S_FC_sv << " x " << GetN_S_FC_sv << "] d1=" << d1 << endl;

    const int E = (int) sz[0]; // assume sz[1] == sz[0]
    const int E2= E*E;
    const int Rd1 = R * d1;

    // long long mult = 0;
        int l;
        int maxL = sz[2];
#pragma omp parallel for
        for (l=0; l < maxL; l++ )
        {
            const int ofs = l * d1;
            const int E2l = l * E2;
            const int Rofs = l * Rd1;

            for (int i=0, m1o = Rofs; i< R; i++, m1o++)
            {
                for (int j=0, Lj = 0, Ej = 0, m2o = ofs; j< R; j++, m2o += L, Ej += E )
                {
                    // ostringstream ostr;
                    double sum= 0.0;
                    for (int k=0, rk=0; k< d1; k++, rk += R )
                    {
                        //                    if (i<=1 && l==0 && j<=1)
                        //                        printf("i=%d j=%d k=%d m1[%d]=%f m2[%d]=%f\n",i,j,k, m1o+rk, M1[m1o+rk], m2o+k, M2[m2o+k]);
                        sum += M1[m1o+rk] * M2[m2o+k];
                        // ++mult;
                        // ostr << "(" << (m1o+rk) << "," << (m2o+k) << ")";
                    }
                    // cout << (i + Ej + E2l) << "=" << ostr.str() << endl;
                    DE[i + Ej + E2l] = sum;
                }
            }
        }

    // cout << "Mult: " << mult << endl;
}


void gamma( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{
// X1=INT(A2/64)*16+MOD(MOD(A2,64),8)
// Y1==INT(A2/64)*2+INT(MOD(A2,64)/8)*128
// X2=X1+8
// Y2=Y1+1

// X1=INT(A2/196)*42+MOD(MOD(A2,196),14)
// Y1==INT(A2/196)*3+INT(MOD(A2,196)/14)*23529
// X2=X1+14
// Y2=Y1+1
// X3=X1+28
// Y3=Y1+2

    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    int maxI = S_sz[0] * S_sz[1] * S_sz[2];

#pragma omp parallel for
        for( int index = 0; index < 1463140/14; ++index )
        {
             //       int      x1 = (( i >> 6 ) << 4 ) + (( i & 63 ) & 7);
             // const int      y1 = (( i >> 6 ) << 1 ) + (( ( i & 63 ) >> 3 ) * GetN_S_zi2E_FCt );
             //       int      x2 = x1 + 8;
             // const int      y2 = y1 + 1;

                      int       i = index * 14;
                      int      x1 = (i/196)*42 + ((i%196)%14);
                const int      y1 = (i/196)*3  + ((i%196)/14) * 23529;
                      int      x2 = x1 + 14;
                const int      y2 = y1 + 1;
                      int      x3 = x1 + 28;
                const int      y3 = y1 + 2;

                const double m2y1 = M2[y1];
                const double m2y2 = M2[y2];
                const double m2y3 = M2[y3];

                DE[i] = M1[x1]*m2y1  +M1[x2]*m2y2  +M1[x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;

                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;

                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;

                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;
                DE[++i] = M1[++x1]*m2y1  +M1[++x2]*m2y2  +M1[++x3]*m2y3  ;

                // cout << i << " " << x1 << "," << y1 << "," << x2 << "," << y2 << "," << x3 << "," << y3 << endl;
        }
}

void delta( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{
// X1=INT(A2/64)*16+MOD(MOD(A2,64),8)
// Y1==INT(A2/64)*2+INT(MOD(A2,64)/8)*128
// X2=X1+8
// Y2=Y1+1

    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    int maxI = (S_sz[0] * S_sz[1] * S_sz[2])/8;

    for( int i = 0; i < maxI; ++i )
    {
        int iter = i << 3;

              int      x1 = (( iter >> 6 ) << 4 ) + (( iter & 63 ) & 7);
        const int      y1 = (( iter >> 6 ) << 1 ) + (( ( iter & 63 ) >> 3 ) * GetN_S_zi2E_FCt );
              int      x2 = x1 + 8;
        const int      y2 = y1 + 1;
        const double m2y1 = M2[y1];
        const double m2y2 = M2[y2];

        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
    }
}

void epsilon( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE )
{
// X1=INT(A2/64)*16+MOD(MOD(A2,64),8)
// Y1==INT(A2/64)*2+INT(MOD(A2,64)/8)*128
// X2=X1+8
// Y2=Y1+1

    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    int maxI = (S_sz[0] * S_sz[1] * S_sz[2])/8;

#pragma omp parallel for
    for( int i = 0; i < maxI; ++i )
    {
        int iter = i << 3;

              int      x1 = (( iter >> 6 ) << 4 ) + (( iter & 63 ) & 7);
        const int      y1 = (( iter >> 6 ) << 1 ) + (( ( iter & 63 ) >> 3 ) * GetN_S_zi2E_FCt );
              int      x2 = x1 + 8;
        const int      y2 = y1 + 1;
        const double m2y1 = M2[y1];
        const double m2y2 = M2[y2];

        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
        DE[iter++] = M1[x1++]*m2y1+M1[x2++]*m2y2;
    }
}

int atoi( char *str )
{
   int ret;

   std::string s( str );
   istringstream istr( s );

   istr >> ret;

   return ret;
}


void runTrial( string name, void (*fn)( double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE ), double *S_sz, double *S_zi2E_FCt, double *S_FC_sv, int GetM_S_zi2E_FCt, int GetN_S_zi2E_FCt, int GetM_S_FC_sv, int GetN_S_FC_sv, double *DE, const int maxrep )
{
    std::cout << std::endl;
    std::cout << "Calling " << name << " (rep=" << maxrep << ")" << std::endl;
    time_t t0 = time( 0 );
    for( int i = 0; i < maxrep; ++i )
    {
        fn( S_sz, S_zi2E_FCt, S_FC_sv, GetM_S_zi2E_FCt, GetN_S_zi2E_FCt, GetM_S_FC_sv, GetN_S_FC_sv, DE );
    }
    cout << "Time: " << ( time(0) - t0 ) << endl;
}

void runCompare( int size, string ref, double *refData, string sample, double *sampleData )
{
    double diff = 0.0;
    for(int i=0; i<size; i++)
    {
        diff += fabs( refData[i] - sampleData[i]); 
    }
    std::cout << "Diff between " << ref << " and " << sample << " is " << diff << std::endl;
}

int main(int argc, char* argv[])
{

    const int rep = (argc==2) ? atoi( argv[1] ) : 1;

// #include "defvars-short.h"
#include "defvars.h"

    struct
    {
        string name;
        void (*fn)( double *S_sz,
                    double *S_zi2E_FCt,
                    double *S_FC_sv,
                    int GetM_S_zi2E_FCt,
                    int GetN_S_zi2E_FCt,
                    int GetM_S_FC_sv,
                    int GetN_S_FC_sv,
                    double *DE );
        double *data;
    } table[] =
    {
        { "Original C", original, (double *)malloc( sizeof(double) * DE_SZ ) },
        { "Alpha",      alpha,    (double *)malloc( sizeof(double) * DE_SZ ) },
        { "Beta",       beta,     (double *)malloc( sizeof(double) * DE_SZ ) },
        { "Gamma",      gamma,    (double *)malloc( sizeof(double) * DE_SZ ) },
        //{ "Delta",      delta,    (double *)malloc( sizeof(double) * DE_SZ ) },
        //{ "Epsilon",    epsilon,  (double *)malloc( sizeof(double) * DE_SZ ) },
    };

    for( int i = 0; i < ((sizeof(table))/(sizeof(table[0]))); ++i )
    {
        runTrial( table[i].name, table[i].fn, S_sz, S_zi2E_FCt, S_FC_sv, GetM_S_zi2E_FCt, GetN_S_zi2E_FCt, GetM_S_FC_sv, GetN_S_FC_sv, table[i].data, rep );
        runCompare( DE_SZ, "Matlab", S_DE, table[i].name, table[i].data );
        if( i != 0 )
        {
            runCompare( DE_SZ, table[0].name, table[0].data, table[i].name, table[i].data );
        }
    }

    for( int i = 0; i < ((sizeof(table))/(sizeof(table[0]))); ++i )
    {
        free( table[i].data );
    }

    return 0;
}
