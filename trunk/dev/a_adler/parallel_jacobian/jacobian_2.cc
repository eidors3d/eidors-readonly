// Jacobian.cpp : Defines the entry point for the console application.
//
#include <omp.h>
#include <iostream>

using std::cout;
using std::endl;

#include <sstream>
using std::istringstream;

#include <stdio.h>
#include <math.h>
#include <time.h>

void original( double *DE, const int maxrep )
{
#include "defvars.h"

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
  printf("DE [%d x %d x %d]\n", dims[0],dims[1],dims[2]);
  
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
printf("M1[%d x %d] M2[%d x %d] d1=%d\n",
       GetM_S_zi2E_FCt, GetN_S_zi2E_FCt,
       GetM_S_FC_sv,    GetN_S_FC_sv,     d1);

  int E = (int) sz[0]; // assume sz[1] == sz[0]
  int E2= E*E;
for( int rep = 0; rep < maxrep; ++rep )
{
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
}

  //double diff = 0.0;
  //for(int i=0; i<DE_SZ; i++) {
    //diff += fabs( DE[i] - S_DE[i]); 
  //}

  //printf("fini. Diff (should be 0) is %g\n",diff);
}

void alpha( double *DE, const int maxrep )
{
#include "defvars.h"

    const double * sz = S_sz;
    const int dims[] = {sz[0],sz[1],sz[2]};

    cout << "DE [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << endl;

    const int d1 = (int) sz[3] - 1;
    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    const int R = GetM_S_zi2E_FCt;
    const int L = GetN_S_zi2E_FCt;

    cout << "M1[" << GetM_S_zi2E_FCt << " x " << GetN_S_zi2E_FCt << "]"
        << " M2[" << GetM_S_FC_sv << " x " << GetN_S_FC_sv << "] d1=" << d1 << endl;

    const int E = (int) sz[0]; // assume sz[1] == sz[0]
    const int E2= E*E;
    const int Rd1 = R * d1;

    long long mult = 0;
    for( int rep = 0; rep < maxrep; ++rep )
    {
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
                        ++mult;
                    }
                    DE[i + Ej + E2l] = sum;
                }
            }
        }
    }

    cout << "Mult: " << mult << endl;
}

void beta( double *DE, const int maxrep )
{
#include "defvars.h"

    const double * sz = S_sz;
    const int dims[] = {sz[0],sz[1],sz[2]};

    cout << "Found " << omp_get_num_procs() << " processors" << endl;
    cout << "DE [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << endl;

    const int d1 = (int) sz[3] - 1;
    const double * M1 = S_zi2E_FCt;
    const double * M2 = S_FC_sv;
    const int R = GetM_S_zi2E_FCt;
    const int L = GetN_S_zi2E_FCt;

    cout << "M1[" << GetM_S_zi2E_FCt << " x " << GetN_S_zi2E_FCt << "]"
        << " M2[" << GetM_S_FC_sv << " x " << GetN_S_FC_sv << "] d1=" << d1 << endl;

    const int E = (int) sz[0]; // assume sz[1] == sz[0]
    const int E2= E*E;
    const int Rd1 = R * d1;

    long long mult = 0;
    for( int rep = 0; rep < maxrep; ++rep )
    {
        int l;
        int maxL = sz[2];
#pragma omp parallel for reduction(+:mult)
        for (l=0; l < maxL; l++ )
        {
            const int ofs = l * d1;
            const int E2l = l * E2;
            const int Rofs = l * Rd1;

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
                        ++mult;
                    }
                    DE[i + Ej + E2l] = sum;
                }
            }
        }
    }

    cout << "Mult: " << mult << endl;
}
int atoi( char *str )
{
   int ret;

   std::string s( str );
   istringstream istr( s );

   istr >> ret;

   return ret;
}

int main(int argc, char* argv[])
{
    #include "defvars.h"

    const int rep = (argc==2) ? atoi( argv[1] ) : 1;

    double DE [DE_SZ];

    std::cout << "Calling original(rep=" << rep << ")" << std::endl;
    clock_t t0 = clock();
    original( DE, rep );
    cout << "Time: " << ((double) clock() - t0 )/CLOCKS_PER_SEC << endl;

    double diff = 0.0;
    for(int i=0; i<DE_SZ; i++)
    {
        diff += fabs( DE[i] - S_DE[i]); 
    }
    std::cout << "Sum fabs(delta): " << diff << std::endl;

    double DEalpha [DE_SZ];

    std::cout << std::endl;
    std::cout << "Calling alpha(rep=" << rep << ")" << std::endl;
    t0 = clock( );
    alpha( DEalpha, rep );
    cout << "Time: " << ((double) clock() - t0 )/CLOCKS_PER_SEC << endl;

    diff = 0.0;
    for(int i=0; i<DE_SZ; i++)
    {
        diff += fabs( DEalpha[i] - S_DE[i]); 
    }
    std::cout << "Sum fabs(delta): " << diff << std::endl;

    double DEbeta [DE_SZ];

    std::cout << std::endl;
    std::cout << "Calling beta(rep=" << rep << ")" << std::endl;
    t0 = clock( );
    beta( DEbeta, rep );
    cout << "Time: " << ((double) clock() - t0 )/CLOCKS_PER_SEC << endl;

    diff = 0.0;
    for(int i=0; i<DE_SZ; i++)
    {
        diff += fabs( DEbeta[i] - S_DE[i]); 
    }
    std::cout << "Sum fabs(delta): " << diff << std::endl;

    std::cout << std::endl;
    std::cout << "Comparing original vs alpha" << std::endl;
    diff = 0.0;
    for(int i=0; i<DE_SZ; i++)
    {
        diff += fabs( DE[i] - DEalpha[i]); 
    }
    std::cout << "Sum fabs(delta alpha): " << diff << std::endl;

    std::cout << std::endl;
    std::cout << "Comparing original vs beta" << std::endl;
    diff = 0.0;
    for(int i=0; i<DE_SZ; i++)
    {
        diff += fabs( DE[i] - DEbeta[i]); 
    }
    std::cout << "Sum fabs(delta beta): " << diff << std::endl;

    return 0;
}
