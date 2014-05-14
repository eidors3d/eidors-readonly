#include <stdio.h>
#include <math.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

int main() {

  double * sz = S_sz;
  
  int dims[] = {sz[0],sz[1],sz[2]};
  double DE [DE_SZ]; // test eq to S_DE
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
  for (int l=0; l< sz[2]; l++) {
     int ofs = d1*l;
     for (int i=0; i< R; i++) {
     for (int j=0; j< R; j++) {
        double sum= 0.0;
        for (int k=0; k< d1; k++) {
// DE(i,j,l) = sum_k M1(i, l*d1 + k) * M2(l*d1+k, j)
          int m1o=  i + R*(ofs + k);
          int m2o=  ofs + k + L*j;
if (i<=1 && l==0 && j<=1)
printf("i=%d j=%d k=%d m1[%d]=%f m2[%d]=%f\n",i,j,k, m1o, M1[m1o], m2o, M2[m2o]);
          sum+= M1[m1o] * M2[m2o];
        }
        DE[i + E*j + E2*l] = sum;
     }
     }
  }

  double diff = 0.0;
  for(int i=0; i<DE_SZ; i++) {
    diff += fabs( DE[i] - S_DE[i]); 
  }

  printf("fini. Diff (should be 0) is %g\n",diff);
}
