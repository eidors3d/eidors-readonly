#include "mex.h"
#include "matrix.h"

/*
function DE = next_make_c(sz, zi2E_FCt, FC_sv);
   DE= zeros(sz(1),sz(2),sz(3) );
   for k= 1: sz(3);
       idx= (sz(4)-1)*(k-1)+1 : (sz(4)-1)*k;
       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
       DE(:,:,k)= dq;
   end
*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x,*y;
  mwSize mrows,ncols;
  
  /* Check for proper number of arguments. */
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

  double * sz = mxGetPr(prhs[0]);
/*
  for (int k=0; k< 4; k++) {
     printf("v%d = %f\n", k, sz[k]);
  }
*/
  
  /* Create matrix for the return argument.
  plhs[0] = mxCreateDoubleMatrix(2,1, mxREAL);
*/
  mwSize dims[] = {sz[0],sz[1],sz[2]};
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double * DE = mxGetPr(plhs[0]);
printf("DE [%d x %d x %d]\n", dims[0],dims[1],dims[2]);
  
/*
   d1 = sz(4)-1;
   for k= 1: sz(3);
       idx= d1*(k-1) + (1:d1);
       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
       DE(:,:,k)= dq;
   end
*/
  int d1 = (int) sz[3] - 1;
  double * M1 = mxGetPr( prhs[1] ); // todo: check it is right size;
  double * M2 = mxGetPr( prhs[2] ); // todo: check it is right size;
  int R = mxGetM(prhs[1]); 
  int L = mxGetN(prhs[1]); 
printf("M1[%d x %d] M2[%d x %d] d1=%d\n",
       mxGetM(prhs[1]), mxGetN(prhs[1]), 
       mxGetM(prhs[2]), mxGetN(prhs[2]), d1);

        
  int E = (int) sz[0]; // assume sz[1] == sz[0]
  int E2= E*E;
  for (int l=0; l< sz[2]; l++) {
     int ofs = d1*l;
     for (int i=0; i< R; i++) {
     for (int j=0; j< R; j++) {
        double sum= 0.0;
        for (int k=0; k< d1; k++) {
/*
   DE(i,j,l) = sum_k M1(i, l*d1 + k) * M2(l*d1+k, j)
*/
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
   printf("fini\n");
}
