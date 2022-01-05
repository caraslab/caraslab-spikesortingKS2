#include "mex.h"

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    double *s1, *s2, *output, *edges, limit, delta;
    int n2,n1,j,i,istop,k,nedges;
    
    s1 = mxGetPr(prhs[0]);
    s2 = mxGetPr(prhs[1]);
    edges = mxGetPr(prhs[2]);
    n1 = mxGetM(prhs[0])*mxGetN(prhs[0]);
    n2 = mxGetM(prhs[1])*mxGetN(prhs[1]);
    nedges = mxGetM(prhs[2])*mxGetN(prhs[2]);
    limit = edges[nedges-1];
    
    plhs[0] = mxCreateDoubleMatrix(nedges-1,1,mxREAL);
    output = mxGetPr(plhs[0]);
    
    j = 0;
    i = 0;
    istop = 0;
    while(j<n2) {
        // Advance i until we are within range of j
        while(i<n1&&s2[j]-s1[i]>limit) {
            i++;
        }
        // Record the earliest position within range of j
        istop = i;
        // Advance i and record ISIs until we are past j
        while(i<n1&&s2[j]-s1[i]>0) {
            delta = s2[j]-s1[i];
            k=0;
            while(k<nedges-1 && edges[k+1]<delta) {
                k++;
            }
            output[k]++;
            i++;
        }
        // Send i back to the earliest position within range of j
        i = istop;
        j++;
    }
    
}