#include <cmath>
#include <matrix.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *cth2 = mxGetPr(prhs[0]);
    double *E = mxGetPr(prhs[1]);
    double dth = mxGetScalar(prhs[2]);

    mwSize m = mxGetDimensions(prhs[0])[0];
    mwSize n = mxGetDimensions(prhs[0])[1];
    mwSize o = mxGetDimensions(prhs[0])[2];
    mwSize p = mxGetDimensions(prhs[0])[3];

    mwSize dims[] = {m, n, o, p};
    plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double *short_array = mxGetPr(plhs[0]);

    // Perform the computation
    for (mwSize tj = 0; tj < p; ++tj) {
        for (mwSize z = 0; z < o; ++z) {
            for (mwSize y = 0; y < n; ++y) {
                for (mwSize x = 0; x < m; ++x) {
                    // Calculate the offset for short_array
                    mwSize offset = x + y * m + z * m * n;
                    double E_value = E[offset];              // Precompute E value
                    double cth2_value = cth2[offset];        // Precompute cth2 value

                    double sum_result = 0.0;

                    for (mwSize i = 0; i < p; ++i) {
                        mwSize idx = (i - tj + p) % p;
                        sum_result += E_value * cth2_value;
                    }

                    short_array[offset + tj * m * n * o] = sum_result * dth;
                }
            }
        }
    }
}
