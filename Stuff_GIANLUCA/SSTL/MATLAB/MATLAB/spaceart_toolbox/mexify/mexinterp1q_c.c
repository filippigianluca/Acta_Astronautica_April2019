/******************************************************************************
 *                     Mex Gateway routine for interp1q_c                     *
 *                                                                            *
 *                                                      Federico Zuiani, 2011 *
 ******************************************************************************/

#include "mexinterp1q_c.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* Declaration */
    double *x_data;		/* Arrays to hold the x data set (inputs) */
	double *y_data;		/* Arrays to hold the y data set (inputs) */
	double x;			/* Scalar to hold the x coordinate at which to interpolate (Input)*/
	int n_data, flag;
	double *out;
    mwSize nr, nc;      /* number of rows and columns for the input vectors*/
    
    
    /* Check the number of inputs and outputs */
    if(nrhs != 3)
        mexErrMsgTxt("interp1q_c accepts 3 inputs only!");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. interp1q_c can return 1 output maximum!");
    
    /* Get values of inputs */
    nr = mxGetM(prhs[0]);
    nc = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
                             || !(nr>1 && nc==1)
                             || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("The first input of interp1q_c must be a column vector of at least 2 reals!");
	n_data=nr;
    x_data= mxGetPr(prhs[0]);
	
	nr = mxGetM(prhs[1]);
    nc = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
                             || !(nr==n_data && nc==1)
                             || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("The second input of interp1q_c must be a column vector with the length as the first!");
    y_data= mxGetPr(prhs[1]);
	
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
                             || mxGetM(prhs[2])!=1
                             || mxGetN(prhs[2])!=1
                             || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("The third input of interp1q_c must be a real scalar!");
    x = mxGetScalar(prhs[2]);
	
    /*Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* Pointer to first output */
    out = mxGetPr(plhs[0]);
    
    /* Compute the outputs */
    flag=interp1q_c(x_data, y_data, x, n_data, out);
	
	if (flag==1)
		mexErrMsgTxt("x must be comprised in the bounds defined by x_data!");

    /* Exit the function */
    return;
}
