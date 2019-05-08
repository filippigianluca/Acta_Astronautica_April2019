/******************************************************************************
 *                     Mex Gateway routine for interp1q_c                     *
 *                                                                            *
 *                                                      Federico Zuiani, 2011 *
 ******************************************************************************/

/* To compile the mex file, use the following command:
mex mexinterp1q_c.c math_utils.c -output interp1q_c -DMEXCOMPILE
*/

#ifndef MEXINTERP1Q_C_H
#define MEXINTERP1Q_C_H

/* Include */
#include "mex.h"
#include "math_utils.h"

/* Prototype */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif /* MEXINTERP1Q_C_H */
