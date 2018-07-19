/*
  Caluclate integral image of sin and cos
  2018/7/17 
  1 D integral by DC cos sin 
  Usage:
  [cosInte sinInte] = scInte(extVect, K, extType, cosPara, sinPara)
*/
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *vect, *cosPara, *sinPara;
  int     type, ncols, K, extType, P;

  double *cosInte, *sinInte;

  int p, pd, pos, K2, cosPos = 0, sinPos = 0, posK2 = 0;
  double cosTmp[7], sinTmp[6], cosTmp0, sinTmp0, extS, extE;
  
  if (nrhs != 5 ) {
    mexErrMsgIdAndTxt("123", "Number of parameters are not correct: [cosInte sinInte] = scInte(vect, K, paraCos, paraSin");
  }

  /* For input */
  vect    = mxGetPr(prhs[0]);
  ncols   = (int) mxGetN(prhs[0]);
  K       = (int) mxGetScalar(prhs[1]);
  extType = (int) mxGetScalar(prhs[2]);
  cosPara = mxGetPr(prhs[3]);
  sinPara = mxGetPr(prhs[4]);
  P       = (int) mxGetN(prhs[4]);

  /* For output */
  plhs[0] = mxCreateDoubleMatrix(P + 1, ncols, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(P, ncols, mxREAL);
  cosInte = mxGetPr(plhs[0]);
  sinInte = mxGetPr(plhs[1]);
 
/* For initial position */
  if (extType == 0) {
    extS = extE = 0.0;
  } else {
    extS = vect[0]; extE = vect[ncols - 1]; 
  }
  K2 = K + K;
  cosTmp[0] = K * extS;
  for (p = 1 ; p <= P ; ++p) {
    pd = p - 1;
    if (p % K2 == 0) {
      cosTmp[p]  = K * extS;
      sinTmp[pd] = 0.0;
    } else if (p % 2 == 0) {
      cosTmp[p]  = 0.0;
      sinTmp[pd] = 0.0;
    } else {
      cosTmp[p]  =    extS;
      sinTmp[pd] =  - extS * sinPara[pd] / (1 - cosPara[pd]);
    }
  }

  for (pos = 0 ; pos < K ; ++pos) {
    cosTmp[0]  = cosTmp[0] + vect[pos];
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1;
      cosTmp0 = cosPara[pd] * cosTmp[p] + sinPara[pd] * sinTmp[pd] + vect[pos];
      sinTmp0 = - sinPara[pd] * cosTmp[p] + cosPara[pd] * sinTmp[pd];
      cosTmp[p] = cosTmp0; sinTmp[pd] = sinTmp0;
    }
  }
  for (pos = K ; pos < K2 ; ++pos) {
    cosInte[cosPos] = cosTmp[0] + vect[pos];
    cosTmp[0]       = cosInte[cosPos++] - extS;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1;
      cosInte[cosPos] =   cosPara[pd] * cosTmp[p] + sinPara[pd] * sinTmp[pd] + vect[pos];
      sinInte[sinPos] = -  sinPara[pd] * cosTmp[p] + cosPara[pd] * sinTmp[pd];
      cosTmp[p] = cosInte[cosPos++] - extS; sinTmp[pd]  = sinInte[sinPos++];
    }
  }
  for (pos = K2 ; pos < ncols ; ++pos) {
    cosInte[cosPos] = cosTmp[0] + vect[pos];
    cosTmp[0]       = cosInte[cosPos++] - vect[posK2];
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1;
      cosInte[cosPos] =   cosPara[pd] * cosTmp[p] + sinPara[pd] * sinTmp[pd] + vect[pos];
      sinInte[sinPos] = - sinPara[pd] * cosTmp[p] + cosPara[pd] * sinTmp[pd];
      cosTmp[p] = cosInte[cosPos++] - vect[posK2]; sinTmp[pd] = sinInte[sinPos++];
    }
    ++posK2;
  }
  for (pos = 0 ; pos < K ; ++pos) {
    cosInte[cosPos] = cosTmp[0] + extE;
    cosTmp[0]       = cosInte[cosPos++] - vect[posK2];
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1;
      cosInte[cosPos] =   cosPara[pd] * cosTmp[p] + sinPara[pd] * sinTmp[pd] + extE;
      sinInte[sinPos] = - sinPara[pd] * cosTmp[p] + cosPara[pd] * sinTmp[pd];
      cosTmp[p] = cosInte[cosPos++] - vect[posK2]; sinTmp[pd] = sinInte[sinPos++];
    }
    ++posK2;
  }
  return;
}
