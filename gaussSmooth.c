/*
  Caluclate Gaussian smoothing of image (with smoothed differential and LOG)
  2018/7/20
  Usage:
    type = 0:
      smoothedImg   = gaussSmooth(inImg, type, sigma, P, extType)
    type = 1:
      [diffX diffY] = gaussSmooth(inImg, type, sigma, P, extType)
    type = 2:
      logImg        = gaussSmooth(inImg, type, sigma, P, extType)
  Paramter
    inImg   : input image
    type    : 0 blur, 1 diff (two output), 2 : LOG
    sigma   : sigma
    P       : 2 or 4 or 6 : order of Fourier series
    extType :  0: zero extension,  1: the value of edge is used for extension

  Output: 
    smoothImg: Gaussian smoothed image
    diffX    : Gaussian smoothed x-directional differential image
    diffY    : Gaussian smoothed x-directional differential image
    lapImg   : Laplacian of Gaussian smoothed (LOG) image
*/

#include "mex.h"
#include <math.h>
#define PI 3.14159265

/* filter coefficients P = 2 sgimaN = 1.1 * [0.9 0.05 1.0 1.05, 1.1] */
double coef02[] = {1.5890479993e-01, -1.9546463415e-01, 4.4433708754e-02, 1.5872508026e-01, -1.8517758499e-01, 3.5193829762e-02, 1.5846315202e-01, -1.7508079293e-01, 2.7323762984e-02, 1.5810609593e-01, -1.6524776412e-01, 2.0681709055e-02, 1.5764419663e-01, -1.5573913229e-01, 1.5122716379e-02};
double coef12[] = {-1.9544398402e-01, 8.8908717758e-02, -1.8515950340e-01, 7.0423822709e-02, -1.7506270969e-01, 5.4683692457e-02, -1.6522793689e-01, 4.1403072576e-02, -1.5571649590e-01, 3.0290705538e-02};
double coef22[] = {-2.6365863313e-03, 1.9017087850e-01, -1.8309033963e-01, -3.8005988197e-03, 1.7755840140e-01, -1.4844846049e-01, -5.0840517439e-03, 1.6489473361e-01, -1.1953497880e-01, -6.4250671203e-03, 1.5237796345e-01, -9.5655636200e-02, -7.7635909137e-03, 1.4018950834e-01, -7.6107815825e-02};

/* filter coefficients P = 4 sgimaN = 0.8 * [0.9 0.05 1.0 1.05, 1.1] */
double coef04[] = {1.5915197146e-01, -2.4563562311e-01, 1.1286311162e-01, -3.0888860214e-02, 5.0271861757e-03, 1.5914871818e-01, -2.3847783114e-01, 1.0025453514e-01, -2.3671129824e-02, 3.1248009011e-03, 1.5914081930e-01, -2.3116780826e-01, 8.8476942331e-02, -1.7890179087e-02, 1.8836364529e-03, 1.5912513408e-01, -2.2373944516e-01, 7.7568327529e-02, -1.3344924015e-02, 1.0890989979e-03, 1.5909736625e-01, -2.1622839417e-01, 6.7544370560e-02, -9.8403987857e-03, 5.8391830816e-04};
double coef14[] = {-2.4563375548e-01, 2.2572995849e-01, -9.2660977760e-02, 2.0116215211e-02, -2.3847674082e-01, 2.0051125093e-01, -7.1010118505e-02, 1.2503564895e-02, -2.3116693143e-01, 1.7695563833e-01, -5.3667906760e-02, 7.5380531476e-03, -2.2373838915e-01, 1.5513876708e-01, -4.0031604007e-02, 4.3606200458e-03, -2.1622684967e-01, 1.3509183012e-01, -2.9516562861e-02, 2.3418512347e-03};

double coef24[] = {-5.8486698162e-05, 2.4551678405e-01, -4.5157688250e-01, 2.7786597760e-01, -8.0581802738e-02, -1.6945020880e-04, 2.3813784485e-01, -4.0136138452e-01, 2.1269149507e-01, -5.0353088930e-02, -3.4905652299e-04, 2.3046882714e-01, -3.5460935466e-01, 1.6030568607e-01, -3.0850185473e-02, -6.2259675783e-04, 2.2249321114e-01, -3.1152266568e-01, 1.1884975799e-01, -1.8687425716e-02, -1.0086574679e-03, 2.1420955983e-01, -2.7220087480e-01, 8.6532599483e-02, -1.1384318381e-02};

/* filter coefficients P = 6 sgimaN = 0.7 * [0.9 0.05 1.0 1.05, 1.1] */
double coef06[] = {1.5915480802e-01, -2.6101461924e-01, 1.4391510624e-01, -5.3355659590e-02, 1.3300581149e-02, -2.2297154933e-03, 2.5107856342e-04, 1.5915455493e-01, -2.5516640193e-01, 1.3144254476e-01, -4.3511641759e-02, 9.2549485793e-03, -1.2657346498e-03, 1.1063870838e-04, 1.5915377148e-01, -2.4914489111e-01, 1.1946305955e-01, -3.5095823000e-02, 6.3138382801e-03, -6.9786114192e-04, 4.5673259016e-05, 1.5915183612e-01, -2.4296997730e-01, 1.0804230050e-01, -2.7999754721e-02, 4.2213838519e-03, -3.7551843901e-04, 1.5748785727e-05, 1.5914766619e-01, -2.3666253739e-01, 9.7231115375e-02, -2.2098290403e-02, 2.7627467159e-03, -2.0080104338e-04, 3.3186313163e-07};

double coef16[] = {-2.6101454448e-01, 2.8783036200e-01, -1.6006675449e-01, 5.3202623640e-02, -1.1148203661e-02, 1.5069199468e-03, -2.5516636049e-01, 2.6288517241e-01, -1.3053480093e-01, 3.7019960107e-02, -6.3284660116e-03, 6.6408093577e-04, -2.4914483525e-01, 2.3892623081e-01, -1.0528730143e-01, 2.5255576544e-02, -3.4890264302e-03, 2.7437469024e-04, -2.4296986453e-01, 2.1608482654e-01, -8.3998925847e-02, 1.6885986494e-02, -1.8770283353e-03, 9.5169347827e-05, -2.3666231416e-01, 1.9446267722e-01, -6.6294201510e-02, 1.1051879798e-02, -1.0028890470e-03, 3.3305861405e-06};

  double coef26[] = {-4.8526928608e-06, 2.6100483926e-01, -5.7567042875e-01, 4.8019055951e-01, -2.1282019739e-01, 5.5731316905e-02, -9.0512193281e-03, -1.9174786423e-05, 2.5512801140e-01, -5.2580869246e-01, 3.9156605759e-01, -1.4811818225e-01, 3.1603992605e-02, -4.0228177347e-03, -4.9872068291e-05, 2.4904509235e-01, -4.7795220081e-01, 3.1576217126e-01, -1.0112203057e-01, 1.7345418863e-02, -1.7459478551e-03, -1.0990253007e-04, 2.4275006219e-01, -4.3238944727e-01, 2.5177699695e-01, -6.7763707535e-02, 9.1654045904e-03, -7.9072326202e-04, -2.1478494543e-04, 2.3623274959e-01, -3.8935490302e-01, 1.9845308256e-01, -4.4637003891e-02, 4.5850084572e-03, -4.4936171980e-04};

void slideFilter(int type, double *vect, int L, int K, int extType, int P, double *cosPara, double *sinPara, double *filter1, double *filter2, double *outVect1, double *outVect2);
void calFilterVal(int type, double *cosTmp, double *sinTmp, double *cosPara, double *sinPara, double add, double subtract, int P, double *filter1, double *filter2, double *outVect1, double *outVect2, int pos);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *inImg, sigma;
  int     type, P, extType, M, N;

  double *outImg1, *outImg2;

  double *xBlurImg, *xTranImg, filter1[7], filter2[7], *outVect1, *outVect2, sigmaPN, sigmaN, sigmaDiffN, interCoefL, interCoefU;
  double *coef0, *coef1, *coef2, coefN; 
  double paraTheta, cosPara[6], sinPara[6];
  int    K, m, n, pos, maxMN, p, pd, P1, interType, interPosC, interPosS;

  if (nrhs != 5 ) {
    mexErrMsgIdAndTxt("In0", "Number of parameters are not correct: gaussSmooth(inImg, type, sigma, P, extType)");
  }
  /* For input */
  inImg   = mxGetPr(prhs[0]);
  M       = (int) mxGetM(prhs[0]);
  N       = (int) mxGetN(prhs[0]);
  type    = (int) mxGetScalar(prhs[1]);
  sigma   = mxGetScalar(prhs[2]);
  P       = (int) mxGetScalar(prhs[3]);
  extType = (int) mxGetScalar(prhs[4]);

  /* For output */
  if (type == 0 || type ==2) {
    if (nlhs != 1 ) {
      mexErrMsgIdAndTxt("Out0", "Number of output is not correct: outImg = gaussSmooth(inImg, type, sigma, P, extType)");
    } else {
      plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
      outImg1 = mxGetPr(plhs[0]);
    }
  } else if (type == 1) {
    if (nlhs != 2 ) {
      mexErrMsgIdAndTxt("Out1", "Number of output is not correct: [diffX diffY] = gaussSmooth(inImg, type, sigma, P, extType)");
    } else {
      plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(M, N, mxREAL);
      outImg1 = mxGetPr(plhs[0]);
      outImg2 = mxGetPr(plhs[1]);
    }
  } else {
    mexErrMsgIdAndTxt("Out2", "type should be 0 or 1 or 2.");
  }
  
  /* For work */
  P1 = P + 1;
  maxMN = (M > N) ? M : N;
  
  xBlurImg = (double *) mxMalloc(sizeof(double) * M * N);
  outVect1 = (double *) mxMalloc(sizeof(double) * maxMN);
  if (type != 0) {
    xTranImg = (double *) mxMalloc(sizeof(double) * M * N);
    outVect2 = (double *) mxMalloc(sizeof(double) * maxMN);
  }

  /* Fix K and filter coefficients */
  switch(P) {
  case 2:
    sigmaPN = 1.1; coef0 = coef02; coef1 = coef12; coef2 = coef22; break;
  case 4:
    sigmaPN = 0.8; coef0 = coef04; coef1 = coef14; coef2 = coef24; break;
  case 6:
    sigmaPN = 0.7; coef0 = coef06; coef1 = coef16; coef2 = coef26; break;
  }

  /* Decide K */ 
  K = (int) (sigma * PI / sigmaPN + 0.5);
  // printf("K = %d\n", K);

  /* Interporation of sigma for Fourier coefficients */
  coefN      = PI / K;
  sigmaN     = sigma * coefN;
  sigmaDiffN = (sigmaN - sigmaPN) / sigmaPN;

  if (sigmaDiffN < -0.5)     interType = 0;
  else if (sigmaDiffN < 0.0) interType = 1;
  else if (sigmaDiffN < 0.5) interType = 2;
  else                       interType = 3;

  interCoefU = (sigmaDiffN - 0.05 * (interType - 2)) / 0.05;
  interCoefL = 1.0 - interCoefU;
  interCoefU *= coefN;
  interCoefL *= coefN;
  interPosC  = P1 * interType;
  interPosS  = P  * interType;
  //printf("sigmaDiffN = %lf type = %d  L %lf U %lf  \n",  sigmaDiffN, interType, interCoefL, interCoefU);
  switch(type) {
    case 0:
      for (p = 0 ; p <= P ; ++p) {
	filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p];
      }
      break;
    case 1:
      p = 0;
      filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p];
      for (p = 1 ; p <= P ; ++p) {
	pd = p - 1;
	filter1[p]  = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p];
	filter2[pd] = interCoefL * coef1[interPosS + pd] + interCoefU * coef1[interPosS + P + pd];
      }
      break;
    case 2:
      for (p = 0 ; p <= P ; ++p) {
	filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p];
	filter2[p] = interCoefL * coef2[interPosC + p] + interCoefU * coef2[interPosC + P1 + p];
      }
      break;
  }
  paraTheta =  PI / K;
  for(p = 1 ; p <= P ; ++p) {
    pd = p - 1;
    cosPara[pd] = cos(p * paraTheta);
    sinPara[pd] = sin(p * paraTheta);
  }
  
  /* Vertical transformation */
  for(n = 0 ; n < N ; ++n) {
    switch(type) {
    case 0:
      slideFilter(1, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = n;
      for(m = 0 ; m < M ; ++m) {
	xBlurImg[pos] = outVect1[m];
	pos += N;
      }
      break;
    case 1:
      slideFilter(4, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = n;
      for(m = 0 ; m < M ; ++m) {
	xBlurImg[pos] = outVect1[m];
	xTranImg[pos] = outVect2[m];
	pos += N;
      }
      break;
    case 2:
      slideFilter(3, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = n;
      for(m = 0 ; m < M ; ++m) {
	xBlurImg[pos] = outVect1[m];
	xTranImg[pos] = outVect2[m];
	pos += N;
      }
      break;
    }
  }

  /* Horizontal transformation */
  for(m = 0 ; m < M ; ++m) {
    switch(type) {
    case 0:
      slideFilter(1, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = m;
      for(n = 0 ; n < N ; ++n) {
	outImg1[pos] = outVect1[n];
	pos += M;
      }
      break;
    case 1:
      slideFilter(2, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, filter2, filter2, outVect1, outVect2);
      pos = m;
      for(n = 0 ; n < N ; ++n) {
	outImg1[pos] = outVect1[n];
	pos += M;
      }
      slideFilter(1, & (xTranImg[m * N]), N, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = m;
      for(n = 0 ; n < N ; ++n) {
	outImg2[pos] = outVect1[n];
	pos += M;
      }
      break;
    case 2:
      slideFilter(1, & (xTranImg[m * N]), N, K, extType, P, cosPara, sinPara, filter1, filter2, outVect1, outVect2);
      pos = m;
      for(n = 0 ; n < N ; ++n) {
	outImg1[pos] = outVect1[n];
	pos += M;
      }
      slideFilter(1, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, filter2, filter2, outVect1, outVect2);
      pos = m;
      for(n = 0 ; n < N ; ++n) {
	outImg1[pos] += outVect1[n];
	pos += M;
      }
      break;
    }
  }

  /* Free memories */
  mxFree(xBlurImg);
  mxFree(outVect1);
  if (type != 0) {
    mxFree(xTranImg);
    mxFree(outVect2);
  }
  return;
}

/* type : 0: cos, 1 sin, 2: cos cos, 3: cos sin */
void slideFilter(int type, double *vect, int L, int K, int extType, int P, double *cosPara, double *sinPara, double *filter1, double *filter2, double *outVect1, double *outVect2) {
  double cosTmp[7], sinTmp[6], extS, extE;
  int    p, pd, pos, K2, posK = 0, posK2 = 0;
  
/* For initial position */
  if (extType == 0) {
    extS = extE = 0.0;
  } else {
    extS = vect[0]; extE = vect[L - 1]; 
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
      cosTmp[p]  =   extS;
      sinTmp[pd] = - extS * sinPara[pd] / (1 - cosPara[pd]);
    }
  }

  if (L >= K2) {
    /*      |     L     |
     | K . K |                         */
    for (pos = 0 ; pos < K ; ++pos) {
      calFilterVal(0, cosTmp, sinTmp, cosPara, sinPara, vect[pos], 0.0, P, filter1, filter2, outVect1, outVect2, 0);
    }

    /*      |     L     |
          | K . K |                         */
    for (pos = K ; pos < K2 ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, vect[pos], extS, P, filter1, filter2, outVect1, outVect2, posK);
      ++posK;
    }
    /*      |     L     |
              | K . K |                         */
    for (pos = K2 ; pos < L ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, vect[pos], vect[posK2], P, filter1, filter2, outVect1, outVect2, posK);
      ++posK; ++posK2;
    }
    /*      |     L     |
                  | K . K |                         */
    for (pos = 0 ; pos < K ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, extE, vect[posK2], P, filter1, filter2, outVect1, outVect2, posK);
      ++posK; ++posK2;
    }
  } else if (L >= K) {
    /*      |  L  |
      | K . K |                         */
    for (pos = 0 ; pos < K ; ++pos) {
      calFilterVal(0, cosTmp, sinTmp, cosPara, sinPara, vect[pos], 0.0, P, filter1, filter2, outVect1, outVect2, 0);
    }
    /*      |  L  |
         | K . K |                         */
    for (pos = K ; pos < L ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, vect[pos], extS, P, filter1, filter2, outVect1, outVect2, posK);
      ++posK;
    }
    /*      |  L  |
           | K . K |                         */
    for (pos = L ; pos < K2 ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, extE, extS, P, filter1, filter2, outVect1, outVect2, posK);
      ++posK; ++posK2;
    }
    /*      |  L  |
             | K . K |                         */
    for (pos = K2 ; pos < L + K ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, extE, vect[posK2], P, filter1, filter2, outVect1, outVect2, posK);
      ++posK; ++posK2;
    }
  } else { /* L < K */
    /*          |  L  |
    |   K   .   K   |                         */
    for (pos = 0 ; pos < L ; ++pos) {
      calFilterVal(0, cosTmp, sinTmp, cosPara, sinPara, vect[pos], 0.0, P, filter1, filter2, outVect1, outVect2, 0);
    }
    /*          |  L  |
       |   K   .   K   |                         */
    for (pos = L ; pos < K ; ++pos) {
      calFilterVal(0, cosTmp, sinTmp, cosPara, sinPara, extE, 0.0, P, filter1, filter2, outVect1, outVect2, 0);
    }
    /*          |  L  |
           |   K   .   K   |                         */
    for (pos = K ; pos < L + K ; ++pos) {
      calFilterVal(type, cosTmp, sinTmp, cosPara, sinPara, extE, extS, P, filter1, filter2, outVect1, outVect2, posK);
      ++posK;
    }
  }
  return;
}

 void calFilterVal(int type, double *cosTmp, double *sinTmp, double *cosPara, double *sinPara, double add, double subtract, int P, double *filter1, double *filter2, double *outVect1, double *outVect2, int pos) {
  int p, pd;
  double cosTmp0, sinTmp0, val1, val2;

  switch (type) {
  case 0: /* no output */
    cosTmp[0] += add;
    cosTmp[0] -= subtract;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1; cosTmp0 = cosTmp[p]; sinTmp0 = sinTmp[pd];
      cosTmp[p]  =   cosPara[pd] * cosTmp0 + sinPara[pd] * sinTmp0 + add;
      sinTmp[pd] = - sinPara[pd] * cosTmp0 + cosPara[pd] * sinTmp0;
      cosTmp[p] -= subtract;
    }
    break;
  case 1: /* cos */
    cosTmp[0] += add;
    val1 = filter1[0] * cosTmp[0];
    cosTmp[0] -= subtract;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1; cosTmp0 = cosTmp[p]; sinTmp0 = sinTmp[pd];
      cosTmp[p]  =   cosPara[pd] * cosTmp0 + sinPara[pd] * sinTmp0 + add;
      sinTmp[pd] = - sinPara[pd] * cosTmp0 + cosPara[pd] * sinTmp0;
      val1 += filter1[p] * cosTmp[p];
      cosTmp[p] -= subtract;
    }
    outVect1[pos] = val1;
    break;
  case 2: /* sin */
    cosTmp[0] += add;
    val1 = 0.0;
    cosTmp[0] -= subtract;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1; cosTmp0 = cosTmp[p]; sinTmp0 = sinTmp[pd];
      cosTmp[p]  =   cosPara[pd] * cosTmp0 + sinPara[pd] * sinTmp0 + add;
      sinTmp[pd] = - sinPara[pd] * cosTmp0 + cosPara[pd] * sinTmp0;
      val1 += filter1[pd] * sinTmp[pd];
      cosTmp[p] -= subtract;
    }
    outVect1[pos] = val1;
    break;
  case 3: /* cos cos */
    cosTmp[0] += add;
    val1 = filter1[0] * cosTmp[0];
    val2 = filter2[0] * cosTmp[0];
    cosTmp[0] -= subtract;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1; cosTmp0 = cosTmp[p]; sinTmp0 = sinTmp[pd];
      cosTmp[p]  =   cosPara[pd] * cosTmp0 + sinPara[pd] * sinTmp0 + add;
      sinTmp[pd] = - sinPara[pd] * cosTmp0 + cosPara[pd] * sinTmp0;
      val1 += filter1[p] * cosTmp[p];
      val2 += filter2[p] * cosTmp[p];
      cosTmp[p] -= subtract;
    }
    outVect1[pos] = val1;
    outVect2[pos] = val2;
    break;
  case 4: /* cos sin */
    cosTmp[0] += add;
    val1 = filter1[0] * cosTmp[0];
    val2 = 0.0;
    cosTmp[0] -= subtract;
    for (p = 1 ; p <= P ; ++p) {
      pd = p - 1; cosTmp0 = cosTmp[p]; sinTmp0 = sinTmp[pd];
      cosTmp[p]  =   cosPara[pd] * cosTmp0 + sinPara[pd] * sinTmp0 + add;
      sinTmp[pd] = - sinPara[pd] * cosTmp0 + cosPara[pd] * sinTmp0;
      val1 += filter1[p]  * cosTmp[p];
      val2 += filter2[pd] * sinTmp[pd];
      cosTmp[p] -= subtract;
    }
    outVect1[pos] = val1;
    outVect2[pos] = val2;
    break;
  }
}
