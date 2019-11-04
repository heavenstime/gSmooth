/*
  Caluclate Gaussian smoothing of image (with smoothed differential and LOG (Laplacian of Gaussian))
  2018/7/20
  Usage:
    type = 0 or 10:
      smoothedImg   = gaussSmooth(inImg, type, sigma, P, extType)
    type = 1 or 11:
      diff          = gaussSmooth(inImg, type, sigma, P, extType)   (1D data)
      [diffX diffY] = gaussSmooth(inImg, type, sigma, P, extType)   (2D data)
    type = 2 or 12:
      logImg        = gaussSmooth(inImg, type, sigma, P, extType)
  Paramter
    inImg   : input image  (1D/2D data)
    type    : 0 blur, 1 diff (two output), 2 : LOG 
              (+10 : smoothed boundary. P = 1 and type = 12 is impossible. )
    sigma   : sigma
    P       : 2 or 4 or 6 : order of Fourier series
    extType :  0: zero extension,  1: the value of edge is used for extension

  Output: 
    smoothImg: Gaussian smoothed image (1D/2D data)
    diff     : Gaussian smoothed differential data for 1D data
    diffX    : Gaussian smoothed x-directional differential image (2D data)
    diffY    : Gaussian smoothed x-directional differential image (2D data)
    lapImg   : Laplacian of Gaussian smoothed (LOG) image  (1D/2D data)
*/

#include "mex.h"
#include <math.h>
#define PI 3.14159265

void decideCoef(int P, int type, double *sigmaPN, double **coef0, double **coef1, double **coef2, double *coefTypeN);

void slideFilterExt(int type, double *vect, int L, int K, int extType, int P, double *cosPara, double *sinPara, double *cosInitCoef, double *sinInitCoef, double *filter1, double *filter2, double *outVect1, double *outVect2, int outStep, double *lineExt);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *inImg, sigma;
  int     type, typeProc, P, extType, M, N, H;

  double *outImg1, *outImg2;

  double *xBlurImg, *xTranImg, filter1[7], filter2[7], sigmaPN, sigmaN, sigmaDiffN, interCoefL, interCoefU;
  double *lineExt, *coef0, *coef1, *coef2, coefN, coefTypeN;  
  double paraTheta, cosPara[6], sinPara[6], cosInitCoef[6], sinInitCoef[6], *nullD = NULL;
  int    K, K2, m, n, pos, maxMN, p, pd, P1, interType, interPosC, interPosS, Dim1 = 0, MDim1;

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

  typeProc = type % 10;
  
  if (P <= 0 || P >= 7) {
    mexErrMsgIdAndTxt("MATLAB", "P is not correct. P should be 1, 2, 3, 4 5, or 6");
  }
  if (M == 1 || N == 1) { /* 1D data */
   Dim1 = 1;
  /* For output */
   if (nlhs != 1 ) {
     mexErrMsgIdAndTxt("MATLAB", "Number of output is not correct: outImg = gaussSmooth(inImg, type, sigma, P, extType)");
   } else {
     plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
     outImg1 = mxGetPr(plhs[0]);
   }
 } else { /* 2D case */
   /* For output */
   if (typeProc == 0 || typeProc == 2) {
     if (nlhs != 1 ) {
       mexErrMsgIdAndTxt("MATLAB", "Number of output is not correct: outImg = gaussSmooth(inImg, type, sigma, P, extType)");
     } else {
       plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
       outImg1 = mxGetPr(plhs[0]);
     }
   } else if (typeProc == 1) {
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
  
   /* For work data (2D data) */
   xBlurImg = (double *) mxMalloc(sizeof(double) * M * N);
   if (typeProc != 0) {
     xTranImg = (double *) mxMalloc(sizeof(double) * M * N);
   }
 }

  /* Decide coefficients */
  decideCoef(P, type, &sigmaPN, &coef0, &coef1, &coef2, & coefTypeN);

  /* Decide K */ 
  K  = (int) (sigma * PI / sigmaPN + 0.5);
  K2 = K + K;

  /* Line extension  */
  lineExt = (double *) mxMalloc(sizeof(double) * (K2 + ((M > N)? M : N)));

  /* Interporation of sigma to interporate Fourier coefficients */
  coefN      = PI / K;
  switch (typeProc) {
  case 0:
    coefTypeN = coefN; break;
  case 1:
    coefTypeN = coefN * coefN; break;
  case 2:
    coefTypeN = coefN * coefN * coefN; break;
  }
  
  sigmaN     = sigma * coefN;
  sigmaDiffN = (sigmaN - sigmaPN) / sigmaPN;

  /* Set Interporation type */
  if (sigmaDiffN < -0.5)     interType = 0; 
  else if (sigmaDiffN < 0.0) interType = 1;
  else if (sigmaDiffN < 0.5) interType = 2;
  else                       interType = 3;

  /* Interporation of Fourier coeficients */
  P1 = P + 1;
  interCoefU  = (sigmaDiffN - 0.05 * (interType - 2)) / 0.05; /* 0.005 step for normalized sigma */
  interCoefL  = 1.0 - interCoefU;
  interCoefU *= coefTypeN;             /* To normize output */
  interCoefL *= coefTypeN;             /* To normize output */
  interPosC   = P1 * interType;
  interPosS   = P  * interType;

  switch(typeProc) {
  case 0: /* Gauss smooth */
    for (p = 0 ; p <= P ; ++p) {
      filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p]; /* Gauss */
    }
    break;
  case 1: /* Gauss differential */
    for (p = 0 ; p <= P ; ++p) {
      filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p]; /* Gauss */
    }
    for (pd = 0 ; pd <P ; ++pd) {
      filter2[pd] = interCoefL * coef1[interPosS + pd] + interCoefU * coef1[interPosS + P + pd]; /* differential of Gauss */
    }
    break;
  case 2: /* Gauss Laplacian */
    for (p = 0 ; p <= P ; ++p) {
      filter1[p] = interCoefL * coef0[interPosC + p] + interCoefU * coef0[interPosC + P1 + p];  /* Gauss */
      filter2[p] = interCoefL * coef2[interPosC + p] + interCoefU * coef2[interPosC + P1 + p];  /* Laplacian */
    }
    break;
  }

  /* Parameter for sliding Fourier transform */
  /* Initial coefficient of integral signal */
  paraTheta =  PI / K;
  for(p = 0 ; p < P ; ++p) {
    cosPara[p] = cos((p + 1) * paraTheta);
    sinPara[p] = sin((p + 1) * paraTheta);
    if ((p + 1) % K2 == 0) {
      cosInitCoef[p] = K;
      sinInitCoef[p] = 0.0;
    } else if ((p + 1) % 2 == 0) {
      cosInitCoef[p] = 0.0;
      sinInitCoef[p] = 0.0;
    } else {
      cosInitCoef[p] = 1.0;
      sinInitCoef[p] = sinPara[p] / (1 - cosPara[p]);
    }
  }

  if (Dim1 == 1) {  /* 1D case */
    MDim1 = (M > 1)? M : N;
    /* Transformation */
    if (extType == 0) {
      for (pos = 0     ; pos < K      ; ++pos) lineExt[pos] = 0.0;
      for (pos = MDim1 + K ; pos < MDim1 + K2 ; ++pos) lineExt[pos] = 0.0;
    }
    switch(typeProc) {
    case 0: /* Gauss smooth */
      slideFilterExt(1, inImg, MDim1, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, nullD,  outImg1, nullD, 1, lineExt);
      break;
    case 1: /* Differential of Gauss */
      slideFilterExt(2, inImg, MDim1, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter2, nullD, outImg1, nullD, 1, lineExt);
      break;
    case 2: /* Laplacian of Gauss */
      slideFilterExt(1, inImg, MDim1, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter2, nullD, outImg1, nullD, 1, lineExt);
      break;
    }
  } else {  /* 2D case */
    /* Vertical transformation */
    if (extType == 0) {
      for (pos = 0     ; pos < K      ; ++pos) lineExt[pos] = 0.0;
      for (pos = M + K ; pos < M + K2 ; ++pos) lineExt[pos] = 0.0;
    }
    for(n = 0 ; n < N ; ++n) {
      switch(typeProc) {
      case 0: /* Gauss smooth */
	slideFilterExt(1, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, nullD, &(xBlurImg[n]), nullD, N, lineExt);
	break;
      case 1: /* Differential of Gauss */
	slideFilterExt(4, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, filter2, &(xBlurImg[n]), &(xTranImg[n]), N, lineExt);
	break;
      case 2: /* Laplacian of Gauss */
	slideFilterExt(3, & (inImg[n * M]), M, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, filter2, &(xBlurImg[n]), &(xTranImg[n]), N, lineExt);
	break;
      }
    }

    /* Horizontal transformation */
    if (extType == 0) {
      for (pos = 0     ; pos < K      ; ++pos) lineExt[pos] = 0.0;
      for (pos = N + K ; pos < N + K2 ; ++pos) lineExt[pos] = 0.0;
    }
    for(m = 0 ; m < M ; ++m) {
      switch(typeProc) {
      case 0: /* Gauss smooth */
	slideFilterExt(1, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, nullD, &(outImg1[m]), nullD, M, lineExt);
	break;
      case 1: /* Differential of Gauss */
	slideFilterExt(2, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter2, nullD, &(outImg1[m]), nullD, M, lineExt);
	slideFilterExt(1, & (xTranImg[m * N]), N, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, nullD, &(outImg2[m]), nullD, M, lineExt);
	break;
      case 2: /* Laplacian of Gauss */
	slideFilterExt(1, & (xTranImg[m * N]), N, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter1, nullD, &(outImg1[m]), nullD, M, lineExt);
	slideFilterExt(5, & (xBlurImg[m * N]), N, K, extType, P, cosPara, sinPara, cosInitCoef, sinInitCoef, filter2, nullD, &(outImg1[m]), nullD, M, lineExt);
	break;
      }
    }
    /* Free work memories */
    mxFree(xBlurImg);
    if (typeProc != 0)  mxFree(xTranImg);
  }
  mxFree(lineExt);
  
  return;
}

/* Type 0: no output, 1: cos, 2: sin, 3: cos cos, 4: cos sin, 5: cos with output addition */
void slideFilterExt(int type, double *vect, int L, int K, int extType, int P, double *cosPara, double *sinPara, double *cosInitCoef, double *sinInitCoef, double *filter1, double *filter2, double *outVect1, double *outVect2, int outStep, double *lineExt) {
  double inte, inteCos[6], inteSin[6], extS, inteCosTmp, inteSinTmp, add, sub, val1, val2;
  int p, ordPos, pos, outPos = 0, K2 = K + K;

  for (pos = 0 ; pos < L ; ++pos) lineExt[pos + K] = vect[pos];
  if (extType == 1) {
    for (pos = 0     ; pos < K      ; ++pos) lineExt[pos] = vect[0];
    for (pos = L + K ; pos < L + K2 ; ++pos) lineExt[pos] = vect[L - 1];
  }

  /* Initial value of integral signal */
  extS = (extType == 0) ? 0.0 : vect[0]; 
  inte = K * extS;
  for (p = 0 ; p < P ; ++p) {
    inteCos[p] = cosInitCoef[p] * extS;
    inteSin[p] = sinInitCoef[p] * extS;
  }
  
  switch (type) {
  case 1: /* cos */
    for (pos = 0 ; pos < K ; ++pos) {
      inte += add = lineExt[pos + K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
      }
    }
    for (pos = K ; pos < L + K ; ++pos) {
      inte += add = lineExt[pos + K];
      val1 = filter1[0] * inte;
      inte -= sub = lineExt[pos - K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
	val1 += filter1[p + 1] * inteCos[p];
	inteCos[p] -= sub;
      }
      outVect1[outPos] = val1; outPos += outStep;
    }
    break;
  case 2: /* sin */
    for (pos = 0 ; pos < K ; ++pos) {
      add = lineExt[pos + K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
      }
    }
    for (pos = K ; pos < L + K ; ++pos) {
      add = lineExt[pos + K];
      val1 = 0.0;
      sub = lineExt[pos - K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
	val1 += filter1[p] * inteSin[p];
	inteCos[p] -= sub;
      }
      outVect1[outPos] = val1; outPos += outStep;
    }
    break;
  case 3: /* cos cos for Laplacian  */
    for (pos = 0 ; pos < K ; ++pos) {
      inte += add = lineExt[pos + K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
      }
    }
    for (pos = K ; pos < L + K ; ++pos) {
      inte += add = lineExt[pos + K];
      val1 = filter1[0] * inte;
      val2 = filter2[0] * inte;
      inte -= sub = lineExt[pos - K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
	val1 += filter1[p + 1] * inteCos[p];
	val2 += filter2[p + 1] * inteCos[p];
	inteCos[p] -= sub;
      }
      outVect1[outPos] = val1;
      outVect2[outPos] = val2; outPos += outStep;
    }
    break;
  case 4: /* cos sin */
    for (pos = 0 ; pos < K ; ++pos) {
      inte += add = lineExt[pos + K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
      }
    }
    for (pos = K ; pos < L + K ; ++pos) {
      inte += add = lineExt[pos + K];
      val1 = filter1[0] * inte;
      val2 = 0.0;
      inte -= sub = lineExt[pos - K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
	val1 += filter1[p + 1] * inteCos[p];
	val2 += filter2[p] * inteSin[p];
	inteCos[p] -= sub;
      }
      outVect1[outPos] = val1;
      outVect2[outPos] = val2; outPos += outStep;
    }
    break;
  case 5: /* cos with output addition */
    for (pos = 0 ; pos < K ; ++pos) {
      inte += add = lineExt[pos + K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
      }
    }
    for (pos = K ; pos < L + K ; ++pos) {
      inte += add = lineExt[pos + K];
      val1 = filter1[0] * inte;
      inte -= sub = lineExt[pos - K];
      for (p = 0 ; p < P ; ++p) {
	inteCosTmp = inteCos[p]; inteSinTmp = inteSin[p]; 
	inteCos[p] = cosPara[p] * inteCosTmp - sinPara[p] * inteSinTmp + add;
	inteSin[p] = sinPara[p] * inteCosTmp + cosPara[p] * inteSinTmp;
	val1 += filter1[p + 1] * inteCos[p];
	inteCos[p] -= sub;
      }
      outVect1[outPos] += val1; outPos += outStep;
    }
    break;
  }
}

void decideCoef(int P, int type, double *sigmaPN, double **coef0, double **coef1, double **coef2, double *coefTypeN) {
  /* Invtype 1*/
static double sigmaI1[] = {1.2400000000e+00, 1.0200000000e+00, 8.8000000000e-01, 7.9000000000e-01, 7.2000000000e-01, 6.7000000000e-01};
static double coefGP1I1[] = {1.5841866036e-01, -1.7209198406e-01, 1.5796153279e-01, -1.6115813726e-01, 1.5736998700e-01, -1.5065819241e-01, 1.5663733553e-01, -1.4065118498e-01, 1.5576238574e-01, -1.3117575039e-01};
static double coefGDP1I1[] = {1.7217168122e-01, 1.6120804118e-01, 1.5068229793e-01, 1.4065337990e-01, 1.3115967709e-01};
/* filter coefficients P = 1 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP1I1[] = {-5.6886688489e-03, 1.6079448054e-01, -7.1572515210e-03, 1.4689371300e-01, -8.6038352226e-03, 1.3347483955e-01, -9.9601722462e-03, 1.2073328229e-01, -1.1177671384e-02, 1.0880461242e-01};
static double coefGP2I1[] = {1.5904136047e-01, -2.0907708286e-01, 5.8810354526e-02, 1.5895487664e-01, -1.9942521467e-01, 4.8350584150e-02, 1.5881608670e-01, -1.8983478272e-01, 3.9215803701e-02, 1.5861236247e-01, -1.8037606324e-01, 3.1311210452e-02, 1.5833236709e-01, -1.7111223318e-01, 2.4527720571e-02};
static double coefGDP2I1[] = {2.0904759351e-01, -1.1767968775e-01, 1.9940270967e-01, -9.6746178304e-02, 1.8981590703e-01, -7.8469358782e-02, 1.8035824721e-01, -6.2658052983e-02, 1.7109362001e-01, -4.9092667502e-02};
/* filter coefficients P = 2 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP2I1[] = {-1.3809848228e-03, 2.0628566092e-01, -2.3812119693e-01, -2.2354669342e-03, 1.9493183322e-01, -1.9796306079e-01, -3.2531388896e-03, 1.8330971144e-01, -1.6344466659e-01, -4.3961643051e-03, 1.7156602894e-01, -1.3410799319e-01, -5.6170399505e-03, 1.5985968076e-01, -1.0941885229e-01};
static double coefGP3I1[] = {1.5914693169e-01, -2.3263233798e-01, 9.0771877184e-02, -1.8933051456e-02, 1.5912938906e-01, -2.2448215788e-01, 7.8623918729e-02, -1.3745044056e-02, 1.5909849134e-01, -2.1622614400e-01, 6.7546620727e-02, -9.8381486188e-03, 1.5904696173e-01, -2.0791382321e-01, 5.7535413269e-02, -6.9710520851e-03, 1.5896656739e-01, -1.9959504616e-01, 4.8561734803e-02, -4.9294004810e-03};
static double coefGDP3I1[] = {2.3263950832e-01, -1.8152941370e-01, 5.6820665377e-02, 2.2448557908e-01, -1.5724099505e-01, 4.1245395775e-02, 2.1622684967e-01, -1.3509183012e-01, 2.9516562861e-02, 2.0791232335e-01, -1.1507382627e-01, 2.0908656648e-02, 1.9959150834e-01, -9.7130545242e-02, 1.4777587982e-02};
/* filter coefficients P = 3 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP3I1[] = {-3.7199908684e-04, 2.3189551787e-01, -3.6380279469e-01, 1.6971806744e-01, -6.2828114792e-04, 2.2322903149e-01, -3.1573849359e-01, 1.2247975735e-01, -1.0305925708e-03, 2.1416568962e-01, -2.7224474500e-01, 8.6488729277e-02, -1.5873721648e-03, 2.0473761822e-01, -2.3332224007e-01, 5.9551578423e-02, -2.2939214255e-03, 1.9500372247e-01, -1.9884870542e-01, 3.9745433905e-02};
static double coefGP4I1[] = {1.5915228801e-01, -2.4722184936e-01, 1.1580869816e-01, -3.2731159398e-02, 5.5736570555e-03, 1.5914981887e-01, -2.4019281495e-01, 1.0317574361e-01, -2.5247744922e-02, 3.5073995798e-03, 1.5914338104e-01, -2.3300776552e-01, 9.1340959054e-02, -1.9211745798e-02, 2.1447320801e-03, 1.5913022647e-01, -2.2569887337e-01, 8.0346474011e-02, -1.4429568242e-02, 1.2641625022e-03, 1.5910649479e-01, -2.1830007879e-01, 7.0212668746e-02, -1.0711084412e-02, 7.0108357984e-04};
static double coefGDP4I1[] = {2.4721970041e-01, -2.3162169421e-01, 9.8187031346e-02, -2.2303224021e-02, 2.4019160005e-01, -2.0635391703e-01, 7.5739590054e-02, -1.4034457935e-02, 2.3300687690e-01, -1.8268369534e-01, 5.7632571545e-02, -8.5824827866e-03, 2.2569789666e-01, -1.6069490143e-01, 4.3285774607e-02, -5.0605568373e-03, 2.1829869671e-01, -1.4042810164e-01, 3.2129107016e-02, -2.8098626194e-03};
/* filter coefficients P = 4 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP4I1[] = {-4.0101325536e-05, 2.4713949937e-01, -4.6332358466e-01, 2.9448090582e-01, -8.9293073068e-02, -1.3791703462e-04, 2.3991576969e-01, -4.1298365328e-01, 2.2694296948e-01, -5.6413606448e-02, -2.9624826565e-04, 2.3241438784e-01, -3.6595985731e-01, 1.7230528537e-01, -3.4922308085e-02, -5.4050841575e-04, 2.2461689330e-01, -3.2247076583e-01, 1.2877642820e-01, -2.1323028693e-02, -8.9055554745e-04, 2.1651760777e-01, -2.8263722576e-01, 9.4606409327e-02, -1.3020207122e-02};
static double coefGP5I1[] = {1.5915503671e-01, -2.5802917112e-01, 1.3744432403e-01, -4.8108193838e-02, 1.1065336182e-02, -1.6720449053e-03, 1.5915435823e-01, -2.5191757166e-01, 1.2487394499e-01, -3.8772003690e-02, 7.5388137004e-03, -9.1907371368e-04, 1.5915291348e-01, -2.4563373907e-01, 1.1286499565e-01, -3.0886976177e-02, 5.0290702127e-03, -4.9079162876e-04, 1.5914970525e-01, -2.3920001418e-01, 1.0147981895e-01, -2.4323942904e-02, 3.2821823431e-03, -2.5751054259e-04, 1.5914316964e-01, -2.3263986209e-01, 9.0764353079e-02, -1.8940575561e-02, 2.0906930624e-03, -1.3818790752e-04};
static double coefGDP5I1[] = {2.5802975475e-01, -2.7488748079e-01, 1.4432633242e-01, -4.4259010183e-02, 8.3631427094e-03, 2.5191779206e-01, -2.4974744916e-01, 1.1631667229e-01, -3.0154373174e-02, 4.5964706025e-03, 2.4563375548e-01, -2.2572995849e-01, 9.2660977760e-02, -2.0116215211e-02, 2.4540401926e-03, 2.3919986390e-01, -2.0295993845e-01, 7.2971377884e-02, -1.3129330479e-02, 1.2868013277e-03, 2.3263950832e-01, -1.8152941370e-01, 5.6820665377e-02, -8.3641873270e-03, 6.8917068732e-04};
/* filter coefficients P = 5 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP5I1[] = {-2.2307113655e-05, 2.5798514082e-01, -5.4981957465e-01, 4.3293438565e-01, -1.7708065033e-01, 4.1771106557e-02, -3.7945743694e-05, 2.5184190140e-01, -4.9957078652e-01, 3.4887413277e-01, -1.2069337105e-01, 2.2906482052e-02, -8.1724321087e-05, 2.4547030881e-01, -4.5162335774e-01, 2.7781950236e-01, -8.0628277984e-02, 1.2106801544e-02, -1.6742793230e-04, 2.3886501216e-01, -4.0625471628e-01, 2.1857931488e-01, -5.2852111838e-02, 6.0992538119e-03, -3.1177236467e-04, 2.3201597131e-01, -3.6368234124e-01, 1.6983852088e-01, -3.4080170519e-02, 2.8225017095e-03};
static double coefGP6I1[] = {1.5915483517e-01, -2.6539564269e-01, 1.5382352551e-01, -6.1978545477e-02, 1.7359527693e-02, -3.3803279141e-03, 4.5732120504e-04, 1.5915478419e-01, -2.5994250638e-01, 1.4156496718e-01, -5.1415305003e-02, 1.2452802114e-02, -2.0117335041e-03, 2.1644494008e-04, 1.5915448582e-01, -2.5431634174e-01, 1.2969931000e-01, -4.2224232376e-02, 8.7735612826e-03, -1.1645467854e-03, 9.7977991487e-05, 1.5915364359e-01, -2.4853396660e-01, 1.1829488612e-01, -3.4328791413e-02, 6.0702760907e-03, -6.5655091178e-04, 4.1558471882e-05, 1.5915167095e-01, -2.4261291886e-01, 1.0740768657e-01, -2.7631589916e-02, 4.1227943589e-03, -3.6227897192e-04, 1.4608360457e-05};
static double coefGDP6I1[] = {2.6539548692e-01, -3.0764736257e-01, 1.8593516910e-01, -6.9438733884e-02, 1.6900860682e-02, -2.7448618962e-03, 2.5994244255e-01, -2.8313006203e-01, 1.5424572351e-01, -4.9811463787e-02, 1.0058348358e-02, -1.2990526361e-03, 2.5431630078e-01, -2.5939870192e-01, 1.2667257423e-01, -3.5094408992e-02, 5.8225290995e-03, -5.8811374224e-04, 2.4853390706e-01, -2.3658989132e-01, 1.0298619562e-01, -2.4281342521e-02, 3.2824568609e-03, -2.4970806984e-04, 2.4261280137e-01, -2.1481560812e-01, 8.2894417278e-02, -1.6491647398e-02, 1.8108074052e-03, -8.8355109887e-05};
/* filter coefficients P = 6 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP6I1[] = {1.3037307758e-06, 2.6539809444e-01, -6.1529211745e-01, 5.5780811528e-01, -2.7775232714e-01, 8.4506912327e-02, -1.6466561820e-02, -6.7631527270e-06, 2.5992891644e-01, -5.6627364957e-01, 4.6272364602e-01, -1.9925937827e-01, 5.0278220457e-02, -7.8078349579e-03, -2.2289828165e-05, 2.5427172168e-01, -5.1884198126e-01, 3.7997314807e-01, -1.4042220667e-01, 2.9068079831e-02, -3.5732419643e-03, -5.4314851895e-05, 2.4842527870e-01, -4.7328840697e-01, 3.0884996924e-01, -9.7233978293e-02, 1.6303688186e-02, -1.6068297585e-03, -1.1453215318e-04, 2.4238373990e-01, -4.2986026920e-01, 2.4845421303e-01, -6.6195608560e-02, 8.8250435642e-03, -7.5909294683e-04};


  /* Invtype 2*/
 static double sigmaI2[] = {1.2400000000e+00, 1.0200000000e+00, 8.8000000000e-01, 7.9000000000e-01, 7.2000000000e-01, 6.7000000000e-01};
static double coefGP1I2[] = {1.6296862415e-01, -1.6299205647e-01, 1.5902524158e-01, -1.5903071968e-01, 1.5513655619e-01, -1.5512505402e-01, 1.5131775064e-01, -1.5129035477e-01, 1.4758088553e-01, -1.4753875080e-01};
static double coefGDP1I2[] = {1.7217168122e-01, 1.6120804118e-01, 1.5068229793e-01, 1.4065337990e-01, 1.3115967709e-01};
/* filter coefficients P = 1 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP1I2[] = {-5.7302003188e-02, 5.7567811860e-02, -5.3656248736e-02, 5.3895718572e-02, -5.0156171820e-02, 5.0370166354e-02, -4.6821264052e-02, 4.7011098675e-02, -4.3664216493e-02, 4.3831522201e-02};
static double coefGP2I2[] = {1.5728745386e-01, -2.1258489608e-01, 5.5302541304e-02, 1.5737974328e-01, -2.0257548138e-01, 4.5200317437e-02, 1.5717761786e-01, -1.9311172040e-01, 3.5938866017e-02, 1.5670397018e-01, -1.8419284783e-01, 2.7494425867e-02, 1.5598416157e-01, -1.7580864423e-01, 1.9831309528e-02};
static double coefGDP2I2[] = {2.0904759351e-01, -1.1767968775e-01, 1.9940270967e-01, -9.6746178304e-02, 1.8981590703e-01, -7.8469358782e-02, 1.8035824721e-01, -6.2658052983e-02, 1.7109362001e-01, -4.9092667502e-02};
/* filter coefficients P = 2 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP2I2[] = {5.2584588115e-03, 2.1956454819e-01, -2.2484230966e-01, -1.1827401464e-03, 1.9703728680e-01, -1.9585760721e-01, -6.5735893904e-03, 1.7666881043e-01, -1.7008556759e-01, -1.1004696036e-02, 1.5834896548e-01, -1.4732505665e-01, -1.4576588088e-02, 1.4194058449e-01, -1.2733794856e-01};
static double coefGP3I2[] = {1.5938118083e-01, -2.3216383971e-01, 9.1240375459e-02, -1.8464553181e-02, 1.5919680714e-01, -2.2434732172e-01, 7.8758754882e-02, -1.3610207903e-02, 1.5901586164e-01, -2.1639140339e-01, 6.7381361332e-02, -1.0003408014e-02, 1.5880546863e-01, -2.0839680941e-01, 5.7052427076e-02, -7.4540382786e-03, 1.5853922688e-01, -2.0044972717e-01, 4.7707053791e-02, -5.7840814938e-03};
static double coefGDP3I2[] = {2.3263950832e-01, -1.8152941370e-01, 5.6820665377e-02, 2.2448557908e-01, -1.5724099505e-01, 4.1245395775e-02, 2.1622684967e-01, -1.3509183012e-01, 2.9516562861e-02, 2.0791232335e-01, -1.1507382627e-01, 2.0908656648e-02, 1.9959150834e-01, -9.7130545242e-02, 1.4777587982e-02};
/* filter coefficients P = 3 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP3I2[] = {-5.6981913249e-03, 2.2124313339e-01, -3.7445517916e-01, 1.5906568296e-01, -4.8025934409e-03, 2.1488040690e-01, -3.2408711818e-01, 1.1413113276e-01, -4.9256501963e-03, 2.0637557437e-01, -2.8003486025e-01, 7.8698614026e-02, -5.7670293603e-03, 1.9637830383e-01, -2.4168155446e-01, 5.1192264032e-02, -7.0749207421e-03, 1.8544172384e-01, -2.0841070405e-01, 3.0183435272e-02};
static double coefGP4I2[] = {1.5908778005e-01, -2.4735086529e-01, 1.1567968223e-01, -3.2860175324e-02, 5.4446411297e-03, 1.5910629829e-01, -2.4027985610e-01, 1.0308870246e-01, -2.5334786070e-02, 3.4203584318e-03, 1.5909795743e-01, -2.3309861274e-01, 9.1250111835e-02, -1.9302593017e-02, 2.0538848608e-03, 1.5906230399e-01, -2.2583471833e-01, 8.0210629049e-02, -1.4565413203e-02, 1.1283175409e-03, 1.5899457923e-01, -2.1852390991e-01, 6.9988837620e-02, -1.0934915538e-02, 4.7725245370e-04};
static double coefGDP4I2[] = {2.4721970041e-01, -2.3162169421e-01, 9.8187031346e-02, -2.2303224021e-02, 2.4019160005e-01, -2.0635391703e-01, 7.5739590054e-02, -1.4034457935e-02, 2.3300687690e-01, -1.8268369534e-01, 5.7632571545e-02, -8.5824827866e-03, 2.2569789666e-01, -1.6069490143e-01, 4.3285774607e-02, -5.0605568373e-03, 2.1829869671e-01, -1.4042810164e-01, 3.2129107016e-02, -2.8098626194e-03};
/* filter coefficients P = 4 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP4I2[] = {1.1839195282e-03, 2.4958754108e-01, -4.6087554295e-01, 2.9692894753e-01, -8.6845031360e-02, 1.5892152928e-04, 2.4050944682e-01, -4.1238997616e-01, 2.2753664661e-01, -5.5819929320e-02, -6.8900257520e-04, 2.3162887922e-01, -3.6674536593e-01, 1.7151977675e-01, -3.5707816704e-02, -1.5452268021e-03, 2.2260745653e-01, -3.2448020260e-01, 1.2676699143e-01, -2.3332465466e-02, -2.5071549508e-03, 2.1328440897e-01, -2.8587042457e-01, 9.1373210521e-02, -1.6253405929e-02};
static double coefGP5I2[] = {1.5916808213e-01, -2.5800308027e-01, 1.3747041488e-01, -4.8082102991e-02, 1.1091427030e-02, -1.6459540582e-03, 1.5915810223e-01, -2.5191008367e-01, 1.2488143297e-01, -3.8764515702e-02, 7.5463016883e-03, -9.1158572573e-04, 1.5914971575e-01, -2.4564013454e-01, 1.1285860019e-01, -3.0893371644e-02, 5.0226747453e-03, -4.9718709615e-04, 1.5913796461e-01, -2.3922349545e-01, 1.0145633768e-01, -2.4347424175e-02, 3.2587010729e-03, -2.8099181283e-04, 1.5911796545e-01, -2.3269027047e-01, 9.0713944693e-02, -1.8990983948e-02, 2.0402846757e-03, -1.8859629419e-04};
static double coefGDP5I2[] = {2.5802975475e-01, -2.7488748079e-01, 1.4432633242e-01, -4.4259010183e-02, 8.3631427094e-03, 2.5191779206e-01, -2.4974744916e-01, 1.1631667229e-01, -3.0154373174e-02, 4.5964706025e-03, 2.4563375548e-01, -2.2572995849e-01, 9.2660977760e-02, -2.0116215211e-02, 2.4540401926e-03, 2.3919986390e-01, -2.0295993845e-01, 7.2971377884e-02, -1.3129330479e-02, 1.2868013277e-03, 2.3263950832e-01, -1.8152941370e-01, 5.6820665377e-02, -8.3641873270e-03, 6.8917068732e-04};
/* filter coefficients P = 5 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP5I2[] = {-5.4228356064e-04, 2.5694518792e-01, -5.5085952755e-01, 4.3189443276e-01, -1.7812060322e-01, 4.0731153663e-02, -3.3727069151e-04, 2.5124325150e-01, -5.0016943642e-01, 3.4827548288e-01, -1.2129202094e-01, 2.2307832156e-02, -3.5786706807e-04, 2.4491802331e-01, -4.5217564324e-01, 2.7726721686e-01, -8.1180563478e-02, 1.1554516050e-02, -5.5229434994e-04, 2.3809527932e-01, -4.0702444912e-01, 2.1780958204e-01, -5.3621844674e-02, 5.3295209767e-03, -9.0698621624e-04, 2.3082554361e-01, -3.6487276894e-01, 1.6864809318e-01, -3.5270598222e-02, 1.6320740064e-03};
static double coefGP6I2[] = {1.5915171749e-01, -2.6540187805e-01, 1.5381729015e-01, -6.1984780833e-02, 1.7353292338e-02, -3.3865632697e-03, 4.5108584939e-04, 1.5915329378e-01, -2.5994548719e-01, 1.4156198637e-01, -5.1418285818e-02, 1.2449821300e-02, -2.0147143185e-03, 2.1346412565e-04, 1.5915293714e-01, -2.5431943911e-01, 1.2969621263e-01, -4.2227329742e-02, 8.7704639167e-03, -1.1676441512e-03, 9.4880625629e-05, 1.5915049819e-01, -2.4854025740e-01, 1.1828859532e-01, -3.4335082215e-02, 6.0639852885e-03, -6.6284171403e-04, 3.5267669632e-05, 1.5914477783e-01, -2.4262670510e-01, 1.0739390032e-01, -2.7645376163e-02, 4.1090081127e-03, -3.7606521817e-04, 8.2211420722e-07};
static double coefGDP6I2[] = {2.6539548692e-01, -3.0764736257e-01, 1.8593516910e-01, -6.9438733884e-02, 1.6900860682e-02, -2.7448618962e-03, 2.5994244255e-01, -2.8313006203e-01, 1.5424572351e-01, -4.9811463787e-02, 1.0058348358e-02, -1.2990526361e-03, 2.5431630078e-01, -2.5939870192e-01, 1.2667257423e-01, -3.5094408992e-02, 5.8225290995e-03, -5.8811374224e-04, 2.4853390706e-01, -2.3658989132e-01, 1.0298619562e-01, -2.4281342521e-02, 3.2824568609e-03, -2.4970806984e-04, 2.4261280137e-01, -2.1481560812e-01, 8.2894417278e-02, -1.6491647398e-02, 1.8108074052e-03, -8.8355109887e-05};
/* filter coefficients P = 6 sgimaN = sigma0 * [0.9 0.05 1.0 1.05, 1.1] */
static double coefGDDP6I2[] = {1.3894637695e-04, 2.6567337973e-01, -6.1501683216e-01, 5.5808340058e-01, -2.7747704185e-01, 8.4782197619e-02, -1.6191276527e-02, 2.5172724129e-05, 2.5999278819e-01, -5.6620977782e-01, 4.6278751778e-01, -1.9919550651e-01, 5.0342092211e-02, -7.7439632042e-03, -5.7013433774e-05, 2.5420227447e-01, -5.1891142847e-01, 3.7990370086e-01, -1.4049165388e-01, 2.8998632619e-02, -3.6426891756e-03, -1.6122205433e-04, 2.4821146430e-01, -4.7350222137e-01, 3.0863615484e-01, -9.7447792698e-02, 1.6089873781e-02, -1.8206441633e-03, -3.2395514705e-04, 2.4196489391e-01, -4.3027911519e-01, 2.4803536704e-01, -6.6614454548e-02, 8.4061975764e-03, -1.1779389346e-03};

  /* Fix K and filter coefficients */
  if (type < 10) {
    switch(P) {
    case 1:
      *sigmaPN = sigmaI1[0]; *coef0 = coefGP1I1; *coef1 = coefGDP1I1; *coef2 = coefGDDP1I1; break;
    case 2:
      *sigmaPN = sigmaI1[1]; *coef0 = coefGP2I1; *coef1 = coefGDP2I1; *coef2 = coefGDDP2I1; break;
    case 3:
      *sigmaPN = sigmaI1[2]; *coef0 = coefGP3I1; *coef1 = coefGDP3I1; *coef2 = coefGDDP3I1; break;
    case 4:
      *sigmaPN = sigmaI1[3]; *coef0 = coefGP4I1; *coef1 = coefGDP4I1; *coef2 = coefGDDP4I1; break;
    case 5:
      *sigmaPN = sigmaI1[4]; *coef0 = coefGP5I1; *coef1 = coefGDP5I1; *coef2 = coefGDDP5I1; break;
    case 6:
      *sigmaPN = sigmaI1[5]; *coef0 = coefGP6I1; *coef1 = coefGDP6I1; *coef2 = coefGDDP6I1; break;
    }
  } else {
    switch(P) {
    case 1:
      *sigmaPN = sigmaI2[0]; *coef0 = coefGP1I2; *coef1 = coefGDP1I2; *coef2 = coefGDDP1I2; break;
    case 2:
      *sigmaPN = sigmaI2[1]; *coef0 = coefGP2I2; *coef1 = coefGDP2I2; *coef2 = coefGDDP2I2; break;
    case 3:
      *sigmaPN = sigmaI2[2]; *coef0 = coefGP3I2; *coef1 = coefGDP3I2; *coef2 = coefGDDP3I2; break;
    case 4:
      *sigmaPN = sigmaI2[3]; *coef0 = coefGP4I2; *coef1 = coefGDP4I2; *coef2 = coefGDDP4I2; break;
    case 5:
      *sigmaPN = sigmaI2[4]; *coef0 = coefGP5I2; *coef1 = coefGDP5I2; *coef2 = coefGDDP5I2; break;
    case 6:
      *sigmaPN = sigmaI2[5]; *coef0 = coefGP6I2; *coef1 = coefGDP6I2; *coef2 = coefGDDP6I2; break;
    }
  }
}
