# gaussSmooth
* The programs are licenced by GPL v3. 
* Gaussian smoothing (with x and y-directional differential and LOG) are approximately calculated by using sliding discrete Fourier transform of order P (= 2, 4, or 6).
* Let M, N, and K be the numbers of pixcels for vertical, horizontal, and window sizes, then the calculation complexity is (M + K) * (N + K) * P. 

## File for the explanation of the algorithm.
* gaussSmooth.pdf

## A mex function for Gaussian smoothing (with x and y-directional differential and LOG)
The latest version. It works only by the mex function.  
1D and 2D data are seperatly handled. 
  
* gaussSmooth.c
  * Usage:
    * type = 0:  (Gauss smoothing)  
      `smoothedImg   = gaussSmooth(inImg, type, sigma, P, extType)`  
    * type = 1:  (x and y-directional differential with Gauss smoothing)  
        * 1D Data  
      `diff = gaussSmooth(inImg, type, sigma, P, extType)`  
        * 2D Data  
      `[diffX diffY] = gaussSmooth(inImg, type, sigma, P, extType)`  
    * type = 2:  (Laplacian of Gaussian)  
      `logImg        = gaussSmooth(inImg, type, sigma, P, extType)`  
  * Paramter:
    * inImg   : input image  (1D/2D data)
    * type    : 0 blur, 1 diff (one or two outputs), 2 : LOG  
    * sigma   : sigma  
    * P       : 2 or 4 or 6 : order of Fourier series  
    * extType :  0: zero extension,  1: the value of edge is used for extension  
  * Output: (The number of output is one or two.)
    * smoothImg: Gaussian smoothed image (1D/2D data)
    * diff     : Gaussian smoothed differential data (1D data)
    * diffX    : Gaussian smoothed x-directional differential image  (2D data)
    * diffY    : Gaussian smoothed x-directional differential image  (2D data)
    * lapImg   : Laplacian of Gaussian smoothed (LOG) image  (1D/2D data)
  
* testSmooth.m : main program for test.  
  
* img1.pgm     : Boat image for test.  

## A matlab program for calculate coefficients and truncation errors
* totalError.m : main program
  * output
    * coefficientsForC.txt : coefficients for C program.
    * dataPlotFunc.txt : data to plot function of approximated and true Gaussian function
    * errorForPlot.txt : data to plot truncation errors 

## A matlab and mex (for Gaussian smoothing and for sinusoidal integral, respectively)  
  (This is an older version.)  
  
* gSmooth.m  : a function to calculate Gaussian smoothing  
  * Usage:
    * `[blurImg blurImgY] = gSmooth(inImg, type, sigma, P, extType)`  
  * Paramter:
     * inImg   : input image  
     * type    : 0 blur, 1 diff (two output), 2 : LOG 
     * sigma   : sigma  
     * P       : 2 or 4 or 6 : order of Fourier series  
     * extType :  0: zero extension,  1: the value of edge is used for extension  
  * Output:
     * blurImg  : output for Gaussian smoothed (type = 0)  
                       or Gaussian smoothed differential of horizontal direction (type = 1)  
	               or LOG image (type = 2).  
     * blurImgY : output for Gaussian smoothed differential of vertical directoin (type = 1)  
  
     * scInte.c : sliding discrete Fourier transform. (mex scInte.c)  
  
* testSmooth.m : main program for test.
  
* img1.pgm     : Boat image for test.
