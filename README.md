# gaussSmooth
The programs are licenced by GPL v3. 

## File for the explanation of the algorithm.
* gaussSmooth.pdf

## A mex function for Gaussian smoothing (with x and y-directional differential and LOG)
The latest version. It works only by the mex function.  
  
* gaussSmooth.c
  * Usage:
    * type = 0:  (Gauss smoothing)  
      `smoothedImg   = gaussSmooth(inImg, type, sigma, P, extType)`  
    * type = 1:  (x and y-directional differential with Gauss smoothing)  
      `[diffX diffY] = gaussSmooth(inImg, type, sigma, P, extType)`  
    * type = 2:  (Laplacian of Gaussian)  
      `logImg        = gaussSmooth(inImg, type, sigma, P, extType)`  
  * Paramter:
    * inImg   : input image  
    * type    : 0 blur, 1 diff (two outputs), 2 : LOG  
    * sigma   : sigma  
    * P       : 2 or 4 or 6 : order of Fourier series  
    * extType :  0: zero extension,  1: the value of edge is used for extension  
  * Output: (The number of output is one or two.)
    * smoothImg: Gaussian smoothed image
    * diffX    : Gaussian smoothed x-directional differential image  
    * diffY    : Gaussian smoothed x-directional differential image  
    * lapImg   : Laplacian of Gaussian smoothed (LOG) image  
  
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
