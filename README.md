# gSmooth

gSmooth.m  : a function to calculate Gaussian smoothing
   [blurImg blurImgY] = gSmooth(inImg, type, sigma, P, extType)
   inImg   : input image
   type    : 0 blur, 1 diff (two output), 2 : LOG
   sigma   : sigma
   P       : 2 or 4 or 6 : order of Fourier series
   extType :  0: zero extension,  1: the value of edge is used for extension

   blurImg  : output for Gaussian smoothed (type = 0)
                      or Gaussian smoothed differential of horizontal direction (type = 1)
	              or LOG image (type = 2).
   blurImgY : output for Gaussian smoothed differential of vertical directoin (type = 1)

scInte.c : sliding discrete Fourier transform. (mex scInte.c)

testSmooth.m : main program for test.
img1.pgm     : Boat image for test.
