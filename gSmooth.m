% inImg   : input image
% type    : 0 blur, 1 diff (two output), 2 : LOG
% sigma   : sigma
% P       : 2 or 4 or 6 : order of Fourier series
% extType :  0: zero extension,  1: the value of edge is used for extension

function [blurImg blurImgY] = gSmooth(inImg, type, sigma, P, extType)
  M        = size(inImg, 1);
  N        = size(inImg, 2);
  xBlurImg = zeros(M, N);
  blurImg  = zeros(M, N);
  switch type
    case 0
      blurImgY = type;
    case 1
      blurImgY = zeros(M, N);
      xTranImg = zeros(M, N);
    case 2
      blurImgY = type;
      xTranImg = zeros(M, N);
  end
  
% Set coefficient
  switch P
    case 2
      K = round(sigma * pi / 1.1);
      coef0 =   [ 1.5846315202e-01 -(1.7508079293e-01)   2.7323762984e-02];
      coef1 = - [                  -(-1.7506270969e-01) -5.4683692457e-02];
      coef2 =   [-5.0840517439e-03 -(-1.6489473361e-01) -1.1953497880e-01];
    case 4
      K = round(sigma * pi / 0.80);
      coef0 =   [ 1.5914081930e-01 -(2.3116780826e-01)   8.8476942331e-02 -(1.7890179087e-02) 1.8836364529e-03];
      coef1 = - [                  -(-2.3116693143e-01) -1.7695563833e-01 -(-5.3667906760e-02) -7.5380531476e-03];
      coef2 =   [-3.4905652299e-04 -(-2.3046882714e-01) -3.5460935466e-01 -(-1.6030568607e-01) -3.0850185473e-02];
    case 6
      K = round(sigma * pi / 0.70);
      coef0 =   [ 1.5915377148e-01  -( 2.4914489111e-01)  1.1946305955e-01 -( 3.5095823000e-02)  6.3138382801e-03	-( 6.9786114192e-04)  4.5673259016e-05];
      coef1 = - [                   -(-2.4914483525e-01) -2.3892623081e-01 -(-1.0528730143e-01) -2.5255576544e-02 -(-3.4890264302e-03) -2.7437469024e-04];
      coef2 =   [-4.9872068291e-05  -(-2.4904509235e-01) -4.7795220081e-01 -(-3.1576217126e-01) -1.0112203057e-01 -(-1.7345418863e-02) -1.7459478551e-03];
  end

  cosPara = zeros(1, P);
  sinPara = zeros(1, P);
  paraTheta =  pi / K;
  for p = 1:P
    cosPara(p) = cos(p * paraTheta);
    sinPara(p) = sin(p * paraTheta);
  end
% Horizontal transformation
  for m=1:M
    [cosInte sinInte] = scInte(inImg(m, :), K, extType, cosPara, sinPara);
    xBlurImg(m, :)     = coef0 * cosInte;
    switch type
      case 1
	xTranImg(m, :) = coef1 * sinInte;
      case 2
	xTranImg(m, :) = coef2 * cosInte;
    end
  end
% Vertical transformation
 for n=1:N
   [cosInte sinInte] = scInte((xBlurImg(:, n))', K, extType, cosPara, sinPara);
    switch type
      case 0
	blurImg(:, n)  = coef0 * cosInte;
      case 1
	blurImgY(:, n) = coef1 * sinInte;
      case 2
	blurImg(:, n)  = coef2 * cosInte;
    end
    if type ~= 0
      [cosInte sinInte] = scInte((xTranImg(:, n))', K, extType, cosPara, sinPara);
      switch type
	case 1
	  blurImg(:, n)  = coef0 * cosInte;
	case 2
	  blurImg(:, n)  = blurImg(:, n) + (coef0 * cosInte)';
      end
    end
 end
end
