type    = 2; % type  : 0 smooth, 1 diff (two outputs), 2 : LOG
sigma   = 10;
P       = 6;
extType = 1; % extType :  0: zero extension,  1: the value of edge is used for extension

inImg = im2double(imread('img1.pgm'));
%inImg = zeros(100);
%inImg(50, 50) = 1;

% [blurImg blurImgY] = gSmooth(inImg, type, sigma, P, extType);

switch type
  case 0
    blurImg = gaussSmooth(inImg, type, sigma, P, extType);
  case 1
    [blurImg blurImgY] = gaussSmooth(inImg, type, sigma, P, extType);
  case 2
    blurImg = gaussSmooth(inImg, type, sigma, P, extType);
end

% blurImg = blurImgY;
maxV = max(max(blurImg));
minV = min(min(blurImg));
sImg = (blurImg - minV) / (maxV - minV);

figure(1)
imshow(inImg);
figure(2)
imshow(sImg);
