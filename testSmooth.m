type    = 0; % type  : 0 smooth, 1 diff (two outputs), 2 : LOG
sigma   = 5;
P       = 4;
extType = 1; % extType :  0: zero extension,  1: the value of edge is used for extension

inImg = im2double(imread('img1.pgm'));

[blurImg blurImgY] = gSmooth(inImg, type, sigma, P, extType);

% blurImg = blurImgY;
maxV = max(max(blurImg));
minV = min(min(blurImg));
sImg = (blurImg - minV) / (maxV - minV);

figure(1)
imshow(inImg);
figure(2)
imshow(sImg);
