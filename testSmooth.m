type    = 0; % type  : 0 smooth, 1 diffX 2 : diffY, 2 : LOG
sigma   = 11.5;
P       = 4;
extType = 1; % extType :  0: zero extension,  1: the value of edge is used for extension

inImg = im2double(imread('img1.pgm'));
%inImg = zeros(200);
%inImg(100, 100) = 1;
%inImg = ones(64, 1) * [zeros(1,32) 128 * ones(1, 32)];

inVec = inImg(:, 100);

switch type
  case 0
    blurVec = gaussSmooth(inVec, 0, sigma, P, extType);
  case 1
    blurVec = gaussSmooth(inVec, 1, sigma, P, extType);
  case 2
    blurVec = gaussSmooth(inVec, 1, sigma, P, extType);
  case 3
    blurVec = gaussSmooth(inVec, 2, sigma, P, extType);
end

figure(1)
plot(inVec);
figure(2)
plot(blurVec);


switch type
  case 0
    blurImg = gaussSmooth(inImg, 0, sigma, P, extType);
  case 1
    [blurImg blurImgY] = gaussSmooth(inImg, 1, sigma, P, extType);
  case 2
    [blurImgX blurImg] = gaussSmooth(inImg, 1, sigma, P, extType);
  case 3
    blurImg = gaussSmooth(inImg, 2, sigma, P, extType);
end

maxV = max(max(blurImg))
minV = min(min(blurImg))
sImg = (blurImg - minV) / (maxV - minV);

figure(3)
imshow(inImg);
figure(4)
imshow(sImg);

x = 1:200;
inLine = zeros(1, 200);
inLine(100) = 1;
blurLine = gaussSmooth(inLine, 0, sigma, P, extType);
m = sum(x .* blurLine)
v = sqrt(sum((x-m) .* (x -m) .* blurLine))
