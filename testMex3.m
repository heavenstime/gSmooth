type    = 11; % type  : 0 smooth, 1 diff (two outputs), 2 : LOG
sigma   = 10.0;
extType = 1; % extType :  0: zero extension,  1: the value of edge is used for extension

PL      = 1:6;

inType = 2;
figC   = 2;
plotY  = 1;

if inType == 1
  inImg = zeros(1,400);
  inImg(1, 400) = 1.0;
  center = 1;
elseif inType == 2
  inImg = im2double(imread('img1.pgm'));
else
  inImg = zeros(500, 500);
  center = 250;
  inImg(:,center) = ones(500, 1);
end

for pN = 1:size(PL, 2)
  P = PL(pN);
  switch rem(type, 10)
    case 0
      blurImg = gaussSmooth(inImg, type, sigma, P, extType);
    case 1
      if size(inImg, 1) > 1 && size(inImg, 2) > 1
	[blurImg blurImgY] = gaussSmooth(inImg, type, sigma, P, extType);
      else
	blurImg = gaussSmooth(inImg, type, sigma, P, extType);
      end
    case 2
      blurImg = gaussSmooth(inImg, type, sigma, P, extType);
  end

  if figC == 1 
    figure(1);
    plot (inImg(center,:));
    figure(pN + 1);
    plot (blurImg(center,:));
  elseif figC == 2
    if rem(type, 10) == 1 && plotY == 1 
      blurImg = blurImgY;
    end
    maxV = max(max(blurImg));
    minV = min(min(blurImg));
    sImg = (blurImg - minV) / (maxV - minV);
    figure(1)
    imshow(inImg);
    figure(pN + 1)
    imshow(sImg);
  end
end
