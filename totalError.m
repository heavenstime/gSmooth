K     = 256;
N     = 2 * K + 1;
P     = 4; 
wL    = 0.8 * [0.9 0.95 1.0 1.05 1.1]; % Sigma list
% wL    = 0.5:0.05:2.0; % Sigma list
Ex    = 5;                             % Extension 

nW      = size(wL, 2);
KEx     = Ex * K;
NEx     = 2 * KEx + 1;
posL    = (K * (Ex - 1) + 1):(K * (Ex + 1) + 1);
posExtL = [1:(K * (Ex - 1))  (K * (Ex + 1) + 2):NEx]; % Position for extension

thetaL = (2 * pi / (2 * K) ) * [-K:K];
cosLL  = ones(P + 1, N);
sinLL  = zeros(P + 1, N);
for p = 1:P
  cosLL(p + 1, :) = cos(p * thetaL);
  sinLL(p + 1, :) = sin(p * thetaL);
end

coefGL   = zeros(1, (P + 1) * nW);
coefDGL  = zeros(1, P * nW);
coefDDGL = zeros(1, (P + 1) * nW);

fpE = fopen("errorForPlot.txt", "w");
fprintf(fpE, "# Truncation error \n");
fprintf(fpE, "# sigma  G:interval     order      total  GD:interval     order      total   GDD:interval   order      total \n");

for posW = 1:nW
  w        = wL(posW);
  beta     = 1.0 / (2.0 * w * w);
  thetaExL = (pi / K) * [-KEx:KEx];

  gaussExL       = exp(- beta * (thetaExL .* thetaExL)) * sqrt(beta / pi);
  dotGaussExL    = - 2 * beta * (gaussExL .* thetaExL) ;
  dotdotGaussExL = - 2 * beta * gaussExL + 2 * beta * 2 * beta * (gaussExL .* thetaExL .* thetaExL);
  
  gaussExtL       = gaussExL(posExtL);
  dotGaussExtL    = dotGaussExL(posExtL);
  dotdotGaussExtL = dotdotGaussExL(posExtL);
  extErrG         = gaussExtL * gaussExtL'; 
  extErrDG        = dotGaussExtL * dotGaussExtL'; 
  extErrDDG       = dotdotGaussExtL * dotdotGaussExtL';

  gaussL       = gaussExL(posL);
  dotGaussL    = dotGaussExL(posL);
  dotdotGaussL = dotdotGaussExL(posL);

  if 1 == 0  % orthogonal transform in [0, 2K - 1]
    cosLL(1, :) = ones(1, N) / sqrt(2);
    coefG   = cosLL(:, 1:(N-1)) * gaussL(1:(N-1))' / K;
    coefDG  = sinLL(2:(P+1), 1:(N-1)) * dotGaussL(1:(N-1))' / K;
    coefDDG = cosLL(:, 1:(N-1)) * dotdotGaussL(1:(N-1))' / K;
  else       % Mean square error in [0, 2K]
    coefG   = inv(cosLL * cosLL') * cosLL * gaussL';
    A       = sinLL(2:(P+1), :);
    coefDG  = inv(A * A') * A * dotGaussL';
    coefDDG = inv(cosLL * cosLL') * cosLL * dotdotGaussL';
  end
  coefGL( ((posW - 1) * (P + 1) + 1) :(posW * (P + 1)) )  = coefG;
  coefDGL( ((posW - 1) * P + 1)       :(posW * P) )       = coefDG;
  coefDDGL( ((posW - 1) * (P + 1) + 1):(posW * (P + 1)) ) = coefDDG;

  gaussA       = coefG'   * cosLL;  
  dotGaussA    = coefDG'  * sinLL(2:(P+1), :);  
  dotdotGaussA = coefDDG' * cosLL;  

  sqGauss       = gaussExL * gaussExL';
  sqDotGauss    = dotGaussExL * dotGaussExL';
  sqDotDotGauss = dotGaussExL * dotGaussExL';
  apErrG   = (gaussA       - gaussL)       * (gaussA       - gaussL)';
  apErrDG  = (dotGaussA    - dotGaussL)    * (dotGaussA    - dotGaussL)';
  apErrDDG = (dotdotGaussA - dotdotGaussL) * (dotdotGaussA - dotdotGaussL)';
  fprintf(fpE, "%6.2f  %10.4f %10.4f %10.4f  ", w, 100 * sqrt(extErrG / sqGauss), 100 * sqrt(apErrG / sqGauss), 100 * sqrt((extErrG + apErrG) / sqGauss) );
  fprintf(fpE, "%10.4f %10.4f %10.4f  ",  100 * sqrt(extErrDG / sqDotGauss), 100 * sqrt(apErrDG / sqDotGauss), 100 * sqrt((extErrDG + apErrDG) / sqDotGauss) );
  fprintf(fpE, "%10.4f %10.4f %10.4f  \n", 100 * sqrt(extErrDDG / sqDotDotGauss), 100 * sqrt(apErrDDG / sqDotDotGauss), 100 * sqrt((extErrDDG + apErrDDG) / sqDotDotGauss) );
%  fprintf(fpE, "%f  %f %f %f   %f %f %f   %f %f %f \n", w, sqrt(extErrG / sqGauss), sqrt(apErrG / sqGauss), sqrt((extErrG + apErrG / sqGauss)));
end
fclose(fpE);

% Print coefficients for C program 
fpC = fopen("coefficientsForC.txt", "w");
pos = 1;
fprintf(fpC, "{");
for wPos = 1:nW
  signF = 1;
  first = 1;
  for p = 1:(P + 1)
    if first == 1
      fprintf(fpC, "%15.10e", signF * coefGL(pos));
      first = 0;
    else
      fprintf(fpC, ", %15.10e", signF * coefGL(pos));
    end 
    signF = - signF;
    pos = pos + 1;
  end
end
fprintf(fpC, "}\n");
pos = 1;
fprintf(fpC, "{");
for wPos = 1:nW
  signF = 1;
  first = 1;
  for p = 1:P
    if first == 1
      fprintf(fpC, "%15.10e", signF * coefDGL(pos));
	first = 0;
    else
      fprintf(fpC, ", %15.10e", signF * coefDGL(pos));
    end 
    signF = - signF;
    pos = pos + 1;
  end
end
fprintf(fpC, "}\n");
pos = 1;
fprintf(fpC, "{");
for wPos = 1:nW
  signF = 1;
  first = 1;
  for p = 1:(P + 1)
    if first == 1
      fprintf(fpC, "%15.10e", signF * coefDDGL(pos));
      first = 0;
    else
      fprintf(fpC, ", %15.10e", signF * coefDDGL(pos));
    end 
    signF = - signF;
    pos = pos + 1;
  end
end
fprintf(fpC, "}\n");
fclose(fpC);

% File of data to plot Gausfunction functions (of the last w). (Approximated and True) 
fp = fopen("dataPlotFunc.txt", "w");
for posE = (K * (Ex - 2) + 1):(K * (Ex - 1));
  fprintf(fp, "%e %e %e\n", thetaExL(posE), 0.0, gaussExL(posE));
end
pos = 1;
for posE = (K * (Ex - 1) + 1):(K * (Ex + 1) + 1);
  fprintf(fp, "%e %e %e\n", thetaExL(posE), gaussA(pos), gaussExL(posE));
  pos = pos + 1;
end
for posE = (K * (Ex + 1) + 2):(K * (Ex + 2) + 2);
  fprintf(fp, "%e %e %e\n", thetaExL(posE), 0.0, gaussExL(posE));
end
fclose(fp);
