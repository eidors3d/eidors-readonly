% function Ro = stdres(Ri,tol)
% based on http://www.rfcafe.com/references/electrical/resistor-values.htm
% Ri = desired resistance
% tol = tolerance
% Ro = best fit std resistor
function Ro = stdres(Ri,tol)
   Ro = Ri;
   assert(all(all(Ri>=0)));
   switch tol
      case 0.01
         N = 96; % 1% -- EIA E96
         n = 3;  % significant digits
         d = [10.^[0:6] (1+0.1*[1:12])*1e6];
      case 0.02 
         N = 48; % 2% -- EIA E48
         n = 3;
         d = [10.^[0:7] (1+0.1*[1:12])*1e7];
      case 0.05
         N = 24; % 5% -- EIA E24
         n = 2;
         d = [10.^[0:7] (1+0.1*[1:12])*1e7];
      case 0.10
         N = 12; % 10% -- EIA E12
         n = 2;
         d = 10.^[0:6];
      otherwise
         error('bad tol')
   end
   n=n-1;
   i = [0:(N-1)];
   R = round(10.^(i/N),n);
   R = R' * d;
   for j = 1:length(Ri(:))
      if Ri(j) == 0; continue; end
      [~,ii] = min(abs(R(:)-Ri(j)));
      Ro(j) = R(ii);
   end
