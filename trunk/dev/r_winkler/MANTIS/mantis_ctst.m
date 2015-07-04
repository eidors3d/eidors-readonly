function s = mantis_ctst(t,type)
%MANTIS_CTST inverse transformation of conductivity
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

switch type
  case 'alpha'
    s=(t+sqrt(4*(1-0.25)*0.25+t.^2))./(2*0.25);
  case 'log'
    s = exp(-t);
  case 'sigma'
    s = t;
  case 'rho'
    s = 1./t;
  otherwise
    warning('Invalid transformation type: Using identity');
    s = t;
end

end