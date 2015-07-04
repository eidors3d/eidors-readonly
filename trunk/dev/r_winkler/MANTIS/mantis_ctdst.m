function s = mantis_ctdst(t,type)
%MANTIS_CTDST derivative of inverse transformation of conductivity
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

switch type
  case 'alpha'
    s=(1+t./sqrt(4*(1-0.25)*0.25+t.^2))./(2*0.25);
  case 'log'
    s = -exp(-t);
  case 'sigma'
    s = ones(size(t));
  case 'rho'
    s = -1./(t.^2);
  otherwise
    warning('Invalid transformation type: Using identity');
    s = ones(size(t));
end

end