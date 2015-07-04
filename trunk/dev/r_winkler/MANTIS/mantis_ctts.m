function t = mantis_ctts(s,type)
%MANTIS_CTTS transformation of conductivity
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

switch type
  case 'alpha'
    t=0.25*s+ (0.25-1)./s;
  case 'log'
    t = -log(s);
  case 'sigma'
    t = s;
  case 'rho'
    t = 1./s;
  otherwise
    warning('Invalid transformation type: Using identity');
    t = s;
end

