function weights = mantis_ctweights(weightstype,Jacw,c,ww)
%MANTIS_CTWEIGHTS compute weights of the scalar product
%
% (C) 2015 Robert Winkler. License: GPL version 2 or version 3

switch weightstype
  case 'Ssigma'
    weights = sqrt(sum(Jacw.^2))'.*(abs(ww)./c);
  case 'S'
    weights = sqrt(sum(Jacw.^2))';
  case 'sigma'
    weights = abs(ww)./c;
  case 'uniform'
    weights = ones(size(Jacw,2),1);
  otherwise
    warning('Weights not set. Using uniform weights.');
    weights = ones(size(Jacw,2),1);
end

end

