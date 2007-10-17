function Reg= calc_covar_prior( inv_model )
% CALC_COVAR_PRIOR image prior with distance-based interelement covar
% This is a simplification of exponential_covar_prior.m
% Reg= calc_covar_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% P_type--prior type
% 1: elements are globally correlated
% 2: elements within/without electrode rings are correlated to elements in same region.
% 3: only elements within electrode rings are correlated.

% (C) 2007, Tao Dai and Andy Adler. Licenced under the GPL Version 2
% $Id: calc_covar_prior.m,v 1.1 2007-10-17 00:37:31 daitao Exp $

% get average x,y,z of each element
ff= inv_model.fwd_model;
nel= size(ff.elems,1);
eta = .1;%attenuation factor. eta is large when elems're spatially highly correlated

dist= zeros(nel);

z1 = ff.nodes(ff.electrode(1).nodes,3);%upper electrode ring
z2 = ff.nodes(ff.electrode(end).nodes,3);%lower electrode ring
%    ff.nodes(:,3)>z1||ff.nodes(:,3)<z2

for dim= 1: size(ff.nodes,2);
    coords= reshape(ff.nodes(ff.elems,dim),[size(ff.elems)]);%coord of four vertex of each elem
    m_coords= mean( coords,2);%calc center coord
    if dim == 3
        temp = double((m_coords<max(z1,z2))&(m_coords>min(z1,z2)));
        H = temp*temp';%(temp1*temp1')+(temp2*temp2');
        switch inv_model.fourD_prior.P_type
            case 1
                H = eta*(ones(size(H)));
            case 2
                temp = double(m_coords>max(z1,z2));
                H = H + temp*temp';
                temp = double(m_coords<min(z1,z2));
                H = H + temp*temp';
                H = eta*(H+1e-6);
            case 3
                H = eta*(H+1e-6);
            otherwise
                error('no such a 3-D prior type');
        end
    end
    difm = m_coords*ones(1,nel);
    difm = difm - difm';
    dist= dist + difm.^2;%L^2=dx^2+dy^2+dz^2
end
dist = sqrt(dist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = dist/max(dist(:));

Reg = exp(-dist ./ H);
