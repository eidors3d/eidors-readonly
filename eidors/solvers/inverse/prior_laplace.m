function Reg= prior_laplace( inv_model )
% PRIOR_LAPLACE calculate image prior
% Reg= prior_laplace( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%  or
% Reg= prior_laplace( fwd_model )
%
% This image prior is intended to be used as
%  R'*R, but may be used as R for as well.
%
% The Laplace prior is a 2nd order high pass filter.
% On a rectangular mesh, it is a convolution with
%   [-1,-1,-1;      [ 0;-1; 0
%    -1, 8,-1    or  -1; 4;-1
%    -1,-1,-1]        0;-1; 0]
%
% On a finite element mesh, we define the it as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself. 
%
% The EIDORS implementation is equivalent to twice the Graph-Laplacian

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

switch inv_model.type
    case 'inv_model'; fwd_model = inv_model.fwd_model;
    case 'fwd_model'; fwd_model = inv_model;
    otherwise; error('PRIOR_LAPLACE requires input type of inv_model or fwd_model');
end

copt.cache_obj = cache_obj(fwd_model);
copt.fstr = 'prior_laplace';

Reg = eidors_cache(@build_laplace, fwd_model, copt);

% new implementation builds by running through the faces provided by
% fix_model.
function Reg = build_laplace(fwd_model)

% obtain face2elem matrix
fmopt.face2elem=true;
fwd_model=fix_model(fwd_model,fmopt);

n_elem=size(fwd_model.elems,1);

% limit to interior faces
interiorindx=fwd_model.face2elem(:,2)~=0;
facei=fwd_model.face2elem(interiorindx,1);
facej=fwd_model.face2elem(interiorindx,2);

% build Idx as Iidx= [Iidx, ii, ii, jj, jj]; for each face
Iidx=kron(facei,uint32([1;1;0;0]))+kron(facej,uint32([0;0;1;1]));

% build Jidx as Jidx= [Jidx, ii, jj, ii, jj]; for each face
Jidx=kron(facei,uint32([1;0;1;0]))+kron(facej,uint32([0;1;0;1]));

% assign values to be consistent with previous implementation
Vidx=kron(ones(size(facei)),[2;-2;-2;2]);

% build
Reg = sparse(Iidx,Jidx, Vidx, n_elem, n_elem );



% Mapping depends only on elems - remove the other stuff
function c_obj = cache_obj(f_mdl)
c_obj = { f_mdl.elems, f_mdl.nodes};


function do_unit_test

imdl = mk_common_model('a2c2',16);
RtR = prior_laplace( imdl );
subplot(221); spy(RtR);
unit_test_cmp('a2c2', nnz(RtR), 240);


fmdl = mk_circ_tank(2,[],4);
RtR = prior_laplace( fmdl );
subplot(222); spy(RtR);
unit_test_cmp('2-4',RtR(1:4,1:8), [ ...
    6, -2,  0, -2, -2,  0,  0,  0; -2,  6, -2,  0,  0, -2,  0,  0;
    0, -2,  6, -2,  0,  0, -2,  0; -2,  0, -2,  6,  0,  0,  0, -2]);

compmdl = mk_common_model('n3r2',[16 2]);
RtRnew = prior_laplace( compmdl);
RtRold = prior_laplace_old(compmdl);
subplot(223);spy(RtRnew);
subplot(224);spy(RtRold);
unit_test_cmp('consistent behaviour',RtRnew,RtRold);

