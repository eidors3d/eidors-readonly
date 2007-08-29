function [fwd_mdl]=set_fwd_model(vtx,simp,srf,elec,zc,gnd_ind,Ib,indH,df,perm_sym);
% SET_FWD_MODEL: create EIDORS v3 fwd_model from v2 parameters
% [fwd_mdl]=set_fwd_model(vtx,simp,srf,elec,zc,gnd_ind,Ib,indH,df,perm_sym);
%
% Optional parameters
%   srf = specify [] to automatically generate
%   perm_sym = default is '{y}'
%
% Unknown parameters
%   specify '[]' to be used later
%
% (C) 2005 Stephen Murphy. Licenced under GPL Version 2
% $Id: set_fwd_model.m,v 1.8 2007-08-29 09:24:24 aadler Exp $

if nargin<10
    perm_sym= '{y}';
end
if isempty( srf)
    srf= find_boundary( simp );
end
fwd_mdl= eidors_obj('fwd_model','FWD_MDL created by set_fwd_model');
fwd_mdl.nodes=vtx;
fwd_mdl.elems=simp;
fwd_mdl.boundary=srf;
fwd_mdl.gnd_node=gnd_ind;

elec_mdl=[];
for loop1=1:size(elec,1);
    elec_mdl(loop1,1).z_contact=zc(loop1);
    elec_mdl(loop1,1).nodes=elec(loop1,:);
end

fwd_mdl.electrode=elec_mdl;

stim_mdl=[];
for loop1=1:size(df,1)
    stim_mdl(loop1,1).stimulation=abs(max(Ib(:,loop1)));
    stim_mdl(loop1,1).stim_pattern=Ib(:,loop1);
    meas=zeros(df(loop1),size(elec,1));
    loop3=1;
    for loop2=sum(df(1:loop1))-df(loop1)+1:sum(df(1:loop1))
        meas(loop3,indH(loop2,1))=1;
        meas(loop3,indH(loop2,2))=-1;
        loop3=loop3+1;
    end
    stim_mdl(loop1,1).meas_pattern=meas;
end

fwd_mdl.stimulation=stim_mdl;
fwd_mdl.misc.perm_sym=perm_sym;
fwd_mdl.solve='np_fwd_solve';
fwd_mdl.jacobian='np_calc_jacobian';
fwd_mdl.system_mat='np_calc_system_mat';
