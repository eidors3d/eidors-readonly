function image= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% $Id: np_inv_solve.m,v 1.3 2004-07-18 02:46:40 aadler Exp $

% calculate parameters from input structures
fwd_model= inv_model.fwd_model;
vtx= fwd_model.nodes;
simp= fwd_model.elems;
gnd_ind= fwd_model.gnd_node;
% calc num electrodes, nodes, stim_patterns
n_elem= size(simp,1);
n_elec=  length(fwd_model.electrode );
n_nodes= size(fwd_model.nodes,1);
n_stim = length(fwd_model.stimulation );
n_meas = 0;

% Recreate 'df' from fwd_model.stimulation
df= zeros(n_stim,1);
for i=1:n_stim;
    df(i) = size(fwd_model.stimulation(i).meas_pattern ,1);
    n_meas = n_meas + df(i);
end

elec= zeros(length(fwd_model.electrode ), ...
            length(fwd_model.electrode(1).nodes) );
zc  = zeros(length(fwd_model.electrode ), 1);

for i=1:length(fwd_model.electrode);
    elec(i,:)= fwd_model.electrode(i).nodes;
    zc(i)    = fwd_model.electrode(i).z_contact;
end

% Recreate 'indH' from fwd_model.stimulation
indH= zeros(n_stim, 2);
idx=0;
for i=1:n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern';
   elpos= rem( find(meas_pat(:)== 1)-1 , n_elec) + 1;
   indH( idx+(1:df(i)),1 ) = elpos;
   elpos= rem( find(meas_pat(:)==-1)-1 , n_elec) + 1;
   indH( idx+(1:df(i)),2 ) = elpos;
   idx= idx+ df(i);
end

% calculate FEM RHS matrix, i.e., the current patterns padded with zeroes 
I = zeros( n_elec + n_nodes, n_stim );
idx=0;
for i=1:n_stim
   I( n_nodes + (1:n_elec), i ) = ...
         fwd_model.stimulation(i).stim_pattern;
end
I(fwd_model.gnd_node,:) = 0;
Ib= I( n_nodes + (1:n_elec), : );


% calc jacobian with homogeneous background
homg_img.type = 'image';
homg_img.elem_data= ones( n_elem ,1);
homg_img.fwd_model= fwd_model;

J = calc_jacobian( fwd_model, homg_img);

% Calculating a smoothing prior
[Reg] = iso_f_smooth(simp,vtx,3,1);

% Calculating a linear inverse solution

tfac= inv_model.hyperparameter;
dva= data1.meas - data2.meas;
sol = (J'*J +  tfac*Reg'*Reg)\J' * dva;


% create a data structure to return
image.name= 'solved by np_inv_solve';
image.elem_data = sol;
image.type = 'real conductivity differences';
image.fwd_model= fwd_model;
image.inv_model= inv_model;
