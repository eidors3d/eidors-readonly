function transfimp = calc_transferimpedance( img)
% CALC_TRANSFERIMPEDANCE: Calculates transfer impedance matrix
% 
%   transfimp = calc_transferimpedance( img)
%
% fwd_model is a fwd_model structure
% img       is an image structure
%
% transfimp calculates electrode voltages based on electrode currents as
%     V = transfimp*I
% For example, for 4 electrodes, the voltage across [1,2] when 3A is across [3,4]
%    is given by [1,-1,0,0] * transfimp * [0;0;3;-3];

% (C) 2006 Bill Lionheart. License: GPL version 2 or version 3
% $Id$

% create new stim patterns
% stimulate with one ref electrode and then each in turn
% make an equiv set of measurements

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

copt.cache_obj = {img.fwd_model, img.elem_data};
copt.fstr = 'calc_transferimpedance';
transfimp = eidors_cache(@calc_T, img, copt);


function transfimp = calc_T( img);
   n_elecs= length( img.fwd_model.electrode );

    [stim_pat, meas_pat]= trigonometric( n_elecs );
   %[stim_pat, meas_pat]= electrode_wise( n_elecs );
   %[stim_pat, meas_pat]= monopolar( n_elecs );
   %[stim_pat, meas_pat]= monopolar_even( n_elecs );
   img.fwd_model.stimulation = stim_pat;

   imeas_pat= pinv(meas_pat);

   data = fwd_solve(img);

   sz= length(img.fwd_model.stimulation);
   transfimp = reshape( data.meas, sz, sz);
   transfimp = imeas_pat * transfimp * imeas_pat';

function [stim_pat, meas_pat] = trigonometric( n_elecs )
    stim_pat = struct;
    idx= linspace(0,2*pi,n_elecs+1)'; idx(end)= [];
    omega= idx*[1:n_elecs/2];
    meas_pat= [cos(omega), sin(omega) ]';
    for i=1:n_elecs
        stim_pat(i).stimulation='Amp';
        stim_pat(i).stim_pattern= meas_pat(i,:)';
        stim_pat(i).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = electrode_wise( n_elecs)
    stim_pat = struct;
    meas_pat= [-ones(n_elecs-1,1), speye(n_elecs-1)];
    for i=2:n_elecs
        stim_pat(i-1).stimulation='Amp';
        stim_pat(i-1).stim_pattern= sparse([1,i],1,[-1,1],n_elecs,1);
        stim_pat(i-1).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = monopolar( n_elecs)
    stim_pat = struct;
    meas_pat= speye(n_elecs);
    for i=1:n_elecs
        stim_pat(i).stimulation='Amp';
        stim_pat(i).stim_pattern= sparse(i,1,1,n_elecs,1);
        stim_pat(i).meas_pattern= meas_pat;
    end

function [stim_pat, meas_pat] = monopolar_even( n_elecs)
    stim_pat = struct;
    meas_pat= eye(n_elecs) - ones(n_elecs)/n_elecs;
    for i=1:n_elecs
        stim_pat(i).stimulation='Amp';
        stim_pat(i).stim_pattern= meas_pat(i,:)';
        stim_pat(i).meas_pattern= meas_pat;
    end

function do_unit_test
   current = 4; measure=1;
   [R,img] = test_2d_resistor(current,measure)
   T = calc_transferimpedance(img);
   unit_test_cmp('2D resistor',[1,-1]*T*[1;-1],sum(R), 1e-10);
   

function [R,img] = test_2d_resistor(current,measure)
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 3; len = 12; 

   fmdl=mk_grid_model([],linspace(0,wid,3), linspace(0,len,4));
   fmdl.electrode(1).nodes = find(fmdl.nodes(:,2) ==   0);
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,2) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R = len / wid / conduc;
   Contact_R = z_contact/wid;
   R = [Block_R, 2*Contact_R];
