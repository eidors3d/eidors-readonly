% DEMO to show really simple application of EIDORS
%
% This code models a resistor (resitance R), with one electrode at each end. 
%
% Stimulation patterns apply current I1, I2, I3, and measure voltages
% Forward Model: V= IR    => [V1;V2;V3] = [I1;I2*I3]*(R+2*Zc)
% Jacobian:      J= dV/dR =I = [I1; I2; I3]
% Inverse Model: R= inv(J'*J)*J'*V
%    This corresponds to the least squares solution

% $Id: resistor_model.m,v 1.1 2004-07-25 21:59:41 aadler Exp $

function resistor_model;

% Set log level to show all messages
eidors_msg('set_level',4);

%
% Step 1: create FEM model structure
%
% Fwd model:
%  Two nodes are in space at [1,1,1] and [2,2,2]
%  The resistor is connected between them

r_mdl.name = 'demo resistor model';
r_mdl.nodes= [1,1,1;  2,2,2];
r_mdl.elems= [1;2];
r_mdl.solve=      @f_solve;
r_mdl.jacobian=   @c_jacobian;

%
% Step 2: create FEM model electrode definitions
%

r_mdl.electrode(1).z_contact= 10; % ohms
r_mdl.electrode(1).nodes=     1;
r_mdl.gnd_node= 2;

%
% Step 3: create stimulation and measurement patterns
% patterns are 0.010,0.020,0.030 mA

for i=1:3
    r_mdl.stimulation(i).stimulation= 'mA';
    r_mdl.stimulation(i).stim_pattern= ( 0.010*i );
    r_mdl.stimulation(i).meas_pattern= 1; % measure electrode 1
end

r_mdl= eidors_obj('fwd_model', r_mdl);

%
% Step 4: simulate data for medium with R=1 kohms
% This medium is called an 'image'
%

img_1k = eidors_obj('image', 'homogeneous image', ...
                     'elem_data', 1e3, ...
                     'fwd_model', r_mdl );

data_1k =fwd_solve( r_mdl, img_1k );

%
% Step 5: add noise to simulated data
%

data_noise= eidors_obj('data', 'noisy data', ...
                       'meas', data_1k.meas + 1e-3*randn(3,1));

%
% Step 7: create inverse model
%

% create an inv_model structure of name 'demo_inv'
r_inv.name=  'Resistor Model inverse';
r_inv.solve= @i_solve;
r_inv.reconst_type= 'static';
r_inv.fwd_model= r_mdl;
r_inv= eidors_obj('inv_model', r_inv);

%
% solve inverse model');
%

R= inv_solve( r_inv, data_1k );
fprintf('R calculated with clean data= %5.3f kOhm\n', R.elem_data / 1000 );

R= inv_solve( r_inv, data_noise );
fprintf('R calculated with noisy data= %5.3f kOhm\n', R.elem_data / 1000 );


% Forward Model:
% For each stimulation there is I1 into Node1
%  Node2 is connected to gnd with Zcontact
%
% each stim has one measurement pattern where
%  Vmeas= Meas_pat * Node1
%       = Meas_pat * I1 * ( Zcontact*2 + R )
%
% Thus
%  V= IR    => [V1;V2;V3] = [I1;I2*I3]*(R + 2*Zcontact)
function data =f_solve( f_mdl, img )
  R= img.elem_data;

  n_stim= length( f_mdl.stimulation );
  V= zeros(n_stim, 1);

  for i=1:n_stim
    if ~strcmp( f_mdl.stimulation(i).stimulation, 'mA' )
       error('f_solve expects current in mA');
    end

    I        = f_mdl.stimulation(i).stim_pattern / 1000;
    meas_pat = f_mdl.stimulation(i).meas_pattern;

    stim_elec= find( I );
    zc       = f_mdl.electrode( stim_elec ).z_contact;

    V(i)= meas_pat * I * ( R + 2*zc);
  end

  data.name='resistor model data';
  data.meas= V;
  
% Jacobian:      J= dV/dR =I = [I1; I2; I3]
function J= c_jacobian( f_mdl, img)
  n_stim= length( f_mdl.stimulation );
  J= zeros(n_stim, 1);
  for i=1:n_stim
    J(i)     = f_mdl.stimulation(i).stim_pattern / 1000; % mA
  end

% Inverse Model: R= inv(J'*J)*J'*V
%    This corresponds to the least squares solution
function img= i_solve( i_mdl, data )
  % Normally the Jacobian depends on an image. Create a dummy one here
  i_img= eidors_obj('image','Unused');
  f_mdl= i_mdl.fwd_model;
  J = calc_jacobian( f_mdl, i_img); 

  img.name= 'solved by i_solve';
  img.elem_data= (J'*J)\J'* data.meas;
  img.inv_model= i_mdl;
  img.fwd_model= f_mdl;


