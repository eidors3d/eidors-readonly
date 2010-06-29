function inv_mdl= select_imdl( mdl, options )
% SELECT_IMDL: select pre-packaged inverse model features
% inv_mdl= select_imdl( mdl, options )
%
%   mdl - inv_model structure - parameters are replaced as specified
%    OR
%   mdl - fwd_model - a basic GN difference model is created, and parameters replaced
%
% OPTIONS
%   options - {'opt1','opt2'} etc

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1; options = {}; end

switch mdl.type
  case 'inv_model'; inv_mdl = mdl;
  case 'fwd_model'; inv_mdl = basic_imdl( mdl );
  otherwise;        error('select_imdl: expects inv_model or fwd_model input');
end

for i=1:length(options);
  switch options{i}
    case 'NOSER dif';       inv_mdl = NOSER_dif( inv_mdl );
    case 'Basic GN dif';    inv_mdl = Basic_GN_Diff( inv_mdl );
    
    otherwise; error('option {%s} not understood', options{i});
  end
end

function imdl = basic_imdl( fmdl );
   imdl.name= 'Basic imdl from select_imdl';
   imdl.type= 'inv_model';

   imdl.solve= 'aa_inv_solve'
   imdl.hyperparameter.value = .01;
   imdl.RtR_prior = @laplace_image_prior;
   imdl.jacobian_bkgnd.value = 1;
   imdl.reconst_type= 'difference'
   imdl.fwd_model = fmdl;




function do_unit_test
   imdl = mk_common_model('a2c2',16); 
   imdl = select_imdl( imdl );
   imdl = select_imdl( imdl.fwd_model );

