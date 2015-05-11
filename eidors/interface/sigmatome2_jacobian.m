function J= sigmatome2_jacobian( fwd_model, img)
% JACOBIAN_FILTERED: J= jacobian_filtered( fwd_model, img)
%
% Filter a jacobian matrix by a specified filter function
% INPUT:
%  fwd_model = forward model
%  fwd_model.jacobian_filtered.jacobian = actual jacobian function (J0)
%  fwd_model.jacobian_filtered.filter   = Filter Matrix (F)
%  img = image background for jacobian calc
% OUTPUT:
%  J         = Jacobian matrix = F*J0

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

J0 = feval(fwd_model.jacobian_filtered.jacobian, ...
           fwd_model, img);
F  = sigmatome2_filter;

J= F*J0;


% Generate the filter function of the sigmatome system
function filter = sigmatome2_filter;
   Filter3= [0,         0.002992185737274, 0.002885352918496,-0.001992209356618, ...
     0.000247868161165, 0.001320969614758,-0.002994506380400, 0.000217287704644, ...
     0.000820760030266,-0.002388560206626, 0.005467093174149,-0.000465034175608, ...
    -0.014049400948852, 0.011968526958771, 0.017269336581277,-0.033327865383464, ...
    -0.000687996077078, 0.065884070897957,-0.055509942269519,-0.084217637095161, ...
     0.414790173689917, 0.947607669916644, 0.723749489957398, 0.088544295521622, ...
    -0.118281839343952, 0.032044694797133, 0.043642225555493,-0.035476390241121, ...
    -0.007655680300629, 0.022891001199580,-0.005076656322158,-0.007794874431280, ...
     0.006724761261718,-0.000064773581282,-0.002313251865051, 0.001234620138549, ...
    -0.002848224928827,-0.000772726894395, 0.002559973165365,-0.001297698809235, ...
     0.000569993353282, 0.002828133181959, 0.000426156082041];
   filter = sparse(416,208);
   filter(1:length(Filter3),3) = Filter3';
   for i=1:208; 
      filter(:,i) = circshift(filter(:,3),+2*(i-3));
   end
