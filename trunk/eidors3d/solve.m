function data = fwd_solve( fwd_model, image)
% FWD_SOLVE: calculate data from a fwd_model object and an image
%
% $Id: solve.m,v 1.1 2004-07-09 18:51:28 aadler Exp $

data= feval( fwd_model.solve, fwd_model, image);
