% These functions are from the EIDORS2D suite, downloaded from 
% http://venda.uku.fi/~mvauhkon/MatlabEIT2d.html on 14 Sept 2005.
% Code is licensed under the GPL, as per the COPYING.TXT file.
%    - Andy Adler
%
% Functions for solving the EIT inverse problem in 2D 
%
% Jacobian computations
%
%  Jacobian          - Computes the Jacobian for 2D EIT when elementwise basis is used
%  JacobianNode      - Computes the Jacobian for 2D EIT when node basis is used
%  RemoveJacob       - Removes the derivatives corresponding to the current carrying electrodes
%
% Regularisation matrix computations
%
%  MakeRegmatrix     - Computes a regularisation matrix which includes smoothness assumptions
%  MakeRegmatrixNode - Computes a regularisation matrix which includes smoothness assumptions. Conductivity in linear basis.
%
% Dynamic reconstructions
%
%  Dyneit            - Estimates resistivity distributions by Kalman filter and Kalman smoother




