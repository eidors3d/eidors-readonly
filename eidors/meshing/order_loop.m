function [p n] = order_loop(pp,clk)
%ORDER_LOOP Order a list of points on a loop
% P = ORDER_LOOP(PP) orders clockwise a matrix PP (N x D) of N points in D 
% dimensions that constitute a continues loop, under the assumption that 
% the distance between two neighbouring points is always smaller than 
% distance between any non-neighbouring points.

% P = ORDER_LOOP(PP, CLK) orders the loop clockwise if CLK == 1 and
% counter-clockwise if CLK == 0.

% (C) 2012 Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if ischar(pp) && strcmp(pp,'UNIT_TEST'); do_unit_test; return, end;

if nargin == 1
   clk = 1;
end

D = distmat(pp) + diag(inf*ones(length(pp),1));
p = zeros(size(pp));
n = zeros(size(pp,1),1);
N = 1:size(pp,1);
p(1,:) = pp(1,:);
idx1 = 1;
n(1) = 1;
for i = 1:length(pp)-1
   [jnk,idx2] = min(D(idx1,:));
   p(i+1,:) = pp(idx2,:);
   D(:,idx1) = [];
   D(idx1,:) = [];
   pp(idx1,:) = [];
   N(idx1)   = [];
   if idx2 > idx1
      idx1 = idx2-1;
   else
      idx1 = idx2;
   end
   n(i+1) = N(idx1);
end
ctr = mean(p);
tmp = p - repmat(ctr,length(p),1);
[th,r]=cart2pol(tmp(:,1),tmp(:,2)); % OCTAVE BUG in 3.7 - needs 2 outputs
th = unwrap(th);

if th(1) < th(end)
   %counter-clockwise
   p = flipud(p);
end

if clk == -1
   p = flipud(p);
end


function do_unit_test
t = linspace(0,2*pi,101); t(end) = [];
p1(:,1) = sin(t);
p1(:,2) = 4*cos(t);
n = randperm(100);
p2 = p1(n,:);
p3 = order_loop(p2);
expect1= circshift(p1,-n(1));
expect2= circshift(p1,-n(1)+1);
%% No idea why it could be either -- AA, May2015
unit_test_cmp('order_loop:t2', ...
    all(all(p3==expect1)) | all(all(p3== expect2)),1);


