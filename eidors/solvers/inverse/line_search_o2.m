function  [alpha, img, dv, opt] = line_search_o2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt)
% function  [alpha, img, dv, opt] = line_search_o2(imgk, dx, data1, img1, N, W, hps2RtR, dv0, opt)
% function  [alpha, img, dv, opt] = line_search_o2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt)
% line search function with a fitted polynomial of O(2) -- a second order polynomial
% (C) 2013 Alistair Boyle
% License: GPL version 2 or version 3

pf_max = 2; % max poly-fit order = 2 --> second-order polynomial
retry = 0; % first any only call through this function, no recursion
if nargin == 9
   hpt2LLt = 1;
end
[alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt, retry, pf_max);
