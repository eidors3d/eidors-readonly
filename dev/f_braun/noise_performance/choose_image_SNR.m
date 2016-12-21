function HP = choose_image_SNR(imdl)
%% CHOOSE_IMAGE_SNR: choose hyperparameter based on image SNR calculation
% as proposed by Braun et al. in:
% F Braun et al., A Versatile Noise Performance Metric for Electrical
% Impedance Tomography Algorithms, IEEE Trans. Biomed. Eng. 2017 (submitted).
%
%   HP = choose_image_SNR(imdl)
%
% Output:
%   HP   - hyperparameter selected to achieve the desired image SNR
%
% Input:
%   imdl            - inverse model (EIDORS struct)
%      imdl.hyperparameter.image_SNR - desired target image SNR
%
% NOTE
%  This script was adapted from the function CHOOSE_NOISE_FIGURE
%
%
% See also: CALC_IMAGE_SNR, CHOOSE_NOISE_FIGURE
%
% Fabian Braun, December 2016
%

% (C) 2016 Fabian Braun. License: GPL version 2 or version 3
% $Id$


if ~isfield(imdl.hyperparameter, 'SNR_func')
    imdl.hyperparameter.SNR_func = @calc_image_SNR;
end

SnrGoal = imdl.hyperparameter.image_SNR;
assert(SnrGoal > 0, 'SNR is expected to be greater than zero!');

hpInitValue = 1E-1;
try % remove the value field
   hpInitValue = imdl.hyperparameter.value;
   imdl.hyperparameter = rmfield(imdl.hyperparameter,'value');
end

copt.boost_priority = 3;
copt.fstr = 'choose_image_SNR';
HP = eidors_cache(@HP_for_SNR_search,{SnrGoal,imdl,hpInitValue},copt);



function HP = HP_for_SNR_search(SnrGoal,imdl,hpInitValue)
   llv = eidors_msg('log_level');
   if llv==3; eidors_msg('log_level',2); end % HIDE A LOT OF MESSAGES

   if true
       % new version first finding maximal possible SNR
       % first of all look for maximal possible SNR
       % we can't go higher, for once "the sky is the limit" does not apply :-)
       fms_opts.MaxIter = 15;
       fms_opts.TolFun = 0.01 * SnrGoal;   % we don't need to be more accurate
       [HpMaxSnr, NegSnrMax] = fminsearch(@(x) -calc_image_SNR(imdl, x), hpInitValue, fms_opts);
       % TODO: how to ensure we won't fall into a local extremum?
       
       assert(SnrGoal <= abs(NegSnrMax), ...
            sprintf('Desired SNR (%0.4d) exceeds maximal possible SNR of %0.4d!', SnrGoal, abs(NegSnrMax)));
   else
       % old version with very dummy initialization
       HpMaxSnr = hpInitValue;
   end
   
   % now go and start from this to find best one
   [hp, x, y] = search1(SnrGoal, imdl, -log10(HpMaxSnr));
   eidors_msg('Crossed target value, now refining search in this region...', 2);

   dx= hp-linspace(-0.7,0.7,5);
   [hp, Snr] = search2(SnrGoal, imdl, dx);

   dx= hp-linspace(-0.2,0.2,5);
   [hp, Snr] = search2(SnrGoal,imdl, dx);

   dx= hp-linspace(-0.1,0.1,5); 
   [hp, Snr] = search2(SnrGoal,imdl, dx);

   dx= hp-linspace(-0.05,0.05,21);
   [hp, Snr] = search2(SnrGoal,imdl, dx);
   
   % TODO: verify that SNR = f(lambda) has positive derivative
   % else we're in the wrong range (after SNR max)

   HP= 10^-hp;
   assert(HP < HpMaxSnr, 'Something went wrong we''re on the wrong side of the curve!');
   eidors_msg('log_level',llv); %% Reset log_level
   
   % TODO: check that we achieved our goal without too much deviation
  
function [hp, x, y] = search1(SnrGoal, imdl, hp)
  [Snr,SE]=imdl.hyperparameter.SNR_func(imdl, 10^(-hp));
  hpInitial = hp;
  SnrInitial = Snr;
  y = Snr;
  x = 10^(-hp);
   if     Snr > SnrGoal; 
       dir = 1;   % we're too high >> decrease hp
   elseif Snr < SnrGoal; 
       dir = -1;  % we're too low >> increase hp
       warning('this should not happen!!!');
   else   dir = 0; end
   
%    slope = dir;
   while  (dir*Snr+3*SE > dir*SnrGoal)% || (sign(slope) ~= dir)
     hp= hp+0.5*dir;
     [Snr,SE]=imdl.hyperparameter.SNR_func( imdl, 10^(-hp));
     y = [y Snr];
     x = [x 10^(-hp)];
%      slope = (Snr - SnrInitial)/(10^(-hp) - 10^(-hpInitial));
     
     if Snr < 0     % TODO: verify that this approach works in all situations
       % NO GOOD: we're getting in the negative t-Value zone 
       % --> fall back to initial hp and inverse direction
       hp = hpInitial;
       dir = -dir;
       Snr = inf;
       eidors_msg('Falling back to initial hp and reversing search direction due to negative SNRs.', 2);
     end
   end

function [hp, Snr] =search2(SnrGoal, imdl, dx)
   hp = [];
   it = 0;
   Snr = zeros(1,length(dx));
   while isempty(hp) && it < 2  % TODO: one could increase the number of iterations and reduce Nnoise in calc_image_SNR
       it = it+1;
       %if it > 1 keyboard, end
       for k=1:length(dx)
           Snr(k)=Snr(k)+imdl.hyperparameter.SNR_func( imdl, 10^-dx(k));
       end
       if 0 % TODO: validate which behaviour is more accurate!
         log_snr = log10(Snr/it);
         p= polyfit( dx, log_snr-log10(SnrGoal), 3);
       else
         p= polyfit( dx, (Snr/it)-SnrGoal, 3);
       end
       hp = roots(p);
       %eliminate complex roots
       hp = hp(imag(hp)==0);
       %eliminate out of range roots
       hp = hp( hp<max(dx) & hp>min(dx) );  %USE if poly>1
       % pick the root with the smallest derivative
       if numel(hp) >1
           p2 = p.*(numel(p)-1:-1:0); p2(end) = [];
           [jnk, pos] = min(abs(polyval(p2,hp)));
           hp = hp(pos);
       end
       
       if 0  % visualize for debug purposes
        figure(); subplot(121);
        plot(dx, Snr-SnrGoal, ':.')
        hold on;
        xs = linspace(min(dx), max(dx), 100);
        plot(xs, polyval(p, xs), 'r')
        
        subplot(122);
        plot(dx, log10(Snr/it)-log10(SnrGoal))
       end
   end
   if isempty(hp)
       %fallback
       [jnk,idx] = min(abs((Snr/it)-SnrGoal));
%        [jnk,idx] = min(abs(log_snr-log10(SnrGoal)));
       hp = dx(idx);
   end
