function [cur_move] = line_optimize_diff_tangential(dx,...
    data0, img_h,img, Jm,prior_e,cur_move, opt)

%% Bracket the local minimum - Get 3 points such that
%    0 < s1 < s2 < s3 with f(s2) \leq min{ f(s1),f(s3) }
%Step 1: Get step_length s.t. f(step_length) < f(0)
step_length=1;

%Now do a coarse2fine on this
vsim = fwd_solve(img_h);
%Linearise movement v = v + Jm*current
vsim.meas=vsim.meas + Jm*cur_move;
   
%calculate the misfit at 0
mlist(1) = feval(opt.objective_func,data0,vsim, ...
    cur_move,prior_e,img_h,img,opt);   

%Iterate through
for i=1:10 %while loop here   
   %Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
   vsim = fwd_solve(img_h);
   
   %We linearise movement v = v + Jm*current + Jm*dx(m)
   vsim.meas=vsim.meas + Jm*cur_move*step_length;
   
   %calculate the misfit  
   mlist(i+1) = feval(opt.objective_func,data0,vsim, ...
       cur_move*step_length,prior_e,img_h,img,opt);   
     
   if(mlist(i+1) < mlist(1))       
       break; %Bracket now [i-1,i,i+1]
   else
       %Decrease step length and try again
       step_length=step_length/2;
   end  
end

%% Quadratic polynomial with 3 points
%Define some stuff
lambda(1)=0; lambda(2)=step_length;
   
%Image at current iterate
cur_movek=cur_move;

%Current iterate
vsim = fwd_solve(img_h);
%Linearise movement v = v + Jm*current
vsim.meas=vsim.meas + Jm*cur_move;

%calculate the misfit at 0
mlist(1) = feval(opt.objective_func,data0,vsim,...
    cur_movek,prior_e,img_h,img,opt);   

%Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
cur_movei = cur_movek + step_length*dx;

%Now do a coarse2fine on this
vsim = fwd_solve(img_h);
   
%We linearise movement v = v + Jm*current + Jm*dx(m)
vsim.meas=vsim.meas + Jm*cur_movei;
   
%calculate the misfit  
mlist(2) =  feval(opt.objective_func,data0,vsim,...
    cur_movei,prior_e,img_h,img,opt);  
   
%0 < s1 < s2 < s3 with f(s2) \leq min{ f(s1),f(s3) }
for i=1:10 %How many
   %Increase step length
   lambda(i+2) = lambda(i+1) + step_length;
   
   %Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
   cur_movei = cur_movek + lambda(i+2)*dx;

   %Now do a coarse2fine on this
   vsim = fwd_solve(img_h);
   %Linearise movement v = v + Jm*dx(m)
   vsim.meas=vsim.meas + Jm*cur_movei;
   
   %calculate the misfit  
   mlist(i+2) = feval(opt.objective_func,data0,vsim,...
       cur_movei,prior_e,img_h,img,opt);   
     
   if(mlist(i+2) > mlist(i+1))      
       break; %We have a bracket now [i-1,i,i+1]
   else
      %Continute iterations 
   end
   
end

%We have our bracket now
opt.perturb = [lambda(i),lambda(i+1),lambda(i+2)]
mlist = [mlist(i),mlist(i+1),mlist(i+2)]

% perform quadratic line fit in log space
pf = polyfit(opt.perturb, mlist, 2);
fmin = -pf(2)/pf(1)/2 % poly minimum for a 2nd order poly

% RETURN VALUES
cur_move = cur_movek + fmin*dx;