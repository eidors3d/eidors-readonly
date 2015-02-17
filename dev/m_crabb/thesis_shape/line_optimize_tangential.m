function [img] = line_optimize_tangential(imgk, dx, data0, ...
    img_h, c2f2, n_elemsc, Jm,img0, opt)
%Copy the image across
img = imgk;

%% Bracket the local minimum - Get 3 points such that
%    0 < s1 < s2 < s3 with f(s2) \leq min{ f(s1),f(s3) }
%Step 1: Get step_length s.t. f(step_length) < f(0)
step_length=1;

%Now do a coarse2fine on this
img_h.elem_data = c2f2*img.elem_data(1:n_elemsc);  
vsim = fwd_solve(img_h);
%Linearise movement v = v + Jm*current
vsim.meas=vsim.meas + Jm*img.elem_data(n_elemsc+1:end);
   
%Misfit at 0
mlist(1) = feval(opt.objective_func,data0,vsim,img,img0,opt);   

%Iterate through
for i=1:10 %while loop here   
   %Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
   img.elem_data = imgk.elem_data + step_length*dx;

   %Now do a coarse2fine on this
   img_h.elem_data = c2f2*img.elem_data(1:n_elemsc);  
   vsim = fwd_solve(img_h);   
   %Linearise movement v = v + Jm*current + Jm*dx(m)
   vsim.meas=vsim.meas + Jm*img.elem_data(n_elemsc+1:end);
   
   %calculate the misfit  
   mlist(i+1) = feval(opt.objective_func,data0,vsim,img,img0,opt);   
     
   if(mlist(i+1) < mlist(1))      
       break; %Bracket now [i-1,i,i+1]
   else
       step_length=step_length/2; %Else try again
   end  
end

%% Qudaratic minimum with 3 points
%Define some stuff
lambda(1)=0; lambda(2)=step_length;
   
%Image at current iterate
img.elem_data = imgk.elem_data;
img_h.elem_data = c2f2*img.elem_data(1:n_elemsc);  
vsim = fwd_solve(img_h);
%Linearise movement v = v + Jm*current
vsim.meas=vsim.meas + Jm*img.elem_data(n_elemsc+1:end);

%calculate the misfit at 0
mlist(1) = feval(opt.objective_func,data0,vsim,img,img0,opt);   

%Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
img.elem_data = imgk.elem_data + step_length*dx;

%Now do a coarse2fine on this
img_h.elem_data = c2f2*img.elem_data(1:n_elemsc);  
vsim = fwd_solve(img_h);
   
%We linearise movement v = v + Jm*current + Jm*dx(m)
vsim.meas=vsim.meas + Jm*img.elem_data(n_elemsc+1:end);
   
%calculate the misfit  
mlist(2) = feval(opt.objective_func,data0,vsim,img,img0,opt);   
   
%0 < s1 < s2 < s3 with f(s2) \leq min{ f(s1),f(s3) }
for i=1:10 %How many
    %Increase step length
   lambda(i+2) = lambda(i+1) + step_length;
   
   %Perturb image (s,v) |-> (s,v) +p*dx(1:n_elemsc)    
   img.elem_data = imgk.elem_data + lambda(i+2)*dx;

   %Now do a coarse2fine on this
   img_h.elem_data = c2f2*img.elem_data(1:n_elemsc);  
   vsim = fwd_solve(img_h);
   %We linearise movement v = v + Jm*dx(m)
   vsim.meas=vsim.meas + Jm*img.elem_data(n_elemsc+1:end);
   
   %calculate the misfit  
   mlist(i+2) = feval(opt.objective_func,data0,vsim,img,img0,opt);   
     
   if(mlist(i+2) > mlist(i+1))
       break; %Found bracket exit
   else %continue
   end
   
end

%We have our bracket now
opt.perturb = [lambda(i),lambda(i+1),lambda(i+2)];
mlist = [mlist(i),mlist(i+1),mlist(i+2)];

%Perform quadratic line fit in log space
pf = polyfit(opt.perturb, mlist, 2);
fmin = -pf(2)/pf(1)/2; % poly minimum for a 2nd order poly

% RETURN VALUES
img.elem_data = imgk.elem_data + fmin*dx;