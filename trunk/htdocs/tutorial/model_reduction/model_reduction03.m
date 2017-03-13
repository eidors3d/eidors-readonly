M = inf; % resistor scaling where inf = unrestricted
tol = 0.05; %  5% tolerance resistors
calc_err = @model_reduction03_calc_err;
stdres = @model_reduction03_stdres;

% rough search, oscillatory function
mm = 10.^[0:tol/sqrt(2)/50:10];
err = calc_err(RR, mm, tol);
if ~isinf(M);
   em = M;
else
   [~,ii] = min(err); em = mm(ii);
end

% now refine the search within  that range using fminbnd
[f_mm,f_err,flag] = ...
   fminbnd(@(x) log10(calc_err(RR, x, tol)), em*(1-2*tol), em*(1+2*tol));
assert(flag == 1); % success?

% final result
! rm model_reduction03.txt
diary model_reduction03.txt
stdres(RR*f_mm,tol)
diary off

% plot fitting error
mm(end+1) = f_mm; err(end+1) = 10^f_err;
[mm, jj]=sort(mm); err=err(jj);
i=find(mm==em); [me,mi]=min(err);
clf; h=loglog(mm,err,'.--'); hold on;
loglog(mm(i),err(i),'or'); loglog(mm(mi),me,'sk');
mid = floor(length(mm)*2/5);
text(mm(mid),tol/3,sprintf('%d%% resistors',tol*100),'color','red');
legend('test points','approx min','fminbnd');
legend('location','best');
loglog([mm(1),mm(end)],[tol tol],'--r');
xlabel('multiplier'); ylabel('std resistor error [%]');
hold off;
print_convert model_reduction03.png
