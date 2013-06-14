function Reg = prior_combine(inv_model)


if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

opt = inv_model.prior_combine;
Reg = sparse(0);
flds = fieldnames(opt);
use  = regexp(flds,'^prior\d+$');
N    = length(cell2mat(use));
try 
   weights = opt.weights;
catch
   weights = ones(N,1);
end


count= 1;
for i = 1:length(flds)
   if ~isempty(use{i})
      tmp = inv_model;
      tmp.RtR_prior = opt.(flds{i}).func;
      tmp = copyfields(opt.(flds{i}), tmp);
      Reg = Reg + weights(count) * calc_RtR_prior(tmp);
      count = count + 1;
   end
end

function to = copyfields(from,to)
fld = fieldnames(from);
for i = 1:length(fld)
   if ~strcmp(fld{i},'func')
      to.(fld{i}) = from.(fld{i});
   end
end

function do_unit_test
imdl = mk_common_model('a2C',16);
imdl.RtR_prior = @prior_laplace;
R1 = calc_RtR_prior(imdl);
imdl.RtR_prior = @prior_tikhonov;
R2 = calc_RtR_prior(imdl);

imdl.RtR_prior = @prior_combine;
imdl.prior_combine.prior1.func = @prior_laplace;
imdl.prior_combine.prior2.func = @prior_tikhonov;
imdl.prior_combine.weights = [1 2];
Reg = calc_RtR_prior(imdl);
imagesc(full(Reg))

assert(all(all(Reg == (R1 + 2*R2))));