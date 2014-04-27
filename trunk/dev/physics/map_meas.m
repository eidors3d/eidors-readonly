% function b = map_meas(b, out, [N])
function b = map_meas(b, out, N)
   % UNIT_TEST?
   if isstr(b) && strcmp(b,'UNIT_TEST'); b = do_unit_test; return; end

   if ~isstruct(b)
     d = b;
     b = struct();
     b.meas = d;
     b.time = nan;
   end

   if ~isfield(b, 'measured_quantity')
     eidors_msg('warning: map_meas is assuming voltage measurements', 1);
     b.measured_quantity = 'voltage';
   end
   in = b.measured_quantity;

   if nargin < 3
     N = 1;
   end

   b.meas = map_meas_core(b.meas, N, in, out);
   if ~isfield(b, 'name')
     b.name = sprintf('converted %s --> %s', in, out);
   end
 
function b = map_meas_core(b, N, in, out);
   if strcmp(in, out) % in == out
      return; % do nothing
   end

   % resistivity to conductivity conversion
   % we can't get here if in == out
   if     strcmp(in, 'voltage') && strcmp(out, 'apparent_resistivity')
      if N == 1
         error('missing apparent resistivity conversion factor N');
      end
      b = N * b; % voltage -> apparent resistivity
   elseif strcmp(in, 'apparent_resistivity') && strcmp(out, 'voltage')
      if N == 1
         error('missing apparent resistivity conversion factor N');
      end
      b = N \ b; % apparent resistivity -> voltage
   % log conversion
   elseif any(strcmp({in(1:3), out(1:3)}, 'log'))
      % log_10 b -> b
      if strcmp(in(1:6), 'log10_')
         if any(b > log10(realmax)-eps) warning('loss of precision -> inf'); end
         b = map_meas_core(10.^b, N, in(7:end), out);
      % ln b -> b
      elseif strcmp(in(1:4), 'log_')
         if any(b > log(realmax)-eps) warning('loss of precision -> inf'); end
         b = map_meas_core(exp(b), N, in(5:end), out);
      % b -> log_10 b
      elseif strcmp(out(1:6), 'log10_')
         if any(b <= 0+eps) warning('loss of precision -> -inf'); end
         b = log10(map_meas_core(b, N, in, out(7:end)));
      % b -> ln b
      elseif strcmp(out(1:4), 'log_')
         if any(b <= 0+eps) warning('loss of precision -> -inf'); end
         b = log(map_meas_core(b, N, in, out(5:end)));
      else
         error(sprintf('unknown conversion (log conversion?) %s - > %s', in, out));
      end
   elseif any(strcmp({in(1:4), out(1:4)}, 'abs_'))
      % log_10 b -> b
      if strcmp(in(1:4), 'abs_')
         error('can not convert from absolute value');
      else
         b = abs(map_meas_core(b, N, in, out(5:end)));
      end
   else
      error('unknown conversion %s -> %s', in, out);
   end

function pass = do_unit_test()
pass = 1;
d = 1;
while d ~= 1 & d ~= 0
  d = rand(1);
end
disp('TEST: map_meas()');
N = 1/15;
Ninv = 1/N;
% function b = map_meas(b, N, in, out)
elem_types = {'voltage', 'log_voltage', 'log10_voltage', ...
              'apparent_resistivity',  'log_apparent_resistivity',  'log10_apparent_resistivity'};
expected = [d         log(d)         log10(d)      N*d      log(N*d)      log10(N*d); ...
            exp(d)    d              log10(exp(d)) N*exp(d) log(N*exp(d)) log10(N*exp(d)); ...
            10.^d     log(10.^d )    d             N*10.^d  log(N*10.^d ) log10(N*10.^d ); ...
            Ninv*d      log(Ninv*d  )    log10(Ninv*d)   d         log(d)         log10(d); ...
            Ninv*exp(d) log(Ninv*exp(d)) log10(Ninv*exp(d)) exp(d) d              log10(exp(d)); ...
            Ninv*10.^d  log(Ninv*10.^d)  log10(Ninv*10.^d)  10.^d  log(10.^d)     d ];
for i = 1:length(elem_types)
  for j = 1:length(elem_types)
    if isnan(expected(i,j))
      fprintf('skipping map_meas(%s -> %s)\n', elem_types{i}, elem_types{j});
      continue;
    end
    pass = test_map_meas(d, N, elem_types{i}, elem_types{j}, expected(i,j), pass);
  end
end
if pass
   disp('TEST: overall PASS');
else
   disp('TEST: overall FAIL');
end

function pass = test_map_meas(data, N, in, out, expected, pass)
fprintf('TEST: map_meas(%s -> %s)\n', in, out);

b = struct;
b.measured_quantity = in;
b.meas = data;
b = map_meas(b, out, N);

if any(b.meas ~= expected)
   pass = 0;
   fprintf('TEST: FAIL map_data(%s -> %s)\n', in, out);
end
