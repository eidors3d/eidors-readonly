function out = mono2stim(mat, stim)
% MONO2SIM transform monopolar stimulation into another stim pattern
%   S = MONO2STIM(M, STIM) is the measurement corresponding to a
%   stimulation and measurement pattern STIM, as created with e.g.
%   MK_STIM_PATTERNS, obtained from a transformation of the monopolar
%   measurments M. The monopolar stimulation pattern with which M was
%   calculated must be
%       MSP = mk_stim_patterns(E,R,'{mono}','{mono}', {'meas_current'},1);
%   where E is the number of electrodes per ring and R is the number of
%   rings. 
%   
%   STIM must use the same set of electrodes as MSP.
%   
%   M can be an EIDORS data struct with a 'meas' field, a Jacobian matrix, 
%   or a matrix of voltages.
%
% See also MK_STIM_PATTERNS, STIM_MEAS_LIST, FWD_SOLVE, CALC_JACOBIAN

% (C) 2022 Bartek Grychtol. License: GPL version 2 or 3
% $ Id: $
    
    if nargin==1 && ischar(mat) && strcmp(mat, 'UNIT_TEST')
        do_unit_test; return
    end
    
    m2s = [];
    for i = 1:numel(stim)
        m2s = [m2s ; kron(stim(i).stim_pattern', stim(i).meas_pattern)];
    end
    if isstruct(mat) && isfield(mat, 'meas')
        out = mat;
        out.meas = m2s * mat.meas;
    else
        out = m2s * mat;
    end
end

function do_unit_test
th= linspace(0,360,8+1); th(end)=[]; zz=0*th;
epos = [th; zz+0.5]';
fmdl = ng_mk_cyl_models(1,epos,0.1);
fmdl.stimulation = mk_stim_patterns(8,1,'{mono}','{mono}', {'meas_current'},1);
Jm= calc_jacobian(mk_image(fmdl,1));
vm = fwd_solve(mk_image(fmdl,1));
for SP=1:3
    for MP=1:3
        fprintf('SP=%d MP=%d\n',SP,MP);
        fmdl.stimulation = mk_stim_patterns(8,1,[0,SP],[0,MP],{'no_meas_current','rotate_meas'},1);
        Jo= calc_jacobian(mk_image(fmdl,1));
        Jx = mono2stim(Jm, fmdl.stimulation);
        unit_test_cmp('Equal',Jo,Jx, 10*eps);
        img = mk_image(fmdl,1);
        vo = fwd_solve(img);
        vx = mono2stim(vm, fmdl.stimulation);
        unit_test_cmp('Equal',vo.meas,vx.meas, 10*eps);
    end
end
end
