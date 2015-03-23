function progress_msg(varargin)
%PROGRESS_MSG Progress messages and timing.
% 
% PROGRESS_MSG('msg',I,N,opt) where I = 0 initilises the messages. 
%   'msg'   [optional] string to print
%   I       [optional] must equal zero for initialisation
%   N       [optional] if present, prints 'I/N' rather than 'I%'
%   opt     [optional] structure with the following fields and defaults:
%      .log_level = 2        controls verbosity, see EIDORD_MSG for details
%      .final_msg = 'Done'   customises the test of the final message
%
% PROGRESS_MSG(I) displays I%
% PROGRESS_MSG(I,N) displays I/N
%
% PROGRESS_MSG(Inf) prints the final message, including the time elapsed
%  since initialisation, measured using tic/toc.
%
% Example:
%   progress_msg('The big loop:', 0,10);         % prints 0/10
%   for i = 1:10
%       progress_msg(i,10);     % e.g.   3/10
%       pause(.2); % think
%   end
%   progress_msg(Inf)           % final message, e.g. Done (2.015 s)
%
% SEE ALSO EIDORS_MSG, TIC, TOC

% (C) 2015 - Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<

t0 = tic;

if nargin > 0 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST')
    varargin(1) = [];
    do_unit_test(varargin{:});
    return
end

persistent pvar;

nargs = nargin;
opt = default_opts;
if nargs > 0 && isstruct(varargin{end})
    tmp = varargin{end};
    fld = fieldnames(tmp);
    for i = 1:numel(fld)
        opt.(fld{i}) = tmp.(fld{i});
    end
    nargs = nargs-1;
    varargin(end) = [];
end

msg = '';
if nargs > 0 && ischar(varargin{1})
    msg = [repmat(' ',1,opt.log_level-1) varargin{1}];
    nargs = nargs-1;
    varargin(1) = [];
end


if nargs == 0 % there was only a message
    nargs = 1;
    varargin(1) = {0}; % starting number
end

first_msg = false;
if nargs == 0 || varargin{1} == 0
    first_msg = true;
    pvar.log_level = opt.log_level;
    pvar.final_msg = opt.final_msg;
    pvar.timerVal = tic;
    pvar.own_time = 0;
    if nargs <= 1
        pvar.nChars = 4;
    else
        pvar.numLength = floor(log10(varargin{2})) + 1;
        pvar.nChars = 2 * pvar.numLength + 1;
    end
end

if pvar.log_level > eidors_msg('log_level');
    return
end

if ~isempty(msg)
    fprintf('%s\t',msg);
end

if nargs > 0 && isinf(varargin{1})
    print_final_msg(pvar.final_msg, pvar.nChars);
    print_time(pvar);
    if eidors_debug('query','progress_msg')
        fprintf('Self-time: %f\n',pvar.own_time);
    end
    return
end
    

if nargs == 1
    percentmsg(varargin{1},first_msg, pvar.nChars);
else
    outofmsg(varargin{1},varargin{2},first_msg,...
        pvar.numLength, pvar.nChars);
end



pvar.own_time = pvar.own_time + toc(t0);


end


function percentmsg(num, first, N)
    str = '%3d%%';
    if ~first
        str = [repmat('\b',1,N) str ];
    end
    fprintf(str, round(100*num));
end

function outofmsg(a,b,first,N,T)
    str = sprintf('%%%dd/%%%dd',N,N);
    if ~first
        str = [repmat('\b',1,T) str ];
    end
    
    fprintf(str,a,b);
end

function print_final_msg(msg, N)
    str = [repmat('\b',1,N) '%s'];
    fprintf(str, msg);
end

function print_time(pvar)
    fprintf(' (%.3f s)\n',...
        toc(pvar.timerVal) -  pvar.own_time);
end

function opt = default_opts
    opt.final_msg = 'Done';
    opt.log_level = 2;
end

function do_unit_test(N)
    eidors_msg('log_level',2);
    eidors_debug on progress_msg
    if nargin == 0
        N = 1:30;
    end
    for n = N
        switch n
            case 1
                progress_msg(0);

            case 2
                progress_msg('Test 2',0);
            case 3
                progress_msg('Test 3');
            case 4
                progress_msg
            case 5
                opt.final_msg = 'YUPPIE!';
                progress_msg(opt);
            case 10
                progress_msg('Test 10',0,10);
            case 11
                progress_msg(0,10);
            case 20
                opt.final_msg = 'Ready !!!';
                progress_msg('Test 20', 0,10,opt);
            case 21
                eidors_msg('log_level',1);
                progress_msg('Test 21', 0,10,opt);
            case 22
                opt.log_level = 1;
                progress_msg('Test 22', 0,10,opt);
            case 23
                opt.log_level = 4;
                eidors_msg('log_level',5);
                progress_msg('Test 23', 0,10,opt);
            case 24
                eidors_msg('log_level',2);
                progress_msg('Test 24', 0,10);
            otherwise
                continue
        end 
        if n < 10
            for i = 1:10
                pause(.2);
                progress_msg(i/10);
            end
        else
            for i = 1:10
                pause(.2);
                progress_msg(i,10);
            end
        end
        progress_msg(Inf);

    end
    eidors_debug off progress_msg
end