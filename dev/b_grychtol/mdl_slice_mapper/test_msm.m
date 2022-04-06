function test_msm

test = localfunctions;
index = @(x,y) x(y);
idx = cellfun(@(x)strcmp(index(func2str(x),1:4), 'test'), test);
test = test(idx);

eidors_cache off elem_ptr
eidors_cache off node_ptr
eidors_cache off nodeinterp

fh = {@mdl_slice_mapper, @mdl_slice_mapper_dev};

N = 200;

for i = 1:numel(test)
    for j = 1:2
        [ptr{i,j}, tt(i,j)] = feval(test{i}, fh{j}, N);
    end
    fname = func2str(test{i}); fname(1:5) = [];
    name = blanks(15); name(1:numel(fname)) = fname;
    str = sprintf('%s ORIG: %0.4fs\t NEW: %0.4fs\t', name, tt(i,1), tt(i,j));
    res = unit_test_cmp(str, ptr{i,1}, ptr{i,2}, N^2 * eps);
    if ~res
        keyboard
    end
    
end

eidors_cache on elem_ptr
eidors_cache on node_ptr
eidors_cache on nodeinterp

function fmdl = mdl_2d(N)
    imdl = mk_common_model('f2t5',16); fmdl = imdl.fwd_model;
    fmdl.mdl_slice_mapper.npx = N;
    fmdl.mdl_slice_mapper.npy = N;

function fmdl = mdl_3d(N)
    imdl = mk_common_model('f3t5r',[16,1]); fmdl = imdl.fwd_model;
    fmdl.mdl_slice_mapper.npx = N;
    fmdl.mdl_slice_mapper.npy = N;
    fmdl.mdl_slice_mapper.level = [inf inf 1.58];

    
function [ptr, time] = test_elem_2d(fh, N)
tic
ptr = feval(fh, mdl_2d(N), 'elem');
time = toc;

function [ptr, time] = test_elem_3d(fh, N)
tic
ptr = feval(fh, mdl_3d(N), 'elem');
time = toc;

function [ptr, time] = test_node_2d(fh, N)
tic
ptr = feval(fh, mdl_2d(N), 'node');
time = toc;

function [ptr, time] = test_node_3d(fh, N)
tic
ptr = feval(fh, mdl_3d(N), 'node');
time = toc;

function [ptr, time] = test_nodeinterp_2d(fh, N)
tic
ptr = feval(fh, mdl_2d(N), 'nodeinterp');
time = toc;

function [ptr, time] = test_nodeinterp_3d(fh, N)
tic
ptr = feval(fh, mdl_3d(N), 'nodeinterp');
time = toc;
