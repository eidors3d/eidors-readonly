function [pts,FE2CE,FE2pts,CE2pts] = find_edge2edge_intersections(FE,FN,CE,CN, epsilon)
%FIND_EDGE2EDGE_INTERSECTIONS intersections between edges of two models
% Edges are considered intersecting if the minimum distance between them is
% less than epsilon and the closest point is not an endpoint.
% [pts,FE2CE,FE2pts,CE2pts] = find_edge2edge_intersections(FE,FN,CE,CN, epsilon)
% Inputs:
%   FE      - Fine model edges [Nx2] as indices into FN
%   FN      - Fine model nodes [Nx3]
%   CE      - Coarse model edges [Mx2]
%   CN      - Coarse model nodes [Mx3]
%   epsilon - the minimum distance recognised as intersection
% Outputs:
%   pts     - List of intersection points [Px3]
%   FE2CE   - Boolean matrix indicating if two edges interesct [NxM]
%   FE2pts  - Map between fine model edges and intersection points [NxP]
%   CE2pts  - Map between coarse model edges and intersection points [MxP]
%
% This function is inspired by the simple Matlab conversion by Cristian Dima 
% of the C code posted by Paul Bourke at
% http://paulbourke.net/geometry/pointlineplane/lineline.c
% http://paulbourke.net/geometry/pointlineplane/linelineintersect.m
%
% See also: FIX_MODEL, MK_GRID_C2F

% (C) 2015 Bartlomiej Grychtol - all rights reserved by Swisstom AG
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<

if size(CE,1) > size(FE,1)
    % flip the inputs if CE is longer than FE
    [pts,FE2CE,CE2pts,FE2pts] = find_edge2edge_intersections(CE,CN,FE,FN, epsilon);
    FE2CE = FE2CE';
    return
end

epsilon = epsilon.^2; % spares us doing sqrt to calculate distance

% TODO: need a good heuristic to choose the faster approach
try 
    % this is usually faster, but can still, exceptionally, run out of
    % memory
    [pts,FE2CE,FE2pts,CE2pts] = edge2edge_intersections_serial(FE,FN,CE,CN,epsilon);
catch err
    if strcmp(err.identifier, 'MATLAB:nomem')
        % this should not run out of memory
        [pts,FE2CE,FE2pts,CE2pts] = edge2edge_intersections_wrapper(FE,FN,CE,CN,epsilon);
    else
        rethrow(err);
    end
end


%-------------------------------------------------------------------------%
% Wrapper to divide the calculations into smaller chunks
function [intpts, FE2CE, FE2pts, CE2pts] = edge2edge_intersections_wrapper(FE,FN,CE,CN, epsilon)
% Doing this in parallel creates potentially huge matrices
% Do it by parts to prevent out-of-memory errors
    sz = size(FE,1) * size(CE,1) * 8; % result in bytes
    desired_mem = 2*(1024^3 ); % 
    % at least 9 variables of FExCE size are needed in the main function
    n_chunks = ceil(10*sz / desired_mem); 

    len_chnk = ceil(size(CE,1) / n_chunks);
    intpts = [];
    FE2CE = sparse(0,0);
    FE2pts = sparse(0,0);
    CE2pts = sparse(0,0);
    for c = 1:n_chunks
        eidors_msg('@@@: chunk %d of %d',c,n_chunks,2);
        start = 1 + (c-1)*len_chnk;
        stop  = min(1 + c*len_chnk, size(CE,1));
        rng   = start:stop;
        [ip, te2ve, te2ip, ve2ip] = edge2edge_intersections(FE,FN,CE(rng,:),CN, epsilon);
        len    = size(intpts,1);
        intpts = [intpts; ip];
        idx = te2ve>0;
        te2ve(idx) = len + te2ve(idx);
        FE2CE  = [FE2CE te2ve];
        FE2pts = [FE2pts te2ip];
        if size(ve2ip,1) < size(CE2pts,1)
            ve2ip(size(CE2pts,1),1) = 0;
        end
        CE2pts = [CE2pts ve2ip];
    end

    
%-------------------------------------------------------------------------%
% Calculate intersection points between two sets of edges
function [intpts, FE2CE, FE2pts, CE2pts] = edge2edge_intersections_serial(FE,FN,CE,CN, epsilon)
% Based on the simple Matlab conversion by Cristian Dima of the C code 
% posted by Paul Bourke at
% http://paulbourke.net/geometry/pointlineplane/lineline.c
% http://paulbourke.net/geometry/pointlineplane/linelineintersect.m

    progress_msg(0);

    P1 = FN(FE(:,1),:);
    P2 = FN(FE(:,2),:);
    P3 = CN(CE(:,1),:);
    P4 = CN(CE(:,2),:);
    
    % these are too big to keep around
%     p13.x = bsxfun(@minus, P1(:,1), P3(:,1)');
%     p13.y = bsxfun(@minus, P1(:,2), P3(:,2)');
%     p13.z = bsxfun(@minus, P1(:,3), P3(:,3)');
    
    p43 = P4 - P3;
    p21 = P2 - P1;
    
    d4343 = sum(p43.^2,2);
    d2121 = sum(p21.^2,2);
    
    % get rid of degenerate cases - slow, and superfluous
%     JNK = bsxfun(@or,sparse(all(abs(p21)<eps,2)),sparse(all(abs(p43)<eps,2)'));
    IN  = false(size(P1,1),size(P3,1));

    intpts = [];
    
    todo = size(CE,1);
    id = 0;
    step = ceil(todo/100);
    
    pa = zeros(size(P1));
    
    for  v = 1:size(CE,1)
        id = id+1;
        if mod(id,step)==0, progress_msg(id/todo); end
        d1343 = 0; d4321 = 0;  d1321 = 0;
        for i = 1:3
            p13 = P1(:,i) - P3(v,i);
            d1343 = d1343 + p13*p43(v,i);
            d4321 = d4321 + p21(:,i)*p43(v,i);
            d1321 = d1321 + p13 .* p21(:,i);
        end
        
        
        
        denom = d2121 * d4343(v) - d4321.^2;
        % this is slow, and should be taken care of by mua < 1
%         JNK(:,v) = JNK(:,v) | abs(denom)<eps;
        numer = d1343 .* d4321 - d1321 * d4343(v);
        
        mua = numer./denom;
        mub = (d1343 + d4321.*mua) / d4343(v);


        D = 0; % distance
        for i = 1:3
            pa(:,i) = P1(:,i) + p21(:,i).*mua;
            pb = P3(v,i) + p43(v,i) * mub;
            D = D + (pa(:,i) - pb).^2;
        end
%         D = sqrt(D);
%         D(JNK(:,v)) = Inf;
        IN(:,v) = ((mua>0) + (mua<1) + (mub>0) + (mub<1) + (D<=epsilon)) == 5;
        if any(IN(:,v))
            intpts = [intpts; pa(IN(:,v),:)];
        end
    end
    [T, V] = find(IN);
    I = (1:size(intpts,1))';
    FE2CE = sparse(size(P1,1),size(P3,1));
    FE2CE(sub2ind(size(FE2CE),T,V)) = I;
    FE2pts = sparse(T,I,ones(size(I)),size(P1,1),size(I,1));
    CE2pts  = sparse(V,I,ones(size(I)),size(P3,1),size(I,1));
    
    if eidors_msg('log_level')>1
        progress_msg(Inf);
    end
    
    
%-------------------------------------------------------------------------%
% Calculate intersection points between two sets of edges
function [intpts, FE2CE, FE2pts, CE2pts] = edge2edge_intersections(FE,FN,CE,CN, epsilon)
% Based on the simple Matlab conversion by Cristian Dima of the C code 
% posted by Paul Bourke at
% http://paulbourke.net/geometry/pointlineplane/lineline.c
% http://paulbourke.net/geometry/pointlineplane/linelineintersect.m



    P1 = FN(FE(:,1),:);
    P2 = FN(FE(:,2),:);
    P3 = CN(CE(:,1),:);
    P4 = CN(CE(:,2),:);
    
    % these are too big to keep around
%     p13.x = bsxfun(@minus, P1(:,1), P3(:,1)');
%     p13.y = bsxfun(@minus, P1(:,2), P3(:,2)');
%     p13.z = bsxfun(@minus, P1(:,3), P3(:,3)');
    
    p43 = P4 - P3;
    p21 = P2 - P1;
    
    % get rid of degenerate cases
    JNK = bsxfun(@or,sparse(all(abs(p21)<eps,2)),sparse(all(abs(p43)<eps,2)'));
    
    f = {'x','y','z'};
    d1343 = 0; d4321 = 0;  d1321 = 0;

    for i = 1:3;
        p13 = bsxfun(@minus, P1(:,i), P3(:,i)');
        d1343 = d1343 + bsxfun(@times,p13,p43(:,i)');
        d4321 = d4321 + bsxfun(@times,p21(:,i),p43(:,i)');
        d1321 = d1321 + bsxfun(@times,p13,p21(:,i));
    end

    clear p13
    d4343 = sum(p43.^2,2);
    d2121 = sum(p21.^2,2);
    
    denom = bsxfun(@times,d2121,d4343') - d4321.^2;
    JNK = JNK | abs(denom)<eps;
    numer = d1343 .* d4321 - bsxfun(@times,d1321,d4343');
    
    mua = numer./denom;
    mub = bsxfun(@rdivide, d1343 + d4321.*mua , d4343');
    
    clear d1343 d4321 d1321 d4343 d2121
    
    pa = struct; pb = struct;
    D = 0; % distance
    for i = 1:3
        pa.(f{i}) = bsxfun(@plus,P1(:,i),bsxfun(@times,p21(:,i),mua));
        pb = bsxfun(@plus,P3(:,i)',bsxfun(@times,p43(:,i)',mub));
        D = D + (pa.(f{i}) - pb).^2;
    end
%     D = sqrt(D);
    D(JNK) = Inf;
    

    
    IN = mua>0 & mua <1 & mub>0 & mub<1 & D<=epsilon;
    nPts = nnz(IN);
    intpts = zeros(nPts,3);
    for i = 1:3
        intpts(:,i) = pa.(f{i})(IN);
    end
 
    [T, V] = find(IN);
    I = (1:nPts)';
    FE2CE = sparse(size(P1,1),size(P3,1));
    FE2CE(IN) = I;
   
    FE2pts = sparse(T,I,ones(size(I)),size(P1,1),size(I,1));
    CE2pts  = sparse(V,I,ones(size(I)),size(P3,1),size(I,1));
    