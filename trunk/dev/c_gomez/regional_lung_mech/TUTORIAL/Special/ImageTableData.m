function z = ImageTableData(ztab,imgdim,index)
% function z = ImageTableData(ztab,imgdim,index)
% Copyright C. Gomez-Laberge, November 2010

[pixels,time] = size(ztab);
z = zeros([imgdim(1),imgdim(2),time]);

maskindex = find(index(:,3) == 1);

for t = 1:time
    ci = ztab(:,t);
    % Convert into image
    for p = 1:pixels
        coords = index(maskindex(p),:);
        z(coords(1),coords(2),t) = ci(p);
    end
end
return
