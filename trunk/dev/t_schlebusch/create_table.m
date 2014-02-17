load('ws-2014-02-09--big-sim-results.mat')

data = cell2mat(results);
trainGI = data((data(:,2)==2),:) % select contrast=2 as reference

[p_gi,ErrorEst] = polyfit(trainGI(:,6),trainGI(:,1),3);
gi_fit = polyval(p_gi,trainGI(:,6));

%%
volumina = [50:100:550]; % volume in ml
r_real = (volumina./((4/3*pi)*1000*1000)).^(1/3);
cond_real = [1:0.5:3];
noiselevels = [1e1 1e2 1e3 1e4];
erg = zeros(100,100);

for cond = 1:size(cond_real,2)
    dset1 = data(data(:,2)==cond_real(cond),:);
    for noise = 1:size(noiselevels,2)
        dset2 = dset1(dset1(:,3)==noiselevels(noise),:);
        for radius = 1:size(r_real,2)
            dset = dset2(dset2(:,1)==r_real(radius),:);
            spalte = (cond-1)*size(noiselevels,2)+noise;
            zeile = (radius-1)*2+1;
            
            gi = (polyval(p_gi,dset(:,6))-dset(:,1))./dset(:,1);
            pr = (dset(:,4)-dset(:,1))/dset(:,1);
            erg(zeile,spalte) = gi;
            erg(zeile+1,spalte) = pr;
        end
    end
end