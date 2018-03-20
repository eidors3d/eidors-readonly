yposns = [45  20 50]; xposns = [50  40 27]; ofs= [0,22,15];

% Show positions on image
hold on; for i = 1:length(xposns)
    plot(xposns(i),yposns(i),'s','LineWidth',10);
end; hold off;
