% Plot residuals $Id$
subplot(211);
r = rec_img.inv_solve_gn.r;
k = size(r,1)-1;

x = 1:(k+1); % k+1 => look at the solve after last iteration
y = r(x, :);
y = y ./ repmat(max(y,[],1),size(y,1),1) * 100;
plot(x-1, y, 'o-', 'linewidth', 2, 'MarkerSize', 10);
title('residuals');
axis tight;
ylabel('residual (% of max)');
xlabel('iteration');
set(gca, 'xtick', x);
set(gca, 'xlim', [0 max(x)-1]);
legend('residual','meas. misfit','prior misfit');
legend('Location', 'East');
print_convert basic_iterative03a.png
