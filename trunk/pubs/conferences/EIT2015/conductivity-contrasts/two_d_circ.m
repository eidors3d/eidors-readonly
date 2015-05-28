subplot(211)
ds = logspace(-4,4,101);
s = (1-ds)./(1+ds);
semilogx(ds, abs(s),'LineWidth',2);
ylabel('| signal |')
ylim([0,1.1]);
xlabel('<- non-conductive    -    homogeneous    -    conductive ->')
print -dpdf two_d_circ_sens.pdf
