% Netgen simulation $Id$

subplot(211)
plot([vh.meas, vi.meas, vi.meas-vh.meas]);
print_convert netgen_sims05a.png '-density 100'
