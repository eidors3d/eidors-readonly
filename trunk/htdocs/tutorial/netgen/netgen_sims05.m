% Netgen simulation $Id$

subplot(211)
plot([vh.meas, vi.meas, vi.meas-vh.meas]);
print -dpng -r100 netgen_sims05a.png
