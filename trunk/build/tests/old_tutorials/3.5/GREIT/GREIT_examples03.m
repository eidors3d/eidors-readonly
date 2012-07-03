% Simulate obj $Id$

% GREIT v1
show_fem( inv_solve( i_gr, vh, vn) ); axis equal;
print_convert 'GREIT_examples03a.png' '-density 50';

% Sheffield Backprojection
show_fem( inv_solve( i_bp, vh, vn) ); axis equal;
print_convert 'GREIT_examples03b.png' '-density 50';

% 2D Gauss Newton Inverse
show_fem( inv_solve( i_gn, vh, vn) ); axis equal;
print_convert 'GREIT_examples03c.png' '-density 50';
