% Simulate obj $Id$

% GREIT v1
show_fem( inv_solve( i_gr, vh, vn) );
print_convert 'GREIT_examples03a.png';

% Sheffield Backprojection
show_fem( inv_solve( i_bp, vh, vn) );
print_convert 'GREIT_examples03b.png';

% 2D Gauss Newton Inverse
show_fem( inv_solve( i_gn, vt, vn) );
print_convert 'GREIT_examples03c.png';
