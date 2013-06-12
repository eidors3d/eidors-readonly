function out = eidors_template_fcn(in)
%TMPL   EIDORS function template
%  OUT=TMPL(IN) calculates this and that.
%  OUT=TMPL(IN1, IN2) calculates this and that differently.
% 
%  Inputs:
%    IN  - a vector
%
%  Outputs:
%    OUT - a single logical value
% 
%  Examples:
%    eidors_template(rand(1))
%
% See also SOME_OTHER_FUNCTION

% (C) 2012 Author Name. License: GPL version 2 or version 3
% $Id$

% if input is 'UNIT_TEST', run tests
if ischar(in) && strcmp(in,'UNIT_TEST'), do_unit_test; return; end

% uncomment this block to use caching (for long calculations)
% cache_obj = get_cache_obj(in);
% out = eidors_obj('get-cache', cache_obj, 'name_of_stored_value');
% if ~isempty(out)
%    eidors_msg('@@@ Using cached value',3);
%    return
% end

% do calculations in separate functions to ease reading
out = do_calculations(in);

% uncomment this block to use caching
% eidors_obj('set-cache', cache_obj, 'name_of_stored_value', out);

%-------------------------------------------------------------------------%
% The main function
function out = do_calculations(in)
out = logical(in(1) - 0.5);

%-------------------------------------------------------------------------%
% Assemble a reference object for caching
function cache_obj = get_cache_obj(in)
% usually, the computed value will not depend on all inputs (or not all
% their fields). Use this function to select only the relevant ones.
cache_obj = in;
% cache_obj = {in1, in3};
% cache_obj = {in.elems in.nodes};

%-------------------------------------------------------------------------%
% Perfom unit tests
function do_unit_test
% Use this function to provide some code to test/showcase the functionality
% We are not very rigorous about the tests, but of course the more cases
% you test, the more useful this function is. Some ideas are below:

% You could just calculate and display some value / figure
display(eidors_template_fcn(rand(1,2)));

% You could display a more informative message
if eidors_template_fcn(rand(1))
   msg = ('OK');
else
   msg = ('ERROR');
end
display(['Test1: ' msg])

% You could assert that the value is correct
assert(eidors_template_fcn(2) == 1);

% You could have do_unit_test return a boolean indicating if all tests pass.
