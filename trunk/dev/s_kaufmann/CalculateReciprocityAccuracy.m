% Function: CalculateReciprocityAccuracy
%
% Calculates the reciprocity accuracy
%
% Input:    Pattern - structure
%                   .stim           - from mk_stim_patterns
%                   .meas_sel       - from mk_stim_patterns
%
%           Voltage / Current / Impedance Vector
%
% Note:     Requires started EIDORS
%
% Return:   relRA - Vector containing the reciprocity accuracy according to:
%           Evaluation of EIT System Performance Evaluation of EIT System Performance, Yasheng Maimaitijiang, Stephan Böhm, Pascal O. Gaggero, Andy Adler (2011)
%
%           AbsoluteReciprocityError - Vector containing the absolute
%           reciprocity error
%
% Written by Steffen Kaufmann <sk@steffen-kaufmann.com>
% February 2014

function [relRA, AbsoluteReciprocityError] = CalculateReciprocityAccuracy(Pattern, v)
    if (isnan(v))
        relRA = NaN;
        AbsoluteReciprocityError = NaN;
        return;
    end;

    % EIDORS
    fmdl.type = 'fwd_model';
    fmdl.stimulation = Pattern.stim;
    fmdl.meas_select = Pattern.meas_sel;
    fmdl.normalize_measurements = 0;

    idx = reciprocity_idx( fmdl ); % Get idx of reciprocity values

    AbsoluteReciprocityError = v - v(idx, :);
    relRA = 1-abs(AbsoluteReciprocityError ./ v);     % See: Evaluation of EIT System Performance Evaluation of EIT System Performance, Yasheng Maimaitijiang, Stephan Böhm, Pascal O. Gaggero, Andy Adler (2011)
end
