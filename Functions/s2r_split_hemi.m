function [base, hemi] = s2r_split_hemi(names)
%S2R_SPLIT_HEMI  Split ROI names such as 'HIP-rh' into base name and hemisphere.
%   [base, hemi] = s2r_split_hemi(names)
%   names  char, cellstr or string array (any shape)
%   base   cellstr, same shape: name without the final -lh/-rh, _lh/_rh or ' lh'/' rh'
%   hemi   cellstr, same shape: 'L', 'R' or '' (no hemisphere token)
% Same rule as strip_hemi_token / parse_hemi of the original figure scripts.
names = cellstr(names);
names = strtrim(names);
base  = regexprep(names, '([\-_\s])(lh|rh)\s*$', '');
hemi  = repmat({''}, size(names));
low   = lower(names);
isL = ~cellfun('isempty', regexp(low, '[\-_\s]lh$', 'once'));
isR = ~cellfun('isempty', regexp(low, '[\-_\s]rh$', 'once'));
hemi(isL) = {'L'};
hemi(isR) = {'R'};
end
