function H = s2r_tian_hierarchy(out, scales)
%S2R_TIAN_HIERARCHY  Nesting of the Tian scales, read from the voxels (no name rules).
%   H = s2r_tian_hierarchy(out, scales)      scales = {'S1','S2','S3','S4'}
%
% For every ROI of the finest scale, its parent at each coarser scale is the ROI
% that contains its voxels. Voxels are matched across scales through
% out.Tian_Sx.vox.row (row in the common subcortical union), so the scales do
% not need to hold the same voxels in the same order.
%
% OUTPUT  H, one row per ROI of the finest scale (n rows, nS = numel(scales) columns)
%   H.scales   {1 x nS}
%   H.id       [n x nS]  numeric ROI label at each scale
%   H.name     {n x nS}  ROI name at each scale ('HIP-rh', 'aHIP-rh', ...)
%   H.base     {n x nS}  the same without hemisphere ('HIP', 'aHIP', ...)
%   H.key      {n x 1}   base name at the coarsest scale: the region family ('HIP', 'pTHA', ...)
%   H.hemi     {n x 1}   'L' / 'R'
%   H.purity   [n x nS]  share of the voxels of the fine ROI that lie in the parent (1 = strictly nested)
%   H.strict   true if every purity is 1
%   H.order    [n x 1]   row order that keeps children next to their parents (sort by id of S1, S2, ...)
%   H.keys     {1 x nKeys} region families in atlas order
nS = numel(scales);
T = cell(1, nS);
for s = 1:nS
    T{s} = out.(['Tian_' scales{s}]);
end
fine  = T{nS};
rowsF = fine.vox.row(:);
labF  = double(fine.vox.roi_lab(:));
n     = numel(fine.roi_list);

H = struct();
H.scales = scales(:)';
H.id     = nan(n, nS);
H.name   = repmat({''}, n, nS);
H.purity = nan(n, nS);
for s = 1:nS
    [tf, loc] = ismember(rowsF, T{s}.vox.row(:));
    names_s = cellstr(T{s}.roi_names);
    list_s  = double(T{s}.roi_list(:));
    lab_s   = double(T{s}.vox.roi_lab(:));
    for q = 1:n
        m = (labF == double(fine.roi_list(q))) & tf;
        if ~any(m), continue; end
        labs = lab_s(loc(m));
        par  = mode(labs);
        H.id(q, s)     = par;
        H.purity(q, s) = mean(labs == par);
        H.name{q, s}   = names_s{find(list_s == par, 1)};
    end
end
[H.base, hemiAll] = s2r_split_hemi(H.name);
H.key    = H.base(:, 1);
H.hemi   = hemiAll(:, 1);
H.strict = all(H.purity(:) == 1);
[~, H.order] = sortrows(H.id);
[~, firstRow] = unique(H.id(:, 1), 'first');          % one row per coarsest-scale ROI, in label order
H.keys = unique(H.key(firstRow), 'stable')';
H.keys = H.keys(:)';
end
