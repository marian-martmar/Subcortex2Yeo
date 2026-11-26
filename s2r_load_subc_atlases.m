function A = s2r_load_subc_atlases(P)
% Modified 2/11/2026 to match s2r_ijk2mm and any reordering in _0_atlas_indexing.m 
% Read the subcortical atlases once and return them in one list.
%   A = s2r_load_subc_atlases(P)
%   P            paths (subc2rsn_paths)
%
% A(k), one element per atlas, in the order AAL, DBS80, DK84, Tian S1..S4:
%   .key          field name used in `out` ('AAL', 'DBS80', 'DK84', 'Tian_S1', ...)
%   .name         name for titles ('Tian S1')
%   .mm2_subc     3-D label volume on the 2-mm grid, 0 outside the subcortical ROIs
%   .tform        voxel -> mm transform of the 2-mm grid (identical for every atlas, TESTS T0.12);
%                 use it with s2r_ijk2mm, which converts MATLAB 1-based subscripts
%   .coords       always 'ijk' (kept for compatibility with older versions of s2r_roi_tables)
%   .label2name   containers.Map: numeric label -> ROI name
%   .keep_labels  labels to keep ([] = all non-zero)
%
% Same definitions as step 0: AAL relabelled 1..116 (subcortex = labels 37 38 41 42 71:78),
% DBS80 positions 32:49, DK84 FreeSurfer codes below 1000 (cerebellum 8 and 47 removed),
% DK84 1 mm -> 2 mm by majority vote. Cerebellum excluded, as EXCLUDE_CEREB = true in step 0.
% Steps 3, 4 and 5 of the original pipeline each repeated these statements; here they exist once.

% --- AAL (relabelled 1..116, as in step 0) ---
addpath(genpath(P.atlas_aal))
mm2_aal   = niftiread('AAL.nii');
info_aal  = niftiinfo('AAL.nii');
tform_aal = info_aal.Transform;
load('ROI_MNI_V4_List.mat');                  % ROI: ID, Nom_C, Nom_L, order L1, R1, ..., LN, RN
labelsAAL_orig = vertcat(ROI.ID);
Names_aal = {ROI.Nom_L}';
clear ROI;
% Relabel original IDs -> 1..116
[isListed, newLab] = ismember(mm2_aal, labelsAAL_orig);
assert(all(isListed(mm2_aal ~= 0)), 'AAL.nii has labels that are not in ROI_MNI_V4_List.mat');
mm2_aal = cast(newLab, class(mm2_aal));
mm2_aal(mm2_aal > 90) = 0;                    % no cerebellum / vermis
subc_aal_labels = [37 38 41 42 71 72 73 74 75 76 77 78]';   % after relabelling, labels = positions
mm2_subc_aal = mm2_aal;
mm2_subc_aal(~ismember(mm2_subc_aal, subc_aal_labels)) = 0;
aal_label2name = containers.Map(num2cell(1:numel(Names_aal)), Names_aal');

% --- DBS80 ---
addpath(genpath(P.atlas_dbs80))
mm2_dbs80   = niftiread('dbs80symm_2mm.nii');
info_dbs80  = niftiinfo('dbs80symm_2mm.nii');
tform_dbs80 = info_dbs80.Transform;

regions_dbs80 = readtable('dbs80_nodes.xlsx');  % expects id + some name column
labelsDBS80 = double(regions_dbs80.id);
% Manually from table: 32-49
% (32) L THA, L NA,  L CAU, L PUT, L gpe, L gpi, L stn, L AMY, L HIP
% (41) R HIP, R AMY, R stn,  R gpi, R gpe, R PUT, R CAU, R NA, R THA (49)
subc_within_labels_dbs80 = 32:49;
subc_dbs80_labels = labelsDBS80(subc_within_labels_dbs80);

mm2_subc_dbs80 = mm2_dbs80;
mm2_subc_dbs80(~ismember(mm2_subc_dbs80, subc_dbs80_labels)) = 0;

dbs80_label2name = make_label2name_from_table(regions_dbs80, 'id');

% --- DK84 1mm->2mm + LUT ---
addpath(genpath(P.atlas_dk84))
mm1_dk84 = niftiread('IIT_GM_Desikan_atlas.nii.gz');

[sx,sy,sz] = size(mm1_dk84);
assert(mod(sx,2)==0 && mod(sy,2)==0 && mod(sz,2)==0, 'DK84 dims must be divisible by 2.');
A = reshape(mm1_dk84, 2, sx/2, 2, sy/2, 2, sz/2);
A = permute(A, [1 3 5 2 4 6]);
A = reshape(A, 8, sx/2, sy/2, sz/2);
mm2_dk84 = squeeze(mode(A,1));
mm2_dk84 = cast(mm2_dk84, class(mm1_dk84));
clear A mm1_dk84;

fname = 'LUT_GM_Desikan_0to1.txt';
opts = detectImportOptions(fname, 'FileType','text');
opts.Delimiter = '\t';
opts.VariableNames = {'Label','R','G','B','Name'};
opts.VariableTypes = {'double','double','double','double','string'};
opts.ConsecutiveDelimitersRule = 'join';
opts.ExtraColumnsRule = 'ignore';
regions_dk84 = readtable(fname, opts);

labDK84_all = double(regions_dk84.Label);
assert(isequal(sort(labDK84_all), double(unique(mm2_dk84(mm2_dk84 ~= 0)))), ...
    'LUT and DK84 volume do not contain the same labels');
mm2_dk84(ismember(mm2_dk84, [8 47])) = 0;     % FreeSurfer Left/Right-Cerebellum-Cortex
subc_dk84_labels = labDK84_all(labDK84_all < 1000 & ~ismember(labDK84_all, [8 47]));   % below 1000 = subcortical

mm2_subc_dk84 = mm2_dk84;
mm2_subc_dk84(~ismember(mm2_subc_dk84, subc_dk84_labels)) = 0;

dk84_label2name = make_label2name_from_table(regions_dk84, 'Label');

% --- Tian S1..S4 ---
addpath(genpath(P.atlas_tian))
scales = {'S1','S2','S3','S4'};
Tian = struct();
for s = 1:numel(scales)
    S = scales{s};
    niiFile = sprintf('Tian_Subcortex_%s_3T.nii', S);
    txtFile = sprintf('Tian_Subcortex_%s_3T_label.txt', S);

    Tian.(S).mm2   = niftiread(niiFile);
    Tian.(S).info  = niftiinfo(niiFile);
    Tian.(S).tform = Tian.(S).info.Transform;

    lines = strtrim(readlines(txtFile));
    lines(lines=="") = [];
    Tian.(S).regions = cellstr(lines);

    Tian.(S).label2name = containers.Map('KeyType','double','ValueType','char');
    for k = 1:numel(Tian.(S).regions)
        Tian.(S).label2name(double(k)) = char(string(Tian.(S).regions{k}));
    end
end

%% Pack
A = struct('key', {}, 'name', {}, 'mm2_subc', {}, 'tform', {}, 'coords', {}, 'label2name', {}, 'keep_labels', {});
A(end+1) = struct('key', 'AAL',   'name', 'AAL',   'mm2_subc', mm2_subc_aal,   'tform', tform_aal,   'coords', 'ijk', ...
                  'label2name', aal_label2name,   'keep_labels', subc_aal_labels);
A(end+1) = struct('key', 'DBS80', 'name', 'DBS80', 'mm2_subc', mm2_subc_dbs80, 'tform', tform_dbs80, 'coords', 'ijk', ...
                  'label2name', dbs80_label2name, 'keep_labels', subc_dbs80_labels);
A(end+1) = struct('key', 'DK84', 'name', 'DK84', 'mm2_subc', mm2_subc_dk84, 'tform', tform_aal, 'coords', 'ijk', ...
                  'label2name', dk84_label2name, 'keep_labels', subc_dk84_labels);   % same 2-mm grid as AAL (T0.12)
for s = 1:numel(scales)
    S = scales{s};
    A(end+1) = struct('key', ['Tian_' S], 'name', ['Tian ' S], 'mm2_subc', Tian.(S).mm2, 'tform', Tian.(S).tform, 'coords', 'ijk', ...
                      'label2name', Tian.(S).label2name, 'keep_labels', []); %#ok<AGROW>
end
end

function m = make_label2name_from_table(T, labelColName)
    m = containers.Map('KeyType','double','ValueType','char');
    vn = T.Properties.VariableNames;

    if ~any(strcmpi(vn, labelColName))
        return;
    end

    nameCandidates = {'Name','name','Region','region','LabelName','ROIName','structure','Structure'};
    nameCol = '';
    for k = 1:numel(nameCandidates)
        j = find(strcmpi(vn, nameCandidates{k}), 1);
        if ~isempty(j), nameCol = vn{j}; break; end
    end
    if isempty(nameCol)
        for j = numel(vn):-1:1
            if iscellstr(T.(vn{j})) || isstring(T.(vn{j}))
                nameCol = vn{j}; break;
            end
        end
    end
    if isempty(nameCol), return; end

    labs = double(T.(labelColName));
    names = string(T.(nameCol));
    for i = 1:numel(labs)
        if ~isnan(labs(i))
            m(labs(i)) = char(names(i));
        end
    end
end
