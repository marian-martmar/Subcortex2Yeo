%% SUBC2RSN | STEP 3 - Group FC, group network label of every subcortical voxel, ROI tables, Tian hierarchy
% -------------------------------------------------------------------------
% PURPOSE
%   Put cluster results into a SINGLE RESULTS FILE to be read by all output steps (4 to 8) with:
%     avgZ  = tanh(cumZ / NSUB)   Fisher average across subjects (an r value, the main estimator)
%     avgFC = cumFC / NSUB        Arithmetic average of r (kept for comparison)
%     as7/17_Z/FC                 Argmax over networks: group label of every voxel
%     out.<ATLAS>                 structures out.AAL, out.DBS80, out.DK84, out.Tian_S1..S4, see s2r_roi_tables.m
%     hier                        Nesting of the Tian scales S1 > S2 > S3 > S4, read from the voxels, see s2r_tian_hierarchy
%
% INPUTS
%   Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat -> idcs, subset_indices, subjectIDs
%   FC7_group.mat -> cumFC7, cumZ7
%   FC17_group.mat -> cumFC17, cumZ17
%   atlas files (read once, by functions/s2r_load_subc_atlases.m)
%
% OUTPUT  C.file.results = <P.results_dir>/subc2rsn_results.mat
%   labels7, labels17, subjectIDs, NSUB, info
%   out.<ATLAS>         Each atlas struct with 
                        % .roi_list: numeric labels
                        % .roi_names: name label
                        % .tables -> nested struct RNS7/17 -> AvgZ/FC -> count/pct (ROI X RSN tables)
                        % .vox -> structure of linear arrays with voxel's in idcs(subset_indices.<ATLAS>)
                            % .vox.row: corresponding row in the full subcortical voxel list 
                            % .vox.X,Y,Z: coordinates real-space
                            % .roi_lab: ROI's numeric label of the atlas, per voxel
                            % .as7/17_FC/Z: assigned RSN with numeric label
%   hier                
%   subset_idcs_subc    [Nsubc x 1] rows of all FC / label arrays (positions within idcs)
%   idcs_subc      [Nsubc x 1] idcs(subset_idcs_subs) = the same rows as linear indices of the 2-mm volume
%   avgZ/FC7/17         [Nsubc x 7|17] group FC
%   as7/17_Z/FC         [Nsubc x 1]    group labels

% FUNCTIONS Section 3 calls s2r_load_subc_atlases and s2r_roi_tables, which hold the atlas statements
%           and analyze_one_subc_atlas of the original step 3 (without its figures)
%            
% Author: Marian Martinez-Marin. Proposed final version (v2), 16/11/2026.
% -------------------------------------------------------------------------

clear; clc;

C = subc2rsn_config();
P = C.P;
addpath(P.results_cluster)                
addpath(genpath(P.fc_dir))

%% 3.0. Load union indexing + group FC
VariablesFrom0 = fullfile(P.outputs_0,'Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat');
load(VariablesFrom0, 'idcs','subset_indices','subjectIDs');
NSUB = numel(subjectIDs);

% From step 2 (cluster), being nValid7 = NSUBS with valid FC7:
% FC7_group  = tanh(cumZ7  ./ nValid7);  -> Back to r units after Fisher
% FC7_group_meanr  = cumFC7  ./ nValid7; -> direct averaging of r units
ClusterResults7=fullfile(P.results_cluster, 'FC7_group.mat');
ClusterResults17=fullfile(P.results_cluster, 'FC17_group.mat');
load(ClusterResults7,'cumZ7','FC7_group','FC7_group_meanr','nValid7','info','defs');   % [Nsubc x 7]
load(ClusterResults17,'cumZ17',  'FC17_group',  'FC17_group_meanr','nValid17');        % [Nsubc x 7]

labels7  = C.rsn7.labels;
labels17 = C.rsn17.labels;

%% 3.1. Initial FC row definition from step 0 must match cluster job
subset_idcs_subc = unique([
    subset_indices.subc_dbs80(:);
    subset_indices.subc_aal(:);
    subset_indices.subc_dk84(:);
    subset_indices.subc_tian.S1(:);
    subset_indices.subc_tian.S2(:);
    subset_indices.subc_tian.S3(:);
    subset_indices.subc_tian.S4(:)
]);

assert(sum(subset_idcs_subc(:)~=info.rows.subset_idcs_subc(:))==0, ...
    'Defined subset_indices_subc (based on step 0)\ndoes not match the one from cluster variables')
idcs_subc = idcs(subset_idcs_subc);  % linear indices in 2mm reference grid
assert(sum(idcs_subc(:)~=info.rows.idcs_subc(:))==0, ...
    'Defined subset_indices_subc (based on step 0)\ndoes not match the one from cluster variables')

Nsubc = numel(subset_idcs_subc);
fprintf('Nsubc (union rows) = %d\n', Nsubc);

assert(size(FC7_group,1)   == Nsubc,  'FC7_group rows != subset_idcs_subc');
assert(size(FC17_group,1)  == Nsubc,  'FC17_group rows != subset_idcs_subc');
assert(size(FC7_group_meanr,1)  == Nsubc,  'FC7_group_meanr rows != subset_idcs_subc');
assert(size(FC17_group_meanr,1) == Nsubc,  'FC17_group_meanr rows != subset_idcs_subc');


%% 3.2. Voxel-wise RSN assignment (group)
% load('Assigned_7.mat', 'assigned7_group') or recompute:
[~, as7_Z]    = max(FC7_group,   [], 2);      
as7_Z = uint8(as7_Z); 
% load('Assigned_17.mat', 'assigned17_group') or recompute:
[~, as17_Z]   = max(FC17_group,  [], 2);
as17_Z = uint8(as17_Z);
% Compute assignments also for direct average of FC matrices (without Fisher transform)
[~, as7_FC]   = max(FC7_group_meanr,  [], 2); 
as7_FC = uint8(as7_FC); 
[~, as17_FC]  = max(FC17_group_meanr, [], 2);
as17_FC = uint8(as17_FC);


%% 3.3. Every subcortical atlas: voxels, labels, ROI x RSN tables
A = s2r_load_subc_atlases(P);

out = struct();
for a = 1:numel(A)
    out.(A(a).key) = s2r_roi_tables(A(a), idcs_subc, as7_FC,as7_Z,labels7, as17_FC,as17_Z,labels17);
    fprintf('%-8s %3d ROIs, %5d voxels\n', A(a).key, numel(out.(A(a).key).roi_list), numel(out.(A(a).key).vox.row));
end

%% 3.4. Tian hierarchy from the voxels
hier = s2r_tian_hierarchy(out, C.tian_scales);
fprintf('Tian scales strictly nested: %d | region families: %s\n', hier.strict, strjoin(hier.keys, ', '));

%% 3.5. Save the single results file
info = struct('created', datestr(now), 'estimator_for_outputs', C.estimator, ...
              'n_nan_rows_FC7_group', sum(any(isnan(FC7_group), 2)), 'n_nan_rows_avgZ17', sum(any(isnan(FC17_group), 2)));
if info.n_nan_rows_FC7_group > 0 || info.n_nan_rows_avgZ17 > 0
    warning('%d (Yeo-7) and %d (Yeo-17) voxels have NaN group FC and received label 1 from max(). See TESTS T22.1.', ...
        info.n_nan_rows_FC7_group, info.n_nan_rows_avgZ17);
end
if ~exist(P.results_dir, 'dir'); mkdir(P.results_dir); end
save(C.file.results, 'out', 'hier', 'subset_idcs_subc', 'idcs_subc', ...
     'FC7_group', 'FC7_group_meanr', 'FC17_group', 'FC17_group_meanr', 'as7_Z', 'as7_FC', 'as17_Z', 'as17_FC', ...
     'labels7', 'labels17', 'subjectIDs', 'NSUB', 'info', '-v7.3');
fprintf('Saved %s\n', C.file.results);
