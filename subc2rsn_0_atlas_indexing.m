%% SUBC2RSN | STEP 0 - Index the voxels of every atlas in one common 2-mm MNI grid
% -------------------------------------------------------------------------
% PURPOSE
%   Convert atlases into both (a) a list of linear voxel indices in the 91x109x91 MNI152 2-mm array
%                              (b) one label per voxel
%   and builds the two variables the whole project relies on:
%                   - idcs: linear indices (into the 3-D array) of the UNION of all
%                           voxels belonging to any atlas (cortical or subcortical)
%                   - subset_indices (per atlas/RSN): POSITIONS WITHIN idcs of its voxels
%   and also maps every Schaefer-1000 parcel to its Yeo-7 and Yeo-17 network (from the parcel names in the centroid CSVs)
%   and stores the RSN voxel sets.
%
% INPUTS (files)
%   - SubjectIDs_100.mat                                               -> subjectIDs
%   - Schaefer2018_1000Parcels_{7,17}Networks_order_FSLMNI152_2mm.nii  + .Centroid_RAS.csv)
%   - AAL.nii, ROI_MNI_V4_List.mat                                     (AAL116 -> first 90 labels)
%   - IIT_GM_Desikan_atlas.nii.gz (1 mm), LUT_GM_Desikan_0to1.txt      (DK84)
%   - dbs80symm_2mm.nii, dbs80_nodes.xlsx                              (DBS80)
%   - Tian_Subcortex_S{1,2,3,4}_3T.nii + _label.txt                    (Tian)
%
% OUTPUTS (files, written in the current folder when SAVE_OUTPUTS = true)
%   - Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat : subjectIDs, idcs, subset_indices
%   - Variables_metainformation.mat                           : meta, labels7, labels17   (NOT read by any script of this version)
%
% KEY VARIABLES
%   mm2_<atlas>        3-D label volume (numerical atlas label; 0 = background) 
%                      WARNING: mm2 are NOT real-space mm but array elements (MNI 3-D array)
%   idcs_<set>         linear indices of a voxel set in the 3-D array
%   lab<Atlas>_<set>   numeric label of atlas regions used for voxels in mm2_<atlas> and idcs_<set>
%   i_,j_,k_ / X,Y,Z   array subscripts / world coordinates (mm) of the same voxels using s2r_ijk2mm, so (0,0,0) aligns with FSL MNI152 2mm convention (45, 63, 36)
%   subset_indices.cort_schaefer | rsn{1..7} | rsn17{1..17} | cort_aal | cort_dk84 |  cort_dbs80 |
%                  subc_tian.S1..S4                         | subc_aal | subc_dk84 |  subc_dbs80 |
%
% SUBCORTICAL DEFINITIONS (positions within label lists, set by hand)
%   Tian  : the whole atlas
%   AAL   : positions [37 38 41 42 71:78] -> Hipp, Amyg, Caud, Put, Pall, Thal (L/R)
%   DK84  : LUT positions [2:8 10:16] (cerebellum = positions [1 9], removed)
%   DBS80 : positions 32:49

% NOTES
%   1. Variables i, j, k, T, info, mm, idcs are re-used along the way

% Author: Marian Martinez-Marin. v1 December 2025. v2 November 2026
% -------------------------------------------------------------------------
clear all; clc;

SAVE_OUTPUTS = false;   % save: meta.mat & Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat

%% 0.0. Load Data and functions + Parameters
%  Config.
C = subc2rsn_config();
EXCLUDE_CEREB = C.EXCLUDE_CEREB %true;
P = subc2rsn_paths();
folder_out = P.outputs_0; %'C:\Users\maria\OneDrive\Escritorio\CNG\00. THESIS\1.Subc2rsn\outputs_0'

%  Functions
addpath(P.repo_functions) % For ijk2mm (3d array indices to xyz coordinates)
addpath(genpath(P.spm12)) % For loading NIFTI (.nii)
addpath(P.functions_tian) % For smoothing time series
addpath(P.functions_wishart) % For denoising

%  Data MNI loading
addpath(P.hcp_local)
load("SubjectIDs_100.mat")
NSUB = size(subjectIDs,1);

%% 0.1. Data Schaefer 1000 <--> Yeo 7
% Volume Array for Schaefer1000 in files of 7Netwoks
addpath(genpath(P.atlas_rsn)) 
mm2_schaefer = niftiread('Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii'); % Yeo Github for schaefer
Nschaefer = length(unique(mm2_schaefer(mm2_schaefer~=0))); % sum(isnan(mm2_schaefer(:))) == 0
% Mask within array
mask_cort = mm2_schaefer ~= 0;                                              % Non-interesting voxels are 0
idcs_cort = find(mask_cort);                                                % flattened indices
[i_cort, j_cort, k_cort] = ind2sub(size(mm2_schaefer), idcs_cort);
% Numeric atlas label per flattened indices
labSchaefer7 = mm2_schaefer(idcs_cort);  % Schaefer: 1..1000
% Schaefer CENTROIDS (CSV ordering: L1, L2, ... L500, R1, R2, ... R500) have a Yeo7 RSN assigned
SchaeferCOG = readtable('Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'); 

% MNI Coordinates (from NIFTI information)
info = niftiinfo('Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii');
T = info.Transform.T;  % voxel -> world (mm)
tform7 = info.Transform; % affine3d object transformation
pts = [i_cort(:), j_cort(:), k_cort(:)];
% [OLD] mm = transformPointsForward(tform7, pts); % Function on affine3d object
% [OLD] X = mm(:,1); Y = mm(:,2); Z = mm(:,3);
[X, Y, Z] = s2r_ijk2mm(tform7, i_cort, j_cort, k_cort); % [NEW]

% Yeo 7 networks
labels7 = C.rsn7.labels;                                                    % {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
Nrsn = length(labels7);
% MNI Dimensions (shared across the other MNI atlas files, both cortical and subcortical)
[Ni, Nj, Nk] = size(mm2_schaefer);                                          % 3d space
Nijk = Ni*Nj*Nk;

% [NEW] Align real-space origin (0,0,0)
tform_ref = tform7;                                    % the 2 mm grid shared by every atlas (T0.12)
[x0, y0, z0] = s2r_ijk2mm(tform_ref, 46, 64, 37);
assert(all(abs([x0 y0 z0]) < 1e-6), 'Reference is not the FSL MNI152 2 mm grid');

% Visual check NIFTI vs centroids: TESTS T0.1

%% 0.2. Data Schaefer 1000 <--> Yeo 17
% Volume Array for Schaefer1000 in files of 17Netwoks
mm2_schaefer_17 = niftiread('Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii');
Nschaefer17 = length(unique(mm2_schaefer_17(mm2_schaefer_17~=0))); % sum(isnan(mm2_schaefer_17(:))) == 0
% Mask within array
mask_cort17 = mm2_schaefer_17 ~= 0;              % Non-interesting voxels are 0
idcs_cort17 = find(mask_cort17);                % flattened indices
[i_cort17, j_cort17, k_cort17] = ind2sub(size(mm2_schaefer_17), idcs_cort17);
% Numeric atlas label per flattened indices
labSchaefer17 = mm2_schaefer_17(idcs_cort17);  % Schaefer: 1..1000
% Schaefer CENTROIDS (CSV ordering: L1, L2, ... L500, R1, R2, ... R500) have a Yeo17 RSN assigned
SchaeferCOG_17 = readtable('Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');

% MNI Coordinates (from NIFTI information)
info = niftiinfo('Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii');
T = info.Transform.T;  % voxel -> world (mm)
tform17 = info.Transform; % affine3d object transformation
pts = [i_cort17(:), j_cort17(:), k_cort17(:)];
% [OLD] mm = transformPointsForward(tform17, pts); % Function on affine3d object
% [OLD] X17 = mm(:,1); Y17 = mm(:,2); Z17 = mm(:,3);
[X17, Y17, Z17] = s2r_ijk2mm(tform17, i_cort17, j_cort17, k_cort17); % [NEW]

% Yeo 17 networks
labels17 = C.rsn17.labels;
Nrsn_17 = length(labels17);

% Visual check NIFTI vs centroids: TESTS T0.2

%% 0.3 AAL (116 includes cerebellum --> we can later subtract AAL 90)
addpath(genpath(P.atlas_aal))
mm2_aal = niftiread('AAL.nii');   % voxel values = original AAL IDs (2001, 2002, ..., 9170)
load('ROI_MNI_V4_List.mat');      % ROI: ID, Nom_C, Nom_L, in atlas order L1, R1, ..., LN, RN
labelsAAL_orig = vertcat(ROI.ID);
NamesBrief_aal = {ROI.Nom_C}';
Names_aal      = {ROI.Nom_L}'; % Order: L1, R1, ..., LN, RN
clear ROI;

% --- Relabel the volume: original ID -> position in the list (1..116) ---
[isListed, newLab] = ismember(mm2_aal, labelsAAL_orig);   % newLab = 0 for background
assert(all(isListed(mm2_aal ~= 0)), 'AAL.nii has labels that are not in ROI_MNI_V4_List.mat');
mm2_aal   = cast(newLab, class(mm2_aal));
labelsAAL = (1:numel(labelsAAL_orig))';

% Correspondence between both numberings (keep it with the outputs)
regions_aal = table(labelsAAL, labelsAAL_orig, NamesBrief_aal, Names_aal, ...
    'VariableNames', {'Label', 'Label_orig', 'Name_short', 'Name'});
aal_label2name_orig = containers.Map(double(labelsAAL_orig), string(Names_aal));
aal_label2name = containers.Map(double(labelsAAL), string(Names_aal));

if EXCLUDE_CEREB == true
    labelsAAL = labelsAAL(1:90);
    NamesBrief_aal = NamesBrief_aal(1:90);
    Names_aal = Names_aal(1:90);
    mm2_aal(~ismember(mm2_aal, labelsAAL)) = 0;
end
Naal = length(unique(mm2_aal(mm2_aal~=0)))

% Manual separation of subcortical/cortical numerical labels
subc_idx_within_labels_aal = [37 38 41 42 71 72 73 74 75 76 77 78]; % Accumbens in caudate/putamen + single thalamus
if EXCLUDE_CEREB == false
    subc_idx_within_labels_aal = [subc_idx_within_labels_aal 91:116];
end
disp(['Number of subcortical AAL regions: ', num2str(length(subc_idx_within_labels_aal))])
labelsAAL_subc = labelsAAL(subc_idx_within_labels_aal);
labelsAAL_cort = labelsAAL; % Labels are not exactly 1..90/117 so we get them from original list
labelsAAL_cort(ismember(labelsAAL_cort, labelsAAL_subc)) = []; % Get rid of manually identified as subcortical

% --- Mask subcorticals in (i,j,k) array
mask_subc_aal = ismember(mm2_aal, labelsAAL_subc);
mask_non_subc = ~ismember(mm2_aal, labelsAAL_subc);
% mm2 array with subc only
mm2_subc_aal = mm2_aal;
mm2_subc_aal(mask_non_subc) = 0; clear mask_non_subc;
% --- Mask corticals in (i,j,k) array
mask_cort_aal = (mm2_aal ~= 0) & ~mask_subc_aal; 
% mm2 array with cort only
mm2_cort_aal = mm2_aal;
mm2_cort_aal(mask_subc_aal) = 0; % mask_cort_aal could also be defined from here with ismember(mm2_cort_aal, labelsAAL_cort);

% --- 3D Array for AAL subcortex ---
idcs_subc_aal = find(mask_subc_aal);
[i_subc_aal, j_subc_aal, k_subc_aal] = ind2sub(size(mm2_subc_aal), idcs_subc_aal);
labAAL_subc=mm2_subc_aal(idcs_subc_aal);
% --- 3D Array for AAL cortex ---
idcs_cort_aal = find(mask_cort_aal);
[i_cort_aal, j_cort_aal, k_cort_aal] = ind2sub(size(mm2_cort_aal), idcs_cort_aal);
labAAL_cort=mm2_cort_aal(idcs_cort_aal);

% Real-space MNI Coordinates (from NIFTI information of AAL.nii)
info_aal = niftiinfo('AAL.nii');
tform_aal = info_aal.Transform;
% --- subcortical voxels -> mm
pts_s = [i_subc_aal(:), j_subc_aal(:), k_subc_aal(:)];
% [OLD] mm_s = transformPointsForward(tform_aal, pts_s);
% [OLD] XsAAL = mm_s(:,1); YsAAL = mm_s(:,2); ZsAAL = mm_s(:,3);
[XsAAL, YsAAL, ZsAAL] = s2r_ijk2mm(tform_aal, i_subc_aal, j_subc_aal, k_subc_aal); % [NEW]
% --- cortical voxels -> mm
pts_c = [i_cort_aal(:), j_cort_aal(:), k_cort_aal(:)];
% [OLD] mm_c = transformPointsForward(tform_aal, pts_c);
% [OLD] XcAAL = mm_c(:,1); YcAAL = mm_c(:,2); ZcAAL = mm_c(:,3);
[XcAAL, YcAAL, ZcAAL] = s2r_ijk2mm(tform_aal, i_cort_aal, j_cort_aal, k_cort_aal); % [NEW]

% (visual check AAL cortex / subcortex: TESTS T0.3)

%% 0.4 DK84 
addpath(genpath(P.atlas_dk84)) 
% nii.gz file from https://www.nitrc.org/frs/?group_id=432
mm1_dk84=niftiread('IIT_GM_Desikan_atlas.nii.gz');  % CAREFUL 1mm! We'll group 2x2x2 bins into 2mm^3 voxels
% Downsample 1mm -> 2mm(majority label every 2^3 block)
[sx,sy,sz] = size(mm1_dk84);
assert(mod(sx,2)==0 && mod(sy,2)==0 && mod(sz,2)==0, ... % Sizes are divisible by 2 (required for clean 2x downsample)
    'DK84 1mm volume size must be divisible by 2 in each dimension.');
A = reshape(mm1_dk84, 2,sx/2, 2,sy/2, 2,sz/2);
A = permute(A, [1 3 5 2 4 6]);  % [2 2 2 sx/2 sy/2 sz/2]
A = reshape(A, 8, sx/2,sy/2,sz/2);   % stack each 2x2x2 block
mm2_dk84 = squeeze(mode(A, 1)); % majority label in each block
clear A
mm2_dk84 = cast(mm2_dk84, class(mm1_dk84)); % Inherited from mm1 order: Lsubc1,...,LsubcNs,Rsubc1,...,RsubcNs,Lcort1,...,LcortNc,Rcort1,...,RcortNc (Ns=#subc., Nc=#cort.)
Ndk84 = length(unique(mm1_dk84(mm1_dk84~=0))); % sum(isnan(mm1_dk84(:))) == 0 
assert(length(unique(mm2_dk84(mm2_dk84~=0)))==Ndk84) % Check no numerical labels are lost after resampling

% Numeric labels dk84 (numbers in regions_dk84.Label)
fname = 'LUT_GM_Desikan_0to1.txt';
opts = detectImportOptions(fname, 'FileType','text');
opts.Delimiter = '\t';
opts.VariableNames = {'Label','R','G','B','Name'};
opts.VariableTypes = {'double','double','double','double','string'};
opts.ConsecutiveDelimitersRule = 'join';
opts.ExtraColumnsRule = 'ignore';
regions_dk84 = readtable(fname, opts);
labelsDK84=regions_dk84.Label;
assert(isequal(sort(labelsDK84), double(unique(mm1_dk84(mm1_dk84 ~= 0)))),'LUT and volume do not contain the same labels');
Names_dk84=regions_dk84.Name;

% Manual separation of subcortical/cortical numerical labels
IsCereb_within_labels_dk84 = ismember(labelsDK84, [8 47]);          % FreeSurfer Left/Right-Cerebellum-Cortex
cereb_idx_within_labels_dk84 = find(IsCereb_within_labels_dk84);
if EXCLUDE_CEREB == true
    labelsDK84 = labelsDK84(~IsCereb_within_labels_dk84);
    Names_dk84 = Names_dk84(~IsCereb_within_labels_dk84);           % same logical index: stays aligned
    mm2_dk84(~ismember(mm2_dk84, labelsDK84)) = 0;
end
Ndk84 = length(labelsDK84); % Update number of ROIs after removing cerebellum
assert(length(unique(mm2_dk84(mm2_dk84~=0)))==Ndk84)

IsSubc_within_labels_dk84 = labelsDK84 < 1000;  % below 1000 subcortical, cortex lh-1000s/ cortex rh-2000s
subc_idx_within_labels_dk84 = find(IsSubc_within_labels_dk84); %1:find(contains(Names_dk84, "ctx"), 1) - 1;   % everything listed before the first cortical name
disp(['Number of subcortical DK84 regions: ', num2str(length(subc_idx_within_labels_dk84))])
labelsDK84_subc = labelsDK84(subc_idx_within_labels_dk84);
labelsDK84_cort = labelsDK84; % Labels are not exactly 1..84/117 so we get them from original list
labelsDK84_cort(ismember(labelsDK84_cort, labelsDK84_subc)) = []; % Get rid of manually identified as subcortical

% --- Mask subcorticals in (i,j,k) array
mask_subc_dk84 = ismember(mm2_dk84, labelsDK84_subc); 
mask_non_subc = ~ismember(mm2_dk84, labelsDK84_subc);
% mm2 array with subc only
mm2_subc_dk84 = mm2_dk84;
mm2_subc_dk84(mask_non_subc) = 0; clear mask_non_subc;
% --- Mask corticals in (i,j,k) array
mask_cort_dk84 = (mm2_dk84 ~= 0) & ~mask_subc_dk84; 
% mm2 array with cort only
mm2_cort_dk84 = mm2_dk84;
mm2_cort_dk84(mask_subc_dk84) = 0; % mask_cort_dk84 could also be defined from here with ismember(mm2_cort_dk84, labelsDK84_cort);

% --- 3D Array for DK84 subcortex ---
idcs_subc_dk84 = find(mask_subc_dk84);
[i_subc_dk84, j_subc_dk84, k_subc_dk84] = ind2sub(size(mm2_dk84), idcs_subc_dk84);
labDK84_subc = double(mm2_dk84(idcs_subc_dk84));
% --- 3D Array for DK84 cortex ---
idcs_cort_dk84 = find(mask_cort_dk84);
[i_cort_dk84, j_cort_dk84, k_cort_dk84] = ind2sub(size(mm2_dk84), idcs_cort_dk84);
labDK84_cort = double(mm2_dk84(idcs_cort_dk84));

dk_label2name = containers.Map(double(regions_dk84.Label), string(regions_dk84.Name));

% --- ADDITIONAL: Map 2mm voxel indices -> 1mm intrinsic centers -> world (mm) ---
% [OLD] % MNI Coordinates from DK84 NIFTI information
% [OLD] info_dk84_1 = niftiinfo('IIT_GM_Desikan_atlas.nii.gz');
% [OLD] tform_dk84_1 = info_dk84_1.Transform;
% [OLD] % Each 2mm voxel (i2,j2,k2) corresponds to the center of its 2x2x2 1mm block
% [OLD] % --- subcortical
% [OLD] pts1_subc = [2*i_subc_dk84(:) - 0.5, 2*j_subc_dk84(:) - 0.5, 2*k_subc_dk84(:) - 0.5];
% [OLD] mm_subc   = transformPointsForward(tform_dk84_1, pts1_subc);
% [OLD] XsDK = mm_subc(:,1); YsDK = mm_subc(:,2); ZsDK = mm_subc(:,3);
% [OLD] % --- cortical
% [OLD] pts1_cort = [2*i_cort_dk84(:) - 0.5, 2*j_cort_dk84(:) - 0.5, 2*k_cort_dk84(:) - 0.5];
% [OLD] mm_cort   = transformPointsForward(tform_dk84_1, pts1_cort);
% [OLD] XcDK = mm_cort(:,1); YcDK = mm_cort(:,2); ZcDK = mm_cort(:,3);

% [NEW]: Real-space origin (0,0,0) should align with FSL MNI152 2mm convention (45, 63, 36)
% --- subcortical
[XsDK, YsDK, ZsDK] = s2r_ijk2mm(tform_ref, i_subc_dk84, j_subc_dk84, k_subc_dk84);
% --- cortical
[XcDK, YcDK, ZcDK] = s2r_ijk2mm(tform_ref, i_cort_dk84, j_cort_dk84, k_cort_dk84);

% (visual check DK84: TESTS T0.4)

%% 0.5. DBS80
addpath(genpath(P.atlas_dbs80)) % DBS80 subcorticals .nii from Jakub Vohryzek
% TODO ask for RSN-ROI MAP 7x80: https://www.biorxiv.org/content/10.1101/2020.11.04.368274v2.full //
mm2_dbs80 = niftiread('dbs80symm_2mm.nii'); % Cortical + subcortical / background sum(isnan(mm2_dbs80(:)))==0)
Ndbs80 = length(unique(mm2_dbs80(mm2_dbs80(:)~=0)));
regions_dbs80 = readtable('dbs80_nodes.xlsx'); % L1, L2,...LN,RN,...,R2,R1
labelsDBS80 = [regions_dbs80.id];

% Manual separation of subcortical/cortical numerical labels
subc_idx_within_labels_dbs80 = 32:49; % Manually
labelsDBS80_subc = labelsDBS80(subc_idx_within_labels_dbs80);
labelsDBS80_cort = labelsDBS80; % Labels are not exactly 1..90
labelsDBS80_cort(ismember(labelsDBS80, labelsDBS80_subc)) = [];

% --- Mask subcorticals in (i,j,k) array
mask_subc_dbs80 = ismember(mm2_dbs80, labelsDBS80_subc); 
mask_non_subc = ~mask_subc_dbs80; %~ismember(mm2_dbs80, labelsDBS80_subc);
% mm2 array with subc only
mm2_subc_dbs80 = mm2_dbs80;
mm2_subc_dbs80(mask_non_subc) = 0; clear mask_non_subc;
% --- Mask corticals in (i,j,k) array
mask_cort_dbs80 = (mm2_dbs80 ~= 0) & ~mask_subc_dbs80; 
% mm2 array with cort only
mm2_cort_dbs80 = mm2_dbs80;
mm2_cort_dbs80(mask_subc_dbs80) = 0; % mask_subc_dbs80 could also be defined from here with ismember(mm2_cort_dbs80, labelsDBS80_cort);

% --- 3D Array for DBS80 cortex ---
idcs_cort_dbs80 = find(mask_cort_dbs80);
[i_cort_dbs80, j_cort_dbs80, k_cort_dbs80] = ind2sub(size(mm2_cort_dbs80), idcs_cort_dbs80);
labDBS80_cort=mm2_cort_dbs80(idcs_cort_dbs80);
% --- 3D Array for DBS80 subcortex ---
idcs_subc_dbs80 = find(mask_subc_dbs80);
[i_subc_dbs80, j_subc_dbs80, k_subc_dbs80] = ind2sub(size(mm2_subc_dbs80), idcs_subc_dbs80);
labDBS80_subc=mm2_subc_dbs80(idcs_subc_dbs80);

% MNI Coordinates (from NIFTI information of dbs80symm_2mm.nii)
info_dbs80 = niftiinfo('dbs80symm_2mm.nii');
tform_dbs80 = info_dbs80.Transform;

% subcortical voxels -> mm
pts_s = [i_subc_dbs80(:), j_subc_dbs80(:), k_subc_dbs80(:)];
% [OLD] mm_s = transformPointsForward(tform_dbs80, pts_s);
% [OLD] XsDBS80 = mm_s(:,1); YsDBS80 = mm_s(:,2); ZsDBS80 = mm_s(:,3);
[XsDBS80, YsDBS80, ZsDBS80] = s2r_ijk2mm(tform_dbs80, i_subc_dbs80, j_subc_dbs80, k_subc_dbs80); % [NEW]

% cortical voxels -> mm
pts_c = [i_cort_dbs80(:), j_cort_dbs80(:), k_cort_dbs80(:)];
% [OLD] mm_c = transformPointsForward(tform_dbs80, pts_c);
% [OLD] XcDBS80 = mm_c(:,1); YcDBS80 = mm_c(:,2); ZcDBS80 = mm_c(:,3);
[XcDBS80, YcDBS80, ZcDBS80] = s2r_ijk2mm(tform_dbs80, i_cort_dbs80, j_cort_dbs80, k_cort_dbs80); %[NEW]


% (visual check DBS80: TESTS T0.5)
% (visual check manual pick of subcorticals: TESTS T0.15)

%%  0.6 Tian: S1(16), S2, S3, S4
addpath(genpath(P.atlas_tian)) % Tian subcorticals .nii
scales = {'S1','S2','S3','S4'};
Tian = struct();

for s = 1:numel(scales)
    S = scales{s};
    % ---------- Load NIfTI ----------
    niiFile = sprintf('Tian_Subcortex_%s_3T.nii', S);
    % Order R1,...,RN,L1,...,LN
    mm2 = niftiread(niiFile);
    Tian.(S).mm2 = mm2;
    % ---------- Load region names ----------
    txtFile = sprintf('Tian_Subcortex_%s_3T_label.txt', S);
    lines = readlines(txtFile);
    lines = strtrim(lines);
    lines(lines=="") = [];
    regions = cellstr(lines);
    Tian.(S).regions = regions;
    Tian.(S).N = length(unique(mm2(mm2~=0)));
    % ---------- Mask + indices ----------
    mask = mm2 ~= 0;
    idcs = find(mask);
    [i,j,k] = ind2sub(size(mm2), idcs);
    lab = mm2(idcs);
    Tian.(S).mask = mask;
    Tian.(S).idcs = idcs;
    Tian.(S).i = i;
    Tian.(S).j = j;
    Tian.(S).k = k;
    Tian.(S).lab = lab;

    % ---------- NIfTI transform ----------
    info = niftiinfo(niiFile);
    Tian.(S).info = info;
    Tian.(S).T = info.Transform.T;
    Tian.(S).tform = info.Transform;
    % ---------- Voxel -> real space (mm) ----------
    % [OLD]pts = [i(:), j(:), k(:)];
    % [OLD]mm = transformPointsForward(info.Transform, pts);
    % [OLD] Tian.(S).X = mm(:,1);
    % [OLD] Tian.(S).Y = mm(:,2);
    % [OLD] Tian.(S).Z = mm(:,3);
    [Tian.(S).X, Tian.(S).Y, Tian.(S).Z] = s2r_ijk2mm(info.Transform, i, j, k); % [NEW]
    fprintf('Tian %s loaded: %d regions, %d voxels\n', ...
        S, Tian.(S).N, numel(idcs));
end
% (visual check Tian S1..S4: TESTS T0.6)

%% 0.7 Idcs (array space) and subset_indices (within indices)
%  0.7.11 Combine indices (array space) only among identical mm2 grids
idcs = unique([
    idcs_cort(:); % from Schaefer in 7RSN file, equal to idcs_cort17
    idcs_cort_aal(:);
    idcs_subc_aal(:);
    idcs_cort_dk84(:);
    idcs_subc_dk84(:);
    idcs_cort_dbs80(:);
    idcs_subc_dbs80(:);
    Tian.S1.idcs(:);
    Tian.S2.idcs(:);
    Tian.S3.idcs(:);
    Tian.S4.idcs(:)
]);

fprintf('Combined number of unique voxels: %d\n', numel(idcs));

% (0.7 visual check of all atlases side by side: TESTS T0.7)

%%  0.8. subset_indices: positions within "idcs"
subset_indices = struct();

subset_indices.cort_schaefer = find(ismember(idcs, idcs_cort)); % same for idcs_cort17

subset_indices.cort_aal  = find(ismember(idcs, idcs_cort_aal));
subset_indices.subc_aal  = find(ismember(idcs, idcs_subc_aal));

subset_indices.cort_dk84 = find(ismember(idcs, idcs_cort_dk84));
subset_indices.subc_dk84 = find(ismember(idcs, idcs_subc_dk84));

subset_indices.cort_dbs80 = find(ismember(idcs, idcs_cort_dbs80));
subset_indices.subc_dbs80 = find(ismember(idcs, idcs_subc_dbs80));

% Tian as a nested struct for clarity
subset_indices.subc_tian = struct();
subset_indices.subc_tian.S1 = find(ismember(idcs, Tian.S1.idcs));
subset_indices.subc_tian.S2 = find(ismember(idcs, Tian.S2.idcs));
subset_indices.subc_tian.S3 = find(ismember(idcs, Tian.S3.idcs));
subset_indices.subc_tian.S4 = find(ismember(idcs, Tian.S4.idcs));

fprintf('Within combined idcs:\n');
fprintf('  Schaefer cortex: %d voxels\n', numel(subset_indices.cort_schaefer));
fprintf('  AAL cortex: %d voxels\n', numel(subset_indices.cort_aal));
fprintf('  AAL subcortex: %d voxels\n', numel(subset_indices.subc_aal));
fprintf('  DK84 cortex: %d voxels\n', numel(subset_indices.cort_dk84));
fprintf('  DK84 subcortex: %d voxels\n', numel(subset_indices.subc_dk84));
fprintf('  DBS80 cortex: %d voxels\n', numel(subset_indices.cort_dbs80));
fprintf('  DBS80 subcortex: %d voxels\n', numel(subset_indices.subc_dbs80));
fprintf('  Tian S1: %d voxels\n', numel(subset_indices.subc_tian.S1));
fprintf('  Tian S2: %d voxels\n', numel(subset_indices.subc_tian.S2));
fprintf('  Tian S3: %d voxels\n', numel(subset_indices.subc_tian.S3));
fprintf('  Tian S4: %d voxels\n', numel(subset_indices.subc_tian.S4));

%% 0.9. Map Schaefer to Yeo 7- and 17-RSN
%  Requires that the following already exist in workspace:
%  - mm2_schaefer, mm2_schaefer_17
%  - labels7, labels17, Nrsn, Nrsn_17
%  - Nschaefer
%  - idcs (combined indices) and subset_indices struct (optional but used below)
%  - DBS80 subcortex world coords: XsDBS80, YsDBS80, ZsDBS80, labDBS80_subc (optional panel)
%  - Schaefer centroid tables: SchaeferCOG (7) and SchaeferCOG_17 (17)
%  =======================
% Mapping vectors for Yeo 7 and 17 (no swap LH and RH!)


% Get Parcel --> RSN
SchaeferCOG.RSN = zeros(height(SchaeferCOG), 1);
SchaeferCOG_17.RSN = zeros(height(SchaeferCOG_17), 1);

for i = 1:height(SchaeferCOG) % WARNING: tables must have same length so this loop ok for both
    roi7  = string(SchaeferCOG.ROIName{i});
    roi17 = string(SchaeferCOG_17.ROIName{i});
    % Assign in table 7 RSN
    m7  = cellfun(@(p) contains(roi7,  p,  'IgnoreCase',true), labels7);
    k7  = find(m7,  1, 'first');
    if ~isempty(k7),  SchaeferCOG.RSN(i) = uint8(k7); else, warning('Unmatched 7-net: %s', roi7); end
    % Assign in table 17 RSN
    m17 = cellfun(@(p) contains(roi17, p, 'IgnoreCase',true), labels17);
    k17 = find(m17, 1, 'first');    % stops as soon as it finds one
    if ~isempty(k17), SchaeferCOG_17.RSN(i) = uint8(k17); else, warning('Unmatched 17-net: %s', roi17); end
end

% Mapping vector from Schaefer parcel number (1:Nschaefer) to RSN
mapSchaefer2Yeo   = zeros(Nschaefer,1);
mapSchaefer2Yeo(:) = SchaeferCOG.RSN(:);

mapSchaefer2Yeo17 = zeros(Nschaefer,1);
mapSchaefer2Yeo17(:) = SchaeferCOG_17.RSN(:);

% (0.9.1 visual check of the centroids of both Schaefer files: TESTS T0.9)
% WARNING: some numerical Schaefer1000 labels 
% of same centroids differ between 7 vs. 17 files!!!

%% 0.9.2. MNI volume arrays -> indices -> subset_indices for Yeo 7 and 17
mm2_yeo7 = zeros(size(mm2_schaefer), 'like', mm2_schaefer);
mask_yeo7 = mm2_schaefer ~= 0;
mm2_yeo7(mask_yeo7) = mapSchaefer2Yeo(mm2_schaefer(mask_yeo7));

mm2_yeo17 = zeros(size(mm2_schaefer_17), 'like', mm2_schaefer_17);
mask_yeo17 = mm2_schaefer_17 ~= 0;
mm2_yeo17(mask_yeo17) = mapSchaefer2Yeo17(mm2_schaefer_17(mask_yeo17));

% Extract voxel indices per RSN (index space)
idcs_rsn7  = cell(Nrsn,1);
idcs_rsn17 = cell(Nrsn_17,1);

for r = 1:Nrsn
    idcs_rsn7{r} = find(mm2_yeo7 == r);
end
for r = 1:Nrsn_17
    idcs_rsn17{r} = find(mm2_yeo17 == r);
end

% Update subset_indices for RSNs in the combined idcs list (if idcs exists)
if exist('idcs', 'var') && ~isempty(idcs)
    subset_indices.rsn  = cell(Nrsn,1);
    for r = 1:Nrsn
        subset_indices.rsn{r} = find(ismember(idcs, idcs_rsn7{r}));
        fprintf('RSN7  %2d: %d voxels in combined idcs\n', r, numel(subset_indices.rsn{r}));
    end

    subset_indices.rsn17 = cell(Nrsn_17,1);
    for r = 1:Nrsn_17
        subset_indices.rsn17{r} = find(ismember(idcs, idcs_rsn17{r}));
        fprintf('RSN17 %2d: %d voxels in combined idcs\n', r, numel(subset_indices.rsn17{r}));
    end
end

% Prepare Schaefer world transform (mm)
info_sch7 = niftiinfo('Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii');
tform_sch7 = info_sch7.Transform;

info_sch17 = niftiinfo('Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii');
tform_sch17 = info_sch17.Transform;

% Common mm limits based on Schaefer cortex voxels (7-network file)
[i_all, j_all, k_all] = ind2sub(size(mm2_schaefer), find(mm2_schaefer ~= 0));
% [OLD] mm_all = transformPointsForward(tform_sch7, [i_all(:), j_all(:), k_all(:)]);
% [OLD] xlim_mm = [min(mm_all(:,1)) max(mm_all(:,1))];
% [OLD] ylim_mm = [min(mm_all(:,2)) max(mm_all(:,2))];
% [OLD ]zlim_mm = [min(mm_all(:,3)) max(mm_all(:,3))];
[xa, ya, za] = s2r_ijk2mm(tform_ref, i_all, j_all, k_all); % [NEW]
xlim_mm = [min(xa) max(xa)]; ylim_mm = [min(ya) max(ya)]; zlim_mm = [min(za) max(za)]; % [NEW]

clear i_all j_all k_all mm_all

% (figures of the Yeo-7 / Yeo-17 RSNs in real space: TESTS T0.10 and T0.11)

%% Save variables for cluster code: compute voxels' FC(subc,cort) ---> assign RSN
if SAVE_OUTPUTS
    if EXCLUDE_CEREB
        baseFileName = 'Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat';
        fullFileName = fullfile(folder_out, baseFileName);
    elseif EXCLUDE_CEREB == false
        baseFileName = 'Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17_wCereb.mat';
        fullFileName = fullfile(folder_out, baseFileName);
    end
    % Save without or with cerebellum voxels
    save(fullFileName,"subjectIDs","idcs","subset_indices",'-v7.3')
end

% =======================
%% Pack and save meta-information and variables loaded/used
% Variable glossary / metadata (meta struct)
% =======================

meta = struct();

% ---- Global ----
meta.created = datestr(now);
meta.subjects.subjectIDs_var = 'subjectIDs';
meta.subjects.NSUB = NSUB;

meta.conventions = struct();
meta.conventions.labels = 'Vector of ROI IDs (per region), not per voxel';
meta.conventions.lab = 'Vector of labels per voxel aligned with the voxel index list';
meta.conventions.idcs = 'Linear indices into the 3D volume (index space)';
meta.conventions.i_j_k = 'Voxel subscripts aligned with idcs';
meta.conventions.X_Y_Z = 'Real-world coordinates (mm) aligned with idcs';
meta.conventions.mm2 = '3D label volume; values are ROI IDs, 0 is background';

meta.union = struct();
meta.union.idcs = idcs;
meta.union.subset_indices = subset_indices;
meta.union.note = 'Union defined in array index space; valid only if all contributing volumes are on the same voxel grid';

% Helper to safely store transform info
packTransform = @(info) struct( ...
    'ImageSize', info.ImageSize, ...
    'PixelDimensions', info.PixelDimensions, ...
    'TransformT', info.Transform.T ...
);

% ---- Schaefer (7) ----
meta.schaefer7 = struct();
meta.schaefer7.name = 'Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm';
meta.schaefer7.mm2 = 'mm2_schaefer';
meta.schaefer7.info = packTransform(niftiinfo('Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii'));

meta.schaefer7.labels7 = labels7;           % RSN names (per RSN)
meta.schaefer7.Nrsn = Nrsn;
meta.schaefer7.Nschaefer = Nschaefer;

meta.schaefer7.cortex = struct();
meta.schaefer7.cortex.mask = 'mask_cort';
meta.schaefer7.cortex.idcs = idcs_cort;      % per voxel
meta.schaefer7.cortex.i = i_cort;            % per voxel
meta.schaefer7.cortex.j = j_cort;
meta.schaefer7.cortex.k = k_cort;
meta.schaefer7.cortex.lab = labSchaefer7;    % per voxel labels 1..1000
meta.schaefer7.cortex.X = X;                 % per voxel mm coords
meta.schaefer7.cortex.Y = Y;
meta.schaefer7.cortex.Z = Z;

meta.schaefer7.centroids = struct();
meta.schaefer7.centroids.table = 'SchaeferCOG';
meta.schaefer7.centroids.order = 'CSV rows are L1..L500 then R1..R500';
meta.schaefer7.centroids.coords_cols = '[3 4 5] = (X,Y,Z) in mm';

% ---- Schaefer (17) ----
meta.schaefer17 = struct();
meta.schaefer17.name = 'Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm';
meta.schaefer17.mm2 = 'mm2_schaefer_17';
meta.schaefer17.info = packTransform(niftiinfo('Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii'));

meta.schaefer17.labels17 = labels17;
meta.schaefer17.Nrsn_17 = Nrsn_17;

meta.schaefer17.cortex = struct();
meta.schaefer17.cortex.mask = 'mask_cort17';
meta.schaefer17.cortex.idcs = idcs_cort17;
meta.schaefer17.cortex.i = i_cort17;
meta.schaefer17.cortex.j = j_cort17;
meta.schaefer17.cortex.k = k_cort17;
meta.schaefer17.cortex.lab = labSchaefer17;  % per voxel labels 1..1000
meta.schaefer17.cortex.X = X17;
meta.schaefer17.cortex.Y = Y17;
meta.schaefer17.cortex.Z = Z17;

meta.schaefer17.centroids = struct();
meta.schaefer17.centroids.table = 'SchaeferCOG_17';
meta.schaefer17.centroids.order = 'CSV rows are L1..L500 then R1..R500';
meta.schaefer17.centroids.coords_cols = '[3 4 5] = (X,Y,Z) in mm';

% ---- AAL ----
meta.aal = struct();
meta.aal.name = 'AAL';
meta.aal.mm2 = 'mm2_aal';
meta.aal.info = packTransform(niftiinfo('AAL.nii'));

meta.aal.labels = labelsAAL;  % per region IDs (1x90)
meta.aal.N = numel(labelsAAL);

meta.aal.subc_region_ids = labelsAAL_subc;  % per region IDs (not per voxel)

meta.aal.cortex = struct();
meta.aal.cortex.mm2 = 'mm2_cort_aal';
meta.aal.cortex.idcs = idcs_cort_aal;
meta.aal.cortex.i = i_cort_aal;
meta.aal.cortex.j = j_cort_aal;
meta.aal.cortex.k = k_cort_aal;
meta.aal.cortex.lab = labAAL_cort;  % per voxel labels
meta.aal.cortex.X = XcAAL;
meta.aal.cortex.Y = YcAAL;
meta.aal.cortex.Z = ZcAAL;

meta.aal.subcortex = struct();
meta.aal.subcortex.mm2 = 'mm2_subc_aal';
meta.aal.subcortex.idcs = idcs_subc_aal;
meta.aal.subcortex.i = i_subc_aal;
meta.aal.subcortex.j = j_subc_aal;
meta.aal.subcortex.k = k_subc_aal;
meta.aal.subcortex.lab = labAAL_subc; % per voxel labels
meta.aal.subcortex.X = XsAAL;
meta.aal.subcortex.Y = YsAAL;
meta.aal.subcortex.Z = ZsAAL;

% ---- DK84 ----
meta.dk84 = struct();
meta.dk84.name = 'IIT_GM_Desikan (DK84)';
meta.dk84.mm1 = 'mm1_dk84';
meta.dk84.mm2 = 'mm2_dk84';
meta.dk84.info_1mm = packTransform(niftiinfo('IIT_GM_Desikan_atlas.nii.gz'));
meta.dk84.downsample = '2x2x2 majority vote in label space (mode over 8 voxels)';

meta.dk84.labels = labelsDK84;  % per region IDs (vector)
meta.dk84.subc_region_ids = labelsDK84_subc;

meta.dk84.cortex = struct();
meta.dk84.cortex.mm2 = 'mm2_cort_dk84';
meta.dk84.cortex.idcs = idcs_cort_dk84;
meta.dk84.cortex.i = i_cort_dk84;
meta.dk84.cortex.j = j_cort_dk84;
meta.dk84.cortex.k = k_cort_dk84;
meta.dk84.cortex.lab = labDK84_cort;   % per voxel labels
meta.dk84.cortex.X = XcDK;
meta.dk84.cortex.Y = YcDK;
meta.dk84.cortex.Z = ZcDK;

meta.dk84.subcortex = struct();
meta.dk84.subcortex.mm2 = 'mm2_subc_dk84';
meta.dk84.subcortex.idcs = idcs_subc_dk84;
meta.dk84.subcortex.i = i_subc_dk84;
meta.dk84.subcortex.j = j_subc_dk84;
meta.dk84.subcortex.k = k_subc_dk84;
meta.dk84.subcortex.lab = labDK84_subc; % per voxel labels
meta.dk84.subcortex.X = XsDK;
meta.dk84.subcortex.Y = YsDK;
meta.dk84.subcortex.Z = ZsDK;

% ---- DBS80 ----
meta.dbs80 = struct();
meta.dbs80.name = 'dbs80symm_2mm';
meta.dbs80.mm2 = 'mm2_dbs80';
meta.dbs80.info = packTransform(niftiinfo('dbs80symm_2mm.nii'));

meta.dbs80.labels = labelsDBS80; % per region IDs
meta.dbs80.subc_region_ids = labelsDBS80_subc;

meta.dbs80.cortex = struct();
meta.dbs80.cortex.mm2 = 'mm2_cort_dbs80';
meta.dbs80.cortex.idcs = idcs_cort_dbs80;
meta.dbs80.cortex.i = i_cort_dbs80;
meta.dbs80.cortex.j = j_cort_dbs80;
meta.dbs80.cortex.k = k_cort_dbs80;
meta.dbs80.cortex.lab = labDBS80_cort;
meta.dbs80.cortex.X = XcDBS80;
meta.dbs80.cortex.Y = YcDBS80;
meta.dbs80.cortex.Z = ZcDBS80;

meta.dbs80.subcortex = struct();
meta.dbs80.subcortex.mm2 = 'mm2_subc_dbs80';
meta.dbs80.subcortex.idcs = idcs_subc_dbs80;
meta.dbs80.subcortex.i = i_subc_dbs80;
meta.dbs80.subcortex.j = j_subc_dbs80;
meta.dbs80.subcortex.k = k_subc_dbs80;
meta.dbs80.subcortex.lab = labDBS80_subc;
meta.dbs80.subcortex.X = XsDBS80;
meta.dbs80.subcortex.Y = YsDBS80;
meta.dbs80.subcortex.Z = ZsDBS80;

% ---- Tian (S1..S4) ----
meta.tian = struct();
meta.tian.name = 'Tian_Subcortex_*_3T';
meta.tian.scales = scales;

for s = 1:numel(scales)
    S = scales{s};
    meta.tian.(S) = struct();

    meta.tian.(S).mm2 = sprintf('Tian.%s.mm2', S);
    meta.tian.(S).labels = Tian.(S).regions;  % per region names
    meta.tian.(S).N = Tian.(S).N;

    meta.tian.(S).idcs = Tian.(S).idcs;
    meta.tian.(S).i = Tian.(S).i;
    meta.tian.(S).j = Tian.(S).j;
    meta.tian.(S).k = Tian.(S).k;
    meta.tian.(S).lab = Tian.(S).lab;         % per voxel labels

    meta.tian.(S).info = packTransform(Tian.(S).info);
    meta.tian.(S).X = Tian.(S).X;
    meta.tian.(S).Y = Tian.(S).Y;
    meta.tian.(S).Z = Tian.(S).Z;
end

% ---- RSN mapping ----
meta.yeo = struct();
meta.yeo.mapSchaefer2Yeo_var = 'mapSchaefer2Yeo';
meta.yeo.mapSchaefer2Yeo17_var = 'mapSchaefer2Yeo17';
meta.yeo.mm2_yeo7_var = 'mm2_yeo7';
meta.yeo.mm2_yeo17_var = 'mm2_yeo17';

meta.yeo.idcs_rsn7 = idcs_rsn7;
meta.yeo.idcs_rsn17 = idcs_rsn17;

meta.yeo.plotting = struct();
meta.yeo.plotting.xlim_mm = xlim_mm;
meta.yeo.plotting.ylim_mm = ylim_mm;
meta.yeo.plotting.zlim_mm = zlim_mm;

% Quick consistency flags
meta.quickchecks = struct();
meta.quickchecks.lab_vs_labels_explainer = ...
    "labels_* are ROI IDs (per region). lab_* are per-voxel label values aligned to idcs_*.";


%% Save meta-information
% Added save (the original script contains none): record of the atlas definitions.
if SAVE_OUTPUTS
    if EXCLUDE_CEREB
        baseFileName = 'Variables_metainformation.mat';
        fullFileName = fullfile(folder_out, baseFileName);
    elseif EXCLUDE_CEREB == false
        baseFileName = 'Variables_metainformation_wCereb.mat';
        fullFileName = fullfile(folder_out, baseFileName);
    end
    % Save without or with cerebellum voxels
    save(fullFileName, 'meta', 'labels7', 'labels17', '-v7.3')
end
