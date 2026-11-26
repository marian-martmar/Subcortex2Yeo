%% SUBC2RSN | STEP 2 (CLUSTER) - FC(subcortical voxel, network) per subject, group FC and network labels
% -------------------------------------------------------------------------
% PURPOSE
%   For every subject, correlate each subcortical voxel with the cortical (Schaefer)
%   voxels of each Yeo network (7 or 17) and average those correlations within the network:
%       r(i,j)        Pearson correlation, subcortical voxel i, cortical voxel j of network n
%       FC_sub(i,n) = tanh( mean_j atanh(r(i,j)) )     Fisher average of the correlations (RSN_AVERAGE = 'fisher')
%                   = mean_j r(i,j)                    plain average                      (RSN_AVERAGE = 'r')
%   FC_sub is a correlation value (r units) in both cases. 
%   All subcortical x cortical correlations per network are computed, discarded after averaging per RSN
%   (the full subcortical x cortical matrix is never stored). Then
%       label_sub(i) = argmax_n FC_sub(i,n)                   winner network, one subject
%       FC_group     = tanh( mean_subjects atanh(FC_sub) )    Fisher average across subjects (r units)
%       label_group  = argmax_n FC_group(i,n)                 winner network, group
%
% DEPENDENCIES on step 1 script's outputs    ---> must be pasted onto P.cluster_ts_dir:
%              Variables_From_1_xRL_dn_f.mat ---> subjectIDs, idcs, subset_indices (incl. .xRL_dn_f), idcs_filt, meta (not loaded)
%              <ID>_xRL_dn_f.mat             ---> xRL_dn_f [Nt x Nf], subset_indices (.xRL_dn_f), idcs_filt             
%
% OUTPUTS into P.cluster_out. Every file also holds 'info' (see METAINFO).
%   <ID>_avgFC7.mat    FC7_sub  [Nsubc x 7]   single, r units    subjectID, valid_rows, Nt_used
%   <ID>_avgFC17.mat   FC17_sub [Nsubc x 17]
%   FC7_group.mat      FC7_group [Nsubc x 7]  Fisher average across subjects (r units)
%                      FC7_group_meanr        arithmetic average of r across subjects (comparison)
%                      cumZ7, cumFC7, nValid7 sums of atanh(FC7_sub), of FC7_sub, and subject counts
%                      defs                   every row / column definition used
%   FC17_group.mat     same for 17 networks
%   Assigned_7.mat     assigned7_group    [Nsubc x 1]     uint8, label from FC7_group (0 = no valid data)
%                      assigned7_subjects {NSUB x 2}      {subjectID, assigned7_sub [1 x Nsubc]}
%                      assigned7_all      [NSUB x Nsubc]  uint8, subjects x voxels
%   Assigned_17.mat    same for 17 networks
%
% ROWS AND COLUMNS (used in later steps)
%   Rows    = subset_idcs_subc: sorted union of the subcortical positions within idcs of
%             AAL, DBS80, DK84 and Tian S1..S4. Linear index of row i: idcs(subset_idcs_subc(i)).
%             info.rows.in_aal / in_dbs80 / in_dk84 / in_tian tell which atlas each row belongs to;
%             info.rows.is_overlap marks rows that are also Schaefer cortex
%   Columns = networks, in the order of C.rsn7.labels / C.rsn17.labels
%
% METAINFO ('info', saved in every file)
%   info.C                     the configuration (subc2rsn_config) used here
%   info.C_step1               the configuration saved by step 1 with the time series, if step 1 saved it
%   info.switches              OVERLAPPING_VXLS, RSN_AVERAGE, N_DISCARD
%   info.EXCLUDE_CEREB_step1   cerebellum setting of the time series (checked for every subject)
%   info.labels7 / labels17
%   info.subjectIDs
%   info.step1_file
%   info.created
%   info.rows                  rows (positions, linear indices, atlas membership, overlap flag)
%   info.n_cort_used7          cortical voxels per network actually averaged (analogous for 17)
%
% SWITCHES
%   OVERLAPPING_VXLS  true :     voxels that are Schaefer cortex AND subcortical in any atlas are
%                                kept in the cortical network averages (original behaviour)
%                     false:     they are left out of the network averages. 
%                                - Both cases keep their FC row as subcortical voxels
%                                  so tables can count or drop them later using info.rows.is_overlap
%   RSN_AVERAGE       'fisher'   Fisher average the voxel-correlations into a rsn-correlation (original)
%                     'r'        Straight average of voxel-correlations into a rsn-correlation
%   N_DISCARD         frames removed at each end before FC (filtfilt edges);
%                     read from C.discard_t_filt in subc2rsn_config; (0 = original)
%
% NOTES
%   1. Flat or non-finite series: a subcortical row with no signal in a subject gets FC = NaN
%      and label 0 in that subject, and is left out of the group average for that subject
%      (nValid counts the subjects used). A flat cortical voxel is left out of its network
%      average in that subject. Previous version would have given: r = 0 everywhere and label 1 (Vis).
%   2. Group FC is saved already back-transformed (FC7_group, r units). The sums are kept.
%
% MEMORY  The largest network (Default, ~32,000 cortical voxels) gives a correlation block of
%   Nsubc x 32,000 doubles (~3 GB for 11,604 rows) plus temporaries of the same size while it is
%   clipped and Fisher-transformed: request about 12 GB per job (more with cerebellum rows).
%
% Author: Marian Martinez-Marin. v1 December 2025. v2 November 2026
% -------------------------------------------------------------------------

% clear; clc;

P = subc2rsn_paths();
C = subc2rsn_config();

% functions
addpath(P.cluster_functions)
%% SWITCHES
OVERLAPPING_VXLS     = C.OVERLAPPING_VXLS   % usage in RSN-averages FC(subc,RSN) of voxels both schafer & subc for some atlas
RSN_AVERAGE          = 'fisher';            % 'fisher' | 'r'
N_DISCARD            = C.discard_t_filt;    % frames dropped at each end before FC (set in subc2rsn_config)
CHECK_VOXEL_IDENTITY = true;                % stop if a time-series file does not hold the expected voxels
assert(any(strcmp(RSN_AVERAGE, {'fisher', 'r'})), 'RSN_AVERAGE must be ''fisher'' or ''r''.');

fprintf('Summary of switches and configuration settings that affect the results that will come out this code...\n');

%% 2.0. Inputs written by step 1
% Files are read with its full path (no addpath), so that a file of another run
% (e.g. with / without cerebellum) can never be picked up from the MATLAB path.
step1File = fullfile(P.cluster_ts_dir, 'Variables_From_1_xRL_dn_f.mat');
V = load(step1File, 'subjectIDs', 'idcs', 'subset_indices', 'idcs_filt', 'EXCLUDE_CEREB');
EXCLUDE_CEREB_step1 = V.EXCLUDE_CEREB;          % cerebellum setting the time series were made with
subjectIDs     = V.subjectIDs;
idcs           = V.idcs;
subset_indices = V.subset_indices;
idcs_filt      = V.idcs_filt;
clear V;
NSUB = numel(subjectIDs);

outputFolder = P.cluster_out;
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

labels7  = C.rsn7.labels;   Nrsn   = numel(labels7);
labels17 = C.rsn17.labels;  Nrsn17 = numel(labels17);

% Verify loaded idcs_filt (linear indices of 3d MNI array) match expected idcs(subset_indices.xRL_dn_f)
if CHECK_VOXEL_IDENTITY && ~isequal(double(idcs(subset_indices.xRL_dn_f(:))), double(idcs_filt(:)))
    error('idcs_filt is not idcs(subset_indices.xRL_dn_f), CHECK %s: ', step1File);
end

%% 2.1. Rows (subcortical union) and cortical voxels, as positions within idcs
% -------------------------------------------------------------------------
% Reminder of the levels of indexing, all set in step 0:
%   linear index  position of a voxel in the 91x109x91 MNI array. The SAME number addresses
%                 that voxel in every 3-D atlas volume (mm2_*) and in every frame of the 4-D
%                 HCP time series (ts_vol(:,:,:,t)); ind2sub gives its [i j k].
%   idcs          sorted list of the linear indices of every voxel used by any atlas.
%   position      a voxel's place WITHIN idcs: position p means voxel idcs(p).
%                 subset_indices.<atlas> lists positions, not linear indices.
%
% Variables built here (all are positions within idcs):
%   subset_idcs_cort       Schaefer cortex: the voxels averaged into each Yeo network.
%   subset_idcs_subc       Union of the subcortical voxels of AAL, DBS80, DK84 and Tian S1..S4
%   Nsubc =                Number of rows in any  FC/label array
%                          Sorted: row n of every FC/label array of this script is voxel idcs(subset_idcs_subc(n))
%
%   pos_overlap            voxels that are Schaefer cortex AND subcortical in some atlas
%   is_overlap_row         [Nsubc x 1] true for the rows that are such voxels (saved in info.rows.is_overlap)
%   subset_idcs_cort_used  cortical voxels that enter the network averages:
%   OVERLAPPING_VXLS       = true  -> all Schaefer voxels (original behaviour)
%                          = false -> Schaefer voxels minus pos_overlap
%                          Overlap voxels ALWAYS keep their row; 
%                          Only decides whether also used as cortical voxels of their network

subset_idcs_cort = unique(subset_indices.cort_schaefer(:));

subc_pos = [subset_indices.subc_aal(:); subset_indices.subc_dbs80(:); subset_indices.subc_dk84(:)];
tianScales = fieldnames(subset_indices.subc_tian); % redundant subset of positions within linear idcs
for s = 1:numel(tianScales) % include subsets of subcortical linear indices from tian scales
    subc_pos = [subc_pos; subset_indices.subc_tian.(tianScales{s})(:)]; %#ok<AGROW>
end
subset_idcs_subc = unique(subc_pos); % retain only unique positions within linear idcs for subcortical voxels
clear subc_pos;
Nsubc = numel(subset_idcs_subc);

% Voxels that are Schaefer cortex and also subcortical in some atlas
pos_overlap    = intersect(subset_idcs_cort, subset_idcs_subc);
is_overlap_row = ismember(subset_idcs_subc, pos_overlap);
if OVERLAPPING_VXLS
    subset_idcs_cort_used = subset_idcs_cort;
else
    subset_idcs_cort_used = setdiff(subset_idcs_cort, pos_overlap);
end
if OVERLAPPING_VXLS, txt = 'kept in'; else, txt = 'left out of'; end
fprintf('Schaefer cortex: %d voxels | subcortical rows: %d | overlap: %d (%s network averages)\n', ...
    numel(subset_idcs_cort), Nsubc, numel(pos_overlap), txt);

%% 2.2. Columns of xRL_dn_f (positions saved by step 1)
% Maps these (all) subcortical and (kept) cortical positions
% onto columns of xRL_dn_f (subset_indices.xRL_dn_f, saved by step 1).
%
% Example with 7 voxels:
%   idcs                     = [1043 2210 3388 4502 5307 8001 9120]'
%                              linear indices (each being a same voxel in every mm2_* and every ts_vol frame)
%
%   cort_schaefer            = [1 2 3 4 5]'   -> voxels 1043, 2210, 3388, 4502, 5307
%   subc_aal                 = [5 6]'         -> voxels 5307, 8001
%   subc_tian.S1 (..S4)      = [6 7]'         -> voxels 8001, 9120
%
%   subset_idcs_cort         = [1 2 3 4 5]'
%   subset_idcs_subc         = [5 6 7]'       -> FC rows: row 1 = 5307, row 2 = 8001, row 3 = 9120 (Nsubc = 3)
%
%   pos_overlap              = 5              -> voxel 5307 is Schaefer cortex AND AAL subcortex
%   is_overlap_row           = [1 0 0]'       -> row 1 is that voxel
%   subset_idcs_cort_used    = [1 2 3 4 5]'   (true)  voxel 5307 is also a cortical voxel of its network,
%                                                     so row 1 includes its correlation with itself (r ~ 1)
%                            = [1 2 3 4]'     (false) voxel 5307 is only a row; network averages would use
%                                                     1043, 2210, 3388, 4502

subset_idcs_f = subset_indices.xRL_dn_f(:); % subsets of linear idcs included as columns of xRL_dn_f(time,voxels)

[okC, cort_cols] = ismember(subset_idcs_cort_used, subset_idcs_f); % select xRL_dn_f column-voxels to be used into RSN-averages (columns of FC)
[okS, subc_cols] = ismember(subset_idcs_subc,      subset_idcs_f); % select xRL_dn_f column-voxels to be kept as subc voxels    (rows of FC)

if any(~okC), error('%d Schaefer cortex voxels are missing from xRL_dn_f.', sum(~okC)); end
if any(~okS), error('%d subcortical voxels are missing from xRL_dn_f.', sum(~okS)); end

%% 2.3. Network membership as columns of xRL_cort (= subset_idcs_cort_used order)
if ~isfield(subset_indices, 'rsn') || numel(subset_indices.rsn) ~= Nrsn
    error('subset_indices.rsn missing or with %d networks instead of %d.', numel(subset_indices.rsn), Nrsn);
end
if ~isfield(subset_indices, 'rsn17') || numel(subset_indices.rsn17) ~= Nrsn17
    error('subset_indices.rsn17 missing or with %d networks instead of %d.', numel(subset_indices.rsn17), Nrsn17);
end

rsn7_cols_in_cort  = cell(Nrsn, 1);
rsn17_cols_in_cort = cell(Nrsn17, 1);
for r = 1:Nrsn
    [~, rsn7_cols_in_cort{r}] = ismember(intersect(subset_idcs_cort_used, subset_indices.rsn{r}(:)), subset_idcs_cort_used);
end
for r = 1:Nrsn17
    [~, rsn17_cols_in_cort{r}] = ismember(intersect(subset_idcs_cort_used, subset_indices.rsn17{r}(:)), subset_idcs_cort_used);
end

n_cort_used7  = cellfun(@numel, rsn7_cols_in_cort);
n_cort_used17 = cellfun(@numel, rsn17_cols_in_cort);
assert(sum(n_cort_used7)  == numel(subset_idcs_cort_used), 'Some cortical voxels are in no Yeo-7 network, or in two.');
assert(sum(n_cort_used17) == numel(subset_idcs_cort_used), 'Some cortical voxels are in no Yeo-17 network, or in two.');
for r = 1:Nrsn,   fprintf('RSN7  %2d %-12s %6d cortical voxels\n', r, labels7{r},  n_cort_used7(r));  end
for r = 1:Nrsn17, fprintf('RSN17 %2d %-12s %6d cortical voxels\n', r, labels17{r}, n_cort_used17(r)); end

%% 2.4. Metainfo saved in every output file
info = struct();
info.created    = datestr(now);
info.script     = mfilename;
info.C          = C;
info.switches   = struct('OVERLAPPING_VXLS', OVERLAPPING_VXLS, 'RSN_AVERAGE', RSN_AVERAGE, 'N_DISCARD', N_DISCARD);
info.labels7    = labels7;
info.labels17   = labels17;
info.subjectIDs = subjectIDs;
info.step1_file = step1File;
info.EXCLUDE_CEREB_step1 = EXCLUDE_CEREB_step1;
info.rows = struct( ...
    'subset_idcs_subc',        subset_idcs_subc, ...                          % positions within idcs
    'idcs_subc',        idcs(subset_idcs_subc), ...                           % linear indices, 91x109x91 volume
    'in_aal',     ismember(subset_idcs_subc, subset_indices.subc_aal), ...
    'in_dbs80',   ismember(subset_idcs_subc, subset_indices.subc_dbs80), ...
    'in_dk84',    ismember(subset_idcs_subc, subset_indices.subc_dk84), ...
    'in_tian',    ismember(subset_idcs_subc, subset_indices.subc_tian.S1), ...
    'is_overlap', is_overlap_row);
info.n_cort_used7  = n_cort_used7;
info.n_cort_used17 = n_cort_used17;

% configuration saved by step 1 next to the time series (if it was saved)
f1 = fullfile(P.cluster_ts_dir, [subjectIDs{1} '_xRL_dn_f.mat']);
w1 = whos('-file', f1);
if any(strcmp({w1.name}, 'C'))
    tmp = load(f1, 'C');  info.C_step1 = tmp.C;  clear tmp;
end

% full definitions, saved once in the group files
defs = struct('subset_idcs_f', subset_idcs_f, 'subset_idcs_cort', subset_idcs_cort, ...
              'subset_idcs_cort_used', subset_idcs_cort_used, 'subset_idcs_subc', subset_idcs_subc, ...
              'pos_overlap', pos_overlap, 'cort_cols', cort_cols, 'subc_cols', subc_cols, ...
              'rsn7_cols_in_cort', {rsn7_cols_in_cort}, 'rsn17_cols_in_cort', {rsn17_cols_in_cort});

%% 2.4b Run summary to be printed at the top of the cluster .out file
% Step-1 settings read from step 1 outputs: they describe how the time series were really made
% Step-2 settings come from this script 
tmp = whos('-file', step1File);
if any(strcmp({tmp.name}, 'C')), tmp = load(step1File, 'C'); C1 = tmp.C; else, C1 = struct(); end
clear tmp;
step1 = struct('FWHM', NaN, 'voxelsize', NaN, 'TR', NaN, 'fhi', NaN, 'flp', NaN, 'DETRND_PreFilt', NaN, 'ZSCORE_PreFilt', NaN);
fn1 = fieldnames(step1);
for q = 1:numel(fn1)                                       % what step 1 used (NaN = not saved by step 1)
    if isfield(C1, fn1{q}), step1.(fn1{q}) = C1.(fn1{q}); end
end
w_ts = whos('-file', f1, 'xRL_dn_f');                      % size of the time series, without loading them
Nt_saved = w_ts.size(1);
yn = @(x) char(string(x));                                 % true / false / NaN as text
[~, host] = system('hostname');
if EXCLUDE_CEREB_step1, cereb = 'excluded'; else, cereb = 'included as subcortex'; end
if OVERLAPPING_VXLS, ovl = 'kept'; else, ovl = 'left out'; end
tian_same = isequal(subset_indices.subc_tian.S1, subset_indices.subc_tian.S2, subset_indices.subc_tian.S3, subset_indices.subc_tian.S4);
line = repmat('=', 1, 92);

fprintf('\n%s\n SUBC2RSN | STEP 2 (cluster) | run summary\n%s\n', line, line);
fprintf(' started  %s | host %s | SLURM job %s | MATLAB %s\n', datestr(now), strtrim(host), getenv('SLURM_JOB_ID'), version);
fprintf(' inputs   %s  (+ %d subject files <ID>_xRL_dn_f.mat)\n', step1File, NSUB);
fprintf(' outputs  %s\n', outputFolder);

fprintf('\n VOXELS (MNI152 2-mm grid, 91x109x91)\n');
fprintf('   project voxels (idcs)                     %7d\n', numel(idcs));
fprintf('   cerebellum (steps 0-1)                    %s\n', cereb);
fprintf('   cortex     Schaefer | AAL | DBS80 | DK84  %7d | %6d | %6d | %6d\n', numel(subset_indices.cort_schaefer), ...
    numel(subset_indices.cort_aal), numel(subset_indices.cort_dbs80), numel(subset_indices.cort_dk84));
fprintf('   subcortex  Tian S1  | AAL | DBS80 | DK84  %7d | %6d | %6d | %6d  (Tian S2-S4 same voxels: %s)\n', ...
    numel(subset_indices.subc_tian.S1), numel(subset_indices.subc_aal), numel(subset_indices.subc_dbs80), numel(subset_indices.subc_dk84),yn(tian_same));
if isfield(subset_indices, 'xRL_dn')
    fprintf('   step 1 denoised voxels (xRL_dn)           %7d  cortex of every atlas + all subcortex\n', numel(subset_indices.xRL_dn));
end
fprintf('   step 1 filtered voxels (xRL_dn_f)         %7d  Schaefer cortex + all subcortex\n', numel(subset_indices.xRL_dn_f));
fprintf('   step 2 rows (subcortical union)           %7d  of which %d are also Schaefer cortex\n', Nsubc, numel(pos_overlap));
fprintf('   step 2 cortical voxels in network averages %6d  (Schaefer %d; overlap voxels %s)\n', ...
    numel(subset_idcs_cort_used), numel(subset_idcs_cort), ovl);

fprintf('\n TIME SERIES\n');
fprintf('   step 1  input: HCP <ID>_rfMRI_REST1_LR_hp2000_clean.nii.gz (FIX-cleaned, hp2000)\n');
fprintf('     1. spatial smoothing        Gaussian FWHM %g mm on the whole volume (sigma %.2f voxels of %g mm)\n', ...
    step1.FWHM, step1.FWHM / step1.voxelsize / 2.355, step1.voxelsize);
fprintf('     2. denoising                Wishart filter (icaDim) on the xRL_dn voxels\n');
fprintf('     3. demean                   %s  (DETRND_PreFilt)\n', yn(step1.DETRND_PreFilt));
fprintf('     4. divide by SD             %s  (ZSCORE_PreFilt)\n', yn(step1.ZSCORE_PreFilt));
fprintf('     5. temporal filter          low-pass Butterworth order 2, cutoff %g Hz, filtfilt (zero-phase), TR %g s\n', ...
    step1.fhi, step1.TR);
fprintf('                                 (flp = %g Hz is defined but not used: no high-pass)\n', step1.flp);
fprintf('     frames saved per subject    %d\n', Nt_saved);
fprintf('   step 2\n');
fprintf('     6. frames dropped each end  %d  (N_DISCARD) -> frames used %d\n', N_DISCARD, Nt_saved - 2*N_DISCARD);
fprintf('     7. demean + unit SD         always, per voxel, on the frames used (makes X''*Y/(T-1) the Pearson r;\n');
fprintf('                                 so steps 3-4 above change the series but not the correlations)\n');
fprintf('     8. network average          %s  (RSN_AVERAGE)\n', RSN_AVERAGE);
fprintf('     9. group average            Fisher: tanh(mean over subjects of atanh(FC_sub))\n');

% Step-1 fields in the config used here must describe the same processing as the files
diffs = {};
for q = 1:numel(fn1)
    if isfield(C, fn1{q}) && ~isequaln(C.(fn1{q}), step1.(fn1{q}))
        diffs{end+1} = sprintf('%s (step 1: %s, config now: %s)', fn1{q}, mat2str(step1.(fn1{q})), mat2str(C.(fn1{q}))); %#ok<SAGROW>
    end
end
fprintf('\n CONFIG CHECK\n');
if isempty(diffs)
    fprintf('   step-1 settings in the current config equal those saved with the time series\n');
else
    fprintf('   WARNING: the current config differs from the settings saved with the time series:\n');
    fprintf('     %s\n', diffs{:});
    fprintf('   The summary above reports what step 1 actually used.\n');
end
fprintf('%s\n\n', line);
clear C1 step1 fn1 w_ts yn host cereb ovl tian_same line diffs q;

%% 2.5. Preallocate group outputs: 
      % Matrix for cumulative sums (of fisher transformed / raw correlation ) and NSUB=nValid7/17
cumZ7  = zeros(Nsubc, Nrsn);    cumFC7  = zeros(Nsubc, Nrsn);    nValid7  = zeros(Nsubc, Nrsn);
cumZ17 = zeros(Nsubc, Nrsn17);  cumFC17 = zeros(Nsubc, Nrsn17);  nValid17 = zeros(Nsubc, Nrsn17);
assigned7_all      = zeros(NSUB, Nsubc, 'uint8');
assigned17_all     = zeros(NSUB, Nsubc, 'uint8');
assigned7_subjects  = cell(NSUB, 2);
assigned17_subjects = cell(NSUB, 2);

%% 2.6. Subject loop
for sub = 1:NSUB
    subjectID = subjectIDs{sub};
    fprintf('Subject %d/%d: %s\n', sub, NSUB, subjectID);

    % ---- 2.6.1 Load the filtered time series and check its columns -----------------------------------
    % Every subject file must hold the same voxels, in the same order, as the general file from step 1
    D = load(fullfile(P.cluster_ts_dir, [subjectID '_xRL_dn_f.mat']), 'xRL_dn_f', 'subset_indices', 'idcs_filt', 'EXCLUDE_CEREB');
    if ~isequal(D.EXCLUDE_CEREB, EXCLUDE_CEREB_step1)
        error('%s: made with EXCLUDE_CEREB = %d, but %s says %d.', subjectID, D.EXCLUDE_CEREB, step1File, EXCLUDE_CEREB_step1);
    end
    if ~isequal(D.subset_indices.xRL_dn_f(:), subset_idcs_f)
        error('%s: columns of xRL_dn_f differ from %s.', subjectID, step1File);
    end
    if CHECK_VOXEL_IDENTITY && ~isequal(double(D.idcs_filt(:)), double(idcs(subset_idcs_f)))
        error('%s: idcs_filt does not match idcs(subset_indices.xRL_dn_f).', subjectID);
    end
    
    xRL_dn_f = D.xRL_dn_f;                                  % [Nt x Nf]
    clear D;

    % ---- 2.6.2 Optionally: drop filtfilt edges -------------------
    if N_DISCARD > 0
        xRL_dn_f = xRL_dn_f(N_DISCARD+1 : end-N_DISCARD, :);
    end
    Nt_used = size(xRL_dn_f, 1);

    % ---- 2.6.3 Split into cortex (network voxels) and subcortex (rows) --------------------
    xRL_cort = xRL_dn_f(:, cort_cols);                      % [Nt_used x Ncort_used], order = subset_idcs_cort_used
    xRL_subc = xRL_dn_f(:, subc_cols);                      % [Nt_used x Nsubc],      order = subset_idcs_subc
    clear xRL_dn_f;

    % ---- 2.6.4 Demean and scale each voxel to unit SD -------------------------------------
    % After this, X' * Y / (Nt_used - 1) is the Pearson correlation between columns of X and Y.
    % A voxel with no signal in this subject (SD < eps) or with NaN/Inf values is set to 0
    % and flagged as not valid, so it is left out instead of contributing r = 0.
    xRL_cort   = detrend(xRL_cort, 'constant');
    sd_cort    = std(xRL_cort, 0, 1);
    valid_cort = isfinite(sd_cort) & sd_cort >= eps;        % [1 x Ncort_used]
    sd_cort(~valid_cort)     = 1;
    xRL_cort(:, ~valid_cort) = 0;
    xRL_cort   = bsxfun(@rdivide, xRL_cort, sd_cort);

    xRL_subc   = detrend(xRL_subc, 'constant');
    sd_subc    = std(xRL_subc, 0, 1);
    valid_rows = (isfinite(sd_subc) & sd_subc >= eps)';     % [Nsubc x 1]
    sd_subc(~valid_rows)     = 1;
    xRL_subc(:, ~valid_rows) = 0;
    xRL_subc   = bsxfun(@rdivide, xRL_subc, sd_subc);
    clear sd_cort sd_subc;

    % ---- 2.6.5 Yeo-7: correlations with every voxel of a network, averaged over the network
    FC7_sub = nan(Nsubc, Nrsn, 'single');
    for r = 1:Nrsn
        cols = rsn7_cols_in_cort{r};                        % columns of xRL_cort in network r
        cols = cols(valid_cort(cols));                      % without flat voxels of this subject
        if isempty(cols), continue; end
        R = (xRL_subc' * xRL_cort(:, cols)) / (Nt_used - 1);   % Pearson r, [Nsubc x n_r]
        R = max(min(R, 0.999999), -0.999999);               % keep |r| < 1 so that atanh is finite
        if strcmp(RSN_AVERAGE, 'fisher')
            FC7_sub(:, r) = tanh(mean(atanh(R), 2));        % Fisher average of r, back to r units
        else
            FC7_sub(:, r) = mean(R, 2);                     % plain average of r
        end
        clear R;
    end

    % ---- 2.6.6 Yeo-17: the same ------------------------------------------------------------
    FC17_sub = nan(Nsubc, Nrsn17, 'single');
    for r = 1:Nrsn17
        cols = rsn17_cols_in_cort{r};
        cols = cols(valid_cort(cols));
        if isempty(cols), continue; end
        R = (xRL_subc' * xRL_cort(:, cols)) / (Nt_used - 1);
        R = max(min(R, 0.999999), -0.999999);
        if strcmp(RSN_AVERAGE, 'fisher')
            FC17_sub(:, r) = tanh(mean(atanh(R), 2));
        else
            FC17_sub(:, r) = mean(R, 2);
        end
        clear R;
    end
    clear xRL_cort xRL_subc;

    % Rows with no signal in this subject have no FC
    FC7_sub(~valid_rows, :)  = NaN;
    FC17_sub(~valid_rows, :) = NaN;

    % ---- 2.6.7 Winner network of each voxel in this subject (0 = no valid data) -----------
    [~, assigned7_sub] = max(FC7_sub, [], 2);               % column of the largest FC in each row
    assigned7_sub = uint8(assigned7_sub);
    assigned7_sub(any(isnan(FC7_sub), 2)) = 0;
    [~, assigned17_sub] = max(FC17_sub, [], 2);
    assigned17_sub = uint8(assigned17_sub);
    assigned17_sub(any(isnan(FC17_sub), 2)) = 0;

    assigned7_all(sub, :)       = assigned7_sub';           % subjects x voxels
    assigned17_all(sub, :)      = assigned17_sub';
    assigned7_subjects(sub, :)  = {subjectID, assigned7_sub'};
    assigned17_subjects(sub, :) = {subjectID, assigned17_sub'};

    % ---- 2.6.8 Save this subject's FC (r units) with the metainfo -------------------------
    save(fullfile(outputFolder, sprintf('%s_avgFC7.mat',  subjectID)), 'FC7_sub',  'subjectID', 'valid_rows', 'Nt_used', 'info', '-v7.3');
    save(fullfile(outputFolder, sprintf('%s_avgFC17.mat', subjectID)), 'FC17_sub', 'subjectID', 'valid_rows', 'Nt_used', 'info', '-v7.3');

    % ---- 2.6.9 Add this subject to the group sums ------------------------------------------
    % cumZ: sum of Fisher z (for the group Fisher average)
    % cumFC: sum of r (for comparison)
    % nValid: number of subjects that contributed to each entry (NaN entries are skipped)
    ok7 = isfinite(FC7_sub);
    z7  = atanh(double(FC7_sub));
    cumZ7(ok7)   = cumZ7(ok7)   + z7(ok7);
    cumFC7(ok7)  = cumFC7(ok7)  + double(FC7_sub(ok7));
    nValid7(ok7) = nValid7(ok7) + 1;

    ok17 = isfinite(FC17_sub);
    z17  = atanh(double(FC17_sub));
    cumZ17(ok17)   = cumZ17(ok17)   + z17(ok17);
    cumFC17(ok17)  = cumFC17(ok17)  + double(FC17_sub(ok17));
    nValid17(ok17) = nValid17(ok17) + 1;

    fprintf('  valid subcortical rows: %d of %d | frames used: %d\n', sum(valid_rows), Nsubc, Nt_used);
    clear FC7_sub FC17_sub assigned7_sub assigned17_sub ok7 ok17 z7 z17;
end

%% 2.7. Group FC (r units) and group labels
% Fisher average across subjects, back to r units (NaN where no subject has data)
FC7_group  = tanh(cumZ7  ./ nValid7);
FC17_group = tanh(cumZ17 ./ nValid17);
% Arithmetic average of r across subjects, kept for comparison
FC7_group_meanr  = cumFC7  ./ nValid7;
FC17_group_meanr = cumFC17 ./ nValid17;

% Winner network from the group FC (0 = no valid data)
[~, assigned7_group] = max(FC7_group, [], 2);
assigned7_group = uint8(assigned7_group);
assigned7_group(any(isnan(FC7_group), 2)) = 0;
[~, assigned17_group] = max(FC17_group, [], 2);
assigned17_group = uint8(assigned17_group);
assigned17_group(any(isnan(FC17_group), 2)) = 0;

%% 2.8. Save group files
save(fullfile(outputFolder, 'FC7_group.mat'),  'FC7_group',  'FC7_group_meanr',  'cumZ7',  'cumFC7',  'nValid7',  'defs', 'info', '-v7.3');
save(fullfile(outputFolder, 'FC17_group.mat'), 'FC17_group', 'FC17_group_meanr', 'cumZ17', 'cumFC17', 'nValid17', 'defs', 'info', '-v7.3');
save(fullfile(outputFolder, 'Assigned_7.mat'),  'assigned7_group',  'assigned7_subjects',  'assigned7_all',  'info', '-v7.3');
save(fullfile(outputFolder, 'Assigned_17.mat'), 'assigned17_group', 'assigned17_subjects', 'assigned17_all', 'info', '-v7.3');

fprintf('All done. Outputs in %s\n', outputFolder);
fprintf('Rows without a group label: %d (7 networks), %d (17 networks)\n', sum(assigned7_group == 0), sum(assigned17_group == 0));
