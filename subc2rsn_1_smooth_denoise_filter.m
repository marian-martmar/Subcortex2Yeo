%% SUBC2RSN | STEP 1 - Spatial smoothing, Wishart denoising and temporal filtering of the voxel time series
% -------------------------------------------------------------------------
% PURPOSE
%   For every subject: read the 4-D MNI volume, smooth every frame (Gaussian,FWHM = 6 mm),
%   extract the time series of the voxels of interest and, for each,
%   denoise them with Wishart filter (icaDim) ---> xRL_dn(time,voxels)
%   then (opt:z-score and) low-pass filter (2nd-order Butterworth, 0.08 Hz, zero-phase)
%   ---> save one file per subject, xRL_dn_f(time,voxels)

% INPUTS (Folder : <P.outputs_0>) ---> Files:
%   Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17(_wCereb).mat  -> subjectIDs, idcs, subset_indices
%   <ID>_rfMRI_REST1_LR_hp2000_clean.nii.gz                              one per subject, 4-D, 2 mm
%
% OUTPUTS (Folder : <P.dn_f_out> ---> Files: <ID>_xRL_dn_f.mat)
%       xRL_dn_f                      [Nt x Nf] smoothed + denoised + (opt: z-scored) + low-passed
%       subset_indices.xRL_dn_f       [Nf x 1]  POSITIONS of the Nf filtered voxels WITHIN the linear idcs 
%       idcs_filt                     [Nf x 1]  = idcs(subset_indices.xRL_dn_f), linear indices of Nf in the volume
%       i_tempFilt, j_tempFilt, k_tempFilt        array subscripts of the same voxels
%
% KEY VARIABLES
%   subset_indices.xRL_dn             positions within idcs of every voxel that is smoothed + denoised
%                                     Cortex of Schaefer/AAL/DBS80/DK84 + all subcortex ---> Columns of xRL and xRL_dn
%   subset_indices.xRL_dn_f           subset of positions within idcs that is temporally filtered and saved
%                                     Schaefer cortex + all subcortex                   ---> Columns of xRL_dn_f
%   pos_filt_within_dn                Position where voxels in xRL_dn_f sit within the columns of xRL_dn
%
% PARAMETERS set in subc2rsn_config.m
%   FWHM = 6 mm, voxelsize = 2 mm (sigma = FWHM/voxelsize/2.355 voxels), TR = 0.72 s
%   Temporal filter in use: butter(2, 0.08 Hz / Nyquist, 'low')  ("Yeo 2011").
%   flp = 0.008 Hz is defined but NOT used (the band-pass line is commented out).
%   Wishart filter: icaDim(xRL', DEMDT=1, VN=1, Iterate=2, NDist=2)
%
% NOTES
%   1. The first/last samples affected by filtfilt are NOT discarded.
%      If we wanted, should be discarded in step 2 (cluster) before computing FC

% Author: Marian Martinez-Marin. v1 December 2025. v2 November 2026
% -------------------------------------------------------------------------
clear all; clc;

P = subc2rsn_paths();
C = subc2rsn_config();

EXCLUDE_CEREB = C.EXCLUDE_CEREB   % = true;
DETRND_PreFilt = C.DETRND_PreFilt % = true;  % switch whether detrend(xRL_dn)
ZSCORE_PreFilt = C.ZSCORE_PreFilt % = false; % switch whether zscore(detrend(xRL_dn))
SAVE_xRL_dn = C.SAVE_xRL_dn       % = false; % switch whether xRL_dn are saved in external disc (or only xRL_dn_f)

%% 1.0. Data and parameters Filter MNI timeseries
addpath(genpath(P.functions_tian))
FWHM = C.FWHM; % rfMRI spatial Gaussian smoothing (6mm for T3 images)
voxelsize = C.voxelsize; % mm
%  Parameters rfMRI butterworth temporal filter
TR = C.TR;   % HCP 100unrelated is 0.72;  % Repetition Time (seconds)
fnq = C.fnq; % 1/(2*TR);                  % Nyquist frequency
flp = C.flp; % 0.008;                     % lowpass frequency of filter (Hz)
fhi = C.fhi; % 0.08;                      % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn(2),'low');   % Yeo 2011 temporal filter. Old Deco Lab: % [bfilt,afilt]=butter(k,Wn);
% discard_t_filt = C.discard_t_filt;   WARNING: keep all points for now, discard for FC

% Data from external drive
addpath(P.hcp_external) % Where timeseries are
% Define output folder (create it if it doesn't exist)
outputFolderdn = P.dn_out;
if ~exist(outputFolderdn, 'dir')
    mkdir(outputFolderdn);
end
outputFolderdn_f = P.dn_f_out;
if ~exist(outputFolderdn_f, 'dir')
    mkdir(outputFolderdn_f);
end

% Load subject IDs +  mm2(:) indices
%                    'idcs' all used voxels (schaefer, DBS80 cort+subc, Tian S1)
%                    'subset_indices': idx within idcs for each parcellation
folder_in = P.outputs_0;
if EXCLUDE_CEREB
    baseFileName = 'Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17.mat';
    fullFileName = fullfile(folder_in, baseFileName);
elseif EXCLUDE_CEREB == false
    baseFileName = 'Variables_for_FC_AAL_DK84_DBS80_Tian_toYeo_7_and_17_wCereb.mat';
    fullFileName = fullfile(folder_in, baseFileName);
end
load(fullFileName);

NSUB = length(subjectIDs)

% Voxels that will be considered (CORT: Schaefer, SUBC: all atlases)
Ncort_all = length(unique([subset_indices.cort_schaefer; subset_indices.cort_aal; subset_indices.cort_dbs80; subset_indices.cort_dk84])); % Num cort VOXELS (all atlases included)
Ncort_S = length(subset_indices.cort_schaefer); % Num cort VOXELS (just schaefer 1000)
Nsubc_all = length(unique([subset_indices.subc_tian.S1; subset_indices.subc_aal; subset_indices.subc_dbs80; subset_indices.subc_dk84])); % Num subc VOXELS (DBS80 U Tian S1 U AAL)
% WARNING: Ncort_ALL+Nsubc_ALL > |idcs| because some voxels are Cort/Subc depending on atlas
% TODO: set belonging to any subcortical as winner for FC(subc,cortical) purposes
disp(['WARNING: Ncort_ALL+Nsubc_ALL > |idcs|\n Overlapping voxels :',num2str(Ncort_all+Nsubc_all-length(idcs)),'\n']);
disp('-------> TODO: should they count for RSN-average FC?\n\n')

% Positions of these voxels within the linear array 'idcs' of atlas-voxels in 3d arrays
% Denoising will use voxels in any cortical atlas
subset_indices.xRL_dn  = unique([subset_indices.cort_schaefer; subset_indices.cort_aal; ...
                   subset_indices.cort_dbs80;  subset_indices.cort_dk84; ...
                   subset_indices.subc_tian.S1; subset_indices.subc_aal; ...
                   subset_indices.subc_dbs80;  subset_indices.subc_dk84]);
% Voxel time series (with temporal filter) will only be for: 
subset_indices.xRL_dn_f = unique([subset_indices.cort_schaefer; ...   % Schaefer cortex
                   subset_indices.subc_tian.S1; ...                   % + all subcortical
                   subset_indices.subc_aal; ...
                   subset_indices.subc_dbs80; ...
                   subset_indices.subc_dk84]);
[tf, pos_filt_within_dn] = ismember(subset_indices.xRL_dn_f, subset_indices.xRL_dn);

% 3d MNI152 array linear 'idcs' for xRL_dn
idcs_dn = idcs(subset_indices.xRL_dn);
% 3d MNI152 array linear 'idcs' for xRL_dn_f
idcs_filt = idcs(subset_indices.xRL_dn_f);
% Get ijk for both sets (example subject to get the 3D size)
info0 = niftiinfo([subjectIDs{1} '_rfMRI_REST1_LR_hp2000_clean.nii.gz']);
vol_size = info0.ImageSize(1:3);   % [Ni x Nj x Nk] = [91 x 109 x 91]
[i_dn,  j_dn,  k_dn]  = ind2sub(vol_size, idcs_dn);
[i_dn_f, j_dn_f, k_dn_f] = ind2sub(vol_size, idcs_filt);

% We save these so one can check which voxels were used at each step
meta = struct('vol_size',vol_size, ...
              'idcs_all',idcs_dn,'idcs_filt',idcs_filt, ...
              'subset_indices',subset_indices,'pos_filt_within_dn',pos_filt_within_dn, ...
              'ijk_all',[i_dn j_dn k_dn], 'ijk_filt',[i_dn_f j_dn_f k_dn_f]);
output_filename = fullfile(outputFolderdn_f,'Variables_From_1_xRL_dn_f.mat');

if ~isfile(output_filename)
    save(output_filename,'idcs','subjectIDs','subset_indices','NSUB','Ncort_all','Ncort_S','Nsubc_all','idcs_dn','idcs_filt','meta','EXCLUDE_CEREB','C')
    fprintf('File saved: %s\n', output_filename);
else
    load(output_filename);
    fprintf('File already exists\nLoaded: %s\n', output_filename);
end

%% 1.1. Filter MNI timeseries
for sub = 1:NSUB
    fprintf(['Processing subject ', num2str(sub),'\n']); % subjectIDs{sub});

    % Build the filename for the current subject's NIfTI file
    filename_load = [subjectIDs{sub}, '_rfMRI_REST1_LR_hp2000_clean.nii.gz']
    % Read 4D volume
    ts_vol = niftiread(filename_load);
    Nt = size(ts_vol, 4); % Number of timepoints
    % Preallocate time series matrix for the selected voxels (using idcs)
    xRL = zeros(Nt, length(idcs_dn));

    fprintf('Spatial smoothing with FWHM=%d mm...\n', FWHM);
    frst = 0;
    for t = 1:Nt
        % Apply 3D Gaussian smoothing
        ts_vol(:,:,:,t) = imgaussfilt3(ts_vol(:,:,:,t), FWHM/voxelsize/2.355);
        % Extract the time series for voxels in 'idcs'
        tmp = ts_vol(:,:,:,t);
        xRL(t,:) = tmp(idcs_dn);
        % Display progress
        show_progress(t, Nt, frst);
        frst = 1;
    end
    clear ts_vol;  % Free memory after extracting time series

    % Apply Wishart filtering (Glasset et al. 2016) via icaDim
    fprintf('Applying Wishart filtering...\n');
    DEMDT = 1; VN = 1; Iterate = 2; NDist = 2;
    Out = icaDim(xRL', DEMDT, VN, Iterate, NDist);
    xRL_dn = Out.data';  % Transpose back: now size is [time x voxels]
    clear Out;

    % Save dn (smoothed, denoised, NO z-score NOR temporal filter)
    if SAVE_xRL_dn == true
        output_filename = fullfile(outputFolderdn, [subjectIDs{sub}, '_xRL_dn.mat']);
        fprintf('Saving xRL_dn data to %s\n', output_filename);
        save(output_filename, 'xRL_dn','subset_indices','idcs_dn','i_dn','j_dn','k_dn','-v7.3');
    end

    % Demean and standardize (z-score) the time series
    if DETRND_PreFilt == true
        xRL_dn = detrend(xRL_dn, 'constant'); fprintf('Demeaning...\n');
    end
    if ZSCORE_PreFilt == true
        xRL_dn = xRL_dn ./ repmat(std(xRL_dn), Nt, 1); fprintf('Standarizing...\n');
    end

    % Temporal filtering using filtfilt
    fprintf('Applying temporal filtering...\n');
    Nf = numel(pos_filt_within_dn);
    xRL_dn_f      = zeros(Nt, Nf, 'like', xRL_dn);  % compact subset

    for seed = 1:Nf
        if mod(100 * seed / Nf, 10) == 0
            fprintf('%.0f%% filtering\n', 100 * seed / Nf);
        end
        col = pos_filt_within_dn(seed);             % which column in xRL_dn
        y   = xRL_dn(:, col);
        y   = filtfilt(bfilt, afilt, y);            % your low-pass or band-pass
        xRL_dn_f(:, seed) = y;
    end
    % We could remove extreme time points that may be distorted filtering
    % TODO: in next step (2) so save the full recording
    % xRL_dn_f = xRL_dn_f(C.discard_t_filt:end-C.discard_t_filt, :);

    output_filename = fullfile(outputFolderdn_f, [subjectIDs{sub}, '_xRL_dn_f.mat']);
    fprintf('Saving xRL_dn_f data to %s\n', output_filename);
    save(output_filename, 'xRL_dn_f','subset_indices','idcs_filt','i_dn_f','j_dn_f','k_dn_f', ...
     'C','EXCLUDE_CEREB','DETRND_PreFilt','ZSCORE_PreFilt','bfilt','afilt','-v7.3');
    % Clear large variables for the next subject.
    clear xRL xRL_dn xRL_dn_f;
end

% Saved filtered time series for the current subject into SEPARATE MAT files
%    xRL_dn  : Denoised, detrended voxel timeseries (from any atlas)
%    xRL_dn_f: Temporally filtered voxel timeseries (from Schaefer + any subc. atlas)
